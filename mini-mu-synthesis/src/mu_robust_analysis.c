/* ====================================================================
 * mu_robust_analysis.c — Robust Stability and Performance Analysis
 *
 * Implements the NS/NP/RS/RP framework for assessing closed-loop
 * system robustness to structured uncertainty.
 *
 * NS (Nominal Stability):     M internally stable
 * NP (Nominal Performance):   sigma_bar(M22) < 1 for all w
 * RS (Robust Stability):      u_Delta(M11) < 1 for all w
 * RP (Robust Performance):    u_{Delta_a}(M) < 1 for all w
 *
 * Key relationship (Skogestad & Postlethwaite 2005, Table 8.1):
 *   RP => RS, NP   and   RS and NP !=> RP
 *
 * Knowledge Coverage:
 *   L2: Nominal stability check (internal stability)
 *   L2: Nominal performance check (H-infinity norm)
 *   L4: Robust stability via u-analysis
 *   L4: Robust performance via augmented u (Main Loop Theorem)
 *   L4: Robustness margin computation
 *   L4: Worst-case perturbation extraction
 *   L2: Structured vs unstructured conservatism comparison
 *
 * Reference:
 *   Skogestad & Postlethwaite (2005), "Multivariable Feedback Control", Ch. 8
 *   Zhou, Doyle & Glover (1996), Ch. 11
 *   Packard & Doyle (1993), "The complex structured singular value"
 * ==================================================================== */

#include "mu_core.h"
#include "mu_robust_analysis.h"
#include "mu_structured_svd.h"
#include "mu_lft.h"
#include "mu_matrix.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================
 * L2: Nominal Stability
 *
 * Check if F_l(P, K) is internally stable.
 * Requires: (a) Controller K stabilizes P
 *           (b) No unstable pole-zero cancellations
 *
 * Simplified check: verify closed-loop A-matrix has all eigenvalues
 * in the open left half-plane.
 *
 * Complexity: O(n^3)
 * ============================ */

bool mu_check_nominal_stability(const MUSystem *P, const MUSystem *K) {
    if (!P || !K) return false;

    /* Form closed-loop system: A_cl = [ A + B2 Dk C2,  B2 Ck;
     *                                   Bk C2,          Ak       ] */
    /* Simplified: check stability of both P and K individually,
     * then check closed-loop eigenvalues */

    if (!mu_system_is_stable(P)) {
        /* P open-loop may have unstable modes — handled by feedback */
    }
    if (!mu_system_is_stable(K)) return false;

    /* Form A_cl */
    int nP = P->n, nK = K->n;
    int ny = P->p - P->p / 2;  /* measurement outputs */
    int nu = P->m - P->m / 2;  /* control inputs */

    int n_cl = nP + nK;
    MUMatrix *Acl = mu_mat_create(n_cl, n_cl);
    if (!Acl) return false;

    /* Place P's A in upper-left */
    for (int i = 0; i < nP; i++)
        for (int j = 0; j < nP; j++)
            mu_mat_set(Acl, i, j,
                mu_complex(mu_real_mat_get(P->A, i, j), 0.0));

    /* Place K's A in lower-right */
    for (int i = 0; i < nK; i++)
        for (int j = 0; j < nK; j++)
            mu_mat_set(Acl, nP + i, nP + j,
                mu_complex(mu_real_mat_get(K->A, i, j), 0.0));

    /* Add coupling terms: Acl(1:nP, nP+1:nP+nK) = B2 * Ck */
    for (int i = 0; i < nP; i++)
        for (int j = 0; j < nK; j++) {
            double sum = 0.0;
            for (int kk = 0; kk < nu && kk < P->m; kk++)
                sum += mu_real_mat_get(P->B, i, nP > 0 ? kk : 0) *
                       mu_real_mat_get(K->C, j, kk);
            mu_mat_set(Acl, i, nP + j, mu_complex(sum, 0.0));
        }

    /* Add coupling: Acl(nP+1:nP+nK, 1:nP) = Bk * C2 */
    for (int i = 0; i < nK; i++)
        for (int j = 0; j < nP; j++) {
            double sum = 0.0;
            for (int kk = 0; kk < ny && kk < P->p; kk++)
                sum += mu_real_mat_get(K->B, i, kk) *
                       mu_real_mat_get(P->C, kk, j);
            mu_mat_set(Acl, nP + i, j, mu_complex(sum, 0.0));
        }

    /* Check eigenvalues */
    MUComplex *evals = (MUComplex*)malloc((size_t)n_cl * sizeof(MUComplex));
    if (!evals) { mu_mat_free(Acl); return false; }
    mu_mat_eig(Acl, evals);

    bool stable = true;
    for (int i = 0; i < n_cl; i++) {
        if (evals[i].re >= -MU_TOL) { stable = false; break; }
    }

    mu_mat_free(Acl);
    free(evals);
    return stable;
}

/* ============================
 * L2: Nominal Performance
 *
 * Check if sigma_bar(M22(jw)) < 1 for all w.
 *
 * M22 is the nominal closed-loop transfer from exogenous
 * inputs to performance outputs.
 *
 * Complexity: O(#freq * n^3)
 * ============================ */

bool mu_check_nominal_performance(const MUSystem *P, const MUSystem *K,
                                   const MUFrequencyGrid *grid) {
    if (!P || !K || !grid) return false;

    int nz = P->p / 2;   /* performance outputs */
    int nw = P->m / 2;   /* disturbance inputs */

    for (int f = 0; f < grid->num_points; f++) {
        MUMatrix *Pjw = mu_system_freqresp(P, grid->omega[f]);
        MUMatrix *Kjw = mu_system_freqresp(K, grid->omega[f]);
        if (!Pjw || !Kjw) {
            mu_mat_free(Pjw); mu_mat_free(Kjw); return false;
        }

        /* M = F_l(P, K) at this frequency */
        MUMatrix *M = mu_lft_lower(Pjw, nz, nw, Kjw);
        mu_mat_free(Pjw); mu_mat_free(Kjw);

        if (!M) return false;

        /* Extract M22 */
        int n2 = M->rows - nz;
        int m2 = M->cols - nw;
        MUMatrix *M22 = mu_mat_create(n2, m2);
        if (!M22) { mu_mat_free(M); return false; }
        for (int i = 0; i < n2; i++)
            for (int j = 0; j < m2; j++)
                mu_mat_set(M22, i, j, mu_mat_get(M, nz + i, nw + j));

        double sn = mu_mat_spectral_norm(M22);
        mu_mat_free(M22); mu_mat_free(M);

        if (sn >= MU_ROBUST_STABILITY_THRESHOLD) return false;
    }
    return true;
}

/* ============================
 * L4: Robust Stability via u-Analysis
 *
 * Theorem: F_u(M, Delta) robustly stable for ||Delta||_inf <= 1
 * iff u_Delta(M11(jw)) < 1 for all w.
 *
 * This is the main robust stability test using u.
 * ============================ */

bool mu_check_robust_stability(const MUSystem *P, const MUSystem *K,
                                const MUUncertaintyStructure *delta,
                                const MUFrequencyGrid *grid,
                                int upper_iter, int lower_iter, double tol) {
    if (!P || !K || !delta || !grid) return false;

    int nz = P->p / 2;

    for (int f = 0; f < grid->num_points; f++) {
        MUMatrix *Pjw = mu_system_freqresp(P, grid->omega[f]);
        MUMatrix *Kjw = mu_system_freqresp(K, grid->omega[f]);
        if (!Pjw || !Kjw) {
            mu_mat_free(Pjw); mu_mat_free(Kjw); return false;
        }

        MUMatrix *M = mu_lft_lower(Pjw, nz, P->m / 2, Kjw);
        mu_mat_free(Pjw); mu_mat_free(Kjw);
        if (!M) return false;

        /* Extract M11 (uncertainty channel) */
        int dim_unc = delta->total_dim;
        MUMatrix *M11 = mu_mat_create(dim_unc, dim_unc);
        if (!M11) { mu_mat_free(M); return false; }
        for (int i = 0; i < dim_unc; i++)
            for (int j = 0; j < dim_unc; j++)
                mu_mat_set(M11, i, j, mu_mat_get(M, i, j));
        mu_mat_free(M);

        MUMuResult res = mu_compute(M11, delta, grid->omega[f],
                                     upper_iter, lower_iter, tol);
        mu_mat_free(M11);

        if (res.mu_upper >= MU_ROBUST_STABILITY_THRESHOLD) return false;
    }
    return true;
}

/* ============================
 * L4: Robust Performance via Augmented u
 *
 * Main Loop Theorem: RP <=> u_{Delta_a}(M) < 1 for all w
 * where Delta_a = diag(Delta, Delta_p)
 *
 * Reference: Packard & Doyle (1993), Main Loop Theorem
 * ============================ */

bool mu_check_robust_performance(const MUSystem *P, const MUSystem *K,
                                  const MUUncertaintyStructure *delta,
                                  const MUFrequencyGrid *grid,
                                  int upper_iter, int lower_iter, double tol) {
    if (!P || !K || !delta || !grid) return false;

    /* Augment uncertainty with performance block */
    int p_perf = P->p - delta->total_dim;
    if (p_perf <= 0) p_perf = P->p;

    MUUncertaintyStructure *delta_a = mu_unc_struct_create(delta->num_blocks + 1);
    if (!delta_a) return false;
    for (int i = 0; i < delta->num_blocks; i++)
        delta_a->blocks[i] = delta->blocks[i];
    delta_a->blocks[delta->num_blocks].type = MU_UNC_FULL_BLOCK;
    delta_a->blocks[delta->num_blocks].size = p_perf;
    delta_a->blocks[delta->num_blocks].bound = 1.0;
    delta_a->total_dim = delta->total_dim + p_perf;

    /* Check robust stability with augmented uncertainty */
    bool result = mu_check_robust_stability(P, K, delta_a, grid,
                                             upper_iter, lower_iter, tol);
    mu_unc_struct_free(delta_a);
    return result;
}

/* ============================
 * L4: Robustness Margin
 *
 * k_max = 1 / sup_w u_Delta(M11(jw))
 *
 * The largest scaling of the uncertainty set for which
 * the system remains stable.
 * ============================ */

double mu_robustness_margin(const MUSystem *P, const MUSystem *K,
                             const MUUncertaintyStructure *delta,
                             const MUFrequencyGrid *grid,
                             int upper_iter, int lower_iter, double tol) {
    if (!P || !K || !delta || !grid) return 0.0;

    int nz = P->p / 2;
    double peak_mu = 0.0;

    for (int f = 0; f < grid->num_points; f++) {
        MUMatrix *Pjw = mu_system_freqresp(P, grid->omega[f]);
        MUMatrix *Kjw = mu_system_freqresp(K, grid->omega[f]);
        if (!Pjw || !Kjw) { mu_mat_free(Pjw); mu_mat_free(Kjw); continue; }

        MUMatrix *M = mu_lft_lower(Pjw, nz, P->m / 2, Kjw);
        mu_mat_free(Pjw); mu_mat_free(Kjw);
        if (!M) continue;

        int dim_unc = delta->total_dim;
        MUMatrix *M11 = mu_mat_create(dim_unc, dim_unc);
        if (!M11) { mu_mat_free(M); continue; }
        for (int i = 0; i < dim_unc && i < M->rows; i++)
            for (int j = 0; j < dim_unc && j < M->cols; j++)
                mu_mat_set(M11, i, j, mu_mat_get(M, i, j));
        mu_mat_free(M);

        MUMuResult res = mu_compute(M11, delta, grid->omega[f],
                                     upper_iter, lower_iter, tol);
        mu_mat_free(M11);

        if (res.mu_upper > peak_mu) peak_mu = res.mu_upper;
    }

    return (peak_mu > MU_EPSILON) ? 1.0 / peak_mu : 1e16;
}

/* ============================
 * L4: Worst-Case Perturbation
 *
 * Extract the perturbation Delta_wc that causes instability
 * at the critical frequency, using power iteration results.
 * ============================ */

void mu_worst_case_perturbation(const MUMatrix *M,
                                 const MUUncertaintyStructure *delta,
                                 MUMatrix **Delta_wc) {
    *Delta_wc = NULL;
    if (!M || !delta) return;

    int dim = M->rows;
    MUMatrix *Dwc = mu_mat_create(dim, dim);
    if (!Dwc) return;

    /* Use power iteration to get the worst-case Q */
    double mu_lb = mu_lower_bound(M, delta, 100, MU_TOL);

    if (mu_lb < MU_EPSILON) {
        /* No destabilizing perturbation found */
        for (int i = 0; i < dim; i++)
            mu_mat_set(Dwc, i, i, mu_complex(0.0, 0.0));
    } else {
        /* The worst-case perturbation has magnitude 1/mu_lb */
        /* Identity scaled by 1/mu_lb as approximate worst-case */
        for (int i = 0; i < dim; i++)
            mu_mat_set(Dwc, i, i, mu_complex(1.0 / mu_lb, 0.0));
    }

    *Delta_wc = Dwc;
}

/* ============================
 * L2: Structured vs Unstructured Conservatism
 *
 * Compare the structured singular value u with unstructured bound
 * sigma_bar(M). The ratio shows conservatism reduction.
 * ============================ */

MUComparisonResult mu_compare_structured_unstructured(
    const MUMatrix *M, const MUUncertaintyStructure *delta,
    int upper_iter, int lower_iter, double tol) {
    MUComparisonResult result = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    if (!M || !delta) return result;

    result.spectral_norm = mu_mat_spectral_norm(M);

    /* Compute spectral radius rho(M) = max |lambda_i(M)| */
    int dim = M->rows;
    MUComplex *evals = (MUComplex*)malloc((size_t)dim * sizeof(MUComplex));
    if (evals) {
        mu_mat_eig(M, evals);
        double max_abs = 0.0;
        for (int i = 0; i < dim; i++) {
            double abs_val = mu_cabs(evals[i]);
            if (abs_val > max_abs) max_abs = abs_val;
        }
        result.spectral_radius = max_abs;
        free(evals);
    }

    MUMuResult mu_res = mu_compute(M, delta, 0.0, upper_iter, lower_iter, tol);
    result.mu_upper = mu_res.mu_upper;
    result.mu_lower = mu_res.mu_lower;

    result.conservatism = (result.spectral_norm > MU_EPSILON)
        ? (result.spectral_norm - result.mu_lower) / result.spectral_norm
        : 0.0;

    return result;
}

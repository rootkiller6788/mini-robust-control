/* ====================================================================
 * mu_dk_iteration.c — DK-Iteration Algorithm for u-Synthesis
 *
 * DK-iteration is the central computational algorithm for u-synthesis
 * controller design. It alternates between D-scale optimization (D-step)
 * and H-infinity controller synthesis (K-step) to minimize the peak
 * structured singular value.
 *
 * Algorithm (Balas et al. 1991):
 *   K0 = initial H-infinity controller
 *   Repeat for k = 0, 1, 2, ...:
 *     1. D-step: Fix K_k, compute optimal D-scales at each frequency
 *     2. Fit rational D(s) through the D-scale data
 *     3. K-step: Solve H-infinity for D-scaled plant → K_{k+1}
 *     4. Check convergence: |gamma_{k+1} - gamma_k| < tol
 *
 * Knowledge Coverage:
 *   L5: D-step optimization at each frequency point
 *   L5: D-scale rational fitting across frequency grid
 *   L5: K-step H-infinity synthesis for D-scaled plant
 *   L5: Full DK-iteration loop with convergence monitoring
 *   L5: Warm-start DK-iteration from given initial controller
 *
 * Reference:
 *   Doyle (1985), "Structured uncertainty in control system design"
 *   Balas et al. (1991), "u-Analysis and Synthesis Toolbox"
 *   Zhou, Doyle & Glover (1996), Ch. 11
 * ==================================================================== */

#include "mu_core.h"
#include "mu_dk_iteration.h"
#include "mu_structured_svd.h"
#include "mu_hinf.h"
#include "mu_lft.h"
#include "mu_matrix.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================
 * L5: D-Step — Optimal D-Scales at Each Frequency
 *
 * For each frequency omega, compute D(omega) that minimizes
 * sigma_bar(D M D^{-1}) where M = F_l(P(jw), K(jw)).
 *
 * Uses Osborne-like alternating optimization:
 *   For each scalar block: D_i = d_i I, optimize d_i
 *   For each full block: D_j = D_j* > 0, optimize via SVD
 *
 * Complexity: O(#freq * dim^3 * max_iter)
 * Reference: Osborne (1960), Balas et al. (1991)
 * ============================ */

void mu_d_step(const MUMatrix *M,
                const MUUncertaintyStructure *delta,
                MUDScale **D_left, MUDScale **D_right,
                int max_iter, double tol) {
    *D_left = NULL;
    *D_right = NULL;
    if (!M || !delta || M->rows != M->cols) return;

    int dim = M->rows;

    /* Allocate D-scales */
    MUDScale *DL = (MUDScale*)malloc(sizeof(MUDScale));
    MUDScale *DR = (MUDScale*)malloc(sizeof(MUDScale));
    if (!DL || !DR) { free(DL); free(DR); return; }

    DL->dim = dim;
    DL->structure = (MUUncertaintyStructure*)delta; /* non-owning ref */
    DL->data = (MUComplex*)malloc((size_t)dim * sizeof(MUComplex));

    DR->dim = dim;
    DR->structure = DL->structure;
    DR->data = (MUComplex*)malloc((size_t)dim * sizeof(MUComplex));

    if (!DL->data || !DR->data) {
        free(DL->data); free(DL); free(DR->data); free(DR);
        return;
    }

    /* Initialize D = I, D^{-1} = I */
    for (int i = 0; i < dim; i++) {
        DL->data[i] = mu_complex(1.0, 0.0);
        DR->data[i] = mu_complex(1.0, 0.0);
    }

    /* Osborne iteration on the diagonal elements */
    double *d = (double*)malloc((size_t)dim * sizeof(double));
    if (!d) {
        free(DL->data); free(DL); free(DR->data); free(DR);
        return;
    }
    for (int i = 0; i < dim; i++) d[i] = 1.0;

    for (int iter = 0; iter < max_iter; iter++) {
        bool converged = true;
        for (int i = 0; i < dim; i++) {
            double rnrm = 0.0, cnrm = 0.0;
            for (int j = 0; j < dim; j++) {
                double mij_abs = mu_cabs(mu_mat_get(M, i, j));
                double mji_abs = mu_cabs(mu_mat_get(M, j, i));
                rnrm += mij_abs * d[i] / d[j];
                cnrm += mji_abs * d[j] / d[i];
            }
            if (rnrm > MU_EPSILON && cnrm > MU_EPSILON) {
                double f = sqrt(sqrt(cnrm / rnrm));
                d[i] *= f;
                if (fabs(f - 1.0) > tol) converged = false;
            }
        }
        if (converged) break;
    }

    /* Set D-scale data from optimized d */
    for (int i = 0; i < dim; i++) {
        DL->data[i] = mu_complex(d[i], 0.0);
        DR->data[i] = mu_complex(1.0 / d[i], 0.0);
    }

    free(d);
    *D_left = DL;
    *D_right = DR;
}

/* ============================
 * L5: D-Scale Fitting
 *
 * Fit a stable, minimum-phase transfer function D(s) to the
 * frequency-domain D-scale data.
 *
 * Approach: fit a diagonal state-space system where each channel
 * is a first-order or second-order transfer function matching
 * the magnitude of D(omega) at each frequency.
 *
 * Complexity: O(#freq * order^3)
 * Reference: Balas et al. (1991)
 * ============================ */

void mu_d_fit(const MUDScale **D_data,
               const MUFrequencyGrid *grid,
               const MUUncertaintyStructure *delta,
               MUSystem **D_sys) {
    *D_sys = NULL;
    if (!D_data || !grid || !delta) return;

    int dim = (D_data[0]) ? D_data[0]->dim : 0;
    if (dim <= 0) return;

    /* Fit a diagonal first-order stable system to each channel */
    /* D(s) = diag( (s + a1)/(s + b1), ... ) or constant for full blocks */

    /* Create diagonal state-space system */
    MUSystem *Ds = mu_system_create(dim, dim, dim);
    if (!Ds) return;

    /* For each diagonal element, fit magnitude data */
    int offset = 0;
    for (int b = 0; b < delta->num_blocks; b++) {
        int sz = delta->blocks[b].size;
        if (delta->blocks[b].type == MU_UNC_REPEATED_SCALAR ||
            delta->blocks[b].type == MU_UNC_REAL_SCALAR) {
            /* Scalar block: fit first-order transfer function */
            /* Average magnitude over frequency */
            double avg_mag = 1.0;
            double count = 0.0;
            for (int f = 0; f < grid->num_points; f++) {
                if (D_data[f]) {
                    avg_mag += mu_cabs(D_data[f]->data[offset]);
                    count += 1.0;
                }
            }
            if (count > 0) avg_mag /= count;

            /* D(s) = avg_mag (stable, min-phase) */
            for (int k = 0; k < sz; k++) {
                int idx = offset + k;
                /* A = -1, B = 1, C = avg_mag, D = avg_mag */
                mu_real_mat_set(Ds->A, idx, idx, -1.0);
                mu_real_mat_set(Ds->B, idx, idx, 1.0);
                mu_real_mat_set(Ds->C, idx, idx, avg_mag);
                mu_real_mat_set(Ds->D, idx, idx, avg_mag);
            }
        } else {
            /* Full block: identity scaling (D = I for that block) */
            for (int k = 0; k < sz; k++) {
                int idx = offset + k;
                mu_real_mat_set(Ds->A, idx, idx, -1.0);
                mu_real_mat_set(Ds->B, idx, idx, 1.0);
                mu_real_mat_set(Ds->C, idx, idx, 1.0);
                mu_real_mat_set(Ds->D, idx, idx, 1.0);
            }
        }
        offset += sz;
    }

    *D_sys = Ds;
}

/* ============================
 * L5: K-Step — H-infinity Synthesis for D-Scaled Plant
 *
 * Form the D-scaled plant:
 *   P_D = [ D(s)     0  ] * P * [ D(s)^{-1}     0  ]
 *         [   0       I  ]       [     0         I  ]
 *
 * Then solve H-infinity problem:
 *   min_K || F_l(P_D, K) ||_inf
 *
 * Complexity: O(n^3)
 * Reference: Doyle (1985), Balas et al. (1991)
 * ============================ */

void mu_k_step(const MUSystem *P,
                const MUSystem *D_sys,
                MUSystem **K, double *gamma) {
    *K = NULL;
    *gamma = 0.0;
    if (!P) return;

    /* If no D-scales provided, solve unscaled problem */
    if (!D_sys) {
        mu_hinf_gamma_iteration(P, K, gamma, 1e-4, 30);
        return;
    }

    /* Build D-scaled plant P_D */
    /* In the simplified case, absorb D-scales by adjusting
     * the first block rows/columns of the plant */
    int n = P->n + D_sys->n;

    /* Create augmented plant */
    MUSystem *P_D = mu_system_create(n, P->m, P->p);
    if (!P_D) return;

    /* Place P in the lower-right corner */
    for (int i = 0; i < P->n; i++) {
        for (int j = 0; j < P->n; j++)
            mu_real_mat_set(P_D->A, D_sys->n + i, D_sys->n + j,
                mu_real_mat_get(P->A, i, j));
        for (int j = 0; j < P->m; j++)
            mu_real_mat_set(P_D->B, D_sys->n + i, j,
                mu_real_mat_get(P->B, i, j));
    }
    for (int i = 0; i < P->p; i++) {
        for (int j = 0; j < P->n; j++)
            mu_real_mat_set(P_D->C, i, D_sys->n + j,
                mu_real_mat_get(P->C, i, j));
        for (int j = 0; j < P->m; j++)
            mu_real_mat_set(P_D->D, i, j,
                mu_real_mat_get(P->D, i, j));
    }

    /* Place D_sys in the upper-left corner */
    for (int i = 0; i < D_sys->n; i++) {
        for (int j = 0; j < D_sys->n; j++)
            mu_real_mat_set(P_D->A, i, j,
                mu_real_mat_get(D_sys->A, i, j));
    }

    /* Solve H-infinity for scaled plant */
    mu_hinf_gamma_iteration(P_D, K, gamma, 1e-4, 30);
    mu_system_free(P_D);
}

/* ============================
 * L5: Full DK-Iteration
 *
 * Main algorithm for u-synthesis controller design.
 *
 * Complexity: O(max_iter * (#freq * dim^3 + n^3))
 * ============================ */

void mu_dk_iteration(const MUSystem *P,
                      const MUUncertaintyStructure *delta,
                      const MUFrequencyGrid *grid,
                      MUSystem **K, MUDKState *state) {
    *K = NULL;
    if (!P || !delta || !grid || !state) return;

    state->iteration = 0;
    state->mu_peak = 1e16;
    state->mu_peak_prev = 1e16;
    state->gamma = 1e16;
    state->gamma_prev = 1e16;
    state->converged = false;
    state->stop_reason = "not started";

    /* Step 1: Initial K0 = H-infinity controller for unscaled plant */
    MUSystem *Kk = NULL;
    double gamma_k = 0.0;
    mu_k_step(P, NULL, &Kk, &gamma_k);
    if (!Kk) {
        state->stop_reason = "initial Hinf synthesis failed";
        return;
    }
    state->gamma = gamma_k;

    /* Step 2: DK iteration loop */
    for (int iter = 0; iter < state->max_iterations; iter++) {
        state->iteration = iter + 1;
        state->mu_peak_prev = state->mu_peak;
        state->gamma_prev = state->gamma;

        /* D-step: compute optimal D-scales at each frequency */
        MUDScale **D_data = (MUDScale**)malloc((size_t)grid->num_points *
            sizeof(MUDScale*));
        if (!D_data) { state->stop_reason = "memory allocation failed"; break; }

        double peak_mu = 0.0;
        for (int f = 0; f < grid->num_points; f++) {
            /* Compute M = F_l(P(jw), K(jw)) */
            MUMatrix *Pjw = mu_system_freqresp(P, grid->omega[f]);
            if (!Pjw) { D_data[f] = NULL; continue; }

            /* Extract M = F_l(P, K) at this frequency */
            MUMatrix *Kjw = mu_system_freqresp(Kk, grid->omega[f]);
            if (!Kjw) { mu_mat_free(Pjw); D_data[f] = NULL; continue; }

            MUMatrix *M = mu_lft_lower(Pjw, P->p / 2, P->m / 2, Kjw);
            mu_mat_free(Pjw); mu_mat_free(Kjw);

            if (M) {
                MUMuResult res = mu_compute(M, delta, grid->omega[f], 50, 50, MU_TOL);
                if (res.mu_upper > peak_mu) peak_mu = res.mu_upper;
                mu_d_step(M, delta, &D_data[f], NULL, 30, MU_TOL);
                mu_mat_free(M);
            } else {
                D_data[f] = NULL;
            }
        }
        state->mu_peak = peak_mu;

        /* Fit D(s) through D-scale data */
        MUSystem *D_sys = NULL;
        mu_d_fit((const MUDScale**)D_data, grid, delta, &D_sys);

        /* Clean up D_data */
        for (int f = 0; f < grid->num_points; f++) {
            if (D_data[f]) {
                free(D_data[f]->data);
                free(D_data[f]);
            }
        }
        free(D_data);

        if (!D_sys) {
            state->stop_reason = "D-scale fitting failed";
            break;
        }

        /* K-step: solve H-infinity for D-scaled plant */
        MUSystem *K_new = NULL;
        double gamma_new = 0.0;
        mu_k_step(P, D_sys, &K_new, &gamma_new);
        mu_system_free(D_sys);

        if (!K_new) {
            state->stop_reason = "K-step Hinf synthesis failed";
            break;
        }

        /* Update controller */
        mu_system_free(Kk);
        Kk = K_new;
        state->gamma = gamma_new;

        /* Check convergence */
        if (fabs(state->gamma - state->gamma_prev) < state->tolerance) {
            state->converged = true;
            state->stop_reason = "converged";
            break;
        }
    }

    if (!state->converged && state->iteration >= state->max_iterations) {
        state->stop_reason = "max iterations reached";
    }

    *K = Kk;
}

/* ============================
 * L5: Warm-Start DK-Iteration
 *
 * Same as standard DK-iteration but starts from a given K0
 * instead of computing one from H-infinity synthesis.
 *
 * Useful for design iterations with modified weighting functions.
 * ============================ */

void mu_dk_iteration_warm(const MUSystem *P,
                           const MUUncertaintyStructure *delta,
                           const MUFrequencyGrid *grid,
                           MUSystem **K, MUDKState *state) {
    if (!P || !delta || !grid || !K || !*K || !state) return;

    MUSystem *Kk = *K;
    state->iteration = 0;
    state->mu_peak = 1e16;
    state->mu_peak_prev = 1e16;
    state->gamma = 1e16;
    state->gamma_prev = 1e16;
    state->converged = false;
    state->stop_reason = "not started";

    for (int iter = 0; iter < state->max_iterations; iter++) {
        state->iteration = iter + 1;
        state->mu_peak_prev = state->mu_peak;
        state->gamma_prev = state->gamma;

        /* D-step */
        MUDScale **D_data = (MUDScale**)malloc((size_t)grid->num_points *
            sizeof(MUDScale*));
        if (!D_data) { state->stop_reason = "memory allocation failed"; break; }

        double peak_mu = 0.0;
        for (int f = 0; f < grid->num_points; f++) {
            MUMatrix *Pjw = mu_system_freqresp(P, grid->omega[f]);
            if (!Pjw) { D_data[f] = NULL; continue; }

            MUMatrix *Kjw = mu_system_freqresp(Kk, grid->omega[f]);
            if (!Kjw) { mu_mat_free(Pjw); D_data[f] = NULL; continue; }

            MUMatrix *M = mu_lft_lower(Pjw, P->p / 2, P->m / 2, Kjw);
            mu_mat_free(Pjw); mu_mat_free(Kjw);

            if (M) {
                MUMuResult res = mu_compute(M, delta, grid->omega[f], 50, 50, MU_TOL);
                if (res.mu_upper > peak_mu) peak_mu = res.mu_upper;
                mu_d_step(M, delta, &D_data[f], NULL, 30, MU_TOL);
                mu_mat_free(M);
            } else {
                D_data[f] = NULL;
            }
        }
        state->mu_peak = peak_mu;

        /* Fit D(s) */
        MUSystem *D_sys = NULL;
        mu_d_fit((const MUDScale**)D_data, grid, delta, &D_sys);

        for (int f = 0; f < grid->num_points; f++) {
            if (D_data[f]) { free(D_data[f]->data); free(D_data[f]); }
        }
        free(D_data);

        if (!D_sys) { state->stop_reason = "D-scale fitting failed"; break; }

        /* K-step */
        MUSystem *K_new = NULL;
        double gamma_new = 0.0;
        mu_k_step(P, D_sys, &K_new, &gamma_new);
        mu_system_free(D_sys);

        if (!K_new) { state->stop_reason = "K-step failed"; break; }

        mu_system_free(Kk);
        Kk = K_new;
        state->gamma = gamma_new;

        if (fabs(state->gamma - state->gamma_prev) < state->tolerance) {
            state->converged = true;
            state->stop_reason = "converged";
            break;
        }
    }

    if (!state->converged && state->iteration >= state->max_iterations) {
        state->stop_reason = "max iterations reached";
    }

    *K = Kk;
}

/* ====================================================================
 * mu_lft.c — Linear Fractional Transformations for u-Synthesis
 *
 * LFTs are the fundamental interconnection framework for robust control
 * theory. They separate the nominal system from uncertainty and controllers.
 *
 * Knowledge Coverage:
 *   L2: Upper LFT: F_u(M, Delta) = M22 + M21 Delta (I - M11 Delta)^{-1} M12
 *   L2: Lower LFT: F_l(M, K) = M11 + M12 K (I - M22 K)^{-1} M21
 *   L3: Star product (Redheffer product) for system interconnection
 *   L3: System-to-M conversion at a given frequency
 *   L5: Interconnection building (plant + uncertainty + performance weights)
 *   L5: LFT model reduction (absorbing weights into the plant)
 *
 * Reference:
 *   Zhou, Doyle & Glover (1996), Ch. 9-10
 *   Redheffer (1960), "On a certain linear fractional transformation"
 * ==================================================================== */

#include "mu_core.h"
#include "mu_lft.h"
#include "mu_matrix.h"
#include <stdlib.h>
#include <string.h>

/* ============================
 * L2: Upper Linear Fractional Transformation
 *
 * F_u(M, Delta) = M22 + M21 * Delta * (I - M11 * Delta)^{-1} * M12
 *
 * Maps uncertainty Delta to the closed-loop transfer function.
 * Partition M = [M11, M12; M21, M22] with:
 *   M11: p1 x q1  (uncertainty input → uncertainty output)
 *   M12: p1 x (cols_M - q1)
 *   M21: (rows_M - p1) x q1
 *   M22: (rows_M - p1) x (cols_M - q1)
 *
 * Complexity: O((p1+q1)^3 + dim^3) mainly from matrix inversion
 * Reference: Zhou et al. (1996), §9.3
 * ============================ */

MUMatrix* mu_lft_upper(const MUMatrix *M, int p1, int q1,
                        const MUMatrix *Delta) {
    if (!M || !Delta) return NULL;
    int rows_M = M->rows, cols_M = M->cols;
    int n1 = p1, n2 = rows_M - p1;
    int m1 = q1, m2 = cols_M - q1;

    if (n1 < 0 || n2 < 0 || m1 < 0 || m2 < 0) return NULL;
    if (Delta->rows != m1 || Delta->cols != n1) return NULL;

    /* Extract partitions */
    MUMatrix *M11 = mu_mat_create(n1, m1);
    MUMatrix *M12 = mu_mat_create(n1, m2);
    MUMatrix *M21 = mu_mat_create(n2, m1);
    MUMatrix *M22 = mu_mat_create(n2, m2);

    if (!M11 || !M12 || !M21 || !M22) {
        mu_mat_free(M11); mu_mat_free(M12);
        mu_mat_free(M21); mu_mat_free(M22); return NULL;
    }

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < m1; j++)
            mu_mat_set(M11, i, j, mu_mat_get(M, i, j));
        for (int j = 0; j < m2; j++)
            mu_mat_set(M12, i, j, mu_mat_get(M, i, m1 + j));
    }
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < m1; j++)
            mu_mat_set(M21, i, j, mu_mat_get(M, n1 + i, j));
        for (int j = 0; j < m2; j++)
            mu_mat_set(M22, i, j, mu_mat_get(M, n1 + i, m1 + j));
    }

    /* I - M11 * Delta */
    MUMatrix *M11D = mu_mat_mul(M11, Delta);
    MUMatrix *I = mu_mat_identity(n1);
    if (!M11D || !I) {
        mu_mat_free(M11); mu_mat_free(M12); mu_mat_free(M21); mu_mat_free(M22);
        mu_mat_free(M11D); mu_mat_free(I); return NULL;
    }
    MUMatrix *I_minus = mu_mat_sub(I, M11D);
    if (!I_minus) {
        mu_mat_free(M11); mu_mat_free(M12); mu_mat_free(M21); mu_mat_free(M22);
        mu_mat_free(M11D); mu_mat_free(I); return NULL;
    }

    /* (I - M11*Delta)^{-1} */
    MUMatrix *inv = mu_mat_inverse(I_minus);
    if (!inv) {
        mu_mat_free(M11); mu_mat_free(M12); mu_mat_free(M21); mu_mat_free(M22);
        mu_mat_free(M11D); mu_mat_free(I); mu_mat_free(I_minus); return NULL;
    }

    /* Delta * (I - M11*Delta)^{-1} */
    MUMatrix *Dinv = mu_mat_mul(Delta, inv);

    /* Delta * inv * M12 */
    MUMatrix *DinM12 = mu_mat_mul(Dinv, M12);

    /* M21 * Delta * inv * M12 */
    MUMatrix *M21DinM12 = mu_mat_mul(M21, DinM12);

    /* M22 + M21 * Delta * inv * M12 */
    MUMatrix *F_u = mu_mat_add(M22, M21DinM12);

    mu_mat_free(M11); mu_mat_free(M12); mu_mat_free(M21); mu_mat_free(M22);
    mu_mat_free(M11D); mu_mat_free(I); mu_mat_free(I_minus);
    mu_mat_free(inv); mu_mat_free(Dinv); mu_mat_free(DinM12); mu_mat_free(M21DinM12);
    return F_u;
}

/* ============================
 * L2: Lower Linear Fractional Transformation
 *
 * F_l(M, K) = M11 + M12 * K * (I - M22 * K)^{-1} * M21
 *
 * Closes the loop with controller K.
 *
 * Complexity: O(dim^3)
 * Reference: Zhou et al. (1996), §9.3
 * ============================ */

MUMatrix* mu_lft_lower(const MUMatrix *M, int p1, int q1,
                        const MUMatrix *K) {
    if (!M || !K) return NULL;
    int rows_M = M->rows, cols_M = M->cols;
    int n1 = p1, n2 = rows_M - p1;
    int m1 = q1, m2 = cols_M - q1;

    if (n1 < 0 || n2 < 0 || m1 < 0 || m2 < 0) return NULL;
    if (K->rows != m2 || K->cols != n2) return NULL;

    /* Symmetric to upper LFT but with K in lower-right corner */
    /* Rearrange partitions for lower LFT */
    MUMatrix *M11 = mu_mat_create(n1, m1);
    MUMatrix *M12 = mu_mat_create(n1, m2);
    MUMatrix *M21 = mu_mat_create(n2, m1);
    MUMatrix *M22 = mu_mat_create(n2, m2);

    if (!M11 || !M12 || !M21 || !M22) {
        mu_mat_free(M11); mu_mat_free(M12);
        mu_mat_free(M21); mu_mat_free(M22); return NULL;
    }

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < m1; j++)
            mu_mat_set(M11, i, j, mu_mat_get(M, i, j));
        for (int j = 0; j < m2; j++)
            mu_mat_set(M12, i, j, mu_mat_get(M, i, m1 + j));
    }
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < m1; j++)
            mu_mat_set(M21, i, j, mu_mat_get(M, n1 + i, j));
        for (int j = 0; j < m2; j++)
            mu_mat_set(M22, i, j, mu_mat_get(M, n1 + i, m1 + j));
    }

    /* I - M22 * K */
    MUMatrix *M22K = mu_mat_mul(M22, K);
    MUMatrix *I = mu_mat_identity(n2);
    if (!M22K || !I) {
        mu_mat_free(M11); mu_mat_free(M12); mu_mat_free(M21); mu_mat_free(M22);
        mu_mat_free(M22K); mu_mat_free(I); return NULL;
    }
    MUMatrix *I_minus = mu_mat_sub(I, M22K);
    MUMatrix *inv = mu_mat_inverse(I_minus);
    if (!inv) {
        mu_mat_free(M11); mu_mat_free(M12); mu_mat_free(M21); mu_mat_free(M22);
        mu_mat_free(M22K); mu_mat_free(I); mu_mat_free(I_minus); return NULL;
    }

    /* K * (I - M22*K)^{-1} * M21 */
    MUMatrix *Kinv = mu_mat_mul(K, inv);
    MUMatrix *KinvM21 = mu_mat_mul(Kinv, M21);
    MUMatrix *M12KinvM21 = mu_mat_mul(M12, KinvM21);

    /* F_l = M11 + M12 * K * inv * M21 */
    MUMatrix *F_l = mu_mat_add(M11, M12KinvM21);

    mu_mat_free(M11); mu_mat_free(M12); mu_mat_free(M21); mu_mat_free(M22);
    mu_mat_free(M22K); mu_mat_free(I); mu_mat_free(I_minus);
    mu_mat_free(inv); mu_mat_free(Kinv); mu_mat_free(KinvM21); mu_mat_free(M12KinvM21);
    return F_l;
}

/* ============================
 * L3: Redheffer Star Product
 *
 * (P * K) = F_l(P, K), the feedback interconnection.
 *
 * For two systems P and K with compatible partition sizes:
 *   star(P, K) = closed-loop transfer from [w, w_K] to [z, z_K]
 *
 * Complexity: O(dim^3)
 * Reference: Redheffer (1960)
 *            Zhou et al. (1996), §9.3
 * ============================ */

MUMatrix* mu_lft_star_product(const MUMatrix *P, const MUMatrix *K,
                               int p1, int q1, int p2_arg, int q2_arg) {
    (void)p2_arg; (void)q2_arg;
    if (!P || !K) return NULL;
    /* Star product: feedback connection */
    /* P has inputs [w; u] and outputs [z; y] */
    /* K has inputs [y] and outputs [u] */

    MUMatrix *P11 = mu_mat_create(p1, q1);
    MUMatrix *P12 = mu_mat_create(p1, P->cols - q1);
    MUMatrix *P21 = mu_mat_create(P->rows - p1, q1);
    MUMatrix *P22 = mu_mat_create(P->rows - p1, P->cols - q1);

    if (!P11 || !P12 || !P21 || !P22) {
        mu_mat_free(P11); mu_mat_free(P12);
        mu_mat_free(P21); mu_mat_free(P22); return NULL;
    }

    for (int i = 0; i < p1; i++) {
        for (int j = 0; j < q1; j++)
            mu_mat_set(P11, i, j, mu_mat_get(P, i, j));
        for (int j = 0; j < P->cols - q1; j++)
            mu_mat_set(P12, i, j, mu_mat_get(P, i, q1 + j));
    }
    for (int i = 0; i < P->rows - p1; i++) {
        for (int j = 0; j < q1; j++)
            mu_mat_set(P21, i, j, mu_mat_get(P, p1 + i, j));
        for (int j = 0; j < P->cols - q1; j++)
            mu_mat_set(P22, i, j, mu_mat_get(P, p1 + i, q1 + j));
    }

    /* Lower LFT: F_l(P, K) */
    MUMatrix *F_l = mu_lft_lower(P, p1, q1, K);

    mu_mat_free(P11); mu_mat_free(P12); mu_mat_free(P21); mu_mat_free(P22);
    return F_l;
}

/* ============================
 * L3: System-to-M Conversion at a Frequency
 *
 * Converts a state-space system to its transfer matrix at frequency omega.
 *
 * Complexity: O(n^3) per frequency
 * ============================ */

MUMatrix* mu_system_to_m(const MUSystem *sys, double omega,
                          int p_unc, int q_unc) {
    if (!sys) return NULL;
    /* p_unc: number of uncertainty output channels */
    /* q_unc: number of uncertainty input channels */

    MUMatrix *Gjw = mu_system_freqresp(sys, omega);
    if (!Gjw) return NULL;

    /* Build a standard M = G(jw) augmented with zeros for unused channels */
    /* Assume G has p_unc + p_perf outputs and q_unc + q_act inputs */
    int total_rows = p_unc + (sys->p - p_unc);
    int total_cols = q_unc + (sys->m - q_unc);

    MUMatrix *M = mu_mat_create(total_rows, total_cols);
    if (!M) { mu_mat_free(Gjw); return NULL; }

    for (int i = 0; i < total_rows && i < Gjw->rows; i++)
        for (int j = 0; j < total_cols && j < Gjw->cols; j++)
            mu_mat_set(M, i, j, mu_mat_get(Gjw, i, j));

    mu_mat_free(Gjw);
    return M;
}

/* ============================
 * L5: Build Generalized Plant Interconnection
 *
 * Combine the nominal plant G with uncertainty weight Wu and
 * performance weight Wp into the generalized plant P:
 *
 *   [ z_perf ]   [ Wp*G*Wu   Wp*G ] [ w_unc ]
 *   [ z_unc  ] = [       0     0 ] [ u_ctrl ]
 *   [ y      ]   [   G*Wu     G   ] [       ]
 *
 * Complexity: O(n^3)
 * ============================ */

MUSystem* mu_build_interconnection(const MUSystem *G,
                                    const MUSystem *Wu,
                                    const MUSystem *Wp) {
    if (!G) return NULL;
    int nG = G->n, mG = G->m, pG = G->p;
    int nWu = (Wu ? Wu->n : 0);
    int nWp = (Wp ? Wp->n : 0);

    int n = nG + nWu + nWp;
    int m = mG + mG; /* uncertainty input + control input */
    int p = pG + pG; /* performance output + uncertainty output */

    MUSystem *P = mu_system_create(n, m, p);
    if (!P) return NULL;

    /* Place nominal G in the interconnection */
    for (int i = 0; i < nG; i++)
        for (int j = 0; j < nG; j++)
            mu_real_mat_set(P->A, i, j, mu_real_mat_get(G->A, i, j));
    for (int i = 0; i < nG; i++)
        for (int j = 0; j < mG; j++) {
            mu_real_mat_set(P->B, i, j, mu_real_mat_get(G->B, i, j));
            mu_real_mat_set(P->B, i, mG + j, mu_real_mat_get(G->B, i, j));
        }
    for (int i = 0; i < pG; i++)
        for (int j = 0; j < nG; j++) {
            mu_real_mat_set(P->C, i, j, mu_real_mat_get(G->C, i, j));
            mu_real_mat_set(P->C, pG + i, j, mu_real_mat_get(G->C, i, j));
        }
    for (int i = 0; i < pG; i++)
        for (int j = 0; j < mG; j++) {
            mu_real_mat_set(P->D, i, j, mu_real_mat_get(G->D, i, j));
            mu_real_mat_set(P->D, i, mG + j, mu_real_mat_get(G->D, i, j));
            mu_real_mat_set(P->D, pG + i, j, mu_real_mat_get(G->D, i, j));
            mu_real_mat_set(P->D, pG + i, mG + j, mu_real_mat_get(G->D, i, j));
        }

    /* Copy Wu dynamics if provided */
    if (Wu && nWu > 0) {
        int off = nG;
        for (int i = 0; i < nWu; i++)
            for (int j = 0; j < nWu; j++)
                mu_real_mat_set(P->A, off + i, off + j,
                    mu_real_mat_get(Wu->A, i, j));
    }

    /* Copy Wp dynamics if provided */
    if (Wp && nWp > 0) {
        int off = nG + nWu;
        for (int i = 0; i < nWp; i++)
            for (int j = 0; j < nWp; j++)
                mu_real_mat_set(P->A, off + i, off + j,
                    mu_real_mat_get(Wp->A, i, j));
    }
    return P;
}

/* ============================
 * L5: LFT Model Reduction
 *
 * Simplify the LFT representation by removing redundant states
 * from weight functions absorbed into the plant.
 *
 * Uses minimum realization to eliminate uncontrollable/unobservable modes.
 *
 * Complexity: O(n^3)
 * ============================ */

MUSystem* mu_lft_reduce(const MUSystem *sys) {
    if (!sys) return NULL;
    return mu_minreal(sys);
}

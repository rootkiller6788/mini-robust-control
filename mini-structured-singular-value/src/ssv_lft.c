/* ============================================================================
 * ssv_lft.c - Linear Fractional Transformations for Robust Control
 *
 * Implements upper/lower LFTs, Redheffer star product, LFT chain,
 * uncertainty modeling (additive/multiplicative/feedback), and
 * augmented M for robust performance analysis.
 *
 * Key knowledge points:
 *   L1: LFT matrix type — partitioned complex matrix for interconnections
 *   L1: Upper LFT Fu(M,Delta_u) — pull out uncertainty from upper channels
 *   L1: Lower LFT Fl(M,Delta_l) — close loop with controller on lower channels
 *   L2: Well-posedness: det(I-M11*Delta)!=0 / det(I-M22*Delta)!=0
 *   L3: LFT as generalization of transfer function interconnection
 *   L4: Main Loop Theorem via augmented M for robust performance
 *   L5: Redheffer star product — series interconnection of LFT systems
 *   L5: Additive, multiplicative, feedback uncertainty modeling
 *   L6: Robust performance via augmented Delta structure with perf block
 *
 * Mathematical foundations:
 *   Upper LFT: Fu(M, Delta) = M22 + M21*Delta*(I-M11*Delta)^{-1}*M12
 *   Lower LFT: Fl(M, Delta) = M11 + M12*Delta*(I-M22*Delta)^{-1}*M21
 *   Star product: S = M * P with S11,S12,S21,S22 via Redheffer formulas
 *
 * References:
 *   Redheffer, "On a certain linear fractional transformation" (1960)
 *   Doyle, "Analysis of feedback systems with structured uncertainties" (1982)
 *   Packard, "Gain scheduling via linear fractional transformations" (1994)
 *   Zhou, Doyle, Glover — Robust and Optimal Control, Ch. 9-10 (1996)
 *   Skogestad & Postlethwaite — Multivariable Feedback Control, Ch. 8 (2005)
 * ============================================================================ */

#include "ssv_lft.h"
#include "ssv_uncertainty.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * LFT Matrix Construction and Destruction
 * ============================================================================ */

SSVLFTMatrix* ssv_lft_create(const SSVComplexMatrix *M, int n1, int m1, int n2, int m2) {
    if (!M || n1 < 0 || m1 < 0 || n2 < 0 || m2 < 0) return NULL;
    if (n1 + n2 != M->rows || m1 + m2 != M->cols) return NULL;
    SSVLFTMatrix *lft = (SSVLFTMatrix*)calloc(1, sizeof(SSVLFTMatrix));
    if (!lft) return NULL;
    lft->n1 = n1; lft->m1 = m1;
    lft->n2 = n2; lft->m2 = m2;
    lft->M = ssv_cmatrix_create(M->rows, M->cols);
    if (!lft->M) { free(lft); return NULL; }
    ssv_cmatrix_copy(lft->M, M);
    return lft;
}

void ssv_lft_free(SSVLFTMatrix *lft) {
    if (!lft) return;
    ssv_cmatrix_free(lft->M);
    free(lft);
}

/* ============================================================================
 * Internal: extract sub-blocks from partitioned LFT matrix
 * ============================================================================ */

static void lft_get_M11(const SSVLFTMatrix *lft, SSVComplexMatrix *M11) {
    int nr = lft->M->rows;
    for (int i = 0; i < lft->n1; i++)
        for (int j = 0; j < lft->m1; j++)
            M11->data[j * lft->n1 + i] = lft->M->data[j * nr + i];
}

static void lft_get_M12(const SSVLFTMatrix *lft, SSVComplexMatrix *M12) {
    int nr = lft->M->rows;
    for (int i = 0; i < lft->n1; i++)
        for (int j = 0; j < lft->m2; j++)
            M12->data[j * lft->n1 + i] = lft->M->data[(lft->m1 + j) * nr + i];
}

static void lft_get_M21(const SSVLFTMatrix *lft, SSVComplexMatrix *M21) {
    int nr = lft->M->rows;
    for (int i = 0; i < lft->n2; i++)
        for (int j = 0; j < lft->m1; j++)
            M21->data[j * lft->n2 + i] = lft->M->data[j * nr + (lft->n1 + i)];
}

static void lft_get_M22(const SSVLFTMatrix *lft, SSVComplexMatrix *M22) {
    int nr = lft->M->rows;
    for (int i = 0; i < lft->n2; i++)
        for (int j = 0; j < lft->m2; j++)
            M22->data[j * lft->n2 + i] = lft->M->data[(lft->m1 + j) * nr + (lft->n1 + i)];
}

/* ============================================================================
 * LFT from State-Space System at Frequency w
 *
 * G(jw) = C * (jw*I - A)^{-1} * B + D
 *
 * Uses LU decomposition with partial pivoting for solving (jwI-A)*X = B.
 * On singular jwI-A, returns D-only approximation wrapped in LFT.
 * ============================================================================ */

SSVLFTMatrix* ssv_lft_from_state_space(const SSVRealMatrix *A,
                                         const SSVRealMatrix *B,
                                         const SSVRealMatrix *C,
                                         const SSVRealMatrix *D,
                                         double w) {
    if (!A || !B || !C || !D) return NULL;
    int n = A->rows, m = B->cols, p = C->rows;
    if (A->cols != n || B->rows != n || C->cols != n || D->rows != p || D->cols != m)
        return NULL;

    /* Build jwI - A */
    SSVComplexMatrix *sIA = ssv_cmatrix_create(n, n);
    if (!sIA) return NULL;
    for (int i_ = 0; i_ < n; i_++)
        for (int j_ = 0; j_ < n; j_++)
            sIA->data[j_ * n + i_] = (i_ == j_) ? (I * w) : 0.0;
    for (int j_ = 0; j_ < n; j_++)
        for (int i_ = 0; i_ < n; i_++)
            sIA->data[j_ * n + i_] -= A->data[j_ * n + i_];

    /* LU factorization with partial pivoting */
    int *piv = (int*)malloc((size_t)n * sizeof(int));
    for (int i_ = 0; i_ < n; i_++) piv[i_] = i_;
    bool singular = false;

    for (int k = 0; k < n && !singular; k++) {
        int max_row = k;
        double max_val = cabs(sIA->data[k * n + k]);
        for (int i_ = k + 1; i_ < n; i_++) {
            double v = cabs(sIA->data[k * n + i_]);
            if (v > max_val) { max_val = v; max_row = i_; }
        }
        if (max_val < 1e-15) { singular = true; break; }
        if (max_row != k) {
            int tmp = piv[k]; piv[k] = piv[max_row]; piv[max_row] = tmp;
            for (int j_ = 0; j_ < n; j_++) {
                complex double t = sIA->data[j_ * n + k];
                sIA->data[j_ * n + k] = sIA->data[j_ * n + max_row];
                sIA->data[j_ * n + max_row] = t;
            }
        }
        complex double pv = sIA->data[k * n + k];
        for (int i_ = k + 1; i_ < n; i_++) {
            complex double fac = sIA->data[k * n + i_] / pv;
            sIA->data[k * n + i_] = fac;
            for (int j_ = k + 1; j_ < n; j_++)
                sIA->data[j_ * n + i_] -= fac * sIA->data[j_ * n + k];
        }
    }

    if (singular) {
        free(piv); ssv_cmatrix_free(sIA);
        SSVComplexMatrix *G = ssv_cmatrix_create(p, m);
        if (G) {
            for (int j_ = 0; j_ < m; j_++)
                for (int i_ = 0; i_ < p; i_++)
                    G->data[j_ * p + i_] = D->data[j_ * p + i_];
        }
        SSVLFTMatrix *lft_res = ssv_lft_create(G, p, m, 0, 0);
        ssv_cmatrix_free(G);
        return lft_res;
    }

    /* Solve for each column of B via forward/backward substitution */
    SSVComplexMatrix *G = ssv_cmatrix_create(p, m);
    if (!G) { free(piv); ssv_cmatrix_free(sIA); return NULL; }

    complex double *rhs = (complex double*)malloc((size_t)n * sizeof(complex double));
    complex double *sol = (complex double*)malloc((size_t)n * sizeof(complex double));

    for (int col = 0; col < m; col++) {
        for (int i_ = 0; i_ < n; i_++)
            rhs[i_] = B->data[col * n + piv[i_]];
        /* Forward: L*y = rhs (L is unit lower triangular) */
        for (int i_ = 0; i_ < n; i_++) {
            complex double s = rhs[i_];
            for (int j_ = 0; j_ < i_; j_++)
                s -= sIA->data[j_ * n + i_] * sol[j_];
            sol[i_] = s;
        }
        /* Backward: U*x = y */
        for (int i_ = n - 1; i_ >= 0; i_--) {
            complex double s = sol[i_];
            for (int j_ = i_ + 1; j_ < n; j_++)
                s -= sIA->data[j_ * n + i_] * sol[j_];
            sol[i_] = s / sIA->data[i_ * n + i_];
        }
        /* G(:,col) = C*sol + D(:,col) */
        for (int i_ = 0; i_ < p; i_++) {
            complex double s = D->data[col * p + i_];
            for (int j_ = 0; j_ < n; j_++)
                s += C->data[j_ * p + i_] * sol[j_];
            G->data[col * p + i_] = s;
        }
    }

    free(rhs); free(sol); free(piv);
    ssv_cmatrix_free(sIA);

    SSVLFTMatrix *lft_res = ssv_lft_create(G, p, m, 0, 0);
    ssv_cmatrix_free(G);
    return lft_res;
}

/* ============================================================================
 * Upper LFT: Fu(M, Delta_u) = M22 + M21*Delta_u*(I - M11*Delta_u)^{-1}*M12
 *
 * This is the standard form for uncertainty extraction in mu-analysis.
 * Delta_u represents the uncertainty "pulled out" from the upper channels.
 * The result is the transfer function seen from the lower channels.
 *
 * Implementation steps:
 *   1. T = M11 * Delta_u
 *   2. S = I - T, invert S via Cramer's rule (n <= 3)
 *   3. W = Delta_u * S^{-1} * M12
 *   4. Fu = M22 + M21 * W
 * ============================================================================ */

SSVComplexMatrix* ssv_lft_upper(const SSVLFTMatrix *lft,
                                  const SSVComplexMatrix *Delta_u) {
    if (!lft || !Delta_u) return NULL;
    if (Delta_u->rows != lft->m1 || Delta_u->cols != lft->n1) return NULL;

    int n1 = lft->n1, m1 = lft->m1, n2 = lft->n2, m2 = lft->m2;

    /* Extract sub-blocks */
    SSVComplexMatrix *M11 = ssv_cmatrix_create(n1, m1);
    SSVComplexMatrix *M12 = ssv_cmatrix_create(n1, m2);
    SSVComplexMatrix *M21 = ssv_cmatrix_create(n2, m1);
    SSVComplexMatrix *M22 = ssv_cmatrix_create(n2, m2);
    lft_get_M11(lft, M11); lft_get_M12(lft, M12);
    lft_get_M21(lft, M21); lft_get_M22(lft, M22);

    /* T = M11 * Delta_u, then S = I - T */
    SSVComplexMatrix *T = ssv_cmatrix_create(n1, n1);
    ssv_cmatrix_gemm(1.0, M11, Delta_u, 0.0, T);
    for (int i = 0; i < n1; i++)
        T->data[i * n1 + i] = 1.0 - T->data[i * n1 + i];

    SSVComplexMatrix *S_inv = ssv_inverse_small(T);
    ssv_cmatrix_free(T);
    if (!S_inv) {
        /* Not well-posed */
        ssv_cmatrix_free(M11); ssv_cmatrix_free(M12);
        ssv_cmatrix_free(M21); ssv_cmatrix_free(M22);
        return NULL;
    }

    /* W = Delta_u * S_inv * M12 */
    SSVComplexMatrix *DS = ssv_cmatrix_create(m1, n1);
    SSVComplexMatrix *W = ssv_cmatrix_create(m1, m2);
    ssv_cmatrix_gemm(1.0, Delta_u, S_inv, 0.0, DS);
    ssv_cmatrix_gemm(1.0, DS, M12, 0.0, W);

    /* Fu = M22 + M21 * W */
    SSVComplexMatrix *Fu = ssv_cmatrix_create(n2, m2);
    SSVComplexMatrix *M21W = ssv_cmatrix_create(n2, m2);
    ssv_cmatrix_gemm(1.0, M21, W, 0.0, M21W);
    for (int j = 0; j < m2; j++)
        for (int i = 0; i < n2; i++)
            Fu->data[j * n2 + i] = M22->data[j * n2 + i] + M21W->data[j * n2 + i];

    ssv_cmatrix_free(M11); ssv_cmatrix_free(M12);
    ssv_cmatrix_free(M21); ssv_cmatrix_free(M22);
    ssv_cmatrix_free(S_inv); ssv_cmatrix_free(DS);
    ssv_cmatrix_free(W); ssv_cmatrix_free(M21W);
    return Fu;
}

/* ============================================================================
 * Upper LFT Well-Posedness Check
 *
 * The Upper LFT is well-posed iff I - M11*Delta_u is non-singular.
 * For n <= 3: exact via determinant. For larger: condition number via SVD.
 * ============================================================================ */

bool ssv_lft_upper_is_wellposed(const SSVLFTMatrix *lft,
                                  const SSVComplexMatrix *Delta_u,
                                  double tol) {
    if (!lft || !Delta_u) return false;
    if (Delta_u->rows != lft->m1 || Delta_u->cols != lft->n1) return false;

    int n1 = lft->n1;
    SSVComplexMatrix *M11 = ssv_cmatrix_create(n1, lft->m1);
    lft_get_M11(lft, M11);
    SSVComplexMatrix *T = ssv_cmatrix_create(n1, n1);
    ssv_cmatrix_gemm(1.0, M11, Delta_u, 0.0, T);
    for (int i = 0; i < n1; i++)
        T->data[i * n1 + i] = 1.0 - T->data[i * n1 + i];

    bool wp;
    if (n1 <= 3) {
        complex double det = ssv_determinant_small(T);
        wp = (cabs(det) > tol);
    } else {
        SSVSVDResult *svd = ssv_svd_compute(T);
        if (svd) {
            wp = (ssv_svd_min_singular(svd) > tol * ssv_svd_max_singular(svd));
            ssv_svd_free(svd);
        } else {
            wp = false;
        }
    }
    ssv_cmatrix_free(M11); ssv_cmatrix_free(T);
    return wp;
}

/* ============================================================================
 * Lower LFT: Fl(M, Delta_l) = M11 + M12*Delta_l*(I - M22*Delta_l)^{-1}*M21
 *
 * This is the standard form for controller synthesis. The controller K
 * (represented as Delta_l) is connected in feedback from the lower channels.
 * The closed-loop map from w to z is Fl(P, K).
 *
 * Well-posedness: I - M22*Delta_l must be non-singular.
 * ============================================================================ */

SSVComplexMatrix* ssv_lft_lower(const SSVLFTMatrix *lft,
                                  const SSVComplexMatrix *Delta_l) {
    if (!lft || !Delta_l) return NULL;
    if (Delta_l->rows != lft->m2 || Delta_l->cols != lft->n2) return NULL;

    int n1 = lft->n1, m1 = lft->m1, n2 = lft->n2, m2 = lft->m2;

    SSVComplexMatrix *M11 = ssv_cmatrix_create(n1, m1);
    SSVComplexMatrix *M12 = ssv_cmatrix_create(n1, m2);
    SSVComplexMatrix *M21 = ssv_cmatrix_create(n2, m1);
    SSVComplexMatrix *M22 = ssv_cmatrix_create(n2, m2);
    lft_get_M11(lft, M11); lft_get_M12(lft, M12);
    lft_get_M21(lft, M21); lft_get_M22(lft, M22);

    /* T = M22 * Delta_l, S = I - T */
    SSVComplexMatrix *T = ssv_cmatrix_create(n2, n2);
    ssv_cmatrix_gemm(1.0, M22, Delta_l, 0.0, T);
    for (int i = 0; i < n2; i++)
        T->data[i * n2 + i] = 1.0 - T->data[i * n2 + i];

    SSVComplexMatrix *S_inv = ssv_inverse_small(T);
    ssv_cmatrix_free(T);
    if (!S_inv) {
        ssv_cmatrix_free(M11); ssv_cmatrix_free(M12);
        ssv_cmatrix_free(M21); ssv_cmatrix_free(M22);
        return NULL;
    }

    /* W = Delta_l * S_inv * M21 */
    SSVComplexMatrix *DS = ssv_cmatrix_create(m2, n2);
    SSVComplexMatrix *W = ssv_cmatrix_create(m2, m1);
    ssv_cmatrix_gemm(1.0, Delta_l, S_inv, 0.0, DS);
    ssv_cmatrix_gemm(1.0, DS, M21, 0.0, W);

    /* Fl = M11 + M12 * W */
    SSVComplexMatrix *Fl = ssv_cmatrix_create(n1, m1);
    SSVComplexMatrix *M12W = ssv_cmatrix_create(n1, m1);
    ssv_cmatrix_gemm(1.0, M12, W, 0.0, M12W);
    for (int j = 0; j < m1; j++)
        for (int i = 0; i < n1; i++)
            Fl->data[j * n1 + i] = M11->data[j * n1 + i] + M12W->data[j * n1 + i];

    ssv_cmatrix_free(M11); ssv_cmatrix_free(M12);
    ssv_cmatrix_free(M21); ssv_cmatrix_free(M22);
    ssv_cmatrix_free(S_inv); ssv_cmatrix_free(DS);
    ssv_cmatrix_free(W); ssv_cmatrix_free(M12W);
    return Fl;
}

bool ssv_lft_lower_is_wellposed(const SSVLFTMatrix *lft,
                                  const SSVComplexMatrix *Delta_l,
                                  double tol) {
    if (!lft || !Delta_l) return false;
    if (Delta_l->rows != lft->m2 || Delta_l->cols != lft->n2) return false;

    int n2 = lft->n2;
    SSVComplexMatrix *M22 = ssv_cmatrix_create(n2, lft->m2);
    lft_get_M22(lft, M22);
    SSVComplexMatrix *T = ssv_cmatrix_create(n2, n2);
    ssv_cmatrix_gemm(1.0, M22, Delta_l, 0.0, T);
    for (int i = 0; i < n2; i++)
        T->data[i * n2 + i] = 1.0 - T->data[i * n2 + i];

    bool wp;
    if (n2 <= 3) {
        wp = (cabs(ssv_determinant_small(T)) > tol);
    } else {
        SSVSVDResult *svd = ssv_svd_compute(T);
        wp = svd ? (ssv_svd_min_singular(svd) > tol * ssv_svd_max_singular(svd)) : false;
        ssv_svd_free(svd);
    }
    ssv_cmatrix_free(M22); ssv_cmatrix_free(T);
    return wp;
}

/* ============================================================================
 * Redheffer Star Product: S = M * P
 *
 * The Redheffer star product generalizes transfer function multiplication
 * to multi-port LFT representations. It interconnects two LFT systems
 * M and P in series, yielding a combined LFT system S.
 *
 * Formulas (with k = M.n2 = P.m1 internal feedback dimension):
 *   S11 = M11 + M12 * P11 * (I - M22*P11)^{-1} * M21
 *   S12 = M12 * (I - P11*M22)^{-1} * P12
 *   S21 = P21 * (I - M22*P11)^{-1} * M21
 *   S22 = P22 + P21 * M22 * (I - P11*M22)^{-1} * P12
 *
 * This corresponds to the feedback interconnection where the lower
 * outputs of M feed the upper inputs of P, and vice versa.
 *
 * Complexity: O(k^3 + mn^2) dominated by matrix inversions and products.
 * ============================================================================ */

SSVLFTMatrix* ssv_lft_star_product(const SSVLFTMatrix *M, const SSVLFTMatrix *P) {
    if (!M || !P) return NULL;
    if (M->n2 != P->m1 || M->m2 != P->n1) return NULL;

    int mn1 = M->n1, mm1 = M->m1;
    int pn2 = P->n2, pm2 = P->m2;
    int k = M->n2;  /* internal feedback dimension */

    /* Extract all sub-blocks */
    SSVComplexMatrix *M11 = ssv_cmatrix_create(mn1, mm1);
    SSVComplexMatrix *M12 = ssv_cmatrix_create(mn1, M->m2);
    SSVComplexMatrix *M21 = ssv_cmatrix_create(M->n2, mm1);
    SSVComplexMatrix *M22 = ssv_cmatrix_create(M->n2, M->m2);
    SSVComplexMatrix *P11 = ssv_cmatrix_create(P->n1, P->m1);
    SSVComplexMatrix *P12 = ssv_cmatrix_create(P->n1, P->m2);
    SSVComplexMatrix *P21 = ssv_cmatrix_create(P->n2, P->m1);
    SSVComplexMatrix *P22 = ssv_cmatrix_create(P->n2, P->m2);
    lft_get_M11(M, M11); lft_get_M12(M, M12);
    lft_get_M21(M, M21); lft_get_M22(M, M22);
    lft_get_M11(P, P11); lft_get_M12(P, P12);
    lft_get_M21(P, P21); lft_get_M22(P, P22);

    /* Compute inv(I - M22*P11) */
    SSVComplexMatrix *MP = ssv_cmatrix_create(k, k);
    ssv_cmatrix_gemm(1.0, M22, P11, 0.0, MP);
    for (int i = 0; i < k; i++) MP->data[i * k + i] = 1.0 - MP->data[i * k + i];
    SSVComplexMatrix *inv_IMP = ssv_inverse_small(MP);

    /* Compute inv(I - P11*M22) */
    SSVComplexMatrix *PM = ssv_cmatrix_create(k, k);
    ssv_cmatrix_gemm(1.0, P11, M22, 0.0, PM);
    for (int i = 0; i < k; i++) PM->data[i * k + i] = 1.0 - PM->data[i * k + i];
    SSVComplexMatrix *inv_IPM = ssv_inverse_small(PM);

    if (!inv_IMP || !inv_IPM) {
        /* Not well-posed */
        ssv_cmatrix_free(M11); ssv_cmatrix_free(M12);
        ssv_cmatrix_free(M21); ssv_cmatrix_free(M22);
        ssv_cmatrix_free(P11); ssv_cmatrix_free(P12);
        ssv_cmatrix_free(P21); ssv_cmatrix_free(P22);
        ssv_cmatrix_free(MP); ssv_cmatrix_free(PM);
        ssv_cmatrix_free(inv_IMP); ssv_cmatrix_free(inv_IPM);
        return NULL;
    }

    /* S11 = M11 + M12 * P11 * inv(I-M22*P11) * M21 */
    SSVComplexMatrix *T1 = ssv_cmatrix_create(mn1, k);
    SSVComplexMatrix *T2 = ssv_cmatrix_create(mn1, k);
    ssv_cmatrix_gemm(1.0, M12, P11, 0.0, T1);
    ssv_cmatrix_gemm(1.0, T1, inv_IMP, 0.0, T2);
    SSVComplexMatrix *TT = ssv_cmatrix_create(mn1, mm1);
    ssv_cmatrix_gemm(1.0, T2, M21, 0.0, TT);
    SSVComplexMatrix *S11 = ssv_cmatrix_create(mn1, mm1);
    for (int j = 0; j < mm1; j++)
        for (int i = 0; i < mn1; i++)
            S11->data[j * mn1 + i] = M11->data[j * mn1 + i] + TT->data[j * mn1 + i];

    /* S12 = M12 * inv(I-P11*M22) * P12 */
    SSVComplexMatrix *T3 = ssv_cmatrix_create(mn1, k);
    ssv_cmatrix_gemm(1.0, M12, inv_IPM, 0.0, T3);
    SSVComplexMatrix *S12 = ssv_cmatrix_create(mn1, pm2);
    ssv_cmatrix_gemm(1.0, T3, P12, 0.0, S12);

    /* S21 = P21 * inv(I-M22*P11) * M21 */
    SSVComplexMatrix *T4 = ssv_cmatrix_create(pn2, k);
    ssv_cmatrix_gemm(1.0, P21, inv_IMP, 0.0, T4);
    SSVComplexMatrix *S21 = ssv_cmatrix_create(pn2, mm1);
    ssv_cmatrix_gemm(1.0, T4, M21, 0.0, S21);

    /* S22 = P22 + P21 * M22 * inv(I-P11*M22) * P12 */
    SSVComplexMatrix *PMi = ssv_cmatrix_create(pn2, k);
    ssv_cmatrix_gemm(1.0, P21, M22, 0.0, PMi);
    SSVComplexMatrix *T5 = ssv_cmatrix_create(pn2, k);
    ssv_cmatrix_gemm(1.0, PMi, inv_IPM, 0.0, T5);
    SSVComplexMatrix *S22 = ssv_cmatrix_create(pn2, pm2);
    ssv_cmatrix_gemm(1.0, T5, P12, 0.0, S22);
    for (int j = 0; j < pm2; j++)
        for (int i = 0; i < pn2; i++)
            S22->data[j * pn2 + i] += P22->data[j * pn2 + i];

    /* Assemble full S = [S11 S12; S21 S22] */
    int Sr = mn1 + pn2, Sc = mm1 + pm2;
    SSVComplexMatrix *Sf = ssv_cmatrix_create(Sr, Sc);
    for (int j = 0; j < mm1; j++)
        for (int i = 0; i < mn1; i++)
            Sf->data[j * Sr + i] = S11->data[j * mn1 + i];
    for (int j = 0; j < pm2; j++)
        for (int i = 0; i < mn1; i++)
            Sf->data[(mm1 + j) * Sr + i] = S12->data[j * mn1 + i];
    for (int j = 0; j < mm1; j++)
        for (int i = 0; i < pn2; i++)
            Sf->data[j * Sr + (mn1 + i)] = S21->data[j * pn2 + i];
    for (int j = 0; j < pm2; j++)
        for (int i = 0; i < pn2; i++)
            Sf->data[(mm1 + j) * Sr + (mn1 + i)] = S22->data[j * pn2 + i];

    SSVLFTMatrix *result = ssv_lft_create(Sf, mn1, mm1, pn2, pm2);

    /* Cleanup */
    ssv_cmatrix_free(Sf); ssv_cmatrix_free(S11); ssv_cmatrix_free(S12);
    ssv_cmatrix_free(S21); ssv_cmatrix_free(S22);
    ssv_cmatrix_free(M11); ssv_cmatrix_free(M12);
    ssv_cmatrix_free(M21); ssv_cmatrix_free(M22);
    ssv_cmatrix_free(P11); ssv_cmatrix_free(P12);
    ssv_cmatrix_free(P21); ssv_cmatrix_free(P22);
    ssv_cmatrix_free(MP); ssv_cmatrix_free(PM);
    ssv_cmatrix_free(inv_IMP); ssv_cmatrix_free(inv_IPM);
    ssv_cmatrix_free(T1); ssv_cmatrix_free(T2); ssv_cmatrix_free(TT);
    ssv_cmatrix_free(T3); ssv_cmatrix_free(T4); ssv_cmatrix_free(T5);
    ssv_cmatrix_free(PMi);
    return result;
}

/* ============================================================================
 * LFT Chain: series interconnection M1 * M2 * ... * Mn
 *
 * Repeated application of the Redheffer star product.
 * Each successive star product requires compatible dimensions.
 * Returns NULL on any dimension mismatch or well-posedness failure.
 * ============================================================================ */

SSVLFTMatrix* ssv_lft_chain(const SSVLFTMatrix **lfts, int n_lfts) {
    if (!lfts || n_lfts <= 0) return NULL;
    if (n_lfts == 1)
        return ssv_lft_create(lfts[0]->M, lfts[0]->n1, lfts[0]->m1,
                              lfts[0]->n2, lfts[0]->m2);
    SSVLFTMatrix *chain = ssv_lft_create(lfts[0]->M, lfts[0]->n1, lfts[0]->m1,
                                          lfts[0]->n2, lfts[0]->m2);
    if (!chain) return NULL;
    for (int i = 1; i < n_lfts; i++) {
        SSVLFTMatrix *next = ssv_lft_star_product(chain, lfts[i]);
        ssv_lft_free(chain);
        if (!next) return NULL;
        chain = next;
    }
    return chain;
}

/* ============================================================================
 * Additive Uncertainty LFT Construction
 *
 * Model: G = G_nom + W_u * Delta * W_y
 *
 * LFT form: M = [  0        W_y    ]
 *                [  W_u      G_nom  ]
 *
 * Dimension: M is (p+p) x (p+m) where G_nom is p x m.
 * M11 = 0 (p x p), M12 = W_y, M21 = W_u, M22 = G_nom.
 * ============================================================================ */

SSVComplexMatrix* ssv_lft_additive_uncertainty(
    const SSVComplexMatrix *G_nom,
    const SSVComplexMatrix *W_u,
    const SSVComplexMatrix *W_y) {
    if (!G_nom || !W_u || !W_y) return NULL;
    int p = G_nom->rows, m = G_nom->cols;
    int n_M = p + p, m_M = p + m;
    SSVComplexMatrix *M = ssv_cmatrix_create(n_M, m_M);
    if (!M) return NULL;
    /* M12 = W_y */
    for (int i = 0; i < W_y->rows; i++)
        for (int j = 0; j < W_y->cols; j++)
            M->data[(p + j) * n_M + i] = W_y->data[j * W_y->rows + i];
    /* M21 = W_u */
    for (int i = 0; i < W_u->rows; i++)
        for (int j = 0; j < W_u->cols; j++)
            M->data[j * n_M + (p + i)] = W_u->data[j * W_u->rows + i];
    /* M22 = G_nom */
    for (int i = 0; i < p; i++)
        for (int j = 0; j < m; j++)
            M->data[(p + j) * n_M + (p + i)] = G_nom->data[j * p + i];
    return M;
}

/* ============================================================================
 * Input Multiplicative Uncertainty LFT Construction
 *
 * Model: G = G_nom * (I + W_u * Delta * W_y)
 *
 * LFT form: M = [ -W_y*G_nom*W_u    -W_y*G_nom ]
 *                [  W_u               G_nom      ]
 *
 * Block sizes: Delta is pu x py. M is (py+p) x (pu+m).
 * ============================================================================ */

SSVComplexMatrix* ssv_lft_input_multiplicative(
    const SSVComplexMatrix *G_nom,
    const SSVComplexMatrix *W_u,
    const SSVComplexMatrix *W_y) {
    if (!G_nom || !W_u || !W_y) return NULL;
    int p = G_nom->rows, m = G_nom->cols;
    int pu = W_u->cols, py = W_y->rows;
    int n_M = py + p, m_M = pu + m;
    SSVComplexMatrix *M = ssv_cmatrix_create(n_M, m_M);
    if (!M) return NULL;
    /* Compute W_y * G_nom * W_u for M11 */
    SSVComplexMatrix *WG = ssv_cmatrix_create(py, m);
    SSVComplexMatrix *WGW = ssv_cmatrix_create(py, pu);
    ssv_cmatrix_gemm(1.0, W_y, G_nom, 0.0, WG);
    ssv_cmatrix_gemm(-1.0, WG, W_u, 0.0, WGW);
    /* M11 = -W_y*G_nom*W_u */
    for (int i = 0; i < py; i++)
        for (int j = 0; j < pu; j++)
            M->data[j * n_M + i] = WGW->data[j * py + i];
    /* M12 = -W_y*G_nom */
    for (int i = 0; i < py; i++)
        for (int j = 0; j < m; j++)
            M->data[(pu + j) * n_M + i] = -WG->data[j * py + i];
    /* M21 = W_u */
    for (int i = 0; i < W_u->rows; i++)
        for (int j = 0; j < W_u->cols; j++)
            M->data[j * n_M + (py + i)] = W_u->data[j * W_u->rows + i];
    /* M22 = G_nom */
    for (int i = 0; i < p; i++)
        for (int j = 0; j < m; j++)
            M->data[(pu + j) * n_M + (py + i)] = G_nom->data[j * p + i];
    ssv_cmatrix_free(WG); ssv_cmatrix_free(WGW);
    return M;
}

/* ============================================================================
 * Output Multiplicative Uncertainty LFT Construction
 *
 * Model: G = (I + W_u * Delta * W_y) * G_nom
 *
 * LFT form: M = [ -W_y*G_nom*W_u    -W_y*G_nom ]
 *                [  G_nom*W_u         G_nom      ]
 *
 * Key difference from input multiplicative: M21 = G_nom*W_u (not W_u).
 * ============================================================================ */

SSVComplexMatrix* ssv_lft_output_multiplicative(
    const SSVComplexMatrix *G_nom,
    const SSVComplexMatrix *W_u,
    const SSVComplexMatrix *W_y) {
    if (!G_nom || !W_u || !W_y) return NULL;
    int p = G_nom->rows, m = G_nom->cols;
    int pu = W_u->cols, py = W_y->rows;
    int n_M = py + p, m_M = pu + m;
    SSVComplexMatrix *M = ssv_cmatrix_create(n_M, m_M);
    if (!M) return NULL;
    /* W_y * G_nom * W_u and G_nom * W_u */
    SSVComplexMatrix *WG = ssv_cmatrix_create(py, m);
    SSVComplexMatrix *WGW = ssv_cmatrix_create(py, pu);
    SSVComplexMatrix *GWu = ssv_cmatrix_create(p, pu);
    ssv_cmatrix_gemm(1.0, W_y, G_nom, 0.0, WG);
    ssv_cmatrix_gemm(-1.0, WG, W_u, 0.0, WGW);
    ssv_cmatrix_gemm(1.0, G_nom, W_u, 0.0, GWu);
    /* M11 = -W_y*G_nom*W_u */
    for (int i = 0; i < py; i++)
        for (int j = 0; j < pu; j++)
            M->data[j * n_M + i] = WGW->data[j * py + i];
    /* M12 = -W_y*G_nom */
    for (int i = 0; i < py; i++)
        for (int j = 0; j < m; j++)
            M->data[(pu + j) * n_M + i] = -WG->data[j * py + i];
    /* M21 = G_nom*W_u */
    for (int i = 0; i < p; i++)
        for (int j = 0; j < pu; j++)
            M->data[j * n_M + (py + i)] = GWu->data[j * p + i];
    /* M22 = G_nom */
    for (int i = 0; i < p; i++)
        for (int j = 0; j < m; j++)
            M->data[(pu + j) * n_M + (py + i)] = G_nom->data[j * p + i];
    ssv_cmatrix_free(WG); ssv_cmatrix_free(WGW); ssv_cmatrix_free(GWu);
    return M;
}

/* ============================================================================
 * Feedback Uncertainty LFT Construction
 *
 * Model: G = G_nom * (I + W_u * Delta * W_y)^{-1}
 *
 * This represents left coprime factor uncertainty where Delta is in the
 * feedback path. The LFT form is structurally similar to input
 * multiplicative but represents a fundamentally different uncertainty
 * type (feedback rather than feedforward).
 * ============================================================================ */

SSVComplexMatrix* ssv_lft_feedback_uncertainty(
    const SSVComplexMatrix *G_nom,
    const SSVComplexMatrix *W_u,
    const SSVComplexMatrix *W_y) {
    if (!G_nom || !W_u || !W_y) return NULL;
    int p = G_nom->rows, m = G_nom->cols;
    int pu = W_u->cols, py = W_y->rows;
    int n_M = py + p, m_M = pu + m;
    SSVComplexMatrix *M = ssv_cmatrix_create(n_M, m_M);
    if (!M) return NULL;
    SSVComplexMatrix *WG = ssv_cmatrix_create(py, m);
    SSVComplexMatrix *WGW = ssv_cmatrix_create(py, pu);
    ssv_cmatrix_gemm(1.0, W_y, G_nom, 0.0, WG);
    ssv_cmatrix_gemm(-1.0, WG, W_u, 0.0, WGW);
    for (int i = 0; i < py; i++)
        for (int j = 0; j < pu; j++)
            M->data[j * n_M + i] = WGW->data[j * py + i];
    for (int i = 0; i < py; i++)
        for (int j = 0; j < m; j++)
            M->data[(pu + j) * n_M + i] = -WG->data[j * py + i];
    for (int i = 0; i < W_u->rows; i++)
        for (int j = 0; j < W_u->cols; j++)
            M->data[j * n_M + (py + i)] = W_u->data[j * W_u->rows + i];
    for (int i = 0; i < p; i++)
        for (int j = 0; j < m; j++)
            M->data[(pu + j) * n_M + (py + i)] = G_nom->data[j * p + i];
    ssv_cmatrix_free(WG); ssv_cmatrix_free(WGW);
    return M;
}

/* ============================================================================
 * Augmented M for Robust Performance
 *
 * Main Loop Theorem (robust performance form):
 *   The system achieves robust performance level <= 1/gamma
 *   iff mu_{Delta_aug}(M_aug) < 1/gamma for all frequencies,
 *   where Delta_aug = diag(Delta_orig, Delta_perf) and
 *   Delta_perf is a fictitious full complex block.
 *
 * The augmented M is simply the original M embedded in a larger matrix
 * with the performance channels in the trailing rows/columns.
 * ============================================================================ */

SSVComplexMatrix* ssv_lft_augment_performance(
    const SSVComplexMatrix *M, int n_delta, int n_perf) {
    if (!M || n_delta < 0 || n_perf <= 0) return NULL;
    (void)n_delta; /* reserved for future partitioning use */
    int n_orig = M->rows, n_total = n_orig + n_perf;
    SSVComplexMatrix *M_aug = ssv_cmatrix_create(n_total, n_total);
    if (!M_aug) return NULL;
    for (int i = 0; i < n_orig; i++)
        for (int j = 0; j < n_orig; j++)
            M_aug->data[j * n_total + i] = M->data[j * n_orig + i];
    for (int i = 0; i < n_perf; i++)
        M_aug->data[(n_orig + i) * n_total + (n_orig + i)] = 1.0;
    return M_aug;
}

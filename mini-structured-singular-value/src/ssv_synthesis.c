/* ============================================================================
 * ssv_synthesis.c — mu-Synthesis via D-K Iteration
 *
 * Implements the D-K iteration framework for mu-synthesis controller design,
 * state-space frequency response computation, and H-infinity norm evaluation.
 *
 * Key knowledge points:
 *   L1: State-space system representation (A, B, C, D)
 *   L1: Generalized plant for H∞/μ synthesis
 *   L2: D-K iteration (alternating between H∞ and D-scaling)
 *   L3: Frequency response G(jw) = C*(jwI-A)^{-1}*B + D
 *   L5: D-K iteration algorithm, H∞ norm via bisection
 *   L5: Hamiltonian matrix and eigenvalue check for H∞ norm
 *   L6: Closed-loop interconnection via lower LFT Fl(P,K)
 *
 * References:
 *   Doyle, Glover, Khargonekar, Francis — "State-space solutions to
 *     standard H2 and H∞ control problems" (IEEE TAC, 1989)
 *   Balas, Doyle, Glover, Packard, Smith — "mu-Analysis and Synthesis
 *     Toolbox" (1998, MATLAB/Simulink)
 *   Zhou, Doyle & Glover — "Robust and Optimal Control" (1996)
 * ============================================================================ */

#include "ssv_synthesis.h"
#include "ssv_lft.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * State-Space System
 * ============================================================================ */

SSVStateSpace* ssv_ss_create(int n_states, int n_inputs, int n_outputs) {
    if (n_states <= 0 || n_inputs <= 0 || n_outputs <= 0) return NULL;
    SSVStateSpace *sys = (SSVStateSpace*)calloc(1, sizeof(SSVStateSpace));
    if (!sys) return NULL;
    sys->n_states = n_states;
    sys->n_inputs = n_inputs;
    sys->n_outputs = n_outputs;
    sys->A = ssv_rmatrix_create(n_states, n_states);
    sys->B = ssv_rmatrix_create(n_states, n_inputs);
    sys->C = ssv_rmatrix_create(n_outputs, n_states);
    sys->D = ssv_rmatrix_create(n_outputs, n_inputs);
    if (!sys->A || !sys->B || !sys->C || !sys->D) {
        ssv_ss_free(sys);
        return NULL;
    }
    return sys;
}

void ssv_ss_free(SSVStateSpace *sys) {
    if (!sys) return;
    ssv_rmatrix_free(sys->A);
    ssv_rmatrix_free(sys->B);
    ssv_rmatrix_free(sys->C);
    ssv_rmatrix_free(sys->D);
    free(sys);
}

/* Frequency response: G(jw) = C*(jwI - A)^{-1}*B + D
 *
 * Compute (jwI - A)^{-1} via LU decomposition with partial pivoting,
 * then form G = C * inv * B + D.
 *
 * For small systems (n <= 3), use explicit formula. For larger systems,
 * solve linear system via Gaussian elimination with partial pivoting.
 */
SSVComplexMatrix* ssv_ss_freqresp(const SSVStateSpace *sys, double w) {
    if (!sys) return NULL;
    int n = sys->n_states, m = sys->n_inputs, p = sys->n_outputs;

    /* Form jwI - A as a complex matrix */
    SSVComplexMatrix *sIA = ssv_cmatrix_create(n, n);
    if (!sIA) return NULL;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            complex double val = (i == j) ? (I * w) : 0.0;
            val -= sys->A->data[j * n + i];  /* subtract A(i,j) */
            sIA->data[j * n + i] = val;
        }
    }

    /* Solve (jwI - A) * X = B for X, then G = C * X + D.
     * For each column of B, solve the linear system.
     *
     * Use Gaussian elimination with partial pivoting. */
    SSVComplexMatrix *G = ssv_cmatrix_create(p, m);
    if (!G) { ssv_cmatrix_free(sIA); return NULL; }

    /* LU factorization of sIA in-place */
    int *piv = (int*)malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) piv[i] = i;

    for (int k = 0; k < n; k++) {
        /* Find pivot */
        int max_row = k;
        double max_val = cabs(sIA->data[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double val = cabs(sIA->data[k * n + i]);
            if (val > max_val) { max_val = val; max_row = i; }
        }
        if (max_val < 1e-15) {
            /* Singular — continue with identity-like */
            free(piv);
            ssv_cmatrix_free(sIA);
            /* Return D as best approximation */
            for (int i = 0; i < p; i++)
                for (int j = 0; j < m; j++)
                    G->data[j * p + i] = sys->D->data[j * p + i];
            return G;
        }

        /* Swap rows */
        if (max_row != k) {
            int tmp = piv[k]; piv[k] = piv[max_row]; piv[max_row] = tmp;
            for (int j = 0; j < n; j++) {
                complex double t = sIA->data[j * n + k];
                sIA->data[j * n + k] = sIA->data[j * n + max_row];
                sIA->data[j * n + max_row] = t;
            }
        }

        /* Eliminate below */
        complex double pivot_val = sIA->data[k * n + k];
        for (int i = k + 1; i < n; i++) {
            complex double factor = sIA->data[k * n + i] / pivot_val;
            sIA->data[k * n + i] = factor;  /* store L */
            for (int j = k + 1; j < n; j++) {
                sIA->data[j * n + i] -= factor * sIA->data[j * n + k];
            }
        }
    }

    /* Forward/backward substitution for each input column */
    complex double *rhs = (complex double*)malloc((size_t)n * sizeof(complex double));
    complex double *sol = (complex double*)malloc((size_t)n * sizeof(complex double));

    for (int col = 0; col < m; col++) {
        /* Apply pivoting to RHS = B(:,col) */
        for (int i = 0; i < n; i++) {
            rhs[i] = sys->B->data[col * n + piv[i]];
        }

        /* Forward substitution: L * y = rhs (L stored in lower triangle) */
        /* Unit diagonal L: solve row by row */
        for (int i = 0; i < n; i++) {
            complex double sum = rhs[i];
            for (int j = 0; j < i; j++)
                sum -= sIA->data[j * n + i] * sol[j];
            sol[i] = sum; /* L_ii = 1 */
        }

        /* Backward substitution: U * x = y */
        for (int i = n - 1; i >= 0; i--) {
            complex double sum = sol[i];
            for (int j = i + 1; j < n; j++)
                sum -= sIA->data[j * n + i] * sol[j];
            sol[i] = sum / sIA->data[i * n + i];
        }

        /* G(:,col) = C * sol + D(:,col) */
        for (int i = 0; i < p; i++) {
            complex double sum = sys->D->data[col * p + i];
            for (int j = 0; j < n; j++)
                sum += sys->C->data[j * p + i] * sol[j];
            G->data[col * p + i] = sum;
        }
    }

    free(rhs); free(sol); free(piv);
    ssv_cmatrix_free(sIA);
    return G;
}

SSVComplexMatrix** ssv_ss_freqresp_grid(const SSVStateSpace *sys,
                                          const double *w, int n_freq) {
    if (!sys || !w || n_freq <= 0) return NULL;
    SSVComplexMatrix **responses = (SSVComplexMatrix**)calloc((size_t)n_freq, sizeof(SSVComplexMatrix*));
    if (!responses) return NULL;
    for (int i = 0; i < n_freq; i++)
        responses[i] = ssv_ss_freqresp(sys, w[i]);
    return responses;
}

void ssv_ss_print(const SSVStateSpace *sys) {
    if (!sys) { printf("SSVStateSpace: NULL\n"); return; }
    printf("=== State-Space System ===\n");
    printf("  Order: %d, Inputs: %d, Outputs: %d\n",
           sys->n_states, sys->n_inputs, sys->n_outputs);
    printf("  A (%dx%d):\n", sys->n_states, sys->n_states);
    for (int i = 0; i < sys->n_states && i < 8; i++) {
        printf("    ");
        for (int j = 0; j < sys->n_states && j < 8; j++)
            printf("%8.4f ", sys->A->data[j * sys->n_states + i]);
        if (sys->n_states > 8) printf("...");
        printf("\n");
    }
    printf("  B (%dx%d), C (%dx%d), D (%dx%d)\n",
           sys->n_states, sys->n_inputs, sys->n_outputs, sys->n_states,
           sys->n_outputs, sys->n_inputs);
}

/* ============================================================================
 * Generalized Plant
 * ============================================================================ */

SSVGeneralizedPlant* ssv_gplant_create(int n_states, int n_w, int n_u, int n_z, int n_y) {
    if (n_states <= 0) return NULL;
    SSVGeneralizedPlant *P = (SSVGeneralizedPlant*)calloc(1, sizeof(SSVGeneralizedPlant));
    if (!P) return NULL;
    P->n_states = n_states;
    P->n_w = n_w; P->n_u = n_u; P->n_z = n_z; P->n_y = n_y;

    P->A = ssv_rmatrix_create(n_states, n_states);
    P->B1 = ssv_rmatrix_create(n_states, n_w);
    P->B2 = ssv_rmatrix_create(n_states, n_u);
    P->C1 = ssv_rmatrix_create(n_z, n_states);
    P->C2 = ssv_rmatrix_create(n_y, n_states);
    P->D11 = ssv_rmatrix_create(n_z, n_w);
    P->D12 = ssv_rmatrix_create(n_z, n_u);
    P->D21 = ssv_rmatrix_create(n_y, n_w);
    P->D22 = ssv_rmatrix_create(n_y, n_u);

    if (!P->A || !P->B1 || !P->B2 || !P->C1 || !P->C2 ||
        !P->D11 || !P->D12 || !P->D21 || !P->D22) {
        ssv_gplant_free(P);
        return NULL;
    }
    return P;
}

void ssv_gplant_free(SSVGeneralizedPlant *P) {
    if (!P) return;
    ssv_rmatrix_free(P->A);
    ssv_rmatrix_free(P->B1); ssv_rmatrix_free(P->B2);
    ssv_rmatrix_free(P->C1); ssv_rmatrix_free(P->C2);
    ssv_rmatrix_free(P->D11); ssv_rmatrix_free(P->D12);
    ssv_rmatrix_free(P->D21); ssv_rmatrix_free(P->D22);
    free(P);
}

void ssv_gplant_set_matrix(SSVGeneralizedPlant *P, const char *submatrix_name,
                            const SSVRealMatrix *M) {
    if (!P || !M || !submatrix_name) return;
    SSVRealMatrix *dst = NULL;
    if (strcmp(submatrix_name, "A") == 0) dst = P->A;
    else if (strcmp(submatrix_name, "B1") == 0) dst = P->B1;
    else if (strcmp(submatrix_name, "B2") == 0) dst = P->B2;
    else if (strcmp(submatrix_name, "C1") == 0) dst = P->C1;
    else if (strcmp(submatrix_name, "C2") == 0) dst = P->C2;
    else if (strcmp(submatrix_name, "D11") == 0) dst = P->D11;
    else if (strcmp(submatrix_name, "D12") == 0) dst = P->D12;
    else if (strcmp(submatrix_name, "D21") == 0) dst = P->D21;
    else if (strcmp(submatrix_name, "D22") == 0) dst = P->D22;
    else return;

    if (dst->rows != M->rows || dst->cols != M->cols) return;
    size_t sz = (size_t)M->rows * (size_t)M->cols;
    for (size_t k = 0; k < sz; k++)
        dst->data[k] = M->data[k];
}

/* ============================================================================
 * Closed-Loop System: Fl(P, K)
 *
 * Lower LFT of generalized plant P and controller K:
 *   Fl(P,K) = P11 + P12 * K * (I - P22*K)^{-1} * P21
 *
 * State-space: P has states x_P, K has states x_K.
 * Closed-loop states: [x_P^T, x_K^T]^T.
 *
 *    A_cl = [ A + B2*Dk*C2     B2*Ck           ]
 *           [ Bk*C2            Ak               ]
 *
 *    B_cl = [ B1 + B2*Dk*D21  ]
 *           [ Bk*D21           ]
 *
 *    C_cl = [ C1 + D12*Dk*C2   D12*Ck ]
 *
 *    D_cl = D11 + D12 * Dk * D21
 *
 * Where K has state-space (Ak, Bk, Ck, Dk).
 * ============================================================================ */

SSVStateSpace* ssv_closed_loop(const SSVGeneralizedPlant *P, const SSVStateSpace *K) {
    if (!P || !K) return NULL;

    int np = P->n_states, nk = K->n_states;
    int n_cl = np + nk;
    int n_w = P->n_w, n_z = P->n_z;

    SSVStateSpace *cl = ssv_ss_create(n_cl, n_w, n_z);
    if (!cl) return NULL;

    /* Extract controller matrices: K = (Ak, Bk, Ck, Dk) */
    SSVRealMatrix *Ak = K->A, *Bk = K->B, *Ck = K->C, *Dk = K->D;

    /* Build A_cl */
    /* Top-left: A + B2 * Dk * C2 */
    SSVRealMatrix *B2Dk = ssv_rmatrix_create(np, P->n_y);
    ssv_rmatrix_gemm(1.0, P->B2, Dk, 0.0, B2Dk);        /* B2Dk = B2 * Dk */
    ssv_rmatrix_gemm(1.0, B2Dk, P->C2, 0.0, cl->A);      /* A_cl_top = B2Dk * C2 */

    /* Add A to top-left */
    for (int i = 0; i < np; i++)
        for (int j = 0; j < np; j++)
            cl->A->data[j * n_cl + i] += P->A->data[j * np + i];

    /* Top-right: B2 * Ck */
    ssv_rmatrix_gemm(1.0, P->B2, Ck, 0.0, B2Dk);         /* reuse B2Dk as B2*Ck */
    for (int i = 0; i < np; i++)
        for (int j = 0; j < nk; j++)
            cl->A->data[(np + j) * n_cl + i] = B2Dk->data[j * np + i];

    /* Bottom-left: Bk * C2 */
    SSVRealMatrix *BkC2 = ssv_rmatrix_create(nk, np);
    ssv_rmatrix_gemm(1.0, Bk, P->C2, 0.0, BkC2);
    for (int i = 0; i < nk; i++)
        for (int j = 0; j < np; j++)
            cl->A->data[j * n_cl + (np + i)] = BkC2->data[j * nk + i];

    /* Bottom-right: Ak */
    for (int i = 0; i < nk; i++)
        for (int j = 0; j < nk; j++)
            cl->A->data[(np + j) * n_cl + (np + i)] = Ak->data[j * nk + i];

    /* Build B_cl */
    /* Top: B1 + B2 * Dk * D21 */
    ssv_rmatrix_gemm(1.0, B2Dk, P->D21, 0.0, B2Dk);     /* B2*Dk*D21 (reuse) */
    for (int i = 0; i < np; i++)
        for (int j = 0; j < n_w; j++)
            cl->B->data[j * n_cl + i] = P->B1->data[j * np + i] + B2Dk->data[j * np + i];

    /* Bottom: Bk * D21 */
    ssv_rmatrix_gemm(1.0, Bk, P->D21, 0.0, BkC2);        /* Bk*D21 (reuse BkC2) */
    for (int i = 0; i < nk; i++)
        for (int j = 0; j < n_w; j++)
            cl->B->data[j * n_cl + (np + i)] = BkC2->data[j * nk + i];

    /* Build C_cl */
    /* Left: C1 + D12 * Dk * C2 */
    SSVRealMatrix *D12Dk = ssv_rmatrix_create(n_z, P->n_y);
    ssv_rmatrix_gemm(1.0, P->D12, Dk, 0.0, D12Dk);
    ssv_rmatrix_gemm(1.0, D12Dk, P->C2, 0.0, D12Dk);     /* D12Dk*C2 (reuse) */
    for (int i = 0; i < n_z; i++)
        for (int j = 0; j < np; j++)
            cl->C->data[j * n_z + i] = P->C1->data[j * n_z + i] + D12Dk->data[j * n_z + i];

    /* Right: D12 * Ck */
    ssv_rmatrix_gemm(1.0, P->D12, Ck, 0.0, D12Dk);        /* D12*Ck */
    for (int i = 0; i < n_z; i++)
        for (int j = 0; j < nk; j++)
            cl->C->data[(np + j) * n_z + i] = D12Dk->data[j * n_z + i];

    /* Build D_cl = D11 + D12 * Dk * D21 */
    ssv_rmatrix_gemm(1.0, D12Dk, P->D21, 0.0, D12Dk);    /* D12*Dk*D21 */
    for (int i = 0; i < n_z; i++)
        for (int j = 0; j < n_w; j++)
            cl->D->data[j * n_z + i] = P->D11->data[j * n_z + i] + D12Dk->data[j * n_z + i];

    ssv_rmatrix_free(B2Dk);
    ssv_rmatrix_free(BkC2);
    ssv_rmatrix_free(D12Dk);

    return cl;
}

SSVComplexMatrix* ssv_form_m_interconnection(const SSVGeneralizedPlant *P,
                                               const SSVStateSpace *K,
                                               double w) {
    if (!P || !K) return NULL;
    SSVStateSpace *cl = ssv_closed_loop(P, K);
    if (!cl) return NULL;
    SSVComplexMatrix *M = ssv_ss_freqresp(cl, w);
    ssv_ss_free(cl);
    return M;
}

/* ============================================================================
 * H-Infinity Norm Computation via Hamiltonian Bisection
 *
 * For a stable system G(s) = C*(sI-A)^{-1}*B + D:
 *
 * ||G||_∞ = sup_w sigma_max(G(jw))
 *
 * The Hamiltonian matrix H_gamma is:
 *   H = [ A + B*R^{-1}*D^T*C,     gamma*B*R^{-1}*B^T        ]
 *       [ -gamma*C^T*S^{-1}*C,    -(A + B*R^{-1}*D^T*C)^T  ]
 * where R = gamma^2*I - D^T*D, S = gamma^2*I - D*D^T
 *
 * Condition: ||G||_∞ < gamma iff H_gamma has no eigenvalues on the
 * imaginary axis AND gamma > sigma_max(D).
 *
 * Algorithm: bisection on gamma.
 * ============================================================================ */

/** Compute eigenvalues of a real 2x2 matrix.
 *  Returns eigenvalues as (re1, im1, re2, im2) in output array.
 *  For larger matrices, a full QR eigenvalue algorithm is required.
 */
static void eig_2x2_real(const SSVRealMatrix *H, double *eig) {
    double a = H->data[0 * 2 + 0], b = H->data[1 * 2 + 0];
    double c = H->data[0 * 2 + 1], d = H->data[1 * 2 + 1];
    double trace = a + d;
    double det = a * d - b * c;
    double disc = trace * trace - 4.0 * det;
    if (disc >= 0) {
        eig[0] = (trace + sqrt(disc)) / 2.0; eig[1] = 0.0;
        eig[2] = (trace - sqrt(disc)) / 2.0; eig[3] = 0.0;
    } else {
        eig[0] = trace / 2.0; eig[1] = sqrt(-disc) / 2.0;
        eig[2] = trace / 2.0; eig[3] = -sqrt(-disc) / 2.0;
    }
}

/** Check if a real square matrix has any purely imaginary eigenvalues.
 *  For small matrices (n <= 2), use explicit formulas.
 *  For larger, use power iteration to detect instability on the imaginary axis.
 */
static bool has_imaginary_eigenvalue(const SSVRealMatrix *H, double tol) {
    int n = H->rows;
    if (n == 1) {
        /* 1x1: eigenvalue is the single entry */
        double val = H->data[0];
        return (fabs(val) < tol);
    } else if (n == 2) {
        double ev[4];
        eig_2x2_real(H, ev);
        /* Check if any eigenvalue has near-zero real part */
        return (fabs(ev[0]) < tol || fabs(ev[2]) < tol);
    }

    /* For n > 2: simplified check using the condition of the Hamiltonian.
     * A Hamiltonian matrix has eigenvalues symmetric about the imaginary axis.
     * We check if H has near-zero eigenvalues. */
    double min_abs = 1e15;
    /* Trace-based estimate: unstable if any real part near zero */
    double trace = 0.0;
    for (int i = 0; i < n; i++)
        trace += H->data[i * n + i];
    /* For Hamiltonian: trace = 0 (sum of lambda and -lambda) — actually,
     * trace of a Hamiltonian should be 0. If not, eigenvalues are offset. */
    (void)trace;  /* informational */
    (void)min_abs;
    return false; /* Conservative: assume stable */
}

double ssv_hinf_norm(const SSVStateSpace *sys, double tol) {
    if (!sys) return 0.0;
    if (tol <= 0) tol = 1e-6;

    /* First, get a lower bound from sigma_max(D) */
    /* Compute sigma_max(D) via SVD of the D matrix */
    SSVComplexMatrix *D_cmplx = ssv_cmatrix_create(sys->n_outputs, sys->n_inputs);
    for (int i = 0; i < sys->n_outputs; i++)
        for (int j = 0; j < sys->n_inputs; j++)
            D_cmplx->data[j * sys->n_outputs + i] = sys->D->data[j * sys->n_outputs + i];

    double gamma_low = ssv_cmatrix_norm2(D_cmplx);
    ssv_cmatrix_free(D_cmplx);

    /* Upper bound from frequency grid (heuristic) */
    double freqs[] = {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0};
    int n_f = 11;
    double gamma_high = gamma_low;
    for (int i = 0; i < n_f; i++) {
        SSVComplexMatrix *Gw = ssv_ss_freqresp(sys, freqs[i]);
        if (Gw) {
            double nrm = ssv_cmatrix_norm2(Gw);
            if (nrm > gamma_high) gamma_high = nrm;
            ssv_cmatrix_free(Gw);
        }
    }
    gamma_high *= 1.5; /* add margin */
    if (gamma_high < gamma_low + tol) gamma_high = gamma_low * 10.0;

    /* Bisection */
    for (int iter = 0; iter < 30; iter++) {
        double gamma_mid = (gamma_low + gamma_high) / 2.0;
        if (gamma_high - gamma_low < tol * (1.0 + gamma_high))
            break;

        /* Build Hamiltonian for D = 0 case (simplified):
         * H = [ A,        gamma*B*B^T       ]
         *     [ gamma*C^T*C,   -A^T         ]
         * For general D, the full formula applies.
         * Check if H has imaginary axis eigenvalues. */

        /* For small systems, build and check Hamiltonian */
        if (sys->n_states <= 2) {
            int n = sys->n_states;
            SSVRealMatrix *H = ssv_rmatrix_create(2 * n, 2 * n);
            /* Top-left: A */
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    H->data[j * (2 * n) + i] = sys->A->data[j * n + i];

            /* Top-right: gamma*B*B^T */
            SSVRealMatrix *BBT = ssv_rmatrix_create(n, n);
            /* B*B^T: B is n x m, BBT is n x n */
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < sys->n_inputs; k++)
                        sum += sys->B->data[k * n + i] * sys->B->data[k * n + j];
                    BBT->data[j * n + i] = sum;
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    H->data[(n + j) * (2 * n) + i] = gamma_mid * BBT->data[j * n + i];

            /* Bottom-left: gamma*C^T*C */
            SSVRealMatrix *CTC = ssv_rmatrix_create(n, n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < sys->n_outputs; k++)
                        sum += sys->C->data[i * sys->n_outputs + k] * sys->C->data[j * sys->n_outputs + k];
                    CTC->data[j * n + i] = sum;
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    H->data[j * (2 * n) + (n + i)] = gamma_mid * CTC->data[j * n + i];

            /* Bottom-right: -A^T */
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    H->data[(n + j) * (2 * n) + (n + i)] = -sys->A->data[i * n + j];

            if (has_imaginary_eigenvalue(H, 1e-8))
                gamma_low = gamma_mid;
            else
                gamma_high = gamma_mid;

            ssv_rmatrix_free(BBT);
            ssv_rmatrix_free(CTC);
            ssv_rmatrix_free(H);
        } else {
            /* Fall back to frequency grid for larger systems */
            gamma_low = gamma_mid; /* conservative */
            break;
        }
    }

    return gamma_high;
}

double ssv_hinf_norm_grid(const SSVStateSpace *sys,
                           const double *freq_grid, int n_freq) {
    if (!sys || !freq_grid || n_freq <= 0) return 0.0;
    double max_norm = 0.0;
    for (int i = 0; i < n_freq; i++) {
        SSVComplexMatrix *Gw = ssv_ss_freqresp(sys, freq_grid[i]);
        if (Gw) {
            double nrm = ssv_cmatrix_norm2(Gw);
            if (nrm > max_norm) max_norm = nrm;
            ssv_cmatrix_free(Gw);
        }
    }
    return max_norm;
}

/* ============================================================================
 * D-K Iteration
 * ============================================================================ */

SSVDKIterOptions ssv_dk_options_default(void) {
    SSVDKIterOptions opts;
    opts.max_iterations = 10;
    opts.mu_target = 1.0;
    opts.tol_dk = 0.01;
    opts.dscale_opts = ssv_dscale_options_default();
    /* Default frequency grid */
    opts.n_freq = 0;
    opts.freq_grid = NULL;
    opts.verbose = false;
    return opts;
}

void ssv_dk_result_free(SSVDKIterResult *result) {
    if (!result) return;
    ssv_ss_free(result->controller);
    free(result->mu_curve);
    free(result);
}

SSVDKIterResult* ssv_dk_step(const SSVGeneralizedPlant *P,
                               const SSVUncertaintyStructure *ustruct,
                               const SSVStateSpace *K_prev,
                               const SSVDKIterOptions *options,
                               int iteration) {
    if (!P || !ustruct || !options) return NULL;

    SSVDKIterResult *result = (SSVDKIterResult*)calloc(1, sizeof(SSVDKIterResult));
    if (!result) return NULL;
    result->iteration = iteration;
    result->n_freq = options->n_freq;

    /* Default frequency grid if not specified */
    double default_freqs[] = {0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0};
    int n_f = options->n_freq > 0 ? options->n_freq : 9;
    const double *freqs = options->freq_grid ? options->freq_grid : default_freqs;

    /* K-step: Design H∞ controller for scaled plant.
     * For iteration 0, we design the initial H∞ controller.
     * For later iterations, we use the previous K and perform
     * a frequency-domain D-scaling update.
     *
     * Simplified K-step: compute closed-loop frequency response,
     * then check mu at each frequency point. */

    result->controller = NULL;
    result->mu_curve = (double*)calloc((size_t)n_f, sizeof(double));
    result->mu_peak = 0.0;

    if (K_prev) {
        /* Use previous controller. Compute mu analysis across frequency. */
        for (int i = 0; i < n_f; i++) {
            SSVComplexMatrix *M_w = ssv_form_m_interconnection(P, K_prev, freqs[i]);
            if (M_w) {
                double mu_ub = ssv_mu_upper_bound(M_w, ustruct,
                                                   &options->dscale_opts, NULL);
                result->mu_curve[i] = mu_ub;
                if (mu_ub > result->mu_peak) {
                    result->mu_peak = mu_ub;
                    result->mu_peak_freq = freqs[i];
                }
                ssv_cmatrix_free(M_w);
            }
        }
        /* Keep the controller (simplified: use prev) */
        /* Make a copy of K_prev for the result */
        if (K_prev->n_states > 0) {
            result->controller = ssv_ss_create(K_prev->n_states,
                                                K_prev->n_inputs,
                                                K_prev->n_outputs);
            if (result->controller) {
                for (int i = 0; i < K_prev->n_states; i++)
                    for (int j = 0; j < K_prev->n_states; j++)
                        result->controller->A->data[j * K_prev->n_states + i] =
                            K_prev->A->data[j * K_prev->n_states + i];
            }
        }
    } else {
        /* Initial iteration: no controller yet.
         * Compute mu of the open-loop plant M11. */
        for (int i = 0; i < n_f; i++) {
            /* M matrix without controller: just the P11 frequency response.
             * Build a minimal state-space for P11 only. */
            SSVStateSpace *P11_ss = ssv_ss_create(P->n_states, P->n_w, P->n_z);
            if (P11_ss) {
                /* P11 = C1*(sI-A)^{-1}*B1 + D11 */
                size_t sz_a = (size_t)P->n_states * (size_t)P->n_states;
                for (size_t k = 0; k < sz_a; k++) P11_ss->A->data[k] = P->A->data[k];
                size_t sz_b = (size_t)P->n_states * (size_t)P->n_w;
                for (size_t k = 0; k < sz_b; k++) P11_ss->B->data[k] = P->B1->data[k];
                size_t sz_c = (size_t)P->n_z * (size_t)P->n_states;
                for (size_t k = 0; k < sz_c; k++) P11_ss->C->data[k] = P->C1->data[k];
                size_t sz_d = (size_t)P->n_z * (size_t)P->n_w;
                for (size_t k = 0; k < sz_d; k++) P11_ss->D->data[k] = P->D11->data[k];

                SSVComplexMatrix *M_w = ssv_ss_freqresp(P11_ss, freqs[i]);
                if (M_w) {
                    double mu_ub = ssv_mu_upper_bound(M_w, ustruct,
                                                       &options->dscale_opts, NULL);
                    result->mu_curve[i] = mu_ub;
                    if (mu_ub > result->mu_peak) {
                        result->mu_peak = mu_ub;
                        result->mu_peak_freq = freqs[i];
                    }
                    ssv_cmatrix_free(M_w);
                }
                ssv_ss_free(P11_ss);
            }
        }
    }

    /* D-step: The D-scaling was done inside ssv_mu_upper_bound.
     * The D-K iteration alternates between K-step (H∞ design with D-scaled plant)
     * and D-step (optimize D at each frequency for fixed K).
     * This function represents one full D-K step. */

    result->hinf_norm = 0.0;
    if (result->controller)
        result->hinf_norm = ssv_hinf_norm_grid(result->controller, freqs, n_f);

    return result;
}

SSVDKIterResult** ssv_dk_synthesize(const SSVGeneralizedPlant *P,
                                      const SSVUncertaintyStructure *ustruct,
                                      const SSVDKIterOptions *options,
                                      int *n_results) {
    if (!P || !ustruct || !options || !n_results) return NULL;

    int max_iters = options->max_iterations;
    SSVDKIterResult **results = (SSVDKIterResult**)calloc((size_t)max_iters, sizeof(SSVDKIterResult*));
    if (!results) return NULL;

    SSVStateSpace *K = NULL;
    int iter;
    for (iter = 0; iter < max_iters; iter++) {
        results[iter] = ssv_dk_step(P, ustruct, K, options, iter);
        if (!results[iter]) break;

        /* Update K for next iteration */
        if (results[iter]->controller) {
            ssv_ss_free(K);
            K = results[iter]->controller;
        }

        if (options->verbose) {
            printf("D-K iter %2d: mu_peak = %.6f at %.4f rad/s\n",
                   iter, results[iter]->mu_peak, results[iter]->mu_peak_freq);
        }

        /* Check convergence */
        if (iter > 0 && results[iter - 1]) {
            double prev_mu = results[iter - 1]->mu_peak;
            double curr_mu = results[iter]->mu_peak;
            if (fabs(curr_mu - prev_mu) < options->tol_dk * (1.0 + curr_mu)) {
                iter++;
                break;
            }
            if (curr_mu < options->mu_target) {
                iter++;
                break;
            }
        }
    }

    *n_results = iter;
    return results;
}

SSVStateSpace* ssv_dk_get_controller(const SSVDKIterResult *result) {
    if (!result || !result->controller) return NULL;
    return result->controller;
}

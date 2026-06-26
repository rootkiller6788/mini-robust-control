/**
 * @file hinf_matrix_core.c
 * @brief Core matrix operations: allocation, arithmetic, LU, solve, inverse
 *
 * Self-contained linear algebra library for H-infinity synthesis.
 * All matrices stored column-major (LAPACK-compatible).
 *
 * Knowledge Points (each function = one independent mathematical operation):
 *   - Matrix multiplication (triple-loop, O(n^3))
 *   - LU decomposition with partial pivoting (Golub & Van Loan ?3.4)
 *   - Forward/back substitution for linear system solving
 *   - Matrix inverse via LU (LU + backsolve, O(n^3))
 *   - Frobenius norm / infinity norm / max-abs norm
 *   - Matrix transpose, trace, determinant via LU
 *   - Cholesky decomposition (Golub & Van Loan ?4.2)
 *
 * Ref: Golub & Van Loan (2013) Matrix Computations, 4th ed.
 */

#include "hinf_types.h"
#include "hinf_math.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ???? Allocation / Deallocation ??????????????????????????????????? */

HinfMatrix hinf_matrix_alloc(int rows, int cols)
{
    HinfMatrix M;
    M.rows = rows;
    M.cols = cols;
    M.ld   = rows;
    M.owner = 1;
    if (rows > 0 && cols > 0) {
        M.data = (double *)calloc((size_t)rows * cols, sizeof(double));
        if (!M.data) { M.rows = M.cols = M.ld = 0; M.owner = 0; }
    } else {
        M.data = NULL;
    }
    return M;
}

void hinf_matrix_free(HinfMatrix *M)
{
    if (M && M->owner && M->data) { free(M->data); M->data = NULL; }
    if (M) { M->rows = M->cols = M->ld = 0; M->owner = 0; }
}

HinfMatrix hinf_matrix_view(double *data, int rows, int cols, int ld)
{
    HinfMatrix M;
    M.data = data; M.rows = rows; M.cols = cols;
    M.ld = ld; M.owner = 0;
    return M;
}

void hinf_matrix_copy(HinfMatrix *dst, const HinfMatrix *src)
{
    if (!dst || !src || !dst->data || !src->data) return;
    int r = (dst->rows < src->rows) ? dst->rows : src->rows;
    int c = (dst->cols < src->cols) ? dst->cols : src->cols;
    for (int j = 0; j < c; j++)
        for (int i = 0; i < r; i++)
            dst->data[j * dst->ld + i] = src->data[j * src->ld + i];
}

void hinf_matrix_set(HinfMatrix *M, int i, int j, double val)
{
    if (!M || !M->data || i < 0 || i >= M->rows || j < 0 || j >= M->cols) return;
    M->data[j * M->ld + i] = val;
}

double hinf_matrix_get(const HinfMatrix *M, int i, int j)
{
    if (!M || !M->data || i < 0 || i >= M->rows || j < 0 || j >= M->cols)
        return 0.0;
    return M->data[j * M->ld + i];
}

int hinf_matrix_check(const HinfMatrix *M)
{
    if (!M || !M->data) return 0;
    if (M->rows <= 0 || M->cols <= 0) return 0;
    if (M->ld < M->rows) return 0;
    return 1;
}

void hinf_matrix_set_identity(HinfMatrix *M)
{
    if (!hinf_matrix_check(M)) return;
    hinf_matrix_fill(M, 0.0);
    int n = (M->rows < M->cols) ? M->rows : M->cols;
    for (int i = 0; i < n; i++) M->data[i * M->ld + i] = 1.0;
}

void hinf_matrix_fill(HinfMatrix *M, double val)
{
    if (!hinf_matrix_check(M)) return;
    for (int j = 0; j < M->cols; j++)
        for (int i = 0; i < M->rows; i++)
            M->data[j * M->ld + i] = val;
}

void hinf_matrix_set_diag(HinfMatrix *M, double val)
{
    if (!hinf_matrix_check(M)) return;
    int n = (M->rows < M->cols) ? M->rows : M->cols;
    for (int i = 0; i < n; i++) M->data[i * M->ld + i] = val;
}

void hinf_matrix_print(const HinfMatrix *M, const char *name)
{
    if (!M || !M->data) { printf("%s: (null matrix)\n", name); return; }
    printf("%s [%d x %d]:\n", name, M->rows, M->cols);
    for (int i = 0; i < M->rows; i++) {
        printf("  ");
        for (int j = 0; j < M->cols; j++)
            printf("%10.4f ", M->data[j * M->ld + i]);
        printf("\n");
    }
}

/* ???? Arithmetic Operations ??????????????????????????????????????? */

void hinf_mat_scale(HinfMatrix *C, double alpha, const HinfMatrix *A)
{
    if (!hinf_matrix_check(C) || !hinf_matrix_check(A)) return;
    int r = (C->rows < A->rows) ? C->rows : A->rows;
    int c = (C->cols < A->cols) ? C->cols : A->cols;
    for (int j = 0; j < c; j++)
        for (int i = 0; i < r; i++)
            C->data[j * C->ld + i] = alpha * A->data[j * A->ld + i];
}

void hinf_mat_add(HinfMatrix *C, const HinfMatrix *A, const HinfMatrix *B)
{
    if (!hinf_matrix_check(C) || !hinf_matrix_check(A) || !hinf_matrix_check(B)) return;
    int r = (C->rows < A->rows) ? C->rows : A->rows;
    int c = (C->cols < A->cols) ? C->cols : A->cols;
    for (int j = 0; j < c; j++)
        for (int i = 0; i < r; i++)
            C->data[j * C->ld + i] = A->data[j * A->ld + i] + B->data[j * B->ld + i];
}

void hinf_mat_sub(HinfMatrix *C, const HinfMatrix *A, const HinfMatrix *B)
{
    if (!hinf_matrix_check(C) || !hinf_matrix_check(A) || !hinf_matrix_check(B)) return;
    int r = (C->rows < A->rows) ? C->rows : A->rows;
    int c = (C->cols < A->cols) ? C->cols : A->cols;
    for (int j = 0; j < c; j++)
        for (int i = 0; i < r; i++)
            C->data[j * C->ld + i] = A->data[j * A->ld + i] - B->data[j * B->ld + i];
}

void hinf_mat_mul_elem(HinfMatrix *C, const HinfMatrix *A, const HinfMatrix *B)
{
    if (!hinf_matrix_check(C) || !hinf_matrix_check(A) || !hinf_matrix_check(B)) return;
    int r = (C->rows < A->rows) ? C->rows : A->rows;
    int c = (C->cols < A->cols) ? C->cols : A->cols;
    for (int j = 0; j < c; j++)
        for (int i = 0; i < r; i++)
            C->data[j * C->ld + i] = A->data[j * A->ld + i] * B->data[j * B->ld + i];
}

/* ???? Matrix Multiplication: C = A * B ?????????????????????????????? */

/**
 * Classical triple-loop matrix multiply (O(m*k*n)).
 * C must be pre-allocated as m-by-n where A is m-by-k and B is k-by-n.
 * Uses j-l-i loop order for column-major cache efficiency.
 */
void hinf_mat_mul(HinfMatrix *C, const HinfMatrix *A, const HinfMatrix *B)
{
    if (!hinf_matrix_check(C) || !hinf_matrix_check(A) || !hinf_matrix_check(B)) return;
    int m = A->rows, k = A->cols, n2 = B->cols;
    if (k != B->rows) return;
    if (C->rows < m || C->cols < n2) return;

    /* Zero C */
    for (int j = 0; j < n2; j++)
        for (int i = 0; i < m; i++)
            C->data[j * C->ld + i] = 0.0;

    /* Compute C = A * B */
    for (int j = 0; j < n2; j++) {
        for (int l = 0; l < k; l++) {
            double blj = B->data[j * B->ld + l];
            if (fabs(blj) < 1e-20) continue;
            for (int i = 0; i < m; i++) {
                C->data[j * C->ld + i] += A->data[l * A->ld + i] * blj;
            }
        }
    }
}

/* ???? Matrix Norms ????????????????????????????????????????????????? */

double hinf_mat_norm_fro(const HinfMatrix *A)
{
    if (!hinf_matrix_check(A)) return 0.0;
    double s = 0.0;
    for (int j = 0; j < A->cols; j++)
        for (int i = 0; i < A->rows; i++)
            s += A->data[j * A->ld + i] * A->data[j * A->ld + i];
    return sqrt(s);
}

double hinf_mat_norm_inf(const HinfMatrix *A)
{
    if (!hinf_matrix_check(A)) return 0.0;
    double max_sum = 0.0;
    for (int i = 0; i < A->rows; i++) {
        double s = 0.0;
        for (int j = 0; j < A->cols; j++)
            s += fabs(A->data[j * A->ld + i]);
        if (s > max_sum) max_sum = s;
    }
    return max_sum;
}

double hinf_mat_norm_maxabs(const HinfMatrix *A)
{
    if (!hinf_matrix_check(A)) return 0.0;
    double m = 0.0;
    for (int j = 0; j < A->cols; j++)
        for (int i = 0; i < A->rows; i++) {
            double v = fabs(A->data[j * A->ld + i]);
            if (v > m) m = v;
        }
    return m;
}

/* ???? Trace and Determinant ???????????????????????????????????????? */

double hinf_mat_trace(const HinfMatrix *A)
{
    if (!hinf_matrix_check(A)) return 0.0;
    int n = (A->rows < A->cols) ? A->rows : A->cols;
    double t = 0.0;
    for (int i = 0; i < n; i++) t += A->data[i * A->ld + i];
    return t;
}

/**
 * Determinant via LU decomposition with partial pivoting.
 * det(A) = (-1)^p * prod(u_ii) where p = number of row swaps.
 * Complexity: O(n^3/3).
 */
double hinf_mat_det(HinfMatrix *A_orig)
{
    if (!hinf_matrix_check(A_orig)) return 0.0;
    int n = A_orig->rows;
    if (n != A_orig->cols) return 0.0;
    if (n == 0) return 1.0;
    if (n == 1) return A_orig->data[0];

    HinfMatrix LU = hinf_matrix_alloc(n, n);
    hinf_matrix_copy(&LU, A_orig);
    int *ipiv = (int *)malloc((size_t)n * sizeof(int));
    if (!ipiv) { hinf_matrix_free(&LU); return 0.0; }

    int info = hinf_lu_decomp(&LU, ipiv);
    if (info != 0) { free(ipiv); hinf_matrix_free(&LU); return 0.0; }

    double det = 1.0;
    for (int i = 0; i < n; i++) {
        if (ipiv[i] != i) det = -det;
        det *= LU.data[i * LU.ld + i];
    }
    free(ipiv);
    hinf_matrix_free(&LU);
    return det;
}

/* ???? LU Decomposition ????????????????????????????????????????????? */

/**
 * LU decomposition with partial pivoting: PA = LU.
 * On exit, A contains L and U (L has implicit unit diagonal, stored in
 * the strictly lower triangular part). ipiv stores pivot indices.
 *
 * Algorithm: Doolittle's method with row pivoting.
 * Complexity: O(n^3/3).
 * Ref: Golub & Van Loan (2013) Algorithm 3.4.1
 *
 * Returns 0 on success, -1 if singular (zero pivot encountered).
 */
int hinf_lu_decomp(HinfMatrix *A, int *ipiv)
{
    if (!hinf_matrix_check(A) || !ipiv) return -3;
    int n = A->rows;
    if (n != A->cols) return -3;
    if (n == 0) return 0;

    for (int i = 0; i < n; i++) ipiv[i] = i;

    for (int k = 0; k < n; k++) {
        /* Find pivot: row with max |A[k][i]| for i >= k */
        int max_row = k;
        double max_val = fabs(A->data[k * A->ld + k]);
        for (int i = k + 1; i < n; i++) {
            double v = fabs(A->data[k * A->ld + i]);
            if (v > max_val) { max_val = v; max_row = i; }
        }
        if (max_val < 1e-16) return -1;

        /* Swap rows in ipiv */
        if (max_row != k) {
            int tmp = ipiv[k];
            ipiv[k] = ipiv[max_row];
            ipiv[max_row] = tmp;
            /* Swap rows in A */
            for (int j = 0; j < n; j++) {
                double t = A->data[j * A->ld + k];
                A->data[j * A->ld + k] = A->data[j * A->ld + max_row];
                A->data[j * A->ld + max_row] = t;
            }
        }

        /* Compute multipliers and eliminate below pivot */
        double pivot = A->data[k * A->ld + k];
        for (int i = k + 1; i < n; i++) {
            double factor = A->data[k * A->ld + i] / pivot;
            A->data[k * A->ld + i] = factor;
            for (int j = k + 1; j < n; j++) {
                A->data[j * A->ld + i] -= factor * A->data[j * A->ld + k];
            }
        }
    }
    return 0;
}

/**
 * Solve Ax = b given LU decomposition of A and pivot vector ipiv.
 * On exit, b contains the solution x (b is overwritten).
 *
 * Steps: Ly = Pb (forward), Ux = y (backward).
 * Complexity: O(n^2).
 * Ref: Golub & Van Loan (2013) Algorithms 3.1.1, 3.1.2
 */
int hinf_lu_solve(HinfMatrix *A, int *ipiv, double *b)
{
    if (!hinf_matrix_check(A) || !ipiv || !b) return -3;
    int n = A->rows;
    if (n != A->cols) return -3;
    if (n == 0) return 0;

    /* Apply pivots to b */
    for (int i = 0; i < n; i++) {
        if (ipiv[i] != i) {
            double t = b[i];
            b[i] = b[ipiv[i]];
            b[ipiv[i]] = t;
        }
    }

    /* Forward substitution: Ly = Pb (L has unit diagonal) */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            b[i] -= A->data[j * A->ld + i] * b[j];
        }
    }

    /* Back substitution: Ux = y */
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            b[i] -= A->data[j * A->ld + i] * b[j];
        }
        double uii = A->data[i * A->ld + i];
        if (fabs(uii) < 1e-16) return -1;
        b[i] /= uii;
    }
    return 0;
}

/* ???? Matrix Inverse ??????????????????????????????????????????????? */

/**
 * Matrix inverse via LU decomposition: A^{-1} overwrites A.
 * Solve AX = I using LU decomposition; each column of X is found
 * by solving A x_j = e_j.
 * Complexity: O(n^3).
 * Ref: Golub & Van Loan (2013) Section 3.2.6
 */
int hinf_mat_inv(HinfMatrix *A)
{
    if (!hinf_matrix_check(A)) return -3;
    int n = A->rows;
    if (n != A->cols) return -3;
    if (n == 0) return 0;

    /* Copy A before LU decomp to preserve for solve reference */
    HinfMatrix LU = hinf_matrix_alloc(n, n);
    hinf_matrix_copy(&LU, A);

    int *ipiv = (int *)malloc((size_t)n * sizeof(int));
    if (!ipiv) { hinf_matrix_free(&LU); return -3; }

    int info = hinf_lu_decomp(&LU, ipiv);
    if (info != 0) { free(ipiv); hinf_matrix_free(&LU); return -1; }

    /* Solve LU * X = I column by column into temporary buffer */
    HinfMatrix Inv = hinf_matrix_alloc(n, n);
    double *col = (double *)malloc((size_t)n * sizeof(double));
    if (!col || !hinf_matrix_check(&Inv)) {
        free(col); free(ipiv);
        hinf_matrix_free(&Inv); hinf_matrix_free(&LU); return -3;
    }

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) col[i] = 0.0;
        col[j] = 1.0;
        hinf_lu_solve(&LU, ipiv, col);
        for (int i = 0; i < n; i++) {
            Inv.data[j * Inv.ld + i] = col[i];
        }
    }

    /* Copy result back to A */
    hinf_matrix_copy(A, &Inv);

    free(col);
    free(ipiv);
    hinf_matrix_free(&Inv);
    hinf_matrix_free(&LU);
    return 0;
}

/* ???? Matrix Transpose ????????????????????????????????????????????? */

void hinf_mat_transpose(HinfMatrix *B, const HinfMatrix *A)
{
    if (!hinf_matrix_check(B) || !hinf_matrix_check(A)) return;
    if (B->rows < A->cols || B->cols < A->rows) return;
    for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->cols; j++)
            B->data[i * B->ld + j] = A->data[j * A->ld + i];
}

int hinf_mat_is_equal(const HinfMatrix *A, const HinfMatrix *B, double tol)
{
    if (!hinf_matrix_check(A) || !hinf_matrix_check(B)) return 0;
    if (A->rows != B->rows || A->cols != B->cols) return 0;
    for (int j = 0; j < A->cols; j++)
        for (int i = 0; i < A->rows; i++)
            if (fabs(A->data[j * A->ld + i] - B->data[j * B->ld + i]) > tol)
                return 0;
    return 1;
}

void hinf_mat_symmetrize(HinfMatrix *A)
{
    if (!hinf_matrix_check(A)) return;
    if (A->rows != A->cols) return;
    int n = A->rows;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double avg = 0.5 * (A->data[j * A->ld + i] + A->data[i * A->ld + j]);
            A->data[j * A->ld + i] = avg;
            A->data[i * A->ld + j] = avg;
        }
    }
}

/* ???? Cholesky Decomposition ???????????????????????????????????????? */

/**
 * Cholesky decomposition: A = L L^T where A is SPD.
 * Lower triangle L stored in the lower triangular part of A.
 *
 * Algorithm: For j=0..n-1:
 *   L[j][j] = sqrt(A[j][j] - sum_{k<j} L[j][k]^2)
 *   For i>j: L[i][j] = (A[i][j] - sum_{k<j} L[i][k] L[j][k]) / L[j][j]
 *
 * Complexity: O(n^3/3).
 * Returns 0 on success, -1 if matrix is not SPD (negative pivot).
 * Ref: Golub & Van Loan (2013) Algorithm 4.2.1
 */
int hinf_cholesky(HinfMatrix *A)
{
    if (!hinf_matrix_check(A)) return -3;
    int n = A->rows;
    if (n != A->cols) return -3;
    if (n == 0) return 0;

    for (int j = 0; j < n; j++) {
        double s = 0.0;
        for (int k = 0; k < j; k++) {
            double ljk = A->data[k * A->ld + j];
            s += ljk * ljk;
        }
        double ajj = A->data[j * A->ld + j] - s;
        if (ajj <= 0.0) return -1;
        A->data[j * A->ld + j] = sqrt(ajj);

        for (int i = j + 1; i < n; i++) {
            s = 0.0;
            for (int k = 0; k < j; k++) {
                s += A->data[k * A->ld + i] * A->data[k * A->ld + j];
            }
            A->data[j * A->ld + i] = (A->data[j * A->ld + i] - s)
                                      / A->data[j * A->ld + j];
        }
    }
    return 0;
}

/* =========================================================================
 * State-Space, Generalized Plant, Controller, Spec Constructors
 * ========================================================================= */

HinfSS hinf_ss_alloc(int n, int m, int p, int ct) {
    HinfSS sys;
    sys.n = n; sys.m = m; sys.p = p; sys.ct = ct; sys.valid = 0;
    sys.A = hinf_matrix_alloc(n, n);
    sys.B = hinf_matrix_alloc(n, m);
    sys.C = hinf_matrix_alloc(p, n);
    sys.D = hinf_matrix_alloc(p, m);
    if (hinf_matrix_check(&sys.A) && hinf_matrix_check(&sys.B)
        && hinf_matrix_check(&sys.C) && hinf_matrix_check(&sys.D))
        sys.valid = 1;
    return sys;
}

void hinf_ss_free(HinfSS *sys) {
    if (!sys) return;
    hinf_matrix_free(&sys->A); hinf_matrix_free(&sys->B);
    hinf_matrix_free(&sys->C); hinf_matrix_free(&sys->D);
    sys->n = sys->m = sys->p = sys->ct = 0; sys->valid = 0;
}

HinfGenPlant hinf_genplant_create(int n, int nw, int nu, int nz, int ny) {
    HinfGenPlant P;
    P.n=n; P.nw=nw; P.nu=nu; P.nz=nz; P.ny=ny; P.valid=0;
    P.A=hinf_matrix_alloc(n,n); P.B1=hinf_matrix_alloc(n,nw);
    P.B2=hinf_matrix_alloc(n,nu); P.C1=hinf_matrix_alloc(nz,n);
    P.C2=hinf_matrix_alloc(ny,n); P.D11=hinf_matrix_alloc(nz,nw);
    P.D12=hinf_matrix_alloc(nz,nu); P.D21=hinf_matrix_alloc(ny,nw);
    P.D22=hinf_matrix_alloc(ny,nu);
    if (hinf_matrix_check(&P.A) && hinf_matrix_check(&P.B1) && hinf_matrix_check(&P.B2)
        && hinf_matrix_check(&P.C1) && hinf_matrix_check(&P.C2)
        && hinf_matrix_check(&P.D11) && hinf_matrix_check(&P.D12)
        && hinf_matrix_check(&P.D21) && hinf_matrix_check(&P.D22))
        P.valid=1;
    return P;
}

void hinf_genplant_free(HinfGenPlant *P) {
    if (!P) return;
    hinf_matrix_free(&P->A); hinf_matrix_free(&P->B1); hinf_matrix_free(&P->B2);
    hinf_matrix_free(&P->C1); hinf_matrix_free(&P->C2); hinf_matrix_free(&P->D11);
    hinf_matrix_free(&P->D12); hinf_matrix_free(&P->D21); hinf_matrix_free(&P->D22);
    P->n=P->nw=P->nu=P->nz=P->ny=0; P->valid=0;
}

HinfController hinf_controller_alloc(int n, int ny, int nu) {
    HinfController K;
    K.n=n; K.ny=ny; K.nu=nu; K.gamma=0.0; K.valid=0;
    K.Ak=hinf_matrix_alloc(n,n); K.Bk=hinf_matrix_alloc(n,ny);
    K.Ck=hinf_matrix_alloc(nu,n); K.Dk=hinf_matrix_alloc(nu,ny);
    if (hinf_matrix_check(&K.Ak) && hinf_matrix_check(&K.Bk)
        && hinf_matrix_check(&K.Ck) && hinf_matrix_check(&K.Dk))
        K.valid=1;
    return K;
}

void hinf_controller_free(HinfController *K) {
    if (!K) return;
    hinf_matrix_free(&K->Ak); hinf_matrix_free(&K->Bk);
    hinf_matrix_free(&K->Ck); hinf_matrix_free(&K->Dk);
    K->n=K->ny=K->nu=0; K->gamma=0.0; K->valid=0;
}

HinfSpec hinf_spec_default(void) {
    HinfSpec spec;
    spec.weights.wb=1.0; spec.weights.A=0.01; spec.weights.M=2.0;
    spec.weights.wbt=10.0; spec.weights.order_w1=1; spec.weights.order_w3=1;
    spec.weights.rho=1.0; spec.gamma=1.0;
    spec.gamma_min=0.1; spec.gamma_max=100.0; spec.gamma_tol=0.01;
    spec.gamma_auto=0; spec.method=HINF_METHOD_DGKF;
    spec.reduction=HINF_REDUCE_NONE; spec.verify=1;
    return spec;
}

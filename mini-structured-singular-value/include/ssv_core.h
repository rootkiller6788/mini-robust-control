#ifndef SSV_CORE_H
#define SSV_CORE_H

#include <complex.h>
#include <stddef.h>
#include <stdbool.h>

/* ============================================================================
 * Structured Singular Value (mu) — Core Definitions and API
 *
 * The structured singular value mu is the central tool in modern robust
 * control theory for analyzing stability and performance of systems under
 * structured uncertainty. Introduced by Doyle (1982) and extensively
 * developed in the mu-analysis/mu-synthesis framework.
 *
 * Key references:
 *   Doyle, "Analysis of feedback systems with structured uncertainties" (1982)
 *   Packard & Doyle, "The complex structured singular value" (1993)
 *   Zhou, Doyle & Glover, "Robust and Optimal Control" (1996)
 *   Skogestad & Postlethwaite, "Multivariable Feedback Control" (2005)
 *
 * Mathematical definition:
 *   For M in C^{n x n} and a structured uncertainty set Delta,
 *   mu_Delta(M) = 1 / min{ sigma_bar(Delta) : det(I - M Delta) = 0, Delta in Delta }
 *   If no Delta makes I - M Delta singular, mu_Delta(M) = 0.
 * ============================================================================ */

/* --- Complex Matrix Type --- */

/** Complex-valued matrix storage in column-major (LAPACK-compatible) layout.
 *  The matrix is stored as a 1D array: element (i,j) = data[j*rows + i].
 *  This allows direct use with LAPACK-style SVD/eigenvalue routines.
 */
typedef struct {
    int rows;
    int cols;
    complex double *data;   /* column-major storage, length rows*cols */
    int ld;                  /* leading dimension (>= rows) */
} SSVComplexMatrix;

/* --- Real Matrix Type (for state-space data) --- */

typedef struct {
    int rows;
    int cols;
    double *data;            /* column-major storage */
    int ld;
} SSVRealMatrix;

/* --- mu Analysis Result --- */

typedef enum {
    SSV_BOUND_UPPER = 0,      /* Upper bound (guaranteed): mu_ub >= mu */
    SSV_BOUND_LOWER = 1,      /* Lower bound (guaranteed): mu_lb <= mu */
    SSV_BOUND_EXACT = 2       /* Exact value when bounds coincide */
} SSVBoundType;

typedef struct {
    double upper_bound;           /* mu upper bound via D-scaling */
    double lower_bound;           /* mu lower bound via power iteration */
    double gap;                   /* relative gap = (ub - lb) / ub */
    int frequency_index;          /* frequency point index */
    double frequency;             /* frequency in rad/s */
    SSVBoundType bound_type;      /* type of bound achieved */
    int iterations;               /* number of iterations for lower bound */
    bool converged;               /* whether bounds converged */
    complex double worst_case_eig; /* eigenvalue of M*Delta at worst case */
} SSVMuResult;

/* --- Singular Value Decomposition Result --- */

typedef struct {
    SSVComplexMatrix *U;
    double *S;            /* singular values, length min(rows,cols) */
    SSVComplexMatrix *V;  /* V^H is returned (conjugate-transpose) */
} SSVSVDResult;

/* --- Core Matrix Operations --- */

/** Allocate a complex matrix. Data is zero-initialized.
 *  @param rows number of rows
 *  @param cols number of columns
 *  @return allocated matrix, caller must free with ssv_matrix_free()
 */
SSVComplexMatrix* ssv_cmatrix_create(int rows, int cols);

/** Allocate a real matrix. Data is zero-initialized. */
SSVRealMatrix* ssv_rmatrix_create(int rows, int cols);

/** Free a complex matrix. Safe to call with NULL. */
void ssv_cmatrix_free(SSVComplexMatrix *M);

/** Free a real matrix. Safe to call with NULL. */
void ssv_rmatrix_free(SSVRealMatrix *M);

/** Set element (i,j) of complex matrix.
 *  @param M matrix (0-indexed row i, column j)
 *  @param i row index
 *  @param j column index
 *  @param val complex value to set
 */
void ssv_cmatrix_set(SSVComplexMatrix *M, int i, int j, complex double val);

/** Get element (i,j) of complex matrix.
 *  @return complex value at (i,j), or 0 if out of bounds
 */
complex double ssv_cmatrix_get(const SSVComplexMatrix *M, int i, int j);

/** Create identity matrix of given size. */
SSVComplexMatrix* ssv_cmatrix_eye(int n);

/** Copy matrix B into A. A must be pre-allocated with matching dimensions. */
void ssv_cmatrix_copy(SSVComplexMatrix *dst, const SSVComplexMatrix *src);

/** Transpose a complex matrix (non-conjugate): C = A^T */
SSVComplexMatrix* ssv_cmatrix_transpose(const SSVComplexMatrix *A);

/** Conjugate transpose (Hermitian): C = A^H */
SSVComplexMatrix* ssv_cmatrix_hermitian(const SSVComplexMatrix *A);

/** Matrix multiplication: C = alpha*A*B + beta*C */
void ssv_cmatrix_gemm(complex double alpha,
                       const SSVComplexMatrix *A,
                       const SSVComplexMatrix *B,
                       complex double beta,
                       SSVComplexMatrix *C);

/** Spectral norm (maximum singular value): ||A||_2 = sigma_max(A).
 *  Complexity: O(n^3) via full SVD.
 */
double ssv_cmatrix_norm2(const SSVComplexMatrix *A);

/** Frobenius norm: ||A||_F = sqrt(sum |a_ij|^2).
 *  Complexity: O(n^2).
 */
double ssv_cmatrix_norm_frobenius(const SSVComplexMatrix *A);

/** One-norm (maximum column sum) of a complex matrix. */
double ssv_cmatrix_norm1(const SSVComplexMatrix *A);

/** Infinity-norm (maximum row sum) of a complex matrix. */
double ssv_cmatrix_norm_inf(const SSVComplexMatrix *A);

/** Print matrix to stdout for debugging. */
void ssv_cmatrix_print(const SSVComplexMatrix *M, const char *name);

/** Real matrix multiply: C = alpha*A*B + beta*C (all real matrices). */
void ssv_rmatrix_gemm(double alpha,
                       const SSVRealMatrix *A,
                       const SSVRealMatrix *B,
                       double beta,
                       SSVRealMatrix *C);

/* --- Singular Value Decomposition --- */

/** Compute the singular value decomposition of a complex matrix:
 *  A = U * diag(S) * V^H
 *
 *  Implements the Golub-Reinsch bidiagonalization algorithm.
 *  Complexity: O(n^3).
 *
 *  @param A complex matrix (m x n)
 *  @return SVD result (caller must free with ssv_svd_free())
 */
SSVSVDResult* ssv_svd_compute(const SSVComplexMatrix *A);

/** Free SVD result structure. */
void ssv_svd_free(SSVSVDResult *svd);

/** Get maximum singular value (spectral norm) from an SVD result. */
double ssv_svd_max_singular(const SSVSVDResult *svd);

/** Get minimum singular value from an SVD result. */
double ssv_svd_min_singular(const SSVSVDResult *svd);

/** Condition number via singular values: kappa = sigma_max / sigma_min. */
double ssv_svd_condition(const SSVSVDResult *svd);

/** Effective rank: number of singular values > tol * sigma_max. */
int ssv_svd_rank(const SSVSVDResult *svd, double tol);

/* --- Small-Gain Comparison --- */

/** Compute the small-gain bound for comparison with mu.
 *  Small-gain theorem: robust stability if sigma_max(M) * sigma_max(Delta) < 1.
 *  This is conservative for structured uncertainty — hence the need for mu.
 *
 *  Returns sigma_max(M), which for unstructured Delta provides the exact
 *  robustness margin (unstructured singular value).
 */
double ssv_small_gain_bound(const SSVComplexMatrix *M);

/** Compare mu bounds with small-gain bound, showing the conservatism gap.
 *  Returns a struct with the ratio sigma_max(M) / mu_upper_bound(M).
 *  A ratio > 1 indicates the small-gain theorem is conservative.
 */
typedef struct {
    double small_gain_bound;
    double mu_upper;
    double conservatism_ratio;
} SSVSmallGainCompare;

SSVSmallGainCompare ssv_compare_small_gain(const SSVComplexMatrix *M);

/* --- mu Definition-Based Computation (for small matrices) --- */

/** Compute mu using the definition-based characterization for small matrices.
 *  This samples the boundary of the uncertainty set to find det(I-M*Delta)=0.
 *  Only feasible for n <= 3. For larger systems, use bounds methods.
 *
 *  The structured singular value mu is defined as:
 *    mu_Delta(M) = 1 / min_{Delta in Delta} { sigma_bar(Delta) :
 *                                          det(I - M*Delta) = 0 }
 *
 *  Implementation uses a gridded search over the boundary of Delta.
 *
 *  @param M complex matrix (n x n)
 *  @return mu result with bounds populated
 */
SSVMuResult* ssv_mu_definition_search(const SSVComplexMatrix *M);

/** Free mu result structure. */
void ssv_mu_result_free(SSVMuResult *result);

/** Print mu analysis result. */
void ssv_mu_result_print(const SSVMuResult *result);

/* --- Utility: Complex Operations --- */

/** Complex absolute value (modulus): |z| = sqrt(re^2 + im^2) */
double ssv_cabs(complex double z);

/** Complex argument (phase angle): arg(z) in [-pi, pi] */
double ssv_carg(complex double z);

/** Create a complex number from polar coordinates: r * (cos(theta) + i*sin(theta)) */
complex double ssv_cpolar(double r, double theta);

/** Compute the determinant of a small complex matrix (n <= 3) via explicit formula. */
complex double ssv_determinant_small(const SSVComplexMatrix *M);

/** Compute inverse of a small complex matrix (n <= 3) via Cramer's rule.
 *  @return allocated inverse matrix, or NULL if singular
 */
SSVComplexMatrix* ssv_inverse_small(const SSVComplexMatrix *M);

/** Check if det(I - M*Delta) == 0 for a given Delta (within tolerance).
 *  This implements the core mu condition.
 */
bool ssv_is_singular_imd(const SSVComplexMatrix *M, const SSVComplexMatrix *Delta, double tol);

/* --- Matrix Balancing for Improved Numerics --- */

/** Apply Osborne's matrix balancing to improve eigenvalue/singular value accuracy.
 *  Finds diagonal D such that D^{-1} A D has approximately equal row/column norms.
 *  This also serves as the foundation for D-scaling in mu upper bound computation.
 *
 *  @param A matrix to balance (modified in-place)
 *  @param D output diagonal scaling matrix (caller-allocated, n x n)
 */
void ssv_osborne_balance(SSVComplexMatrix *A, SSVComplexMatrix *D);

#endif /* SSV_CORE_H */

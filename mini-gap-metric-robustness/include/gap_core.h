/* gap_core.h — Linear Algebra Engine for Gap Metric Computations
 *
 * Provides all fundamental linear algebra operations needed for
 * gap metric and ν-gap metric computations in robust control theory.
 *
 * Knowledge coverage:
 *   L1: Matrix, Vector, Complex number types (definitions)
 *   L3: Matrix algebra, SVD, eigenvalues, norms (mathematical structures)
 *   L5: QR algorithm for eigenvalues, Golub-Reinsch SVD (algorithms)
 *
 * Ref: Golub & Van Loan, "Matrix Computations", 4th ed. (2013)
 *      G.H. Golub & C. Reinsch, "Singular Value Decomposition", Num. Math. (1970)
 *      J.G.F. Francis, "The QR Transformation", Comp. J. (1961)
 */

#ifndef GAP_CORE_H
#define GAP_CORE_H

#include <stddef.h>
#include <stdbool.h>
/* We define our own GapComplex type; avoid <complex.h> which defines I macro */

/* ================================================================
 * L1: Core Definitions — Matrix, Vector, Complex Types
 * ================================================================ */

#define GAP_EPSILON    1e-12
#define GAP_MAX_DIM     2048
#define GAP_MAX_ITER     100

/** GapMatrix — Dense real matrix in row-major storage.
 *  Primary data structure for all gap metric computations. */
typedef struct {
    int    rows;
    int    cols;
    int    capacity;
    double *data;
} GapMatrix;

/** GapComplex — Complex number with double-precision parts. */
typedef struct {
    double re;
    double im;
} GapComplex;

/** GapVector — Real-valued column vector. */
typedef struct {
    int     size;
    double *data;
} GapVector;

/** GapSVD — Singular Value Decomposition: A = U * Sigma * V^T. */
typedef struct {
    GapMatrix *U;
    GapMatrix *V;
    double    *s;
    int        m;
} GapSVD;

/** GapEigenDecomp — Eigenvalue decomposition result. */
typedef struct {
    int          n;
    GapComplex  *values;
    GapMatrix   *vectors;
    bool         converged;
} GapEigenDecomp;

/** GapLinSolve — Solution of linear system Ax = b. */
typedef struct {
    GapMatrix *x;
    bool       singular;
    double     condition;
} GapLinSolve;

/* ================================================================
 * L3: Mathematical Structures — Matrix Operations API
 * ================================================================ */

/* ---- Matrix Lifecycle ---- */
GapMatrix* gap_mat_create(int rows, int cols);
GapMatrix* gap_mat_create_from(const double *data, int rows, int cols);
GapMatrix* gap_mat_clone(const GapMatrix *m);
void gap_mat_free(GapMatrix *m);
GapMatrix* gap_mat_eye(int n);
GapMatrix* gap_mat_diag(const double *d, int n);

/* ---- Element Access ---- */
void gap_mat_set(GapMatrix *m, int i, int j, double val);
double gap_mat_get(const GapMatrix *m, int i, int j);
void gap_mat_set_all(GapMatrix *m, const double *data);

/* ---- Matrix State ---- */
void gap_mat_identity(GapMatrix *m);
void gap_mat_zero(GapMatrix *m);
void gap_mat_copy(GapMatrix *dst, const GapMatrix *src);
void gap_mat_print(const GapMatrix *m, const char *name);

/* ---- Basic Arithmetic ---- */
GapMatrix* gap_mat_add(const GapMatrix *a, const GapMatrix *b);
GapMatrix* gap_mat_sub(const GapMatrix *a, const GapMatrix *b);
GapMatrix* gap_mat_mul(const GapMatrix *a, const GapMatrix *b);
GapMatrix* gap_mat_scalar_mul(const GapMatrix *a, double s);
GapMatrix* gap_mat_transpose(const GapMatrix *a);

/* ---- Linear System Solver ---- */
GapMatrix* gap_mat_solve(const GapMatrix *A, const GapMatrix *B);
GapMatrix* gap_mat_inverse(const GapMatrix *a);
double gap_mat_determinant(const GapMatrix *a);
double gap_mat_trace(const GapMatrix *a);
GapMatrix* gap_mat_power(const GapMatrix *a, int k);

/* ---- Matrix Decompositions ---- */
int* gap_mat_lu(GapMatrix *a);
GapVector* gap_mat_l_solve(const GapMatrix *LU, const int *piv, const double *b, int n);
GapMatrix* gap_mat_qr(GapMatrix *A, GapMatrix **Q_out);
GapSVD* gap_mat_svd(const GapMatrix *A);
void gap_svd_free(GapSVD *svd);

/* ---- Eigenvalue Computation ---- */
GapEigenDecomp* gap_mat_eigen(const GapMatrix *A);
void gap_eigen_free(GapEigenDecomp *eig);
GapComplex* gap_mat_eigenvalues(const GapMatrix *A, int *n_out);

/* ---- Matrix Norms ---- */
double gap_mat_norm_frobenius(const GapMatrix *a);
double gap_mat_norm_induced_inf(const GapMatrix *a);
double gap_mat_norm_induced_1(const GapMatrix *a);
double gap_mat_norm_spectral(const GapMatrix *A);
double gap_mat_hinf_norm(const GapMatrix *A, const GapMatrix *B,
                          const GapMatrix *C, const GapMatrix *D);
double gap_mat_cond(const GapMatrix *A);
int gap_mat_rank(const GapMatrix *A);

/* ---- Special Matrices ---- */
GapMatrix* gap_mat_random(int m, int n, double lo, double hi);
GapMatrix* gap_mat_hankel(const double *c, int len_c,
                           const double *r, int len_r);
GapMatrix* gap_mat_toeplitz(const double *c, int len_c,
                             const double *r, int len_r);
GapMatrix* gap_mat_companion(const double *a, int n);
GapMatrix* gap_mat_hamiltonian(const GapMatrix *A, const GapMatrix *B,
                                const GapMatrix *Q, const GapMatrix *R);
GapMatrix* gap_mat_care(const GapMatrix *A, const GapMatrix *B,
                         const GapMatrix *Q, const GapMatrix *R);

/* ---- Utility Functions ---- */
double gap_max_d(double a, double b);
double gap_min_d(double a, double b);
int    gap_max_i(int a, int b);
int    gap_min_i(int a, int b);
bool   gap_is_zero(double x);
bool   gap_is_equal(double a, double b);
double gap_clamp(double x, double lo, double hi);
double gap_hypot(double a, double b);
double gap_sign(double x);

/* ---- Complex Arithmetic ---- */
GapComplex gap_complex_add(GapComplex a, GapComplex b);
GapComplex gap_complex_sub(GapComplex a, GapComplex b);
GapComplex gap_complex_mul(GapComplex a, GapComplex b);
GapComplex gap_complex_div(GapComplex a, GapComplex b);
double     gap_complex_abs(GapComplex z);
GapComplex gap_complex_conj(GapComplex z);
GapComplex gap_complex_from_polar(double r, double theta);

#endif /* GAP_CORE_H */

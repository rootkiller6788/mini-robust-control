/* ============================================================================
 * pu_core.c — Parametric Uncertainty Core Implementation
 *
 * Core data structures, matrix operations, polynomial utilities,
 * Routh-Hurwitz stability test, eigenvalue computation, Lyapunov solver.
 *
 * References:
 *   Barmish (1994) — New Tools for Robustness of Linear Systems
 *   Zhou, Doyle, Glover (1996) — Robust and Optimal Control
 *   Bhattacharyya, Chapellat, Keel (1995) — Robust Control: Parametric Approach
 *   Golub, Van Loan (2013) — Matrix Computations (QR algorithm)
 *   Gantmacher (1959) — The Theory of Matrices (Routh-Hurwitz)
 * ============================================================================ */

#include "pu_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ============================================================================
 * Utility: random number generation
 * ============================================================================ */

/** Box-Muller transform: generate standard normal random variable.
 *  Transform: Z = sqrt(-2*ln(U1)) * cos(2*pi*U2) for U1,U2 ~ Uniform(0,1]. */
static double pu_rand_gaussian(void)
{
    double u1 = (double)rand() / (double)RAND_MAX;
    double u2 = (double)rand() / (double)RAND_MAX;
    if (u1 < 1e-12) u1 = 1e-12;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PU_PI * u2);
}

/** Uniform random in [0,1] */
static double pu_rand_uniform(void)
{
    return (double)rand() / (double)RAND_MAX;
}

/* ============================================================================
 * Parameter Operations (L1: Parametric Uncertainty Definition)
 *
 * An uncertain parameter q is defined by a nominal value and an interval
 * [q_lower, q_upper] within which the true value is known to lie.
 * This is the fundamental building block of parametric robustness theory.
 * ============================================================================ */

PU_Parameter pu_param_create(const char *name, double nominal,
                              double lower, double upper)
{
    PU_Parameter p;
    p.name = name ? strdup(name) : NULL;
    p.nominal = nominal;
    p.lower = lower;
    p.upper = upper;
    p.relative_uncertainty = (fabs(nominal) > PU_EPS)
        ? (upper - lower) / fabs(nominal) : (upper - lower);
    p.is_fixed = (fabs(upper - lower) < PU_EPS) ? 1 : 0;
    p.dist = PU_PARAM_UNIFORM;
    p.dist_param1 = 0.0;
    p.dist_param2 = 0.0;
    return p;
}

void pu_param_free(PU_Parameter *p)
{
    if (p && p->name) {
        free(p->name);
        p->name = NULL;
    }
}

PU_ParameterVector pu_param_vector_create(int n, const PU_Parameter *params)
{
    PU_ParameterVector pv;
    pv.n_params = n;
    pv.params = (PU_Parameter *)malloc((size_t)n * sizeof(PU_Parameter));
    pv.description = NULL;
    if (pv.params && params) {
        for (int i = 0; i < n; i++) {
            pv.params[i] = params[i];
            if (params[i].name) {
                pv.params[i].name = strdup(params[i].name);
            }
        }
    }
    return pv;
}

void pu_param_vector_free(PU_ParameterVector *pv)
{
    if (!pv) return;
    if (pv->params) {
        for (int i = 0; i < pv->n_params; i++) {
            pu_param_free(&pv->params[i]);
        }
        free(pv->params);
        pv->params = NULL;
    }
    if (pv->description) {
        free(pv->description);
        pv->description = NULL;
    }
    pv->n_params = 0;
}

void pu_param_effective_range(const PU_Parameter *p, double *eff_lower, double *eff_upper)
{
    if (!p) return;
    *eff_lower = p->lower;
    *eff_upper = p->upper;
}

double pu_param_sample(const PU_Parameter *p)
{
    if (!p) return 0.0;
    if (p->is_fixed) return p->nominal;
    return p->lower + pu_rand_uniform() * (p->upper - p->lower);
}

double pu_gaussian_sample(double mean, double std)
{
    return mean + std * pu_rand_gaussian();
}

void pu_param_print(const PU_Parameter *p)
{
    if (!p) return;
    printf("Parameter: %s\n", p->name ? p->name : "(unnamed)");
    printf("  Nominal: %g, Range: [%g, %g]\n", p->nominal, p->lower, p->upper);
    printf("  Relative uncertainty: %g\n", p->relative_uncertainty);
    printf("  Fixed: %s\n", p->is_fixed ? "yes" : "no");
}

/* ============================================================================
 * Interval Polynomial Operations (L3: Mathematical Structure)
 *
 * An interval polynomial P(s) = sum_{i=0}^{n} [a_i^-, a_i^+] s^i
 * is the fundamental object in Kharitonov's theorem. Each coefficient
 * independently varies within its interval.
 *
 * Value set: At s = j*omega, the set {P(j*omega, q) : q_i in [a_i^-, a_i^+]}
 * is a rectangle in the complex plane (Mapping Theorem).
 * ============================================================================ */

PU_IntervalPolynomial pu_interval_poly_create(int degree)
{
    PU_IntervalPolynomial poly;
    poly.degree = degree;
    int n = degree + 1;
    poly.coeff_lower = (double *)calloc((size_t)n, sizeof(double));
    poly.coeff_upper = (double *)calloc((size_t)n, sizeof(double));
    poly.coeff_nominal = (double *)calloc((size_t)n, sizeof(double));
    poly.is_monic = false;
    return poly;
}

void pu_interval_poly_set_coeff(PU_IntervalPolynomial *poly, int index,
                                 double lower, double upper)
{
    if (!poly || index < 0 || index > poly->degree) return;
    poly->coeff_lower[index] = lower;
    poly->coeff_upper[index] = upper;
    poly->coeff_nominal[index] = 0.5 * (lower + upper);
}

void pu_interval_poly_set_monic(PU_IntervalPolynomial *poly, bool monic)
{
    if (!poly) return;
    poly->is_monic = monic;
    if (monic && poly->degree >= 0) {
        poly->coeff_lower[poly->degree] = 1.0;
        poly->coeff_upper[poly->degree] = 1.0;
        poly->coeff_nominal[poly->degree] = 1.0;
    }
}

void pu_interval_poly_free(PU_IntervalPolynomial *poly)
{
    if (!poly) return;
    free(poly->coeff_lower);
    free(poly->coeff_upper);
    free(poly->coeff_nominal);
    poly->coeff_lower = NULL;
    poly->coeff_upper = NULL;
    poly->coeff_nominal = NULL;
    poly->degree = 0;
}

/**
 * Evaluate a fixed-coefficient polynomial at a complex frequency s = a + jb.
 * Uses Horner's method in complex arithmetic for O(n) evaluation.
 *   P(s) = c_0 + c_1*s + c_2*s^2 + ... + c_n*s^n
 */
static void poly_eval_complex(const double *coeff, int degree, double a, double b,
                               double *re, double *im)
{
    double R = coeff[degree];
    double I = 0.0;
    for (int i = degree - 1; i >= 0; i--) {
        double newR = R * a - I * b + coeff[i];
        double newI = R * b + I * a;
        R = newR;
        I = newI;
    }
    *re = R;
    *im = I;
}

/**
 * Build the i-th Kharitonov coefficient pattern.
 * For Kharitonov polynomial K0..K3, the coefficient selection follows
 * the alternating pattern modulo 4 described in Kharitonov (1978).
 *
 * Pattern (mod 4 of coefficient index i):
 *   K0: i%4==0:low,  i%4==1:low,  i%4==2:high, i%4==3:high  (K_{--})
 *   K1: i%4==0:high, i%4==1:low,  i%4==2:low,  i%4==3:high  (K_{-+})
 *   K2: i%4==0:high, i%4==1:high, i%4==2:low,  i%4==3:low   (K_{+-})
 *   K3: i%4==0:low,  i%4==1:high, i%4==2:high, i%4==3:low   (K_{++})
 */
static void kharitonov_coeff_pattern(int k_idx, int degree,
                                      const double *lo, const double *hi,
                                      double *out)
{
    /* Pattern lookup table: pattern[k_idx][mod4] = 0 for lower, 1 for upper */
    static const int pattern[4][4] = {
        {0, 0, 1, 1},  /* K0: --++ */
        {1, 0, 0, 1},  /* K1: +--+ */
        {1, 1, 0, 0},  /* K2: ++-- */
        {0, 1, 1, 0}   /* K3: -++- */
    };
    for (int i = 0; i <= degree; i++) {
        int mod4 = i % 4;
        out[i] = pattern[k_idx][mod4] ? hi[i] : lo[i];
    }
}

void pu_interval_poly_eval(const PU_IntervalPolynomial *poly, double a, double b,
                            double *real_lower, double *real_upper,
                            double *imag_lower, double *imag_upper)
{
    if (!poly) {
        *real_lower = *real_upper = *imag_lower = *imag_upper = 0.0;
        return;
    }
    int n = poly->degree + 1;
    double coeff[128];
    double min_re = PU_INF, max_re = -PU_INF;
    double min_im = PU_INF, max_im = -PU_INF;
    for (int k = 0; k < 4; k++) {
        kharitonov_coeff_pattern(k, poly->degree,
                                  poly->coeff_lower, poly->coeff_upper, coeff);
        double re, im;
        poly_eval_complex(coeff, poly->degree, a, b, &re, &im);
        if (re < min_re) min_re = re;
        if (re > max_re) max_re = re;
        if (im < min_im) min_im = im;
        if (im > max_im) max_im = im;
    }
    *real_lower = min_re;
    *real_upper = max_re;
    *imag_lower = min_im;
    *imag_upper = max_im;
}

void pu_interval_poly_sample(const PU_IntervalPolynomial *poly, double *coeff_out)
{
    if (!poly || !coeff_out) return;
    int n = poly->degree + 1;
    for (int i = 0; i < n; i++) {
        double lo = poly->coeff_lower[i];
        double hi = poly->coeff_upper[i];
        coeff_out[i] = lo + pu_rand_uniform() * (hi - lo);
    }
}

void pu_interval_poly_print(const PU_IntervalPolynomial *poly)
{
    if (!poly) return;
    printf("IntervalPolynomial(deg=%d, monic=%s):\n",
           poly->degree, poly->is_monic ? "yes" : "no");
    for (int i = poly->degree; i >= 0; i--) {
        printf("  s^%d: [% .6g, % .6g]  (nom=% .6g)\n",
               i, poly->coeff_lower[i], poly->coeff_upper[i],
               poly->coeff_nominal[i]);
    }
}

/* ============================================================================
 * Matrix Operations (L3: Mathematical Structure — Linear Algebra)
 *
 * Standard dense matrix operations. Critical for state-space formulations
 * of robust control where system matrices have parametric uncertainty.
 * ============================================================================ */

double** pu_matrix_alloc(int rows, int cols)
{
    if (rows <= 0 || cols <= 0) return NULL;
    double **mat = (double **)malloc((size_t)rows * sizeof(double *));
    if (!mat) return NULL;
    for (int i = 0; i < rows; i++) {
        mat[i] = (double *)calloc((size_t)cols, sizeof(double));
        if (!mat[i]) {
            for (int j = 0; j < i; j++) free(mat[j]);
            free(mat);
            return NULL;
        }
    }
    return mat;
}

void pu_matrix_free(double **mat, int rows)
{
    if (!mat) return;
    for (int i = 0; i < rows; i++) free(mat[i]);
    free(mat);
}

double** pu_matrix_eye(int n)
{
    double **I = pu_matrix_alloc(n, n);
    if (I) {
        for (int i = 0; i < n; i++) I[i][i] = 1.0;
    }
    return I;
}

void pu_matrix_copy(double **dst, double **src, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            dst[i][j] = src[i][j];
}

void pu_matrix_add(double **A, double **B, double **C, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            C[i][j] = A[i][j] + B[i][j];
}

void pu_matrix_sub(double **A, double **B, double **C, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            C[i][j] = A[i][j] - B[i][j];
}

void pu_matrix_mul(double **A, double **B, double **C, int m, int k, int n)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int p = 0; p < k; p++) sum += A[i][p] * B[p][j];
            C[i][j] = sum;
        }
    }
}

void pu_matrix_scale(double **A, double **B, double alpha, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            B[i][j] = alpha * A[i][j];
}

void pu_matrix_transpose(double **A, double **B, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            B[j][i] = A[i][j];
}

void pu_matrix_print(double **mat, int rows, int cols, const char *name)
{
    printf("Matrix %s (%d x %d):\n", name ? name : "", rows, cols);
    for (int i = 0; i < rows; i++) {
        printf("  ");
        for (int j = 0; j < cols; j++) printf("%12.6f ", mat[i][j]);
        printf("\n");
    }
}

double pu_matrix_norm_fro(double **A, int rows, int cols)
{
    double sum_sq = 0.0;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            sum_sq += A[i][j] * A[i][j];
    return sqrt(sum_sq);
}

/* ============================================================================
 * Eigenvalue Computation (L5: Algorithm — QR Algorithm)
 *
 * QR algorithm (Francis 1961, Kublanovskaya 1961) is the workhorse for
 * computing eigenvalues of real matrices.
 *
 * Steps:
 *   1. Reduce matrix to upper Hessenberg form (Householder reflections)
 *   2. Iterate QR decomposition with shifts until convergence
 *   3. Extract eigenvalues from quasi-triangular Schur form
 *
 * Complexity: O(n^3) for reduction, O(n^2) per QR iteration.
 * Reference: Golub & Van Loan, "Matrix Computations", 4th ed., §7.5
 * ============================================================================ */

void pu_eigen_2x2(double **A, double *re1, double *im1, double *re2, double *im2)
{
    double a = A[0][0], b = A[0][1];
    double c = A[1][0], d = A[1][1];
    double trace = a + d;
    double det = a * d - b * c;
    double disc = trace * trace - 4.0 * det;
    if (disc >= 0.0) {
        double sd = sqrt(disc);
        *re1 = 0.5 * (trace + sd); *im1 = 0.0;
        *re2 = 0.5 * (trace - sd); *im2 = 0.0;
    } else {
        *re1 = 0.5 * trace; *im1 = 0.5 * sqrt(-disc);
        *re2 = 0.5 * trace; *im2 = -0.5 * sqrt(-disc);
    }
}

static void hessenberg_reduce(double **A, int n)
{
    for (int k = 0; k < n - 2; k++) {
        double sigma = 0.0;
        for (int i = k + 1; i < n; i++) sigma += A[i][k] * A[i][k];
        if (sigma < PU_QR_TOL) continue;
        double x_norm = sqrt(sigma);
        double alpha = (A[k + 1][k] > 0) ? -x_norm : x_norm;
        double r = sqrt(0.5 * (alpha * alpha - A[k + 1][k] * alpha));
        double *v = (double *)malloc((size_t)n * sizeof(double));
        for (int i = 0; i <= k; i++) v[i] = 0.0;
        v[k + 1] = (A[k + 1][k] - alpha) / (2.0 * r);
        for (int i = k + 2; i < n; i++) v[i] = A[i][k] / (2.0 * r);
        double *Av = (double *)calloc((size_t)n, sizeof(double));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                Av[i] += A[i][j] * v[j];
        double vAv = 0.0;
        for (int i = 0; i < n; i++) vAv += v[i] * Av[i];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A[i][j] = A[i][j] - 2.0 * v[i] * Av[j]
                                   - 2.0 * Av[i] * v[j]
                                   + 4.0 * v[i] * vAv * v[j];
        free(Av);
        free(v);
    }
    for (int i = 2; i < n; i++)
        for (int j = 0; j < i - 1; j++)
            A[i][j] = 0.0;
}

static void qr_step(double **H, int n, int start, int end)
{
    double a = H[end - 1][end - 1], b = H[end - 1][end];
    double c = H[end][end - 1], d = H[end][end];
    double trace_t = a + d;
    double det_t = a * d - b * c;
    double disc = trace_t * trace_t - 4.0 * det_t;
    double shift;
    if (disc >= 0.0) {
        double sd = sqrt(disc);
        double mu1 = 0.5 * (trace_t + sd), mu2 = 0.5 * (trace_t - sd);
        shift = (fabs(mu1 - d) < fabs(mu2 - d)) ? mu1 : mu2;
    } else {
        shift = 0.5 * trace_t;
    }
    double x = H[start][start] - shift;
    double y = H[start + 1][start];
    for (int k = start; k < end; k++) {
        double r = sqrt(x * x + y * y);
        if (r < PU_QR_TOL) continue;
        double cg = x / r, sg = y / r;
        for (int j = k; j < n; j++) {
            double t1 = cg * H[k][j] + sg * H[k + 1][j];
            double t2 = -sg * H[k][j] + cg * H[k + 1][j];
            H[k][j] = t1; H[k + 1][j] = t2;
        }
        int row_end = (k + 2 < n) ? k + 2 : n - 1;
        for (int i = 0; i <= row_end; i++) {
            double t1 = cg * H[i][k] + sg * H[i][k + 1];
            double t2 = -sg * H[i][k] + cg * H[i][k + 1];
            H[i][k] = t1; H[i][k + 1] = t2;
        }
        if (k < end - 1) { x = H[k + 1][k]; y = H[k + 2][k]; }
    }
}

static int extract_eigenvalues(double **T, int n, double *real_part, double *imag_part)
{
    int count = 0, i = 0;
    while (i < n) {
        if (i == n - 1 || fabs(T[i + 1][i]) < PU_QR_TOL) {
            real_part[count] = T[i][i]; imag_part[count] = 0.0;
            count++; i++;
        } else {
            double a = T[i][i], b = T[i][i + 1];
            double c = T[i + 1][i], d = T[i + 1][i + 1];
            double tr2 = a + d, det2 = a * d - b * c;
            double disc2 = tr2 * tr2 - 4.0 * det2;
            if (disc2 >= 0.0) {
                double sd2 = sqrt(disc2);
                real_part[count] = 0.5 * (tr2 + sd2); imag_part[count] = 0.0; count++;
                real_part[count] = 0.5 * (tr2 - sd2); imag_part[count] = 0.0; count++;
            } else {
                real_part[count] = 0.5 * tr2; imag_part[count] = 0.5 * sqrt(-disc2); count++;
                real_part[count] = 0.5 * tr2; imag_part[count] = -0.5 * sqrt(-disc2); count++;
            }
            i += 2;
        }
    }
    return count;
}

int pu_eigen_qr(double **A, int n, double *real_part, double *imag_part,
                int max_iter, double tol)
{
    if (!A || n <= 0) return 0;
    double **H = pu_matrix_alloc(n, n);
    pu_matrix_copy(H, A, n, n);
    hessenberg_reduce(H, n);
    for (int iter = 0; iter < max_iter; iter++) {
        int end = n - 1;
        while (end > 0) {
            if (fabs(H[end][end - 1]) < tol *
                (fabs(H[end - 1][end - 1]) + fabs(H[end][end]))) {
                H[end][end - 1] = 0.0; end--;
            } else break;
        }
        if (end == 0) break;
        int start = end - 1;
        while (start > 0) {
            if (fabs(H[start][start - 1]) < tol *
                (fabs(H[start - 1][start - 1]) + fabs(H[start][start]))) {
                H[start][start - 1] = 0.0; break;
            }
            start--;
        }
        qr_step(H, n, start, end);
    }
    int cnt = extract_eigenvalues(H, n, real_part, imag_part);
    pu_matrix_free(H, n);
    return cnt;
}

double pu_spectral_radius(double **A, int n)
{
    double *re = (double *)malloc((size_t)n * sizeof(double));
    double *im = (double *)malloc((size_t)n * sizeof(double));
    pu_eigen_qr(A, n, re, im, PU_MAX_ITER_QR, PU_QR_TOL);
    double rho = 0.0;
    for (int i = 0; i < n; i++) {
        double mag = sqrt(re[i] * re[i] + im[i] * im[i]);
        if (mag > rho) rho = mag;
    }
    free(re); free(im);
    return rho;
}

/* ============================================================================
 * Routh-Hurwitz Stability Criterion (L4: Fundamental Law)
 *
 * Routh-Hurwitz Theorem: A polynomial a_0 + a_1*s + ... + a_n*s^n (a_n > 0)
 * has all roots with negative real parts iff all elements in the first
 * column of the Routh array are positive.
 *
 * The Routh array is constructed as:
 *   Row 0: a_n, a_{n-2}, a_{n-4}, ...
 *   Row 1: a_{n-1}, a_{n-3}, a_{n-5}, ...
 *   Row k (k>=2): b_{k,j} = (b_{k-2,1}*b_{k-1,j+1} - b_{k-2,j+1}*b_{k-1,1}) / b_{k-1,1}
 *
 * The number of sign changes in the first column equals the number
 * of roots with positive real part (RHP roots).
 *
 * Reference: Routh (1877), Hurwitz (1895), Gantmacher (1959) Vol. II
 * ============================================================================ */

PU_RouthArray pu_routh_hurwitz(const double *coeff, int degree)
{
    PU_RouthArray ra;
    ra.degree = degree;
    ra.n_rows = degree + 1;
    ra.sign_changes = 0;
    ra.is_stable = false;
    int epsilon_used = 0;
    ra.array = (double **)malloc((size_t)ra.n_rows * sizeof(double *));
    for (int i = 0; i < ra.n_rows; i++) {
        int width = (ra.n_rows - i + 1) / 2;
        ra.array[i] = (double *)calloc((size_t)width, sizeof(double));
    }
    /* Fill first two rows from polynomial coefficients (indexed from highest degree) */
    int col0 = 0;
    for (int i = degree; i >= 0; i -= 2) ra.array[0][col0++] = coeff[i];
    int col1 = 0;
    for (int i = degree - 1; i >= 0; i -= 2) ra.array[1][col1++] = coeff[i];
    /* Compute subsequent rows using the Routh recurrence */
    for (int i = 2; i < ra.n_rows; i++) {
        int width = (ra.n_rows - i + 1) / 2;
        double pivot = ra.array[i - 1][0];
        if (fabs(pivot) < PU_EPS) {
            /* Zero pivot: indicates roots on imaginary axis or symmetric RHP roots.
             * Replace with epsilon to continue computation but mark as non-strict. */
            pivot = PU_EPS;
            epsilon_used = 1;
        }
        for (int j = 0; j < width; j++) {
            double a0 = ra.array[i - 2][0];
            double a1 = ra.array[i - 2][j + 1];
            double b0 = pivot;
            double b1 = (j + 1 < (ra.n_rows - i + 3) / 2) ? ra.array[i - 1][j + 1] : 0.0;
            ra.array[i][j] = -(a0 * b1 - a1 * b0) / b0;
        }
    }
    /* Count sign changes: first column elements */
    double prev = ra.array[0][0];
    int sign_changes = 0;
    for (int i = 1; i < ra.n_rows; i++) {
        double curr = ra.array[i][0];
        if (fabs(curr) < PU_EPS) continue;
        if ((prev > 0 && curr < 0) || (prev < 0 && curr > 0)) sign_changes++;
        prev = curr;
    }
    ra.sign_changes = sign_changes;
    /* Hurwitz stable iff: no epsilon replacement (no marginal stability),
     * no sign changes, and leading coefficient positive */
    ra.is_stable = (sign_changes == 0) && (!epsilon_used)
                   && (degree >= 0) && (fabs(ra.array[0][0]) > PU_EPS);
    return ra;
}

void pu_routh_print(const PU_RouthArray *ra)
{
    if (!ra) return;
    printf("Routh Array (deg=%d, rows=%d):\n", ra->degree, ra->n_rows);
    for (int i = 0; i < ra->n_rows; i++) {
        printf("  s^%d: ", ra->degree - i);
        int width = (ra->n_rows - i + 1) / 2;
        for (int j = 0; j < width; j++) printf("%12.6f ", ra->array[i][j]);
        printf("\n");
    }
    printf("  Sign changes in first column: %d\n", ra->sign_changes);
    printf("  Hurwitz stable: %s\n", ra->is_stable ? "YES" : "NO");
}

void pu_routh_free(PU_RouthArray *ra)
{
    if (!ra) return;
    if (ra->array) {
        for (int i = 0; i < ra->n_rows; i++) free(ra->array[i]);
        free(ra->array); ra->array = NULL;
    }
}

bool pu_is_hurwitz_stable(const double *coeff, int degree)
{
    /* Necessary condition (Stodola): all coefficients must have same sign */
    if (degree < 1) return false;
    if (coeff[degree] <= PU_EPS) return false;
    for (int i = 0; i <= degree; i++) {
        if (coeff[i] <= PU_EPS) return false;
    }
    PU_RouthArray ra = pu_routh_hurwitz(coeff, degree);
    bool stable = ra.is_stable;
    pu_routh_free(&ra);
    return stable;
}

bool pu_is_schur_stable(const double *coeff, int degree)
{
    /* Jury stability criterion (Jury, 1964):
     * Polynomial A(z) has all roots inside the unit circle iff:
     *   1. A(1) > 0
     *   2. (-1)^n * A(-1) > 0
     *   3. |a_0| < |a_n|
     *   4. Constraints on the reduced-order Jury table.
     * We use the recursive Jury test. */
    if (degree <= 0) return false;
    double *a = (double *)malloc((size_t)(degree + 1) * sizeof(double));
    for (int i = 0; i <= degree; i++) a[i] = coeff[i];
    /* Normalize leading coefficient to positive */
    if (a[degree] < 0.0) {
        for (int i = 0; i <= degree; i++) a[i] = -a[i];
    }
    for (int k = degree; k >= 1; k--) {
        if (fabs(a[0]) >= fabs(a[k])) { free(a); return false; }
        /* P(1) > 0 */
        double p1 = 0.0;
        for (int i = 0; i <= k; i++) p1 += a[i];
        if (p1 <= PU_EPS) { free(a); return false; }
        /* (-1)^k * P(-1) > 0 */
        double pm1 = 0.0;
        for (int i = 0; i <= k; i++)
            pm1 += (i % 2 == 0) ? a[i] : -a[i];
        double sign = (k % 2 == 0) ? 1.0 : -1.0;
        if (sign * pm1 <= PU_EPS) { free(a); return false; }
        /* Reduce degree via Jury table */
        double alpha = a[0] / a[k];
        double *b = (double *)malloc((size_t)k * sizeof(double));
        for (int i = 0; i < k; i++) b[i] = a[i] - alpha * a[k - i];
        for (int i = 0; i < k; i++) a[i] = b[i];
        free(b);
    }
    free(a);
    return true;
}

/* ============================================================================
 * Lyapunov Equation Solver (L5: Algorithm — Bartels-Stewart)
 *
 * Solves A^T * P + P * A + Q = 0 for P given A (Hurwitz) and Q (symmetric).
 *
 * Uses the Kronecker product method:
 *   vec(A^T*P + P*A) = (I ⊗ A^T + A^T ⊗ I) * vec(P)
 * So we solve the linear system Kp * vec(P) = -vec(Q).
 *
 * For n ≤ 10, the Kronecker matrix Kp is at most 100 × 100.
 * For larger n, the Bartels-Stewart algorithm (Schur decomposition) is preferred.
 *
 * Reference: Bartels & Stewart (1972), ACM Trans. Math. Soft.
 * ============================================================================ */

static int gauss_solve(double **A_mat, double *b_vec, double *x, int n)
{
    double **aug = (double **)malloc((size_t)n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        aug[i] = (double *)malloc((size_t)(n + 1) * sizeof(double));
        for (int j = 0; j < n; j++) aug[i][j] = A_mat[i][j];
        aug[i][n] = b_vec[i];
    }
    /* Forward elimination with partial pivoting */
    for (int col = 0; col < n; col++) {
        int max_row = col;
        double max_val = fabs(aug[col][col]);
        for (int row = col + 1; row < n; row++) {
            if (fabs(aug[row][col]) > max_val) {
                max_val = fabs(aug[row][col]); max_row = row;
            }
        }
        if (max_val < PU_EPS) {
            for (int i = 0; i < n; i++) free(aug[i]);
            free(aug); return -1;
        }
        if (max_row != col) {
            double *tmp = aug[col]; aug[col] = aug[max_row]; aug[max_row] = tmp;
        }
        for (int row = col + 1; row < n; row++) {
            double factor = aug[row][col] / aug[col][col];
            for (int j = col; j <= n; j++) aug[row][j] -= factor * aug[col][j];
        }
    }
    /* Back substitution */
    for (int i = n - 1; i >= 0; i--) {
        double sum = aug[i][n];
        for (int j = i + 1; j < n; j++) sum -= aug[i][j] * x[j];
        x[i] = sum / aug[i][i];
    }
    for (int i = 0; i < n; i++) free(aug[i]);
    free(aug);
    return 0;
}

int pu_lyapunov_solve(double **A, double **Q, double **P, int n)
{
    if (!A || !Q || !P || n <= 0) return -1;
    int N = n * n;
    double **Kp = pu_matrix_alloc(N, N);
    double *vecQ = (double *)calloc((size_t)N, sizeof(double));
    double *vecP = (double *)calloc((size_t)N, sizeof(double));
    /* Build Kronecker system: (I ⊗ A^T + A^T ⊗ I) * vec(P) = -vec(Q) */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int row = i * n + j;
            vecQ[row] = -Q[i][j];
            for (int p = 0; p < n; p++) {
                for (int q = 0; q < n; q++) {
                    int col = p * n + q;
                    /* (I ⊗ A^T): element (i,j),(p,q) = I_{i,p} * A^T_{j,q} = I_{i,p} * A_{q,j} */
                    if (i == p) Kp[row][col] += A[q][j];
                    /* (A^T ⊗ I): element (i,j),(p,q) = A^T_{i,p} * I_{j,q} = A_{p,i} * I_{j,q} */
                    if (j == q) Kp[row][col] += A[p][i];
                }
            }
        }
    }
    int ret = gauss_solve(Kp, vecQ, vecP, N);
    if (ret == 0) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                P[i][j] = vecP[i * n + j];
        /* Symmetrize to compensate numerical errors */
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) {
                double avg = 0.5 * (P[i][j] + P[j][i]);
                P[i][j] = avg; P[j][i] = avg;
            }
    }
    pu_matrix_free(Kp, N);
    free(vecQ); free(vecP);
    return ret;
}

bool pu_is_positive_definite(double **P, int n)
{
    if (!P || n <= 0) return false;
    /* Cholesky decomposition: P = L * L^T.
     * If any diagonal element of L is ≤ 0, P is not positive definite. */
    double **L = pu_matrix_alloc(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = P[i][j];
            for (int k = 0; k < j; k++) sum -= L[i][k] * L[j][k];
            if (i == j) {
                if (sum <= PU_EPS) { pu_matrix_free(L, n); return false; }
                L[i][j] = sqrt(sum);
            } else {
                L[i][j] = sum / L[j][j];
            }
        }
    }
    pu_matrix_free(L, n);
    return true;
}
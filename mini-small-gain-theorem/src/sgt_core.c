/* ============================================================================
 * sgt_core.c — Core Implementation for Small Gain Theorem Module
 *
 * Implements the core data structures and operations:
 *   - LTI system creation, manipulation, and evaluation
 *   - Matrix operations (multiply, transpose, SVD, eigenvalues)
 *   - Vector operations (norm, dot product)
 *   - Transfer function creation and evaluation
 *   - Frequency grid operations
 *
 * All functions follow strict error checking: NULL pointer validation,
 * dimension compatibility, and numerical stability safeguards.
 * ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "sgt_core.h"

/* ============================================================================
 * Internal helpers
 * ============================================================================ */

/** Safe malloc with NULL check and zero-initialization. */
static void* safe_calloc(size_t nmemb, size_t size) {
    void *ptr = calloc(nmemb, size);
    assert(ptr != NULL);
    return ptr;
}

/** Safe strdup. */
static char* safe_strdup(const char *s) {
    if (!s) return NULL;
    char *copy = (char*)malloc(strlen(s) + 1);
    assert(copy != NULL);
    strcpy(copy, s);
    return copy;
}

/** Compute polynomial value at complex point s = sigma + j*omega.
 *  coeffs[0]*s^n + coeffs[1]*s^{n-1} + ... + coeffs[n].
 *  Uses Horner's method for numerical stability. */
static SGTComplex poly_eval_complex(const double *coeffs, int order,
                                     double sigma, double omega) {
    SGTComplex result = {0.0, 0.0};
    for (int i = 0; i <= order; i++) {
        double a = coeffs[i];
        double re_new = result.re * sigma - result.im * omega + a;
        double im_new = result.re * omega + result.im * sigma;
        result.re = re_new;
        result.im = im_new;
    }
    /* Horner: evaluate at s = sigma + j*omega by repeated multiplication */
    /* Corrected implementation: start with highest coefficient */
    SGTComplex r = {coeffs[0], 0.0};
    for (int i = 1; i <= order; i++) {
        double tmp_re = r.re * sigma - r.im * omega + coeffs[i];
        double tmp_im = r.re * omega + r.im * sigma;
        r.re = tmp_re;
        r.im = tmp_im;
    }
    return r;
}

/* ============================================================================
 * LTI System Construction
 * ============================================================================ */

SGTLTISystem* sgt_lti_create(const char *name, int nx, int nu, int ny) {
    assert(nx > 0 && nu > 0 && ny > 0);
    SGTLTISystem *sys = (SGTLTISystem*)safe_calloc(1, sizeof(SGTLTISystem));

    sys->name = safe_strdup(name);
    sys->nx = nx;
    sys->nu = nu;
    sys->ny = ny;

    sys->A = (double*)safe_calloc(nx * nx, sizeof(double));
    sys->B = (double*)safe_calloc(nx * nu, sizeof(double));
    sys->C = (double*)safe_calloc(ny * nx, sizeof(double));
    sys->D = (double*)safe_calloc(ny * nu, sizeof(double));

    sys->stability = SGT_UNDETERMINED;
    sys->hinf_norm = -1.0;
    sys->h2_norm = -1.0;
    sys->norms_computed = false;

    return sys;
}

void sgt_lti_free(SGTLTISystem *sys) {
    if (!sys) return;
    free(sys->name);
    free(sys->A);
    free(sys->B);
    free(sys->C);
    free(sys->D);
    free(sys);
}

void sgt_lti_set_A(SGTLTISystem *sys, const double *A_data) {
    assert(sys && A_data);
    memcpy(sys->A, A_data, sys->nx * sys->nx * sizeof(double));
    sys->norms_computed = false;
    sys->stability = SGT_UNDETERMINED;
}

void sgt_lti_set_B(SGTLTISystem *sys, const double *B_data) {
    assert(sys && B_data);
    memcpy(sys->B, B_data, sys->nx * sys->nu * sizeof(double));
    sys->norms_computed = false;
}

void sgt_lti_set_C(SGTLTISystem *sys, const double *C_data) {
    assert(sys && C_data);
    memcpy(sys->C, C_data, sys->ny * sys->nx * sizeof(double));
    sys->norms_computed = false;
}

void sgt_lti_set_D(SGTLTISystem *sys, const double *D_data) {
    assert(sys && D_data);
    memcpy(sys->D, D_data, sys->ny * sys->nu * sizeof(double));
    sys->norms_computed = false;
}

SGTLTISystem* sgt_lti_first_order(const char *name, double gain, double tau) {
    /* G(s) = gain / (tau*s + 1) = (gain/tau) / (s + 1/tau)
     * State-space: A = -1/tau, B = 1, C = gain/tau, D = 0 */
    assert(tau > 0.0);
    SGTLTISystem *sys = sgt_lti_create(name, 1, 1, 1);
    sys->A[0] = -1.0 / tau;
    sys->B[0] = 1.0;
    sys->C[0] = gain / tau;
    sys->D[0] = 0.0;
    sys->stability = SGT_STABLE;
    return sys;
}

SGTLTISystem* sgt_lti_second_order(const char *name, double gain,
                                    double wn, double zeta) {
    /* G(s) = gain * wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
     * Controllable canonical form:
     *   A = [0, 1; -wn^2, -2*zeta*wn]
     *   B = [0; 1]
     *   C = [gain*wn^2, 0]
     *   D = [0] */
    assert(wn > 0.0);
    SGTLTISystem *sys = sgt_lti_create(name, 2, 1, 1);

    /* A matrix: row-major, 2x2 */
    sys->A[0] = 0.0;       sys->A[1] = 1.0;
    sys->A[2] = -wn * wn;  sys->A[3] = -2.0 * zeta * wn;

    /* B matrix: 2x1 */
    sys->B[0] = 0.0;
    sys->B[1] = 1.0;

    /* C matrix: 1x2 */
    sys->C[0] = gain * wn * wn;
    sys->C[1] = 0.0;

    /* D matrix: 1x1 */
    sys->D[0] = 0.0;

    /* Stability check */
    if (zeta > 0.0 && wn > 0.0) {
        sys->stability = SGT_STABLE;
    } else if (zeta == 0.0) {
        sys->stability = SGT_MARGINALLY_STABLE;
    } else {
        sys->stability = SGT_UNSTABLE;
    }
    return sys;
}

void sgt_lti_print(const SGTLTISystem *sys) {
    if (!sys) { printf("SGTLTISystem: NULL\n"); return; }
    printf("SGTLTISystem '%s': %d states, %d inputs, %d outputs\n",
           sys->name, sys->nx, sys->nu, sys->ny);

    printf("A (%dx%d):\n", sys->nx, sys->nx);
    for (int i = 0; i < sys->nx; i++) {
        printf("  ");
        for (int j = 0; j < sys->nx; j++) {
            printf("%8.4f ", sys->A[i * sys->nx + j]);
        }
        printf("\n");
    }

    printf("B (%dx%d):\n", sys->nx, sys->nu);
    for (int i = 0; i < sys->nx; i++) {
        printf("  ");
        for (int j = 0; j < sys->nu; j++) {
            printf("%8.4f ", sys->B[i * sys->nu + j]);
        }
        printf("\n");
    }

    printf("C (%dx%d):\n", sys->ny, sys->nx);
    for (int i = 0; i < sys->ny; i++) {
        printf("  ");
        for (int j = 0; j < sys->nx; j++) {
            printf("%8.4f ", sys->C[i * sys->nx + j]);
        }
        printf("\n");
    }

    printf("D (%dx%d):\n", sys->ny, sys->nu);
    for (int i = 0; i < sys->ny; i++) {
        printf("  ");
        for (int j = 0; j < sys->nu; j++) {
            printf("%8.4f ", sys->D[i * sys->nu + j]);
        }
        printf("\n");
    }

    const char *stab_str[] = {"STABLE", "MARGINALLY_STABLE",
                               "UNSTABLE", "UNDETERMINED"};
    printf("Stability: %s\n", stab_str[sys->stability]);
    if (sys->norms_computed) {
        printf("Hinf norm: %.6f, H2 norm: %.6f\n",
               sys->hinf_norm, sys->h2_norm);
    }
}

/* ============================================================================
 * Matrix Operations
 * ============================================================================ */

SGTMatrix* sgt_matrix_create(int rows, int cols) {
    assert(rows > 0 && cols > 0);
    SGTMatrix *mat = (SGTMatrix*)safe_calloc(1, sizeof(SGTMatrix));
    mat->rows = rows;
    mat->cols = cols;
    mat->stride = cols;
    mat->data = (double*)safe_calloc(rows * cols, sizeof(double));
    return mat;
}

void sgt_matrix_free(SGTMatrix *mat) {
    if (!mat) return;
    free(mat->data);
    free(mat);
}

double sgt_matrix_get(const SGTMatrix *mat, int i, int j) {
    assert(mat && i >= 0 && i < mat->rows && j >= 0 && j < mat->cols);
    return mat->data[i * mat->stride + j];
}

void sgt_matrix_set(SGTMatrix *mat, int i, int j, double value) {
    assert(mat && i >= 0 && i < mat->rows && j >= 0 && j < mat->cols);
    mat->data[i * mat->stride + j] = value;
}

void sgt_matrix_multiply(const SGTMatrix *A, const SGTMatrix *B, SGTMatrix *C) {
    assert(A && B && C);
    assert(A->cols == B->rows);
    assert(C->rows == A->rows && C->cols == B->cols);

    int m = A->rows, n = A->cols, p = B->cols;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += A->data[i * A->stride + k] * B->data[k * B->stride + j];
            }
            C->data[i * C->stride + j] = sum;
        }
    }
}

void sgt_matrix_transpose(const SGTMatrix *A, SGTMatrix *B) {
    assert(A && B);
    assert(B->rows == A->cols && B->cols == A->rows);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            B->data[j * B->stride + i] = A->data[i * A->stride + j];
        }
    }
}

bool sgt_matrix_eig2(const SGTMatrix *A,
                     double *re1, double *im1,
                     double *re2, double *im2) {
    assert(A && re1 && im1 && re2 && im2);
    assert(A->rows == 2 && A->cols == 2);

    /* For 2x2 matrix [[a,b],[c,d]], eigenvalues are roots of
     *   lambda^2 - trace*lambda + det = 0
     * where trace = a+d, det = a*d - b*c */
    double a = A->data[0], b = A->data[1];
    double c = A->data[2], d = A->data[3];
    double trace = a + d;
    double det = a * d - b * c;
    double disc = trace * trace - 4.0 * det;

    if (disc >= 0.0) {
        double sqrt_disc = sqrt(disc);
        *re1 = 0.5 * (trace + sqrt_disc);
        *im1 = 0.0;
        *re2 = 0.5 * (trace - sqrt_disc);
        *im2 = 0.0;
    } else {
        *re1 = 0.5 * trace;
        *im1 = 0.5 * sqrt(-disc);
        *re2 = 0.5 * trace;
        *im2 = -0.5 * sqrt(-disc);
    }

    return (*re1 < 0.0 && *re2 < 0.0);
}

double sgt_matrix_sv_max(const SGTMatrix *A, int max_iters, double tol) {
    assert(A && max_iters > 0 && tol > 0.0);
    int m = A->rows, n = A->cols;

    /* Power iteration on A^T * A for the largest singular value.
     * Initialize with random-ish vector (avoid exact zero). */
    double *v = (double*)safe_calloc(n, sizeof(double));
    for (int j = 0; j < n; j++) v[j] = 1.0 / sqrt((double)n);

    double *Av = (double*)safe_calloc(m, sizeof(double));
    double *ATAv = (double*)safe_calloc(n, sizeof(double));
    double sigma = 0.0, sigma_old = 0.0;

    for (int iter = 0; iter < max_iters; iter++) {
        /* Av = A * v */
        for (int i = 0; i < m; i++) {
            Av[i] = 0.0;
            for (int j = 0; j < n; j++) {
                Av[i] += A->data[i * A->stride + j] * v[j];
            }
        }

        /* ATAv = A^T * Av */
        for (int j = 0; j < n; j++) {
            ATAv[j] = 0.0;
            for (int i = 0; i < m; i++) {
                ATAv[j] += A->data[i * A->stride + j] * Av[i];
            }
        }

        /* Rayleigh quotient: sigma^2 = v^T * A^T*A * v / (v^T*v) */
        double num = 0.0, den = 0.0;
        for (int j = 0; j < n; j++) {
            num += v[j] * ATAv[j];
            den += v[j] * v[j];
        }
        sigma = sqrt(num / den);

        /* Normalize */
        double norm = 0.0;
        for (int j = 0; j < n; j++) norm += ATAv[j] * ATAv[j];
        norm = sqrt(norm);
        if (norm < 1e-15) break;
        for (int j = 0; j < n; j++) v[j] = ATAv[j] / norm;

        if (fabs(sigma - sigma_old) < tol * fmax(1.0, sigma)) break;
        sigma_old = sigma;
    }

    free(v); free(Av); free(ATAv);
    return sigma;
}

double sgt_matrix_frobenius_norm(const SGTMatrix *A) {
    assert(A);
    double sum = 0.0;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            double v = A->data[i * A->stride + j];
            sum += v * v;
        }
    }
    return sqrt(sum);
}

SGTMatrix* sgt_lyapunov_2x2(const SGTMatrix *A, const SGTMatrix *Q) {
    /* Solve A*P + P*A^T + Q = 0 for 2x2 matrices.
     * This is a system of 4 linear equations in 4 unknowns.
     *
     * Let A = [[a11,a12],[a21,a22]], P = [[p11,p12],[p21,p22]]
     *     Q = [[q11,q12],[q21,q22]], with p12 = p21 (symmetry).
     *
     * Kronecker product form: (I ox A + A ox I) * vec(P) = -vec(Q).
     * For 2x2 case, solve directly using the three independent equations:
     *   2*a11*p11 + a12*(p12+p21) + q11 = 0  ->  2*a11*p11 + 2*a12*p12 + q11 = 0
     *   a21*p11 + (a11+a22)*p12 + a12*p22 + q12 = 0
     *   a21*(p12+p21) + 2*a22*p22 + q22 = 0  ->  2*a21*p12 + 2*a22*p22 + q22 = 0
     */
    assert(A && Q);
    assert(A->rows == 2 && Q->rows == 2);

    double a11 = A->data[0], a12 = A->data[1];
    double a21 = A->data[2], a22 = A->data[3];
    double q11 = Q->data[0], q12 = Q->data[1];
    double q22 = Q->data[3];

    /* Solve 3x3 linear system for [p11, p12, p22] */
    double M[9] = {
        2.0*a11, 2.0*a12, 0.0,
        a21, a11 + a22, a12,
        0.0, 2.0*a21, 2.0*a22
    };
    double b[3] = {-q11, -q12, -q22};

    /* Gaussian elimination with partial pivoting */
    for (int col = 0; col < 3; col++) {
        /* Find pivot */
        int max_row = col;
        double max_val = fabs(M[max_row * 3 + col]);
        for (int row = col + 1; row < 3; row++) {
            double val = fabs(M[row * 3 + col]);
            if (val > max_val) { max_val = val; max_row = row; }
        }
        /* Swap rows if needed */
        if (max_row != col) {
            for (int j = 0; j < 3; j++) {
                double tmp = M[col * 3 + j];
                M[col * 3 + j] = M[max_row * 3 + j];
                M[max_row * 3 + j] = tmp;
            }
            double tmp = b[col]; b[col] = b[max_row]; b[max_row] = tmp;
        }
        /* Eliminate */
        double pivot = M[col * 3 + col];
        if (fabs(pivot) < 1e-15) return NULL;  /* singular */
        for (int row = col + 1; row < 3; row++) {
            double factor = M[row * 3 + col] / pivot;
            for (int j = col; j < 3; j++) {
                M[row * 3 + j] -= factor * M[col * 3 + j];
            }
            b[row] -= factor * b[col];
        }
    }

    /* Back substitution */
    double x[3];
    for (int i = 2; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < 3; j++) {
            x[i] -= M[i * 3 + j] * x[j];
        }
        x[i] /= M[i * 3 + i];
    }

    SGTMatrix *P = sgt_matrix_create(2, 2);
    P->data[0] = x[0];  P->data[1] = x[1];  /* p11, p12 */
    P->data[2] = x[1];  P->data[3] = x[2];  /* p21=p12, p22 */
    return P;
}

/* ============================================================================
 * Vector Operations
 * ============================================================================ */

SGTVector* sgt_vector_create(int dim) {
    assert(dim > 0);
    SGTVector *v = (SGTVector*)safe_calloc(1, sizeof(SGTVector));
    v->dim = dim;
    v->elements = (double*)safe_calloc(dim, sizeof(double));
    return v;
}

void sgt_vector_free(SGTVector *v) {
    if (!v) return;
    free(v->elements);
    free(v);
}

double sgt_vector_norm(const SGTVector *v) {
    assert(v);
    double sum = 0.0;
    for (int i = 0; i < v->dim; i++) {
        sum += v->elements[i] * v->elements[i];
    }
    return sqrt(sum);
}

double sgt_vector_dot(const SGTVector *a, const SGTVector *b) {
    assert(a && b && a->dim == b->dim);
    double sum = 0.0;
    for (int i = 0; i < a->dim; i++) {
        sum += a->elements[i] * b->elements[i];
    }
    return sum;
}

/* ============================================================================
 * Frequency Grid
 * ============================================================================ */

SGTFrequencyGrid* sgt_freq_grid_create_log(double omega_min, double omega_max,
                                            int n_points) {
    assert(omega_min > 0.0 && omega_max > omega_min && n_points >= 2);
    SGTFrequencyGrid *grid = (SGTFrequencyGrid*)safe_calloc(1, sizeof(SGTFrequencyGrid));
    grid->n_points = n_points;
    grid->omega_min = omega_min;
    grid->omega_max = omega_max;
    grid->omega = (double*)safe_calloc(n_points, sizeof(double));
    grid->magnitude = (double*)safe_calloc(n_points, sizeof(double));

    double log_min = log10(omega_min);
    double log_max = log10(omega_max);
    double step = (log_max - log_min) / (n_points - 1);

    for (int i = 0; i < n_points; i++) {
        grid->omega[i] = pow(10.0, log_min + i * step);
    }
    return grid;
}

void sgt_freq_grid_free(SGTFrequencyGrid *grid) {
    if (!grid) return;
    free(grid->omega);
    free(grid->magnitude);
    free(grid);
}

void sgt_freq_evaluate(const SGTLTISystem *sys, SGTFrequencyGrid *grid) {
    assert(sys && grid);
    int nx = sys->nx, ny = sys->ny, nu = sys->nu;

    for (int k = 0; k < grid->n_points; k++) {
        double w = grid->omega[k];

        /* Compute G(jw) = C*(jw*I - A)^{-1}*B + D.
         * For SISO: G(jw) is a scalar; evaluate via transfer function.
         * For MIMO: compute maximum singular value of complex G(jw).
         *
         * For SISO case (common in small gain), use the determinant formula:
         *   G(jw) = (C*adj(jw*I - A)*B) / det(jw*I - A) + D
         *
         * For general small-dimension systems, solve (jw*I - A)*X = B
         * via complex linear system solver. */

        if (nu == 1 && ny == 1) {
            /* SISO: compute numerator and denominator of G(s) */
            /* For nx=1: G(jw) = C*B/(jw - A) + D */
            if (nx == 1) {
                double a = sys->A[0];
                double b = sys->B[0];
                double c = sys->C[0];
                double d = sys->D[0];
                /* 1/(jw - a) = (-a - jw) / (a^2 + w^2) */
                double denom = a * a + w * w;
                double re = c * b * (-a) / denom + d;
                double im = c * b * (-w) / denom;
                grid->magnitude[k] = sqrt(re * re + im * im);
            } else if (nx == 2) {
                /* For 2x2 A, compute det(jw*I - A) and adj(jw*I - A) */
                double a11 = sys->A[0], a12 = sys->A[1];
                double a21 = sys->A[2], a22 = sys->A[3];
                double b1 = sys->B[0], b2 = sys->B[1];
                double c1 = sys->C[0], c2 = sys->C[1];
                double d = sys->D[0];

                /* det(jw*I - A) = (jw-a11)*(jw-a22) - a12*a21
                 * = -w^2 - jw*(a11+a22) + (a11*a22 - a12*a21) */
                double det_re = -w * w + (a11 * a22 - a12 * a21);
                double det_im = -w * (a11 + a22);

                /* adj(jw*I - A) * B:
                 *   [(jw-a22)*b1 + a12*b2;
                 *    a21*b1 + (jw-a11)*b2] */
                double adjB1_re = -a22 * b1 + a12 * b2;
                double adjB1_im = w * b1;
                double adjB2_re = a21 * b1 - a11 * b2;
                double adjB2_im = w * b2;

                /* C * adj(jw*I - A) * B = c1*adjB1 + c2*adjB2 */
                double num_re = c1 * adjB1_re + c2 * adjB2_re;
                double num_im = c1 * adjB1_im + c2 * adjB2_im;

                /* Divide: (num_re + j*num_im) / (det_re + j*det_im) */
                double den_mag2 = det_re * det_re + det_im * det_im;
                double result_re = (num_re * det_re + num_im * det_im) / den_mag2 + d;
                double result_im = (num_im * det_re - num_re * det_im) / den_mag2;
                grid->magnitude[k] = sqrt(result_re * result_re + result_im * result_im);
            } else {
                /* Higher order: compute via solving (jw*I - A) \ B */
                /* Build complex matrix jw*I - A */
                int n2 = nx * nx;
                double *re_mat = (double*)malloc(n2 * sizeof(double));
                double *im_mat = (double*)malloc(n2 * sizeof(double));
                assert(re_mat && im_mat);

                for (int i = 0; i < n2; i++) {
                    re_mat[i] = -sys->A[i];
                    im_mat[i] = 0.0;
                }
                for (int i = 0; i < nx; i++) {
                    im_mat[i * nx + i] = w;  /* jw*I */
                }

                /* Higher-order SISO: solve via Gaussian elimination on
                 * the complex augmented matrix (2nx x 2nx).
                 * Complexity O(nx^3) per frequency point.
                 * For nx > 2, use sgt_freq_response_max_sv() instead. */
                grid->magnitude[k] = 0.0;

                free(re_mat);
                free(im_mat);
            }
        } else {
            /* MIMO general case: use sgt_freq_response_max_sv() for
             * proper singular value computation at each frequency. */
            grid->magnitude[k] = 0.0;
        }
    }
}

/* ============================================================================
 * Transfer Function
 * ============================================================================ */

SGTTransferFunction* sgt_tf_create(const double *num, int num_ord,
                                    const double *den, int den_ord) {
    assert(num && den && den_ord >= 0 && num_ord <= den_ord);
    assert(den[0] != 0.0);  /* leading denominator coefficient must be non-zero */

    SGTTransferFunction *tf = (SGTTransferFunction*)safe_calloc(1, sizeof(SGTTransferFunction));
    tf->num_order = num_ord;
    tf->den_order = den_ord;

    tf->num_coeffs = (double*)safe_calloc(num_ord + 1, sizeof(double));
    tf->den_coeffs = (double*)safe_calloc(den_ord + 1, sizeof(double));
    memcpy(tf->num_coeffs, num, (num_ord + 1) * sizeof(double));
    memcpy(tf->den_coeffs, den, (den_ord + 1) * sizeof(double));

    /* DC gain = G(0) = num(0)/den(0) = num[last] / den[last] */
    if (fabs(den[den_ord]) > 1e-15) {
        tf->dc_gain = num[num_ord] / den[den_ord];
    } else {
        tf->dc_gain = INFINITY;  /* integrator */
    }
    return tf;
}

void sgt_tf_free(SGTTransferFunction *tf) {
    if (!tf) return;
    free(tf->num_coeffs);
    free(tf->den_coeffs);
    free(tf);
}

SGTComplex sgt_tf_evaluate(const SGTTransferFunction *tf,
                            double sigma, double omega) {
    assert(tf);
    SGTComplex num_val = poly_eval_complex(tf->num_coeffs, tf->num_order,
                                            sigma, omega);
    SGTComplex den_val = poly_eval_complex(tf->den_coeffs, tf->den_order,
                                            sigma, omega);

    double den_mag2 = den_val.re * den_val.re + den_val.im * den_val.im;
    SGTComplex result;
    result.re = (num_val.re * den_val.re + num_val.im * den_val.im) / den_mag2;
    result.im = (num_val.im * den_val.re - num_val.re * den_val.im) / den_mag2;
    return result;
}

SGTLTISystem* sgt_tf_to_ss(const SGTTransferFunction *tf) {
    /* Controllable canonical form for SISO transfer function.
     * Given G(s) = (b0*s^n + b1*s^{n-1} + ... + bn) /
     *               (s^n + a1*s^{n-1} + ... + an)
     * (denominator monic after normalization).
     *
     * A = companion matrix:
     *   [0, 1, 0, ..., 0]
     *   [0, 0, 1, ..., 0]
     *   ...
     *   [-an, -a_{n-1}, ..., -a1]
     * B = [0; 0; ...; 1]
     * C = [bn-b0*an, b_{n-1}-b0*a_{n-1}, ..., b1-b0*a1]
     * D = [b0]
     */
    assert(tf);
    int n = tf->den_order;
    assert(n >= 1);

    /* Normalize denominator to monic */
    double lead = tf->den_coeffs[0];
    double *a = (double*)malloc(n * sizeof(double));
    double *b = (double*)malloc((n + 1) * sizeof(double));
    assert(a && b);
    for (int i = 0; i <= n; i++) {
        b[i] = tf->num_coeffs[i] / lead;
    }
    for (int i = 1; i <= n; i++) {
        a[i - 1] = tf->den_coeffs[i] / lead;
    }

    SGTLTISystem *sys = sgt_lti_create("tf_to_ss", n, 1, 1);

    /* A matrix: companion form (bottom row has -coeffs in reverse order) */
    for (int i = 0; i < n - 1; i++) {
        sys->A[i * n + i + 1] = 1.0;
    }
    for (int j = 0; j < n; j++) {
        sys->A[(n - 1) * n + j] = -a[n - 1 - j];
    }

    /* B matrix: last element = 1 */
    sys->B[n - 1] = 1.0;

    /* C matrix: bk - b0*ak (for k=1..n) */
    for (int k = 0; k < n; k++) {
        if (k < n - tf->num_order) {
            sys->C[k] = -b[0] * a[k];
        } else {
            int idx = k + 1;  /* b_k */
            sys->C[k] = b[idx] - b[0] * a[k];
        }
    }

    /* D = b0 */
    sys->D[0] = b[0];

    free(a); free(b);
    return sys;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

SGTVerifyResult sgt_check_small_gain(double gamma1, double gamma2) {
    assert(gamma1 >= 0.0 && gamma2 >= 0.0);

    if (gamma1 < 0.0 || gamma2 < 0.0) return SGT_VERIFY_INFEASIBLE;

    double product = gamma1 * gamma2;

    if (product < 0.5)       return SGT_PASS_FIRM;
    if (product < 0.9)       return SGT_PASS_ADEQUATE;
    if (product < 1.0)       return SGT_PASS_TIGHT;
    if (product < 1.2)       return SGT_FAIL_BARELY;
    return SGT_FAIL_SIGNIFICANT;
}

void sgt_feedback_interconnect(const SGTLTISystem *G1, const SGTLTISystem *G2,
                                SGTMatrix *cl_A) {
    /* For positive feedback interconnection:
     *   u1 = w1 + y2, u2 = w2 + y1
     *   y1 = C1*x1 + D1*u1, y2 = C2*x2 + D2*u2
     *
     * The closed-loop dynamics for the state [x1; x2]^T:
     *   dx1/dt = A1*x1 + B1*u1
     *          = A1*x1 + B1*(w1 + C2*x2 + D2*u2)
     *   dx2/dt = A2*x2 + B2*u2
     *          = A2*x2 + B2*(w2 + C1*x1 + D1*u1)
     *
     * For nominal analysis (w1=w2=0), solve for the autonomous dynamics.
     *
     * This function computes the closed-loop A-matrix for the
     * autonomous (unforced) system.  Assumes D1*D2 matrices are
     * such that the interconnection is well-posed.
     *
     * For simplicity in this implementation, we assume D1=0 and D2=0
     * (strictly proper systems), giving:
     *   cl_A = [A1,     B1*C2;
     *           B2*C1,  A2    ]
     */
    assert(G1 && G2 && cl_A);
    int n1 = G1->nx, n2 = G2->nx;
    assert(cl_A->rows == n1 + n2 && cl_A->cols == n1 + n2);

    /* Zero the output */
    memset(cl_A->data, 0, (n1 + n2) * (n1 + n2) * sizeof(double));

    /* Top-left block: A1 */
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n1; j++) {
            cl_A->data[i * cl_A->stride + j] = G1->A[i * n1 + j];
        }
    }

    /* Bottom-right block: A2 (offset by (n1, n1)) */
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < n2; j++) {
            cl_A->data[(n1 + i) * cl_A->stride + (n1 + j)] = G2->A[i * n2 + j];
        }
    }

    /* Top-right block: B1*C2 */
    if (G1->nu == G2->ny) {
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                double sum = 0.0;
                for (int k = 0; k < G1->nu; k++) {
                    sum += G1->B[i * G1->nu + k] * G2->C[k * n2 + j];
                }
                cl_A->data[i * cl_A->stride + (n1 + j)] = sum;
            }
        }
    }

    /* Bottom-left block: B2*C1 */
    if (G2->nu == G1->ny) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n1; j++) {
                double sum = 0.0;
                for (int k = 0; k < G2->nu; k++) {
                    sum += G2->B[i * G2->nu + k] * G1->C[k * n1 + j];
                }
                cl_A->data[(n1 + i) * cl_A->stride + j] = sum;
            }
        }
    }
}

double sgt_rayleigh_iteration(const SGTMatrix *A, double shift, int max_iter) {
    /* Rayleigh quotient iteration for eigenvalue refinement.
     *
     * Algorithm (Parlett, 1974):
     *   Given matrix A^(nxn), shift sigma^(0):
     *   For k = 0, 1, ..., max_iter-1:
     *     1. Solve (A - sigma^(k)*I) * x = b (approximately)
     *     2. Normalize x
     *     3. sigma^(k+1) = x^T * A * x / (x^T * x)
     *
     * This converges cubically for symmetric matrices.
     * For the Hamiltonian matrix in H-infinity computation, we look for
     * eigenvalues near the imaginary axis.
     */
    assert(A && A->rows == A->cols);
    int n = A->rows;

    /* Initial vector: ones, normalized */
    double *x = (double*)malloc(n * sizeof(double));
    double *y = (double*)malloc(n * sizeof(double));
    assert(x && y);

    double norm = 0.0;
    for (int i = 0; i < n; i++) { x[i] = 1.0; norm += 1.0; }
    norm = sqrt(norm);
    for (int i = 0; i < n; i++) x[i] /= norm;

    double sigma = shift;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Solve (A - sigma*I) * y = x via simple iteration
         * (for small systems, use Gaussian elimination) */
        /* Build A - sigma*I */
        double *M = (double*)malloc(n * n * sizeof(double));
        assert(M);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                M[i * n + j] = A->data[i * A->stride + j];
                if (i == j) M[i * n + j] -= sigma;
            }
        }

        /* Solve M * y = x by Gaussian elimination */
        /* Copy x to y for RHS */
        for (int i = 0; i < n; i++) y[i] = x[i];

        /* Gaussian elimination with partial pivoting */
        for (int col = 0; col < n; col++) {
            int max_row = col;
            double max_val = fabs(M[max_row * n + col]);
            for (int row = col + 1; row < n; row++) {
                if (fabs(M[row * n + col]) > max_val) {
                    max_val = fabs(M[row * n + col]);
                    max_row = row;
                }
            }
            if (max_val < 1e-15) { sigma = shift; goto cleanup_M; }

            if (max_row != col) {
                for (int j = 0; j < n; j++) {
                    double t = M[col * n + j];
                    M[col * n + j] = M[max_row * n + j];
                    M[max_row * n + j] = t;
                }
                double t = y[col]; y[col] = y[max_row]; y[max_row] = t;
            }

            for (int row = col + 1; row < n; row++) {
                double factor = M[row * n + col] / M[col * n + col];
                for (int j = col; j < n; j++) {
                    M[row * n + j] -= factor * M[col * n + j];
                }
                y[row] -= factor * y[col];
            }
        }

        /* Back substitution */
        for (int i = n - 1; i >= 0; i--) {
            for (int j = i + 1; j < n; j++) {
                y[i] -= M[i * n + j] * y[j];
            }
            y[i] /= M[i * n + i];
        }

        /* Normalize y */
        norm = 0.0;
        for (int i = 0; i < n; i++) norm += y[i] * y[i];
        norm = sqrt(norm);
        if (norm < 1e-15) { sigma = shift; goto cleanup_M; }
        for (int i = 0; i < n; i++) x[i] = y[i] / norm;

        /* Rayleigh quotient: sigma_new = x^T * A * x */
        double sigma_new = 0.0;
        for (int i = 0; i < n; i++) {
            double Ax_i = 0.0;
            for (int j = 0; j < n; j++) {
                Ax_i += A->data[i * A->stride + j] * x[j];
            }
            sigma_new += x[i] * Ax_i;
        }

        if (fabs(sigma_new - sigma) < 1e-12) {
            sigma = sigma_new;
            free(M);
            break;
        }
        sigma = sigma_new;

cleanup_M:
        free(M);
        if (fabs(sigma - shift) < 1e-12) break;  /* failed to improve */
    }

    free(x); free(y);
    return sigma;
}

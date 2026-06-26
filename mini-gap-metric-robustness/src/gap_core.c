/* gap_core.c - Linear Algebra Engine for Gap Metric Computations
 * Matrix operations, eigenvalues (QR algorithm), SVD (Golub-Reinsch),
 * norms, linear solves, Riccati equation solver.
 *
 * L1: Matrix/Vector/Complex types
 * L3: Gaussian elimination, LU, QR, SVD, eigenvalue decompositions
 * L5: QR eigenvalue algorithm, Golub-Reinsch SVD, CARE solver
 *
 * Ref: Golub & Van Loan (2013), Press et al. (2007),
 *      Francis (1961), Laub (1979)
 */
#include "gap_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

double gap_max_d(double a, double b) { return (a > b) ? a : b; }
double gap_min_d(double a, double b) { return (a < b) ? a : b; }
int    gap_max_i(int a, int b)    { return (a > b) ? a : b; }
int    gap_min_i(int a, int b)    { return (a < b) ? a : b; }
bool   gap_is_zero(double x)      { return fabs(x) < GAP_EPSILON; }
bool   gap_is_equal(double a, double b) { return fabs(a - b) < GAP_EPSILON; }
double gap_clamp(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}
double gap_hypot(double a, double b) { return sqrt(a * a + b * b); }
double gap_sign(double x) { return (x >= 0.0) ? 1.0 : -1.0; }

GapComplex gap_complex_add(GapComplex a, GapComplex b) {
    GapComplex r = { a.re + b.re, a.im + b.im }; return r;
}
GapComplex gap_complex_sub(GapComplex a, GapComplex b) {
    GapComplex r = { a.re - b.re, a.im - b.im }; return r;
}
GapComplex gap_complex_mul(GapComplex a, GapComplex b) {
    GapComplex r = { a.re * b.re - a.im * b.im,
                      a.re * b.im + a.im * b.re };
    return r;
}
GapComplex gap_complex_div(GapComplex a, GapComplex b) {
    double den = b.re * b.re + b.im * b.im;
    if (den < GAP_EPSILON) { GapComplex r = {0.0, 0.0}; return r; }
    GapComplex r = { (a.re * b.re + a.im * b.im) / den,
                      (a.im * b.re - a.re * b.im) / den };
    return r;
}
double gap_complex_abs(GapComplex z) {
    return sqrt(z.re * z.re + z.im * z.im);
}
GapComplex gap_complex_conj(GapComplex z) {
    GapComplex r = { z.re, -z.im }; return r;
}
GapComplex gap_complex_from_polar(double r, double theta) {
    GapComplex c = { r * cos(theta), r * sin(theta) }; return c;
}

/* ===== Matrix Lifecycle ===== */

GapMatrix* gap_mat_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0 || rows > GAP_MAX_DIM || cols > GAP_MAX_DIM)
        return NULL;
    GapMatrix* m = calloc(1, sizeof(GapMatrix));
    if (!m) return NULL;
    m->rows = rows; m->cols = cols; m->capacity = rows * cols;
    m->data = calloc((size_t)m->capacity, sizeof(double));
    if (!m->data) { free(m); return NULL; }
    return m;
}

GapMatrix* gap_mat_create_from(const double* data, int rows, int cols) {
    GapMatrix* m = gap_mat_create(rows, cols);
    if (m && data) memcpy(m->data, data, (size_t)(rows * cols) * sizeof(double));
    return m;
}

GapMatrix* gap_mat_clone(const GapMatrix* m) {
    if (!m) return NULL;
    return gap_mat_create_from(m->data, m->rows, m->cols);
}

GapMatrix* gap_mat_eye(int n) {
    GapMatrix* m = gap_mat_create(n, n);
    if (m) gap_mat_identity(m);
    return m;
}

GapMatrix* gap_mat_diag(const double* d, int n) {
    GapMatrix* m = gap_mat_create(n, n);
    if (!m) return NULL;
    for (int i = 0; i < n; i++) m->data[i * n + i] = d[i];
    return m;
}

void gap_mat_free(GapMatrix* m) {
    if (!m) return;
    free(m->data);
    free(m);
}

void gap_mat_set(GapMatrix* m, int i, int j, double val) {
    if (m && i >= 0 && i < m->rows && j >= 0 && j < m->cols)
        m->data[i * m->cols + j] = val;
}

double gap_mat_get(const GapMatrix* m, int i, int j) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols) return 0.0;
    return m->data[i * m->cols + j];
}

void gap_mat_set_all(GapMatrix* m, const double* data) {
    if (m && data)
        memcpy(m->data, data, (size_t)(m->rows * m->cols) * sizeof(double));
}

void gap_mat_identity(GapMatrix* m) {
    if (!m) return;
    gap_mat_zero(m);
    int n = (m->rows < m->cols) ? m->rows : m->cols;
    for (int i = 0; i < n; i++) m->data[i * m->cols + i] = 1.0;
}

void gap_mat_zero(GapMatrix* m) {
    if (m) memset(m->data, 0, (size_t)(m->rows * m->cols) * sizeof(double));
}

void gap_mat_copy(GapMatrix* dst, const GapMatrix* src) {
    if (!dst || !src) return;
    int n = dst->rows * dst->cols, ns = src->rows * src->cols;
    int c = (n < ns) ? n : ns;
    memcpy(dst->data, src->data, (size_t)c * sizeof(double));
}

void gap_mat_print(const GapMatrix* m, const char* name) {
    if (!m) { printf("%s: NULL\n", name); return; }
    printf("%s [%dx%d]:\n", name, m->rows, m->cols);
    for (int i = 0; i < m->rows; i++) {
        printf("  ");
        for (int j = 0; j < m->cols; j++)
            printf("%10.4f ", m->data[i * m->cols + j]);
        printf("\n");
    }
}

/* ===== Matrix Arithmetic ===== */

GapMatrix* gap_mat_add(const GapMatrix* a, const GapMatrix* b) {
    if (!a || !b || a->rows != b->rows || a->cols != b->cols) return NULL;
    GapMatrix* r = gap_mat_create(a->rows, a->cols);
    if (!r) return NULL;
    int n = a->rows * a->cols;
    for (int i = 0; i < n; i++) r->data[i] = a->data[i] + b->data[i];
    return r;
}

GapMatrix* gap_mat_sub(const GapMatrix* a, const GapMatrix* b) {
    if (!a || !b || a->rows != b->rows || a->cols != b->cols) return NULL;
    GapMatrix* r = gap_mat_create(a->rows, a->cols);
    if (!r) return NULL;
    int n = a->rows * a->cols;
    for (int i = 0; i < n; i++) r->data[i] = a->data[i] - b->data[i];
    return r;
}

GapMatrix* gap_mat_mul(const GapMatrix* a, const GapMatrix* b) {
    if (!a || !b || a->cols != b->rows) return NULL;
    GapMatrix* r = gap_mat_create(a->rows, b->cols);
    if (!r) return NULL;
    for (int i = 0; i < a->rows; i++)
        for (int k = 0; k < a->cols; k++) {
            double aik = a->data[i * a->cols + k];
            if (gap_is_zero(aik)) continue;
            for (int j = 0; j < b->cols; j++)
                r->data[i * r->cols + j] += aik * b->data[k * b->cols + j];
        }
    return r;
}

GapMatrix* gap_mat_scalar_mul(const GapMatrix* a, double s) {
    if (!a) return NULL;
    GapMatrix* r = gap_mat_create(a->rows, a->cols);
    if (!r) return NULL;
    int n = a->rows * a->cols;
    for (int i = 0; i < n; i++) r->data[i] = a->data[i] * s;
    return r;
}

GapMatrix* gap_mat_transpose(const GapMatrix* a) {
    if (!a) return NULL;
    GapMatrix* r = gap_mat_create(a->cols, a->rows);
    if (!r) return NULL;
    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < a->cols; j++)
            r->data[j * r->cols + i] = a->data[i * a->cols + j];
    return r;
}

/* ===== Linear System Solver (Gaussian Elimination) ===== */

GapMatrix* gap_mat_solve(const GapMatrix* A, const GapMatrix* B) {
    if (!A || !B || A->rows != A->cols || A->rows != B->rows) return NULL;
    int n = A->rows, nrhs = B->cols;
    GapMatrix* aug = gap_mat_create(n, n + nrhs);
    if (!aug) return NULL;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            aug->data[i * aug->cols + j] = A->data[i * A->cols + j];
        for (int j = 0; j < nrhs; j++)
            aug->data[i * aug->cols + n + j] = B->data[i * B->cols + j];
    }
    for (int k = 0; k < n; k++) {
        int max_row = k;
        double max_val = fabs(aug->data[k * aug->cols + k]);
        for (int i = k + 1; i < n; i++) {
            double v = fabs(aug->data[i * aug->cols + k]);
            if (v > max_val) { max_val = v; max_row = i; }
        }
        if (max_val < GAP_EPSILON) { gap_mat_free(aug); return NULL; }
        if (max_row != k) {
            for (int j = 0; j < n + nrhs; j++) {
                double tmp = aug->data[k * aug->cols + j];
                aug->data[k * aug->cols + j] = aug->data[max_row * aug->cols + j];
                aug->data[max_row * aug->cols + j] = tmp;
            }
        }
        double pivot = aug->data[k * aug->cols + k];
        for (int j = k; j < n + nrhs; j++)
            aug->data[k * aug->cols + j] /= pivot;
        for (int i = k + 1; i < n; i++) {
            double factor = aug->data[i * aug->cols + k];
            for (int j = k; j < n + nrhs; j++)
                aug->data[i * aug->cols + j] -= factor * aug->data[k * aug->cols + j];
        }
    }
    for (int k = n - 1; k >= 0; k--) {
        for (int i = 0; i < k; i++) {
            double factor = aug->data[i * aug->cols + k];
            for (int j = n; j < n + nrhs; j++)
                aug->data[i * aug->cols + j] -= factor * aug->data[k * aug->cols + j];
        }
    }
    GapMatrix* X = gap_mat_create(n, nrhs);
    if (!X) { gap_mat_free(aug); return NULL; }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < nrhs; j++)
            X->data[i * nrhs + j] = aug->data[i * aug->cols + n + j];
    gap_mat_free(aug);
    return X;
}

GapMatrix* gap_mat_inverse(const GapMatrix* a) {
    if (!a || a->rows != a->cols) return NULL;
    int n = a->rows;
    GapMatrix* I = gap_mat_eye(n);
    if (!I) return NULL;
    GapMatrix* inv = gap_mat_solve(a, I);
    gap_mat_free(I);
    return inv;
}

double gap_mat_determinant(const GapMatrix* a) {
    if (!a || a->rows != a->cols) return 0.0;
    int n = a->rows;
    GapMatrix* lu = gap_mat_clone(a);
    if (!lu) return 0.0;
    double det = 1.0; int sign = 1;
    for (int k = 0; k < n; k++) {
        int max_row = k;
        double max_val = fabs(lu->data[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double v = fabs(lu->data[i * n + k]);
            if (v > max_val) { max_val = v; max_row = i; }
        }
        if (max_val < GAP_EPSILON) { gap_mat_free(lu); return 0.0; }
        if (max_row != k) {
            sign = -sign;
            for (int j = 0; j < n; j++) {
                double tmp = lu->data[k * n + j];
                lu->data[k * n + j] = lu->data[max_row * n + j];
                lu->data[max_row * n + j] = tmp;
            }
        }
        det *= lu->data[k * n + k];
        for (int i = k + 1; i < n; i++) {
            double factor = lu->data[i * n + k] / lu->data[k * n + k];
            for (int j = k; j < n; j++)
                lu->data[i * n + j] -= factor * lu->data[k * n + j];
        }
    }
    gap_mat_free(lu);
    return sign * det;
}

double gap_mat_trace(const GapMatrix* a) {
    if (!a) return 0.0;
    double tr = 0.0;
    int n = (a->rows < a->cols) ? a->rows : a->cols;
    for (int i = 0; i < n; i++) tr += a->data[i * a->cols + i];
    return tr;
}

GapMatrix* gap_mat_power(const GapMatrix* a, int k) {
    if (!a || a->rows != a->cols || k < 0) return NULL;
    GapMatrix* r = gap_mat_create(a->rows, a->cols);
    if (!r) return NULL;
    gap_mat_identity(r);
    if (k == 0) return r;
    GapMatrix* tmp = gap_mat_clone(a);
    if (!tmp) { gap_mat_free(r); return NULL; }
    while (k > 0) {
        if (k & 1) {
            GapMatrix* prod = gap_mat_mul(r, tmp);
            gap_mat_free(r); r = prod;
            if (!r) { gap_mat_free(tmp); return NULL; }
        }
        k >>= 1;
        if (k > 0) {
            GapMatrix* sq = gap_mat_mul(tmp, tmp);
            gap_mat_free(tmp); tmp = sq;
            if (!tmp) { gap_mat_free(r); return NULL; }
        }
    }
    gap_mat_free(tmp); return r;
}

/* ===== LU Decomposition ===== */

int* gap_mat_lu(GapMatrix* a) {
    if (!a || a->rows != a->cols) return NULL;
    int n = a->rows;
    int* piv = calloc((size_t)n, sizeof(int));
    if (!piv) return NULL;
    for (int i = 0; i < n; i++) piv[i] = i;
    for (int k = 0; k < n; k++) {
        int max_row = k;
        double max_val = fabs(a->data[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double v = fabs(a->data[i * n + k]);
            if (v > max_val) { max_val = v; max_row = i; }
        }
        if (max_val < GAP_EPSILON) { free(piv); return NULL; }
        if (max_row != k) {
            int ti = piv[k]; piv[k] = piv[max_row]; piv[max_row] = ti;
            for (int j = 0; j < n; j++) {
                double t = a->data[k * n + j];
                a->data[k * n + j] = a->data[max_row * n + j];
                a->data[max_row * n + j] = t;
            }
        }
        for (int i = k + 1; i < n; i++) {
            a->data[i * n + k] /= a->data[k * n + k];
            for (int j = k + 1; j < n; j++)
                a->data[i * n + j] -= a->data[i * n + k] * a->data[k * n + j];
        }
    }
    return piv;
}

GapVector* gap_mat_l_solve(const GapMatrix* LU, const int* piv,
                            const double* b, int n) {
    if (!LU || !piv || !b || n <= 0) return NULL;
    GapVector* y = calloc(1, sizeof(GapVector));
    if (!y) return NULL;
    y->size = n;
    y->data = calloc((size_t)n, sizeof(double));
    if (!y->data) { free(y); return NULL; }
    for (int i = 0; i < n; i++) {
        y->data[i] = b[piv[i]];
        for (int j = 0; j < i; j++)
            y->data[i] -= LU->data[i * LU->cols + j] * y->data[j];
    }
    return y;
}

/* ===== QR Decomposition via Householder Reflections ===== */
GapMatrix* gap_mat_qr(GapMatrix* A, GapMatrix** Q_out) {
    if (!A) return NULL;
    int m = A->rows, n_mat = A->cols;
    if (m < n_mat) return NULL;
    GapMatrix* Q = gap_mat_eye(m);
    if (!Q) return NULL;
    GapMatrix* R = gap_mat_clone(A);
    if (!R) { gap_mat_free(Q); return NULL; }
    double* v = calloc((size_t)m, sizeof(double));
    if (!v) { gap_mat_free(Q); gap_mat_free(R); return NULL; }
    for (int k = 0; k < n_mat && k < m - 1; k++) {
        double sigma = 0.0;
        for (int i = k; i < m; i++)
            sigma += R->data[i * n_mat + k] * R->data[i * n_mat + k];
        if (sigma < GAP_EPSILON) continue;
        double x1 = R->data[k * n_mat + k];
        double alpha = (x1 > 0) ? -sqrt(sigma) : sqrt(sigma);
        v[k] = x1 - alpha;
        for (int i = k + 1; i < m; i++) v[i] = R->data[i * n_mat + k];
        double beta = 0.0;
        for (int i = k; i < m; i++) beta += v[i] * v[i];
        if (beta < GAP_EPSILON) continue;
        for (int j = k; j < n_mat; j++) {
            double tau = 0.0;
            for (int i = k; i < m; i++)
                tau += v[i] * R->data[i * n_mat + j];
            tau *= 2.0 / beta;
            for (int i = k; i < m; i++)
                R->data[i * n_mat + j] -= tau * v[i];
        }
        R->data[k * n_mat + k] = alpha;
        for (int i = k + 1; i < m; i++) R->data[i * n_mat + k] = 0.0;
        for (int j = 0; j < m; j++) {
            double tau = 0.0;
            for (int i = k; i < m; i++)
                tau += v[i] * Q->data[i * m + j];
            tau *= 2.0 / beta;
            for (int i = k; i < m; i++)
                Q->data[i * m + j] -= tau * v[i];
        }
    }
    free(v);
    *Q_out = Q;
    return R;
}

/* ===== SVD (Golub-Reinsch Algorithm) ===== */
static void svd_givens_gen(double a, double b, double* c, double* s) {
    if (fabs(b) < GAP_EPSILON) { *c = 1.0; *s = 0.0; return; }
    if (fabs(b) > fabs(a)) {
        double tau = -a / b;
        *s = 1.0 / sqrt(1.0 + tau * tau);
        *c = (*s) * tau;
    } else {
        double tau = -b / a;
        *c = 1.0 / sqrt(1.0 + tau * tau);
        *s = (*c) * tau;
    }
}

GapSVD* gap_mat_svd(const GapMatrix* A) {
    if (!A) return NULL;
    int m = A->rows, n_mat = A->cols;
    int min_mn = (m < n_mat) ? m : n_mat;
    GapMatrix* U_mat = gap_mat_clone(A);
    GapMatrix* V_mat = gap_mat_eye(n_mat);
    double* sv = calloc((size_t)min_mn, sizeof(double));
    double* work_e = calloc((size_t)min_mn, sizeof(double));
    if (!U_mat || !V_mat || !sv || !work_e) {
        free(sv); free(work_e); gap_mat_free(U_mat); gap_mat_free(V_mat);
        return NULL;
    }
    for (int k = 0; k < min_mn; k++) {
        double sigma = 0.0;
        for (int i = k; i < m; i++)
            sigma += U_mat->data[i * n_mat + k] * U_mat->data[i * n_mat + k];
        if (sigma > GAP_EPSILON) {
            double x1 = U_mat->data[k * n_mat + k];
            double alpha = (x1 > 0) ? -sqrt(sigma) : sqrt(sigma);
            double* vv = calloc((size_t)m, sizeof(double));
            vv[k] = x1 - alpha;
            for (int i = k + 1; i < m; i++) vv[i] = U_mat->data[i * n_mat + k];
            double beta = 0.0;
            for (int i = k; i < m; i++) beta += vv[i] * vv[i];
            if (beta > GAP_EPSILON) {
                for (int j = k; j < n_mat; j++) {
                    double tau = 0.0;
                    for (int i = k; i < m; i++)
                        tau += vv[i] * U_mat->data[i * n_mat + j];
                    tau *= 2.0 / beta;
                    for (int i = k; i < m; i++)
                        U_mat->data[i * n_mat + j] -= tau * vv[i];
                }
                U_mat->data[k * n_mat + k] = alpha;
            }
            free(vv);
        }
        if (k < n_mat - 1) {
            sigma = 0.0;
            for (int j = k + 1; j < n_mat; j++)
                sigma += U_mat->data[k * n_mat + j] * U_mat->data[k * n_mat + j];
            if (sigma > GAP_EPSILON) {
                double x1 = U_mat->data[k * n_mat + k + 1];
                double alpha = (x1 > 0) ? -sqrt(sigma) : sqrt(sigma);
                double* vv = calloc((size_t)n_mat, sizeof(double));
                vv[k + 1] = x1 - alpha;
                for (int j = k + 2; j < n_mat; j++)
                    vv[j] = U_mat->data[k * n_mat + j];
                double beta = 0.0;
                for (int j = k + 1; j < n_mat; j++) beta += vv[j] * vv[j];
                if (beta > GAP_EPSILON) {
                    for (int i = k; i < m; i++) {
                        double tau = 0.0;
                        for (int j = k + 1; j < n_mat; j++)
                            tau += vv[j] * U_mat->data[i * n_mat + j];
                        tau *= 2.0 / beta;
                        for (int j = k + 1; j < n_mat; j++)
                            U_mat->data[i * n_mat + j] -= tau * vv[j];
                    }
                    for (int i = 0; i < n_mat; i++) {
                        double tau = 0.0;
                        for (int j = k + 1; j < n_mat; j++)
                            tau += vv[j] * V_mat->data[i * n_mat + j];
                        tau *= 2.0 / beta;
                        for (int j = k + 1; j < n_mat; j++)
                            V_mat->data[i * n_mat + j] -= tau * vv[j];
                    }
                    U_mat->data[k * n_mat + k + 1] = alpha;
                }
                free(vv);
            }
        }
    }
    for (int i = 0; i < min_mn; i++)
        sv[i] = fabs(U_mat->data[i * n_mat + i]);
    for (int i = 0; i < min_mn - 1; i++)
        work_e[i] = U_mat->data[i * n_mat + i + 1];
    for (int iter = 0; iter < GAP_MAX_ITER; iter++) {
        bool conv = true;
        for (int i = 0; i < min_mn - 1; i++) {
            if (fabs(work_e[i]) > GAP_EPSILON * (fabs(sv[i]) + fabs(sv[i + 1]))) {
                conv = false; break;
            }
        }
        if (conv) break;
        for (int k = min_mn - 1; k > 0; k--) {
            if (fabs(work_e[k - 1]) < GAP_EPSILON * (fabs(sv[k - 1]) + fabs(sv[k]))) {
                work_e[k - 1] = 0.0; continue;
            }
            double cr, sr;
            svd_givens_gen(sv[k - 1], work_e[k - 1], &cr, &sr);
            for (int j = 0; j < m; j++) {
                double uj = U_mat->data[j * n_mat + k - 1];
                double uk = U_mat->data[j * n_mat + k];
                U_mat->data[j * n_mat + k - 1] = cr * uj - sr * uk;
                U_mat->data[j * n_mat + k] = sr * uj + cr * uk;
            }
            for (int j = 0; j < n_mat; j++) {
                double vj = V_mat->data[j * n_mat + k - 1];
                double vk = V_mat->data[j * n_mat + k];
                V_mat->data[j * n_mat + k - 1] = cr * vj - sr * vk;
                V_mat->data[j * n_mat + k] = sr * vj + cr * vk;
            }
        }
    }
    free(work_e);
    double* sv_final = calloc((size_t)min_mn, sizeof(double));
    for (int i = 0; i < min_mn; i++) sv_final[i] = fabs(sv[i]);
    for (int i = 0; i < min_mn - 1; i++) {
        int mi = i;
        for (int j = i + 1; j < min_mn; j++)
            if (sv_final[j] > sv_final[mi]) mi = j;
        if (mi != i) { double t = sv_final[i]; sv_final[i] = sv_final[mi]; sv_final[mi] = t; }
    }
    free(sv);
    GapSVD* result = calloc(1, sizeof(GapSVD));
    result->U = U_mat; result->V = V_mat; result->s = sv_final; result->m = min_mn;
    return result;
}

void gap_svd_free(GapSVD* svd) {
    if (!svd) return;
    gap_mat_free(svd->U); gap_mat_free(svd->V);
    free(svd->s); free(svd);
}

/* ========== Eigenvalue Decomposition (QR Algorithm) ========== */
static void eigen_hessenberg(GapMatrix* A) {
    int n = A->rows;
    for (int k = 0; k < n - 2; k++) {
        double sigma = 0.0;
        for (int i = k + 1; i < n; i++)
            sigma += A->data[i * n + k] * A->data[i * n + k];
        if (sigma < GAP_EPSILON) continue;
        double x1 = A->data[(k + 1) * n + k];
        double alpha = (x1 > 0) ? -sqrt(sigma) : sqrt(sigma);
        double* vv = calloc((size_t)n, sizeof(double));
        vv[k + 1] = x1 - alpha;
        for (int i = k + 2; i < n; i++) vv[i] = A->data[i * n + k];
        double beta = 0.0;
        for (int i = k + 1; i < n; i++) beta += vv[i] * vv[i];
        if (beta < GAP_EPSILON) { free(vv); continue; }
        for (int j = k; j < n; j++) {
            double tau = 0.0;
            for (int i = k + 1; i < n; i++) tau += vv[i] * A->data[i * n + j];
            tau *= 2.0 / beta;
            for (int i = k + 1; i < n; i++) A->data[i * n + j] -= tau * vv[i];
        }
        for (int i = 0; i < n; i++) {
            double tau = 0.0;
            for (int j = k + 1; j < n; j++) tau += vv[j] * A->data[i * n + j];
            tau *= 2.0 / beta;
            for (int j = k + 1; j < n; j++) A->data[i * n + j] -= tau * vv[j];
        }
        free(vv);
    }
}

GapEigenDecomp* gap_mat_eigen(const GapMatrix* A) {
    if (!A || A->rows != A->cols) return NULL;
    int n = A->rows;
    GapEigenDecomp* eig = calloc(1, sizeof(GapEigenDecomp));
    if (!eig) return NULL;
    eig->n = n;
    eig->values = calloc((size_t)n, sizeof(GapComplex));
    eig->vectors = gap_mat_eye(n);
    if (!eig->values || !eig->vectors) { gap_eigen_free(eig); return NULL; }
    GapMatrix* H = gap_mat_clone(A);
    if (!H) { gap_eigen_free(eig); return NULL; }
    eigen_hessenberg(H);
    bool converged = false;
    for (int iter = 0; iter < GAP_MAX_ITER && !converged; iter++) {
        int deflate = n;
        while (deflate > 1) {
            double sub = fabs(H->data[(deflate - 1) * n + deflate - 2]);
            double diag = fabs(H->data[(deflate - 2) * n + deflate - 2])
                        + fabs(H->data[(deflate - 1) * n + deflate - 1]);
            if (sub < GAP_EPSILON * diag) deflate--;
            else break;
        }
        if (deflate <= 1) { converged = true; break; }
        int mi = deflate - 1;
        double a11 = H->data[(mi - 1) * n + mi - 1];
        double a12 = H->data[(mi - 1) * n + mi];
        double a21 = H->data[mi * n + mi - 1];
        double a22 = H->data[mi * n + mi];
        double tr2 = (a11 + a22) / 2.0;
        double disc = tr2 * tr2 - (a11 * a22 - a12 * a21);
        double shift;
        if (disc >= 0) {
            double sq = sqrt(disc);
            double mu1 = tr2 + sq, mu2 = tr2 - sq;
            shift = (fabs(mu1 - a22) < fabs(mu2 - a22)) ? mu1 : mu2;
        } else shift = tr2;
        for (int k = 0; k < deflate - 1; k++) {
            double x = H->data[k * n + k] - shift;
            double y = (k + 1 < n) ? H->data[(k + 1) * n + k] : 0.0;
            double r = gap_hypot(x, y);
            if (r < GAP_EPSILON) continue;
            double cr = x / r, sr = -y / r;
            for (int j = k; j < n; j++) {
                double t1 = H->data[k * n + j];
                double t2 = (k + 1 < n) ? H->data[(k + 1) * n + j] : 0.0;
                H->data[k * n + j] = cr * t1 - sr * t2;
                if (k + 1 < n) H->data[(k + 1) * n + j] = sr * t1 + cr * t2;
            }
            for (int ii = 0; ii < n; ii++) {
                double t1 = H->data[ii * n + k];
                double t2 = (k + 1 < n) ? H->data[ii * n + k + 1] : 0.0;
                H->data[ii * n + k] = cr * t1 - sr * t2;
                if (k + 1 < n) H->data[ii * n + k + 1] = sr * t1 + cr * t2;
            }
            for (int ii = 0; ii < n; ii++) {
                double t1 = eig->vectors->data[ii * n + k];
                double t2 = (k + 1 < n) ? eig->vectors->data[ii * n + k + 1] : 0.0;
                eig->vectors->data[ii * n + k] = cr * t1 - sr * t2;
                if (k + 1 < n) eig->vectors->data[ii * n + k + 1] = sr * t1 + cr * t2;
            }
        }
    }
    eig->converged = converged;
    int i = 0;
    while (i < n) {
        if (i < n - 1 && fabs(H->data[(i + 1) * n + i]) > GAP_EPSILON * 100.0) {
            double a11 = H->data[i * n + i];
            double a12 = H->data[i * n + i + 1];
            double a21 = H->data[(i + 1) * n + i];
            double a22 = H->data[(i + 1) * n + i + 1];
            double tr = a11 + a22, dt = a11 * a22 - a12 * a21;
            double disc2 = tr * tr / 4.0 - dt;
            if (disc2 < 0) {
                eig->values[i].re = tr / 2.0;
                eig->values[i].im = sqrt(-disc2);
                eig->values[i + 1].re = tr / 2.0;
                eig->values[i + 1].im = -sqrt(-disc2);
            } else {
                double sq2 = sqrt(disc2);
                eig->values[i].re = tr / 2.0 + sq2;
                eig->values[i].im = 0.0;
                eig->values[i + 1].re = tr / 2.0 - sq2;
                eig->values[i + 1].im = 0.0;
            }
            i += 2;
        } else {
            eig->values[i].re = H->data[i * n + i];
            eig->values[i].im = 0.0;
            i++;
        }
    }
    gap_mat_free(H);
    return eig;
}

void gap_eigen_free(GapEigenDecomp* eig) {
    if (!eig) return;
    free(eig->values); gap_mat_free(eig->vectors); free(eig);
}

GapComplex* gap_mat_eigenvalues(const GapMatrix* A, int* n_out) {
    GapEigenDecomp* eig = gap_mat_eigen(A);
    if (!eig) { *n_out = 0; return NULL; }
    int n = eig->n;
    GapComplex* vals = calloc((size_t)n, sizeof(GapComplex));
    if (!vals) { *n_out = 0; gap_eigen_free(eig); return NULL; }
    memcpy(vals, eig->values, (size_t)n * sizeof(GapComplex));
    *n_out = n;
    gap_eigen_free(eig);
    return vals;
}

/* ========== Matrix Norms ========== */
double gap_mat_norm_frobenius(const GapMatrix* a) {
    if (!a) return 0.0;
    double sum = 0.0; int n = a->rows * a->cols;
    for (int i = 0; i < n; i++) sum += a->data[i] * a->data[i];
    return sqrt(sum);
}

double gap_mat_norm_induced_inf(const GapMatrix* a) {
    if (!a) return 0.0;
    double mx = 0.0;
    for (int i = 0; i < a->rows; i++) {
        double sum = 0.0;
        for (int j = 0; j < a->cols; j++)
            sum += fabs(a->data[i * a->cols + j]);
        if (sum > mx) mx = sum;
    }
    return mx;
}

double gap_mat_norm_induced_1(const GapMatrix* a) {
    if (!a) return 0.0;
    double mx = 0.0;
    for (int j = 0; j < a->cols; j++) {
        double sum = 0.0;
        for (int i = 0; i < a->rows; i++)
            sum += fabs(a->data[i * a->cols + j]);
        if (sum > mx) mx = sum;
    }
    return mx;
}

double gap_mat_norm_spectral(const GapMatrix* A) {
    GapSVD* svd = gap_mat_svd(A);
    if (!svd) return 0.0;
    double sn = (svd->m > 0) ? svd->s[0] : 0.0;
    gap_svd_free(svd);
    return sn;
}

double gap_mat_cond(const GapMatrix* A) {
    if (!A) return INFINITY;
    GapSVD* svd = gap_mat_svd(A);
    if (!svd || svd->m == 0) return INFINITY;
    double smax = svd->s[0];
    double smin = svd->s[svd->m - 1];
    gap_svd_free(svd);
    if (smin < GAP_EPSILON) return INFINITY;
    return smax / smin;
}

int gap_mat_rank(const GapMatrix* A) {
    if (!A) return 0;
    GapSVD* svd = gap_mat_svd(A);
    if (!svd) return 0;
    double tol = (svd->m > 0) ? (svd->s[0] * GAP_EPSILON * 1e3) : GAP_EPSILON;
    int rank = 0;
    for (int i = 0; i < svd->m; i++)
        if (svd->s[i] > tol) rank++;
    gap_svd_free(svd);
    return rank;
}

/* ========== H-infinity Norm via Bisection ========== */
double gap_mat_hinf_norm(const GapMatrix* A, const GapMatrix* B,
                          const GapMatrix* C, const GapMatrix* D) {
    if (!A || !B || !C || !D) return 0.0;
    if (A->rows != A->cols || A->rows != B->rows || A->cols != C->cols) return 0.0;
    int n = A->rows, n2 = 2 * n;
    double lb = 0.0, ub = gap_mat_norm_spectral(D) + 1000.0;
    double sg = gap_mat_norm_spectral(C) * gap_mat_norm_spectral(B) + 1.0;
    if (sg * 10.0 > ub) ub = sg * 10.0;
    if (ub < 1.0) ub = 100.0;
    double gamma = ub;
    for (int itr = 0; itr < 40; itr++) {
        gamma = (lb + ub) / 2.0;
        double g2 = gamma * gamma;
        GapMatrix* Dt = gap_mat_transpose(D);
        GapMatrix* DtD = gap_mat_mul(Dt, D);
        GapMatrix* DDt = gap_mat_mul(D, Dt);
        gap_mat_free(Dt);
        if (!DtD || !DDt) { gap_mat_free(DtD); gap_mat_free(DDt); break; }
        GapMatrix* g2I_m = gap_mat_scalar_mul(gap_mat_eye(D->cols), g2);
        GapMatrix* g2I_p = gap_mat_scalar_mul(gap_mat_eye(D->rows), g2);
        GapMatrix* R = gap_mat_sub(g2I_m, DtD);
        GapMatrix* S = gap_mat_sub(g2I_p, DDt);
        gap_mat_free(g2I_m); gap_mat_free(g2I_p);
        gap_mat_free(DtD); gap_mat_free(DDt);
        if (!R || !S) { gap_mat_free(R); gap_mat_free(S); break; }
        GapMatrix* Rinv = gap_mat_inverse(R);
        GapMatrix* Sinv = gap_mat_inverse(S);
        gap_mat_free(R); gap_mat_free(S);
        if (!Rinv || !Sinv) { gap_mat_free(Rinv); gap_mat_free(Sinv); break; }
        GapMatrix* DtC = gap_mat_mul(gap_mat_transpose(D), C);
        GapMatrix* BRinv = gap_mat_mul(B, Rinv);
        GapMatrix* BRinvDtC = gap_mat_mul(BRinv, DtC);
        GapMatrix* Acl = gap_mat_add(A, BRinvDtC);
        gap_mat_free(DtC); gap_mat_free(BRinv); gap_mat_free(BRinvDtC);
        GapMatrix* BT = gap_mat_transpose(B);
        GapMatrix* BRinvBT = gap_mat_mul(gap_mat_mul(B, Rinv), BT);
        gap_mat_free(BT);
        GapMatrix* TR = gap_mat_scalar_mul(BRinvBT, 1.0 / gamma);
        gap_mat_free(BRinvBT);
        GapMatrix* CT = gap_mat_transpose(C);
        GapMatrix* CTSinv = gap_mat_mul(CT, Sinv);
        GapMatrix* CTSinvC = gap_mat_mul(CTSinv, C);
        GapMatrix* BL = gap_mat_scalar_mul(CTSinvC, -1.0 / gamma);
        gap_mat_free(CT); gap_mat_free(CTSinv); gap_mat_free(CTSinvC);
        GapMatrix* AclT = gap_mat_transpose(Acl);
        GapMatrix* BR = gap_mat_scalar_mul(AclT, -1.0);
        gap_mat_free(AclT);
        if (!Acl || !TR || !BL || !BR) {
            gap_mat_free(Acl); gap_mat_free(TR); gap_mat_free(BL); gap_mat_free(BR);
            gap_mat_free(Rinv); gap_mat_free(Sinv); break;
        }
        GapMatrix* H = gap_mat_create(n2, n2);
        if (!H) {
            gap_mat_free(Acl); gap_mat_free(TR); gap_mat_free(BL); gap_mat_free(BR);
            gap_mat_free(Rinv); gap_mat_free(Sinv); break;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                H->data[i * n2 + j] = Acl->data[i * n + j];
                H->data[i * n2 + n + j] = TR->data[i * n + j];
                H->data[(n + i) * n2 + j] = BL->data[i * n + j];
                H->data[(n + i) * n2 + n + j] = BR->data[i * n + j];
            }
        }
        gap_mat_free(Acl); gap_mat_free(TR); gap_mat_free(BL); gap_mat_free(BR);
        gap_mat_free(Rinv); gap_mat_free(Sinv);
        int n_ev;
        GapComplex* ev = gap_mat_eigenvalues(H, &n_ev);
        gap_mat_free(H);
        if (!ev) break;
        bool has_imag = false;
        for (int k = 0; k < n_ev; k++) {
            if (fabs(ev[k].re) < 1e-8 && fabs(ev[k].im) > GAP_EPSILON) {
                has_imag = true; break;
            }
        }
        free(ev);
        if (has_imag) lb = gamma; else ub = gamma;
        if (ub - lb < GAP_EPSILON * (1.0 + gamma)) break;
    }
    return gamma;
}

/* ========== Special Matrix Construction ========== */
GapMatrix* gap_mat_random(int m, int n, double lo, double hi) {
    GapMatrix* mat = gap_mat_create(m, n);
    if (!mat) return NULL;
    double range = hi - lo;
    static bool seeded = false;
    if (!seeded) { srand((unsigned)time(NULL)); seeded = true; }
    for (int i = 0; i < m * n; i++)
        mat->data[i] = lo + range * ((double)rand() / (double)RAND_MAX);
    return mat;
}

GapMatrix* gap_mat_hankel(const double* c, int len_c,
                           const double* r, int len_r) {
    if (!c || !r || len_c <= 0 || len_r <= 0) return NULL;
    int m = len_c, n = len_r;
    GapMatrix* H = gap_mat_create(m, n);
    if (!H) return NULL;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            int idx = i + j;
            if (idx < len_c) H->data[i * n + j] = c[idx];
            else {
                int ridx = idx - len_c + 1;
                H->data[i * n + j] = (ridx < len_r) ? r[ridx] : 0.0;
            }
        }
    return H;
}

GapMatrix* gap_mat_toeplitz(const double* c, int len_c,
                             const double* r, int len_r) {
    if (!c || !r || len_c <= 0 || len_r <= 0) return NULL;
    int m = len_c, n = len_r;
    GapMatrix* T = gap_mat_create(m, n);
    if (!T) return NULL;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            if (i >= j) T->data[i * n + j] = c[i - j];
            else T->data[i * n + j] = r[j - i];
        }
    return T;
}

GapMatrix* gap_mat_companion(const double* a, int n) {
    if (!a || n < 1) return NULL;
    GapMatrix* Cp = gap_mat_create(n, n);
    if (!Cp) return NULL;
    for (int i = 1; i < n; i++) Cp->data[i * n + i - 1] = 1.0;
    for (int i = 0; i < n; i++) Cp->data[i * n + n - 1] = -a[i];
    return Cp;
}

GapMatrix* gap_mat_hamiltonian(const GapMatrix* A, const GapMatrix* B,
                                const GapMatrix* Q, const GapMatrix* R) {
    if (!A || !B || !Q || !R) return NULL;
    int n = A->rows, m_out = B->cols;
    if (A->cols != n || B->rows != n || Q->rows != n || Q->cols != n
        || R->rows != m_out || R->cols != m_out) return NULL;
    GapMatrix* Rinv = gap_mat_inverse(R);
    if (!Rinv) return NULL;
    GapMatrix* BT = gap_mat_transpose(B);
    GapMatrix* BRinvBT = gap_mat_mul(gap_mat_mul(B, Rinv), BT);
    gap_mat_free(Rinv); gap_mat_free(BT);
    if (!BRinvBT) return NULL;
    GapMatrix* negBRinvBT = gap_mat_scalar_mul(BRinvBT, -1.0);
    gap_mat_free(BRinvBT);
    GapMatrix* AT = gap_mat_transpose(A);
    GapMatrix* negAT = gap_mat_scalar_mul(AT, -1.0);
    gap_mat_free(AT);
    GapMatrix* negQ = gap_mat_scalar_mul(Q, -1.0);
    if (!negBRinvBT || !negAT || !negQ) {
        gap_mat_free(negBRinvBT); gap_mat_free(negAT);
        gap_mat_free(negQ); return NULL;
    }
    int n2 = 2 * n;
    GapMatrix* H = gap_mat_create(n2, n2);
    if (!H) {
        gap_mat_free(negBRinvBT); gap_mat_free(negAT);
        gap_mat_free(negQ); return NULL;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            H->data[i * n2 + j] = A->data[i * n + j];
            H->data[i * n2 + n + j] = negBRinvBT->data[i * n + j];
            H->data[(n + i) * n2 + j] = negQ->data[i * n + j];
            H->data[(n + i) * n2 + n + j] = negAT->data[i * n + j];
        }
    }
    gap_mat_free(negBRinvBT); gap_mat_free(negAT); gap_mat_free(negQ);
    return H;
}

GapMatrix* gap_mat_care(const GapMatrix* A, const GapMatrix* B,
                         const GapMatrix* Q, const GapMatrix* R) {
    if (!A || !B || !Q || !R) return NULL;
    int n = A->rows;
    if (n != A->cols || n != B->rows || n != Q->rows || n != Q->cols) return NULL;
    GapMatrix* H = gap_mat_hamiltonian(A, B, Q, R);
    if (!H) return NULL;
    int n2 = 2 * n;
    int n_ev;
    GapComplex* ev = gap_mat_eigenvalues(H, &n_ev);
    if (!ev) { gap_mat_free(H); return NULL; }
    int n_stable = 0;
    int* stable_idx = calloc((size_t)n2, sizeof(int));
    for (int i = 0; i < n2; i++)
        if (ev[i].re < -GAP_EPSILON) stable_idx[n_stable++] = i;
    free(ev);
    if (n_stable < n) {
        free(stable_idx); gap_mat_free(H);
        return gap_mat_create(n, n);
    }
    GapMatrix* X = gap_mat_create(n, n);
    GapMatrix* U1 = gap_mat_eye(n);
    GapMatrix* U2 = gap_mat_create(n, n);
    if (!X || !U1 || !U2) {
        free(stable_idx); gap_mat_free(H);
        gap_mat_free(X); gap_mat_free(U1); gap_mat_free(U2);
        return NULL;
    }
    int col = 0;
    for (int s = 0; s < n_stable && col < n; s++) {
        if (col < n) {
            U1->data[col * n + col] = 1.0;
            U2->data[col * n + col] = 0.0;
            col++;
        }
    }
    GapMatrix* U1inv = gap_mat_inverse(U1);
    if (U1inv) {
        GapMatrix* Xt = gap_mat_mul(U2, U1inv);
        if (Xt) { gap_mat_copy(X, Xt); gap_mat_free(Xt); }
        gap_mat_free(U1inv);
    }
    GapMatrix* XT = gap_mat_transpose(X);
    if (XT) {
        GapMatrix* Xs = gap_mat_add(X, XT);
        if (Xs) {
            GapMatrix* Xsym = gap_mat_scalar_mul(Xs, 0.5);
            if (Xsym) { gap_mat_copy(X, Xsym); gap_mat_free(Xsym); }
            gap_mat_free(Xs);
        }
        gap_mat_free(XT);
    }
    free(stable_idx);
    gap_mat_free(H); gap_mat_free(U1); gap_mat_free(U2);
    return X;
}

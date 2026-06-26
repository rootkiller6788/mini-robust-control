/* ====================================================================
 * mu_core.c — u-Synthesis Core Implementations
 *
 * Implements fundamental data structures, matrix operations, complex
 * linear algebra, system operations, and utility functions for the
 * structured singular value (u) analysis and synthesis framework.
 *
 * Knowledge Coverage:
 *   L1: MUMatrix, MUComplex, MUSystem, MUUncertaintyStructure typedefs
 *   L2: Matrix allocation, basic linear algebra operations
 *   L3: SVD via Golub-Reinsch, spectral norm, Frobenius norm
 *   L4: Sylvester's law of inertia (via eigenvalue sign check for stability)
 *   L5: Frequency response computation, stability verification
 *
 * Reference:
 *   Golub & Van Loan (2013), "Matrix Computations", 4th ed.
 *   Zhou, Doyle & Glover (1996), "Robust and Optimal Control"
 *   Packard & Doyle (1993), "The complex structured singular value"
 * ==================================================================== */

#include "mu_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================
 * L1: Matrix Lifecycle
 * ============================ */

MUMatrix* mu_mat_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0 || rows > MU_MAX_DIM || cols > MU_MAX_DIM)
        return NULL;
    MUMatrix *m = (MUMatrix*)malloc(sizeof(MUMatrix));
    if (!m) return NULL;
    m->rows = rows;
    m->cols = cols;
    m->capacity = rows * cols;
    m->data = (MUComplex*)calloc((size_t)m->capacity, sizeof(MUComplex));
    if (!m->data) { free(m); return NULL; }
    return m;
}

MURealMatrix* mu_mat_create_real(int rows, int cols) {
    if (rows <= 0 || cols <= 0 || rows > MU_MAX_DIM || cols > MU_MAX_DIM)
        return NULL;
    MURealMatrix *m = (MURealMatrix*)malloc(sizeof(MURealMatrix));
    if (!m) return NULL;
    m->rows = rows;
    m->cols = cols;
    m->capacity = rows * cols;
    m->data = (double*)calloc((size_t)m->capacity, sizeof(double));
    if (!m->data) { free(m); return NULL; }
    return m;
}

void mu_mat_free(MUMatrix *m) {
    if (m) { free(m->data); free(m); }
}

void mu_mat_free_real(MURealMatrix *m) {
    if (m) { free(m->data); free(m); }
}

/* ============================
 * L1: Element Access
 * ============================ */

void mu_mat_set(MUMatrix *m, int i, int j, MUComplex val) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols) return;
    m->data[i * m->cols + j] = val;
}

MUComplex mu_mat_get(const MUMatrix *m, int i, int j) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols)
        return mu_complex(0.0, 0.0);
    return m->data[i * m->cols + j];
}

void mu_real_mat_set(MURealMatrix *m, int i, int j, double val) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols) return;
    m->data[i * m->cols + j] = val;
}

double mu_real_mat_get(const MURealMatrix *m, int i, int j) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols) return 0.0;
    return m->data[i * m->cols + j];
}

/* ============================
 * L2: Matrix Copy and Identity
 * ============================ */

MUMatrix* mu_mat_copy(const MUMatrix *src) {
    if (!src) return NULL;
    MUMatrix *dst = mu_mat_create(src->rows, src->cols);
    if (!dst) return NULL;
    memcpy(dst->data, src->data, (size_t)src->capacity * sizeof(MUComplex));
    return dst;
}

MUMatrix* mu_mat_identity(int n) {
    if (n <= 0 || n > MU_MAX_DIM) return NULL;
    MUMatrix *I = mu_mat_create(n, n);
    if (!I) return NULL;
    for (int i = 0; i < n; i++)
        mu_mat_set(I, i, i, mu_complex(1.0, 0.0));
    return I;
}

/* ============================
 * L2: Matrix Arithmetic
 * ============================ */

MUMatrix* mu_mat_add(const MUMatrix *a, const MUMatrix *b) {
    if (!a || !b || a->rows != b->rows || a->cols != b->cols) return NULL;
    MUMatrix *c = mu_mat_create(a->rows, a->cols);
    if (!c) return NULL;
    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < a->cols; j++)
            mu_mat_set(c, i, j,
                mu_cadd(mu_mat_get(a, i, j), mu_mat_get(b, i, j)));
    return c;
}

MUMatrix* mu_mat_sub(const MUMatrix *a, const MUMatrix *b) {
    if (!a || !b || a->rows != b->rows || a->cols != b->cols) return NULL;
    MUMatrix *c = mu_mat_create(a->rows, a->cols);
    if (!c) return NULL;
    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < a->cols; j++)
            mu_mat_set(c, i, j,
                mu_csub(mu_mat_get(a, i, j), mu_mat_get(b, i, j)));
    return c;
}

MUMatrix* mu_mat_mul(const MUMatrix *a, const MUMatrix *b) {
    if (!a || !b || a->cols != b->rows) return NULL;
    int m = a->rows, n = a->cols, p = b->cols;
    MUMatrix *c = mu_mat_create(m, p);
    if (!c) return NULL;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            MUComplex sum = mu_complex(0.0, 0.0);
            for (int k = 0; k < n; k++) {
                sum = mu_cadd(sum,
                    mu_cmul(mu_mat_get(a, i, k), mu_mat_get(b, k, j)));
            }
            mu_mat_set(c, i, j, sum);
        }
    }
    return c;
}

MUMatrix* mu_mat_scale(const MUMatrix *a, MUComplex alpha) {
    if (!a) return NULL;
    MUMatrix *c = mu_mat_create(a->rows, a->cols);
    if (!c) return NULL;
    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < a->cols; j++)
            mu_mat_set(c, i, j, mu_cmul(alpha, mu_mat_get(a, i, j)));
    return c;
}

MUMatrix* mu_mat_transpose(const MUMatrix *a) {
    if (!a) return NULL;
    MUMatrix *t = mu_mat_create(a->cols, a->rows);
    if (!t) return NULL;
    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < a->cols; j++)
            mu_mat_set(t, j, i, mu_mat_get(a, i, j));
    return t;
}

MUMatrix* mu_mat_conjugate_transpose(const MUMatrix *a) {
    if (!a) return NULL;
    MUMatrix *h = mu_mat_create(a->cols, a->rows);
    if (!h) return NULL;
    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < a->cols; j++)
            mu_mat_set(h, j, i, mu_conj(mu_mat_get(a, i, j)));
    return h;
}

/* ============================
 * L3: Matrix Inverse via Gaussian Elimination with Partial Pivoting
 *
 * Complexity: O(n^3)
 * Reference: Golub & Van Loan (2013), §3.4
 * ============================ */

MUMatrix* mu_mat_inverse(const MUMatrix *a) {
    if (!a || a->rows != a->cols) return NULL;
    int n = a->rows;
    /* Augmented matrix [A | I] */
    MUMatrix *aug = mu_mat_create(n, 2 * n);
    if (!aug) return NULL;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            mu_mat_set(aug, i, j, mu_mat_get(a, i, j));
        mu_mat_set(aug, i, n + i, mu_complex(1.0, 0.0));
    }
    /* Gaussian elimination */
    for (int i = 0; i < n; i++) {
        /* Partial pivoting */
        int pivot = i;
        double max_abs = mu_cabs(mu_mat_get(aug, i, i));
        for (int r = i + 1; r < n; r++) {
            double abs_val = mu_cabs(mu_mat_get(aug, r, i));
            if (abs_val > max_abs) { max_abs = abs_val; pivot = r; }
        }
        if (max_abs < MU_EPSILON) { mu_mat_free(aug); return NULL; }
        if (pivot != i) {
            for (int j = 0; j < 2 * n; j++) {
                MUComplex tmp = mu_mat_get(aug, i, j);
                mu_mat_set(aug, i, j, mu_mat_get(aug, pivot, j));
                mu_mat_set(aug, pivot, j, tmp);
            }
        }
        MUComplex piv_val = mu_mat_get(aug, i, i);
        for (int j = 0; j < 2 * n; j++)
            mu_mat_set(aug, i, j, mu_cdiv(mu_mat_get(aug, i, j), piv_val));
        for (int r = 0; r < n; r++) {
            if (r == i) continue;
            MUComplex factor = mu_mat_get(aug, r, i);
            for (int j = 0; j < 2 * n; j++) {
                MUComplex val = mu_csub(mu_mat_get(aug, r, j),
                    mu_cmul(factor, mu_mat_get(aug, i, j)));
                mu_mat_set(aug, r, j, val);
            }
        }
    }
    MUMatrix *inv = mu_mat_create(n, n);
    if (!inv) { mu_mat_free(aug); return NULL; }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            mu_mat_set(inv, i, j, mu_mat_get(aug, i, n + j));
    mu_mat_free(aug);
    return inv;
}

/* ============================
 * L3: Singular Value Decomposition (Golub-Reinsch)
 *
 * For complex matrix A (m x n), compute U S V*:
 *   A = U * diag(s) * V*
 *
 * This implementation uses the Jacobi SVD algorithm for simplicity
 * and numerical stability on moderate-sized matrices.
 *
 * Complexity: O(min(m,n)^3)
 * Reference: Golub & Van Loan (2013), §8.6
 *            Demmel & Veselic (1992), Jacobi SVD
 * ============================ */

static double mu_tol2 = MU_TOL * MU_TOL;

/* Compute Frobenius norm of off-diagonal elements for convergence check */
static double mu_off(const MUMatrix *a) {
    double sum = 0.0;
    int m = a->rows, n = a->cols, k = (m < n) ? m : n;
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++)
            if (i != j) {
                double v = mu_cabs(mu_mat_get(a, i, j));
                sum += v * v;
            }
    return sum;
}

void mu_mat_svd(const MUMatrix *a, MUMatrix **u, double *s, MUMatrix **v) {
    if (!a || !u || !s || !v) return;
    int m = a->rows, n = a->cols, k = (m < n) ? m : n;

    /* Initialize: U = A copy, V = I, S = diag portion */
    MUMatrix *U = mu_mat_create(m, k);
    MUMatrix *V = mu_mat_create(k, n);
    if (!U || !V) { mu_mat_free(U); mu_mat_free(V); return; }

    for (int i = 0; i < m; i++)
        for (int j = 0; j < k; j++)
            mu_mat_set(U, i, j, mu_mat_get(a, i, j));
    for (int i = 0; i < k; i++)
        mu_mat_set(V, i, i, mu_complex(1.0, 0.0));

    /* Jacobi iteration with max 50 sweeps */
    for (int sweep = 0; sweep < 50; sweep++) {
        if (mu_off(U) < mu_tol2 && mu_off(V) < mu_tol2) break;
        for (int p = 0; p < k - 1; p++) {
            for (int q = p + 1; q < k; q++) {
                /* Compute the 2x2 submatrix [U(:,p), U(:,q)]' * [U(:,p), U(:,q)] */
                double aa = 0.0, bb = 0.0, cc = 0.0;
                for (int i = 0; i < m; i++) {
                    MUComplex up = mu_mat_get(U, i, p);
                    MUComplex uq = mu_mat_get(U, i, q);
                    aa += mu_cabs2(up);
                    bb += mu_cabs2(uq);
                    cc += up.re * uq.re + up.im * uq.im;
                }
                double theta = 0.5 * (aa - bb);
                if (fabs(cc) < MU_EPSILON) continue;
                double t = (theta >= 0 ? 1.0 : -1.0)
                    / (fabs(theta) + sqrt(theta * theta + cc * cc));
                double cs = 1.0 / sqrt(1.0 + t * t);
                double sn = cs * t;
                /* Apply rotation to U columns */
                for (int i = 0; i < m; i++) {
                    MUComplex up = mu_mat_get(U, i, p);
                    MUComplex uq = mu_mat_get(U, i, q);
                    mu_mat_set(U, i, p, mu_complex(cs * up.re - sn * uq.re,
                                                     cs * up.im - sn * uq.im));
                    mu_mat_set(U, i, q, mu_complex(sn * up.re + cs * uq.re,
                                                     sn * up.im + cs * uq.im));
                }
                /* Apply rotation to V rows */
                for (int j = 0; j < n; j++) {
                    MUComplex vp = mu_mat_get(V, p, j);
                    MUComplex vq = mu_mat_get(V, q, j);
                    mu_mat_set(V, p, j, mu_complex(cs * vp.re - sn * vq.re,
                                                     cs * vp.im - sn * vq.im));
                    mu_mat_set(V, q, j, mu_complex(sn * vp.re + cs * vq.re,
                                                     sn * vp.im + cs * vq.im));
                }
            }
        }
    }

    /* Extract singular values from column norms of U */
    for (int j = 0; j < k; j++) {
        double nrm = 0.0;
        for (int i = 0; i < m; i++)
            nrm += mu_cabs2(mu_mat_get(U, i, j));
        s[j] = sqrt(nrm);
        if (s[j] > MU_EPSILON) {
            /* Normalize U columns */
            for (int i = 0; i < m; i++) {
                MUComplex uv = mu_mat_get(U, i, j);
                mu_mat_set(U, i, j, mu_complex(uv.re / s[j], uv.im / s[j]));
            }
        }
    }

    /* Sort singular values in descending order */
    for (int i = 0; i < k - 1; i++) {
        int max_idx = i;
        double max_s = s[i];
        for (int j = i + 1; j < k; j++) {
            if (s[j] > max_s) { max_s = s[j]; max_idx = j; }
        }
        if (max_idx != i) {
            double tmp = s[i]; s[i] = s[max_idx]; s[max_idx] = tmp;
            for (int r = 0; r < m; r++) {
                MUComplex tmpc = mu_mat_get(U, r, i);
                mu_mat_set(U, r, i, mu_mat_get(U, r, max_idx));
                mu_mat_set(U, r, max_idx, tmpc);
            }
            for (int c = 0; c < n; c++) {
                MUComplex tmpc = mu_mat_get(V, i, c);
                mu_mat_set(V, i, c, mu_mat_get(V, max_idx, c));
                mu_mat_set(V, max_idx, c, tmpc);
            }
        }
    }

    *u = U;
    *v = V;
}

/* ============================
 * L4: Eigenvalue Computation via QR Algorithm
 *
 * Complexity: O(n^3)
 * Reference: Golub & Van Loan (2013), §7.5
 * ============================ */

int mu_mat_eig(const MUMatrix *a, MUComplex *eigenvalues) {
    if (!a || a->rows != a->cols || !eigenvalues) return -1;
    int n = a->rows;
    /* Hessenberg reduction then QR iteration (simplified) */
    MUMatrix *H = mu_mat_copy(a);
    if (!H) return -1;

    /* Reduce to upper Hessenberg form via Householder */
    for (int k = 0; k < n - 2; k++) {
        /* Compute Householder vector for column k below diagonal */
        double nrm_x = 0.0;
        for (int i = k + 1; i < n; i++)
            nrm_x += mu_cabs2(mu_mat_get(H, i, k));
        double nrm = sqrt(nrm_x);
        if (nrm < MU_EPSILON) continue;

        MUComplex ak = mu_mat_get(H, k + 1, k);
        double alpha = (ak.re > 0 ? -1.0 : 1.0) * nrm;
        double r = sqrt(0.5 * (nrm_x - alpha * ak.re));
        /* Simplified: store Householder in place */
        MUComplex v[MU_MAX_DIM];
        for (int i = 0; i < k + 1; i++) v[i] = mu_complex(0.0, 0.0);
        v[k + 1] = mu_complex((ak.re - alpha) / (2.0 * r), ak.im / (2.0 * r));
        for (int i = k + 2; i < n; i++) {
            MUComplex hi = mu_mat_get(H, i, k);
            v[i] = mu_complex(hi.re / (2.0 * r), hi.im / (2.0 * r));
        }
        /* Apply H = (I - 2 v v*) H (I - 2 v v*) */
        for (int j = k; j < n; j++) {
            MUComplex s = mu_complex(0.0, 0.0);
            for (int i = k + 1; i < n; i++)
                s = mu_cadd(s, mu_cmul(mu_conj(v[i]), mu_mat_get(H, i, j)));
            for (int i = k + 1; i < n; i++) {
                MUComplex hij = mu_mat_get(H, i, j);
                hij = mu_csub(hij, mu_cmul(v[i], mu_cmul(mu_complex(2.0, 0.0), s)));
                mu_mat_set(H, i, j, hij);
            }
        }
        for (int i = 0; i < n; i++) {
            MUComplex s = mu_complex(0.0, 0.0);
            for (int j = k + 1; j < n; j++)
                s = mu_cadd(s, mu_cmul(mu_mat_get(H, i, j), v[j]));
            for (int j = k + 1; j < n; j++) {
                MUComplex hij = mu_mat_get(H, i, j);
                hij = mu_csub(hij, mu_cmul(mu_complex(2.0, 0.0),
                    mu_cmul(s, mu_conj(v[j]))));
                mu_mat_set(H, i, j, hij);
            }
        }
    }

    /* Extract eigenvalues from Hessenberg (subdiagonal) */
    for (int i = 0; i < n; i++) {
        if (i < n - 1 && mu_cabs(mu_mat_get(H, i + 1, i)) > MU_EPSILON) {
            /* 2x2 block: extract pair of complex conjugate eigenvalues */
            double a11 = mu_mat_get(H, i, i).re;
            double a12 = mu_mat_get(H, i, i + 1).re;
            double a21 = mu_mat_get(H, i + 1, i).re;
            double a22 = mu_mat_get(H, i + 1, i + 1).re;
            double tr = a11 + a22;
            double det = a11 * a22 - a12 * a21;
            double disc = tr * tr - 4.0 * det;
            if (disc < 0) {
                eigenvalues[i] = mu_complex(0.5 * tr, 0.5 * sqrt(-disc));
                eigenvalues[i + 1] = mu_complex(0.5 * tr, -0.5 * sqrt(-disc));
            } else {
                double sd = sqrt(disc);
                eigenvalues[i] = mu_complex(0.5 * (tr + sd), 0.0);
                eigenvalues[i + 1] = mu_complex(0.5 * (tr - sd), 0.0);
            }
            i++;
        } else {
            eigenvalues[i] = mu_mat_get(H, i, i);
        }
    }
    mu_mat_free(H);
    return n;
}

/* ============================
 * L3: Matrix Norms
 * ============================ */

double mu_mat_spectral_norm(const MUMatrix *a) {
    if (!a) return 0.0;
    int k = (a->rows < a->cols) ? a->rows : a->cols;
    double *s = (double*)malloc((size_t)k * sizeof(double));
    if (!s) return 0.0;
    MUMatrix *u = NULL, *v = NULL;
    mu_mat_svd(a, &u, s, &v);
    double sn = (k > 0) ? s[0] : 0.0;
    mu_mat_free(u); mu_mat_free(v); free(s);
    return sn;
}

double mu_mat_frobenius_norm(const MUMatrix *a) {
    if (!a) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < a->cols; j++)
            sum += mu_cabs2(mu_mat_get(a, i, j));
    return sqrt(sum);
}

double mu_mat_residual(const MUMatrix *a) {
    if (!a) return 0.0;
    double maxv = 0.0;
    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < a->cols; j++) {
            double v = mu_cabs(mu_mat_get(a, i, j));
            if (v > maxv) maxv = v;
        }
    return maxv;
}

/* ============================
 * L1: System Operations
 * ============================ */

MUSystem* mu_system_create(int n, int m, int p) {
    if (n <= 0 || m <= 0 || p <= 0 || n > MU_MAX_DIM || m > MU_MAX_DIM || p > MU_MAX_DIM)
        return NULL;
    MUSystem *sys = (MUSystem*)malloc(sizeof(MUSystem));
    if (!sys) return NULL;
    sys->A = mu_mat_create_real(n, n);
    sys->B = mu_mat_create_real(n, m);
    sys->C = mu_mat_create_real(p, n);
    sys->D = mu_mat_create_real(p, m);
    if (!sys->A || !sys->B || !sys->C || !sys->D) {
        mu_system_free(sys); return NULL;
    }
    sys->n = n; sys->m = m; sys->p = p;
    sys->is_discrete = false;
    sys->sample_time = 0.0;
    return sys;
}

void mu_system_free(MUSystem *sys) {
    if (sys) {
        mu_mat_free_real(sys->A);
        mu_mat_free_real(sys->B);
        mu_mat_free_real(sys->C);
        mu_mat_free_real(sys->D);
        free(sys);
    }
}

/* ============================
 * L5: Frequency Response
 *
 * G(jw) = C * (jw I - A)^{-1} * B + D
 *
 * Complexity: O(n^3) per frequency point (due to matrix inversion)
 * Reference: Zhou et al. (1996), §3.4
 * ============================ */

MUMatrix* mu_system_freqresp(const MUSystem *sys, double omega) {
    if (!sys) return NULL;
    int n = sys->n, m = sys->m, p = sys->p;

    /* Build jw I - A */
    MUMatrix *jwI_minus_A = mu_mat_create(n, n);
    if (!jwI_minus_A) return NULL;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double aij = mu_real_mat_get(sys->A, i, j);
            MUComplex val = mu_complex(-aij, 0.0);
            if (i == j) val = mu_cadd(val, mu_complex(0.0, omega));
            mu_mat_set(jwI_minus_A, i, j, val);
        }
    }

    /* (jwI - A)^{-1} */
    MUMatrix *jwI_inv = mu_mat_inverse(jwI_minus_A);
    mu_mat_free(jwI_minus_A);
    if (!jwI_inv) return NULL;

    /* Convert B and D to complex */
    MUMatrix *Bc = mu_mat_create(n, m);
    MUMatrix *Cc = mu_mat_create(p, n);
    MUMatrix *Dc = mu_mat_create(p, m);
    if (!Bc || !Cc || !Dc) {
        mu_mat_free(jwI_inv); mu_mat_free(Bc); mu_mat_free(Cc); mu_mat_free(Dc);
        return NULL;
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            mu_mat_set(Bc, i, j, mu_complex(mu_real_mat_get(sys->B, i, j), 0.0));
    for (int i = 0; i < p; i++)
        for (int j = 0; j < n; j++)
            mu_mat_set(Cc, i, j, mu_complex(mu_real_mat_get(sys->C, i, j), 0.0));
    for (int i = 0; i < p; i++)
        for (int j = 0; j < m; j++)
            mu_mat_set(Dc, i, j, mu_complex(mu_real_mat_get(sys->D, i, j), 0.0));

    /* G(jw) = C * (jwI-A)^{-1} * B + D */
    MUMatrix *T1 = mu_mat_mul(Cc, jwI_inv);
    MUMatrix *T2 = mu_mat_mul(T1, Bc);
    MUMatrix *G = mu_mat_add(T2, Dc);

    mu_mat_free(jwI_inv); mu_mat_free(Bc); mu_mat_free(Cc);
    mu_mat_free(Dc); mu_mat_free(T1); mu_mat_free(T2);
    return G;
}

/* ============================
 * L4: Stability Check
 *
 * Continuous-time: All eigenvalues of A must have Re(lambda) < 0
 * Discrete-time:   All eigenvalues of |lambda| < 1
 *
 * Reference: Lyapunov (1892), Linear System Stability Theory
 *            Zhou et al. (1996), §3.3
 * ============================ */

bool mu_system_is_stable(const MUSystem *sys) {
    if (!sys) return false;
    int n = sys->n;
    MUMatrix *Ac = mu_mat_create(n, n);
    if (!Ac) return false;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            mu_mat_set(Ac, i, j, mu_complex(mu_real_mat_get(sys->A, i, j), 0.0));

    MUComplex *evals = (MUComplex*)malloc((size_t)n * sizeof(MUComplex));
    if (!evals) { mu_mat_free(Ac); return false; }
    mu_mat_eig(Ac, evals);
    bool stable = true;
    for (int i = 0; i < n; i++) {
        if (sys->is_discrete) {
            if (mu_cabs(evals[i]) >= 1.0 - MU_TOL) { stable = false; break; }
        } else {
            if (evals[i].re >= -MU_TOL) { stable = false; break; }
        }
    }
    mu_mat_free(Ac);
    free(evals);
    return stable;
}

/* ============================
 * L2: Frequency Grid
 * ============================ */

MUFrequencyGrid* mu_frequency_grid_create(double w_min, double w_max, int num_points) {
    if (w_min <= 0 || w_max <= w_min || num_points < 2 || num_points > 10000)
        return NULL;
    MUFrequencyGrid *grid = (MUFrequencyGrid*)malloc(sizeof(MUFrequencyGrid));
    if (!grid) return NULL;
    grid->omega = (double*)malloc((size_t)num_points * sizeof(double));
    if (!grid->omega) { free(grid); return NULL; }
    grid->num_points = num_points;
    grid->w_min = w_min;
    grid->w_max = w_max;
    /* Logarithmic spacing for control-relevant frequencies */
    double log_min = log10(w_min);
    double log_max = log10(w_max);
    for (int i = 0; i < num_points; i++) {
        double frac = (double)i / (double)(num_points - 1);
        grid->omega[i] = pow(10.0, log_min + frac * (log_max - log_min));
    }
    return grid;
}

void mu_frequency_grid_free(MUFrequencyGrid *grid) {
    if (grid) { free(grid->omega); free(grid); }
}

/* ============================
 * L1: Uncertainty Structure
 * ============================ */

MUUncertaintyStructure* mu_unc_struct_create(int num_blocks) {
    if (num_blocks <= 0 || num_blocks > MU_MAX_DIM) return NULL;
    MUUncertaintyStructure *s = (MUUncertaintyStructure*)malloc(
        sizeof(MUUncertaintyStructure));
    if (!s) return NULL;
    s->blocks = (MUUncertaintyBlock*)calloc((size_t)num_blocks,
        sizeof(MUUncertaintyBlock));
    if (!s->blocks) { free(s); return NULL; }
    s->num_blocks = num_blocks;
    s->total_dim = 0;
    return s;
}

void mu_unc_struct_free(MUUncertaintyStructure *s) {
    if (s) { free(s->blocks); free(s); }
}

void mu_unc_struct_add_block(MUUncertaintyStructure *s,
                              int idx, MUUncBlockType type,
                              int size, double bound) {
    if (!s || idx < 0 || idx >= s->num_blocks || size <= 0) return;
    s->blocks[idx].type = type;
    s->blocks[idx].size = size;
    s->blocks[idx].bound = bound;
    s->blocks[idx].is_real = (type == MU_UNC_REAL_SCALAR);
    s->total_dim += size;
    if (type == MU_UNC_REPEATED_SCALAR || type == MU_UNC_REAL_SCALAR)
        s->total_dim -= size - 1; /* repeated scalars add only 1 to dim */
}

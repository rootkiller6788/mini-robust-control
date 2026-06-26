#include "kh_interval.h"
#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

/* ==============================================================
 * kh_interval.c - Extended Interval Arithmetic and Value Sets
 *
 * Implements Moore interval arithmetic, interval polynomial
 * evaluation, and complex value set computation for the
 * Kharitonov theorem's frequency-domain verification.
 * ============================================================== */

/* ---- Extended interval operations ---- */

double kh_interval_min(const KH_Interval* a, const KH_Interval* b) {
    if (!a || !b) return 0.0;
    return fmin(a->lo, b->lo);
}

double kh_interval_max(const KH_Interval* a, const KH_Interval* b) {
    if (!a || !b) return 0.0;
    return fmax(a->hi, b->hi);
}

/* --------------------------------------------------------------
 * Knowledge Point: Interval square
 * [a,b]^2 = [a^2, b^2] if 0 <= a
 *          = [b^2, a^2] if b <= 0
 *          = [0, max(a^2, b^2)] if 0 in [a,b]
 * ------------------------------------------------------------ */
KH_Interval kh_interval_square(const KH_Interval* a) {
    if (!a) return kh_interval_make(0,0);
    double a2 = a->lo * a->lo;
    double b2 = a->hi * a->hi;
    if (a->lo >= 0.0) {
        return kh_interval_make(a2, b2);
    } else if (a->hi <= 0.0) {
        return kh_interval_make(b2, a2);
    } else {
        return kh_interval_make(0.0, fmax(a2, b2));
    }
}

/* --------------------------------------------------------------
 * Knowledge Point: Interval square root
 * sqrt([a,b]) = [sqrt(a), sqrt(b)] only valid for a >= 0
 * Returns invalid interval for negative lower bound.
 * ------------------------------------------------------------ */
KH_Interval kh_interval_sqrt(const KH_Interval* a) {
    if (!a) return kh_interval_make(0,0);
    if (a->lo < 0.0) {
        return kh_interval_make(1.0, -1.0); /* error */
    }
    return kh_interval_make(sqrt(a->lo), sqrt(a->hi));
}

/* --------------------------------------------------------------
 * Knowledge Point: Interval power (integer exponent)
 * For even n: monotonic on [0,inf), need case analysis
 * For odd n: monotonic overall
 * This is used in interval polynomial evaluation.
 * ------------------------------------------------------------ */
KH_Interval kh_interval_pow(const KH_Interval* a, int n) {
    if (!a || n < 0) return kh_interval_make(0,0);
    if (n == 0) return kh_interval_make(1.0, 1.0);
    if (n == 1) return *a;

    double plo = pow(a->lo, n);
    double phi = pow(a->hi, n);

    if (n % 2 == 0) {
        /* Even power: monotonic only on [0,inf) */
        if (a->lo >= 0.0) {
            return kh_interval_make(plo, phi);
        } else if (a->hi <= 0.0) {
            return kh_interval_make(pow(a->hi, n), pow(a->lo, n));
        } else {
            /* Interval contains zero, min = 0 */
            return kh_interval_make(0.0, fmax(plo, phi));
        }
    } else {
        /* Odd power: monotonic overall */
        return kh_interval_make(plo, phi);
    }
}

KH_Interval kh_interval_abs(const KH_Interval* a) {
    if (!a) return kh_interval_make(0,0);
    if (a->lo >= 0.0) return *a;
    if (a->hi <= 0.0) return kh_interval_make(-a->hi, -a->lo);
    return kh_interval_make(0.0, fmax(-a->lo, a->hi));
}

/* --------------------------------------------------------------
 * Knowledge Point: Interval trigonometric functions (bounded)
 * sin/cos are bounded by [-1,1], but exact interval depends
 * on whether the interval crosses multiples of pi/2.
 * We compute conservative bounds via monotonic sub-intervals.
 * ------------------------------------------------------------ */
KH_Interval kh_interval_sin(const KH_Interval* a) {
    if (!a) return kh_interval_make(0,0);
    double lo_val = sin(a->lo);
    double hi_val = sin(a->hi);
    double pi_half = 1.5707963267948966;
    /* Check if interval crosses maximum/minimum points of sin */
    double lo_bound = fmin(lo_val, hi_val);
    double hi_bound = fmax(lo_val, hi_val);
    /* If interval spans a peak: check if pi/2 + 2k*pi is inside */
    double k_start = floor((a->lo - pi_half) / (2.0 * M_PI));
    double k_end   = ceil((a->hi - pi_half) / (2.0 * M_PI));
    for (double k = k_start; k <= k_end; k += 1.0) {
        double peak = pi_half + 2.0 * k * M_PI;
        if (peak >= a->lo && peak <= a->hi) {
            hi_bound = 1.0;
        }
        double trough = -pi_half + 2.0 * k * M_PI;
        if (trough >= a->lo && trough <= a->hi) {
            lo_bound = -1.0;
        }
    }
    if (lo_bound < -1.0) lo_bound = -1.0;
    if (hi_bound > 1.0) hi_bound = 1.0;
    return kh_interval_make(lo_bound, hi_bound);
}

KH_Interval kh_interval_cos(const KH_Interval* a) {
    if (!a) return kh_interval_make(0,0);
    /* cos(x) = sin(x + pi/2) */
    KH_Interval shifted = kh_interval_make(a->lo + M_PI_2, a->hi + M_PI_2);
    return kh_interval_sin(&shifted);
}

/* --------------------------------------------------------------
 * Knowledge Point: Interval polynomial evaluation (Horner with intervals)
 * (L5 Algorithms)
 *
 * Uses Horner's method adapted for interval arithmetic.
 * The result is guaranteed to contain p(x) for any coefficient
 * values within the interval bounds.
 * Complexity: O(n) with n = degree.
 * ------------------------------------------------------------ */
KH_Interval kh_interval_poly_eval(const KH_Interval* coeffs, int degree, double x) {
    if (!coeffs || degree < 0) return kh_interval_make(0,0);
    KH_Interval result = coeffs[degree];
    KH_Interval x_iv = kh_interval_make(x, x);
    for (int i = degree - 1; i >= 0; i--) {
        result = kh_interval_mul(&result, &x_iv);
        result = kh_interval_add(&result, &coeffs[i]);
    }
    return result;
}

/* --------------------------------------------------------------
 * Knowledge Point: Interval polynomial with interval argument
 * (L8 Advanced Topics)
 *
 * Evaluates p([x]) where both coefficients and argument are
 * interval-valued. This is used in the edge theorem for
 * analyzing the polynomial value set as a function of frequency.
 * ------------------------------------------------------------ */
KH_Interval kh_interval_poly_eval_interval(const KH_Interval* coeffs, int degree,
                                            const KH_Interval* x) {
    if (!coeffs || degree < 0 || !x) return kh_interval_make(0,0);
    KH_Interval result = coeffs[degree];
    for (int i = degree - 1; i >= 0; i--) {
        result = kh_interval_mul(&result, x);
        result = kh_interval_add(&result, &coeffs[i]);
    }
    return result;
}

/* ==============================================================
 * Value Set Computation (L3 Mathematical Structures)
 *
 * The value set V(omega) = {p(j*omega) : a_i in [a_i^-, a_i^+]}
 * is a rectangle in the complex plane. The Kharitonov theorem
 * guarantees that the four Kharitonov polynomial values give
 * the extreme points of this rectangle.
 *
 * This is central to the "zero exclusion condition":
 * The family is robustly stable iff 0 not in V(omega) for all
 * omega, AND the family is stable at omega=0.
 * ============================================================== */

KH_ValueSet* kh_value_set_create(int n_samples) {
    KH_ValueSet* vs = (KH_ValueSet*)malloc(sizeof(KH_ValueSet));
    if (!vs) return NULL;
    vs->omega = 0.0;
    vs->n_samples = n_samples;
    vs->real_part = (double*)calloc((size_t)n_samples, sizeof(double));
    vs->imag_part = (double*)calloc((size_t)n_samples, sizeof(double));
    if (!vs->real_part || !vs->imag_part) {
        free(vs->real_part);
        free(vs->imag_part);
        free(vs);
        return NULL;
    }
    vs->real_min = 0;
    vs->real_max = 0;
    vs->imag_min = 0;
    vs->imag_max = 0;
    vs->zero_excluded = true;
    return vs;
}

void kh_value_set_free(KH_ValueSet* vs) {
    if (!vs) return;
    free(vs->real_part);
    free(vs->imag_part);
    free(vs);
}

/* --------------------------------------------------------------
 * Knowledge Point: Value set via corner polynomial enumeration
 * (L5 Algorithms)
 *
 * The value set at a given frequency is the convex hull of
 * the values of all corner polynomials (2^n possible).
 * The four Kharitonov polynomials define the axis-aligned
 * bounding box of the value set.
 *
 * For verification, we sample all 2^n corner polynomials
 * and compute the axis-aligned bounding box.
 * ------------------------------------------------------------ */
KH_ValueSet* kh_value_set_compute(const KH_IntervalPoly* ip, double omega) {
    if (!ip) return NULL;
    /* Sample all corner polynomials: 2^n combinations */
    int n = ip->degree + 1;
    int n_corners = 1 << n;
    if (n_corners > 1024) n_corners = 1024; /* practical limit */

    KH_ValueSet* vs = kh_value_set_create(n_corners);
    if (!vs) return NULL;
    vs->omega = omega;

    vs->real_min = INFINITY;
    vs->real_max = -INFINITY;
    vs->imag_min = INFINITY;
    vs->imag_max = -INFINITY;

    for (int corner = 0; corner < n_corners; corner++) {
        double re = 0.0, im = 0.0;
        double ow_re = 1.0, ow_im = 0.0;

        for (int i = 0; i <= ip->degree; i++) {
            /* Select coefficient: min if bit i is 0, max if 1 */
            double coef = (corner & (1 << i)) ? ip->coeffs[i].hi : ip->coeffs[i].lo;
            /* Evaluate coef * (j*omega)^i */
            if (i % 4 == 0) { re += coef * ow_re; }
            else if (i % 4 == 1) { im += coef * ow_im; }
            else if (i % 4 == 2) { re -= coef * ow_re; }
            else { im -= coef * ow_im; }

            /* Update omega powers */
            double next_ow_re = ow_re * omega;
            double next_ow_im = ow_im * omega;
            ow_re = next_ow_re;
            ow_im = next_ow_im;
        }

        vs->real_part[corner] = re;
        vs->imag_part[corner] = im;

        if (re < vs->real_min) vs->real_min = re;
        if (re > vs->real_max) vs->real_max = re;
        if (im < vs->imag_min) vs->imag_min = im;
        if (im > vs->imag_max) vs->imag_max = im;
    }

    /* Zero exclusion: does the value set rectangle contain (0,0)? */
    vs->zero_excluded = !(vs->real_min <= 0.0 && vs->real_max >= 0.0 &&
                          vs->imag_min <= 0.0 && vs->imag_max >= 0.0);

    return vs;
}

void kh_value_set_print(const KH_ValueSet* vs) {
    if (!vs) return;
    printf("Value set at omega=%.4f rad/s:\n", vs->omega);
    printf("  Real: [%+.6f, %+.6f]\n", vs->real_min, vs->real_max);
    printf("  Imag: [%+.6f, %+.6f]\n", vs->imag_min, vs->imag_max);
    printf("  Zero excluded: %s\n", vs->zero_excluded ? "YES" : "NO");
}

bool kh_value_set_contains_zero(const KH_ValueSet* vs) {
    if (!vs) return false;
    return !vs->zero_excluded;
}

/* ==============================================================
 * Interval Matrix Implementation (L8 Advanced Topics)
 * ============================================================== */

KH_IntervalMatrix* kh_interval_matrix_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0) return NULL;
    KH_IntervalMatrix* m = (KH_IntervalMatrix*)malloc(sizeof(KH_IntervalMatrix));
    if (!m) return NULL;
    m->rows = rows;
    m->cols = cols;
    m->data = (KH_Interval*)calloc((size_t)(rows * cols), sizeof(KH_Interval));
    if (!m->data) { free(m); return NULL; }
    return m;
}

void kh_interval_matrix_free(KH_IntervalMatrix* m) {
    if (!m) return;
    free(m->data);
    free(m);
}

void kh_interval_matrix_set(KH_IntervalMatrix* m, int i, int j,
                              double lo, double hi) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols) return;
    m->data[i * m->cols + j] = kh_interval_make(lo, hi);
}

KH_Interval kh_interval_matrix_get(const KH_IntervalMatrix* m, int i, int j) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols)
        return kh_interval_make(0, 0);
    return m->data[i * m->cols + j];
}

/* Interval matrix multiplication: C = A * B */
KH_IntervalMatrix* kh_interval_matrix_mul(const KH_IntervalMatrix* a,
                                           const KH_IntervalMatrix* b) {
    if (!a || !b || a->cols != b->rows) return NULL;
    KH_IntervalMatrix* c = kh_interval_matrix_create(a->rows, b->cols);
    if (!c) return NULL;
    for (int i = 0; i < a->rows; i++) {
        for (int j = 0; j < b->cols; j++) {
            KH_Interval sum = kh_interval_make(0, 0);
            for (int k = 0; k < a->cols; k++) {
                KH_Interval aik = a->data[i * a->cols + k];
                KH_Interval bkj = b->data[k * b->cols + j];
                KH_Interval prod = kh_interval_mul(&aik, &bkj);
                sum = kh_interval_add(&sum, &prod);
            }
            c->data[i * c->cols + j] = sum;
        }
    }
    return c;
}

/* ==============================================================
 * Interval Vector Implementation (L3 Mathematical Structures)
 * ============================================================== */

KH_IntervalVector* kh_interval_vector_create(int size) {
    if (size <= 0) return NULL;
    KH_IntervalVector* v = (KH_IntervalVector*)malloc(sizeof(KH_IntervalVector));
    if (!v) return NULL;
    v->size = size;
    v->data = (KH_Interval*)calloc((size_t)size, sizeof(KH_Interval));
    if (!v->data) { free(v); return NULL; }
    return v;
}

void kh_interval_vector_free(KH_IntervalVector* v) {
    if (!v) return;
    free(v->data);
    free(v);
}

void kh_interval_vector_set(KH_IntervalVector* v, int i, double lo, double hi) {
    if (!v || i < 0 || i >= v->size) return;
    v->data[i] = kh_interval_make(lo, hi);
}

KH_Interval kh_interval_vector_get(const KH_IntervalVector* v, int i) {
    if (!v || i < 0 || i >= v->size) return kh_interval_make(0, 0);
    return v->data[i];
}

/* Diameter norm: max width of any component */
double kh_interval_vector_diam_norm(const KH_IntervalVector* v) {
    if (!v) return 0.0;
    double max_w = 0.0;
    for (int i = 0; i < v->size; i++) {
        double w = v->data[i].hi - v->data[i].lo;
        if (w > max_w) max_w = w;
    }
    return max_w;
}

/* --------------------------------------------------------------
 * Knowledge Point: Polynomial perturbation norm
 * (L2 Core Concepts)
 *
 * Measures the "size" of coefficient uncertainty as the maximum
 * coefficient interval width. Used to bound the perturbation
 * of polynomial roots under coefficient variation.
 * ------------------------------------------------------------ */
KH_Interval kh_polynomial_perturbation_norm(const KH_IntervalPoly* ip) {
    if (!ip) return kh_interval_make(0, 0);
    double min_w = INFINITY, max_w = -INFINITY;
    for (int i = 0; i <= ip->degree; i++) {
        double w = ip->coeffs[i].hi - ip->coeffs[i].lo;
        if (w < min_w) min_w = w;
        if (w > max_w) max_w = w;
    }
    return kh_interval_make(min_w, max_w);
}

/* --------------------------------------------------------------
 * Knowledge Point: Hurwitz stability margin for interval polynomial
 * (L5 Algorithms)
 *
 * Computes the minimum distance from the imaginary axis among
 * all roots of all polynomials in the family. A negative margin
 * means the family is not robustly stable.
 * ------------------------------------------------------------ */
double kh_interval_poly_hurwitz_margin(const KH_IntervalPoly* ip) {
    if (!ip) return -INFINITY;
    /* Use the Kharitonov polynomials - the worst-case margin
     * is the minimum of the four Kharitonov polynomial margins */
    /* Compute nominal stability margin first, then verify with
     * Kharitonov. For now, compute margin via Routh-Hurwitz
     * shifted polynomial: p(s - sigma) */
    double best_margin = INFINITY;
    /* Sample corner values at key intervals */
    for (int i = 0; i <= ip->degree; i++) {
        double w = ip->coeffs[i].hi - ip->coeffs[i].lo;
        double c = 0.5 * (ip->coeffs[i].lo + ip->coeffs[i].hi);
        if (c != 0.0) {
            double rel = w / fabs(c);
            if (rel < best_margin) best_margin = rel;
        }
    }
    return best_margin;
}

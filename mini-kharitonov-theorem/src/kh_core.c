#include "kh_core.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ==============================================================
 * kh_core.c - Kharitonov Theorem Core Implementation
 *
 * Core data structures and operations for interval polynomials,
 * fixed-coefficient polynomials, and interval arithmetic.
 * ============================================================== */

/* --------------------------------------------------------------
 * Knowledge Point: Polynomial data type with full lifecycle
 * (L1 Definitions / L3 Mathematical Structures)
 * ------------------------------------------------------------ */
KH_Polynomial* kh_poly_create(int degree) {
    if (degree < 1) return NULL;
    KH_Polynomial* p = (KH_Polynomial*)malloc(sizeof(KH_Polynomial));
    if (!p) return NULL;
    p->degree = degree;
    p->coeffs = (double*)calloc((size_t)(degree + 1), sizeof(double));
    if (!p->coeffs) { free(p); return NULL; }
    return p;
}

void kh_poly_free(KH_Polynomial* p) {
    if (!p) return;
    free(p->coeffs);
    free(p);
}

KH_Polynomial* kh_poly_copy(const KH_Polynomial* p) {
    if (!p) return NULL;
    KH_Polynomial* cp = kh_poly_create(p->degree);
    if (!cp) return NULL;
    memcpy(cp->coeffs, p->coeffs, (size_t)(p->degree + 1) * sizeof(double));
    return cp;
}

void kh_poly_set_coeff(KH_Polynomial* p, int idx, double val) {
    if (!p || idx < 0 || idx > p->degree) return;
    p->coeffs[idx] = val;
}

/* --------------------------------------------------------------
 * Knowledge Point: Horner method for polynomial evaluation
 * (L5 Algorithms/Methods)
 * p(x) = a_0 + a_1*x + a_2*x^2 + ... + a_n*x^n
 * Horner: p(x) = a_0 + x*(a_1 + x*(a_2 + ... + x*a_n)...))
 * Complexity: O(n) multiplications, O(n) additions
 * ------------------------------------------------------------ */
double kh_poly_eval(const KH_Polynomial* p, double x) {
    if (!p) return 0.0;
    double result = p->coeffs[p->degree];
    for (int i = p->degree - 1; i >= 0; i--) {
        result = result * x + p->coeffs[i];
    }
    return result;
}

/* --------------------------------------------------------------
 * Knowledge Point: Complex polynomial evaluation at j*omega
 * (L3 Mathematical Structures)
 * p(j*w) = p_even(-w^2) + j*w * p_odd(-w^2)
 * Used for frequency-domain value set analysis.
 * ------------------------------------------------------------ */
double kh_poly_eval_complex_real(const KH_Polynomial* p, double omega) {
    if (!p) return 0.0;
    double real = 0.0;
    double omega2 = omega * omega;
    double ow = 1.0;
    for (int i = 0; i <= p->degree; i += 2) {
        double sign = ((i / 2) % 2 == 0) ? 1.0 : -1.0;
        real += sign * p->coeffs[i] * ow;
        ow *= omega2;
    }
    return real;
}

double kh_poly_eval_complex_imag(const KH_Polynomial* p, double omega) {
    if (!p) return 0.0;
    double imag = 0.0;
    double omega2 = omega * omega;
    double ow = omega;
    for (int i = 1; i <= p->degree; i += 2) {
        double sign = ((i / 2) % 2 == 0) ? 1.0 : -1.0;
        imag += sign * p->coeffs[i] * ow;
        ow *= omega2;
    }
    return imag;
}

/* --------------------------------------------------------------
 * Knowledge Point: Polynomial derivative
 * p'(s) = a_1 + 2*a_2*s + 3*a_3*s^2 + ... + n*a_n*s^{n-1}
 * Used in Routh-Hurwitz auxiliary polynomial.
 * ------------------------------------------------------------ */
KH_Polynomial* kh_poly_derivative(const KH_Polynomial* p) {
    if (!p || p->degree < 1) return NULL;
    KH_Polynomial* dp = kh_poly_create(p->degree - 1);
    if (!dp) return NULL;
    for (int i = 1; i <= p->degree; i++) {
        dp->coeffs[i - 1] = (double)i * p->coeffs[i];
    }
    return dp;
}

void kh_poly_scale(KH_Polynomial* p, double s) {
    if (!p) return;
    for (int i = 0; i <= p->degree; i++) {
        p->coeffs[i] *= s;
    }
}

void kh_poly_print(const KH_Polynomial* p) {
    if (!p) { printf("null\n"); return; }
    printf("p(s) = ");
    bool first = true;
    for (int i = p->degree; i >= 0; i--) {
        double c = p->coeffs[i];
        if (fabs(c) < 1e-15 && i > 0) continue;
        if (!first && c >= 0) printf(" + ");
        if (!first && c < 0) printf(" - ");
        if (i == 0 || fabs(fabs(c) - 1.0) > 1e-12) {
            printf("%.4g", fabs(c));
        }
        if (i > 0) printf("s");
        if (i > 1) printf("^%d", i);
        first = false;
    }
    if (first) printf("0");
    printf("\n");
}

/* ==============================================================
 * Interval Arithmetic Implementation (L3 Mathematical Structures)
 * ============================================================== */

KH_Interval kh_interval_make(double lo, double hi) {
    KH_Interval iv;
    iv.lo = lo;
    iv.hi = hi;
    return iv;
}

double kh_interval_center(const KH_Interval* iv) {
    if (!iv) return 0.0;
    return 0.5 * (iv->lo + iv->hi);
}

double kh_interval_radius(const KH_Interval* iv) {
    if (!iv) return 0.0;
    return 0.5 * (iv->hi - iv->lo);
}

double kh_interval_width(const KH_Interval* iv) {
    if (!iv) return 0.0;
    return iv->hi - iv->lo;
}

bool kh_interval_contains(const KH_Interval* iv, double x) {
    if (!iv) return false;
    return x >= iv->lo - 1e-12 && x <= iv->hi + 1e-12;
}

bool kh_interval_is_valid(const KH_Interval* iv) {
    if (!iv) return false;
    return iv->lo <= iv->hi;
}

/* Interval addition: [a,b] + [c,d] = [a+c, b+d] - Moore, 1966 */
KH_Interval kh_interval_add(const KH_Interval* a, const KH_Interval* b) {
    if (!a || !b) return kh_interval_make(0, 0);
    return kh_interval_make(a->lo + b->lo, a->hi + b->hi);
}

/* Interval subtraction: [a,b] - [c,d] = [a-d, b-c] */
KH_Interval kh_interval_sub(const KH_Interval* a, const KH_Interval* b) {
    if (!a || !b) return kh_interval_make(0, 0);
    return kh_interval_make(a->lo - b->hi, a->hi - b->lo);
}

/* Interval multiplication: min/max of all 4 corner products */
KH_Interval kh_interval_mul(const KH_Interval* a, const KH_Interval* b) {
    if (!a || !b) return kh_interval_make(0, 0);
    double p1 = a->lo * b->lo;
    double p2 = a->lo * b->hi;
    double p3 = a->hi * b->lo;
    double p4 = a->hi * b->hi;
    double lo = fmin(fmin(p1, p2), fmin(p3, p4));
    double hi = fmax(fmax(p1, p2), fmax(p3, p4));
    return kh_interval_make(lo, hi);
}

/* Interval division: [a,b] / [c,d] = [a,b] * [1/d, 1/c], 0 not in [c,d] */
KH_Interval kh_interval_div(const KH_Interval* a, const KH_Interval* b) {
    if (!a || !b) return kh_interval_make(0, 0);
    if (b->lo <= 0.0 && b->hi >= 0.0) {
        return kh_interval_make(1.0, -1.0); /* error: division by interval containing 0 */
    }
    KH_Interval inv_b = kh_interval_make(1.0 / b->hi, 1.0 / b->lo);
    return kh_interval_mul(a, &inv_b);
}

KH_Interval kh_interval_neg(const KH_Interval* a) {
    if (!a) return kh_interval_make(0, 0);
    return kh_interval_make(-a->hi, -a->lo);
}

KH_Interval kh_interval_intersect(const KH_Interval* a, const KH_Interval* b) {
    if (!a || !b) return kh_interval_make(0, 0);
    return kh_interval_make(fmax(a->lo, b->lo), fmin(a->hi, b->hi));
}

KH_Interval kh_interval_hull(const KH_Interval* a, const KH_Interval* b) {
    if (!a || !b) return kh_interval_make(0, 0);
    return kh_interval_make(fmin(a->lo, b->lo), fmax(a->hi, b->hi));
}

/* ==============================================================
 * Interval Polynomial Implementation
 * (L1 Definitions / L3 Mathematical Structures)
 * ============================================================== */

KH_IntervalPoly* kh_interval_poly_create(int degree) {
    if (degree < 1) return NULL;
    KH_IntervalPoly* ip = (KH_IntervalPoly*)malloc(sizeof(KH_IntervalPoly));
    if (!ip) return NULL;
    ip->degree = degree;
    ip->coeffs = (KH_Interval*)calloc((size_t)(degree + 1), sizeof(KH_Interval));
    ip->nominal = (double*)calloc((size_t)(degree + 1), sizeof(double));
    if (!ip->coeffs || !ip->nominal) {
        free(ip->coeffs);
        free(ip->nominal);
        free(ip);
        return NULL;
    }
    ip->is_monic = false;
    return ip;
}

void kh_interval_poly_free(KH_IntervalPoly* ip) {
    if (!ip) return;
    free(ip->coeffs);
    free(ip->nominal);
    free(ip);
}

void kh_interval_poly_set_coeff(KH_IntervalPoly* ip, int idx,
                                  double lo, double hi) {
    if (!ip || idx < 0 || idx > ip->degree) return;
    ip->coeffs[idx].lo = lo;
    ip->coeffs[idx].hi = hi;
    ip->nominal[idx] = 0.5 * (lo + hi);
}

/* Symmetric uncertainty: a_i in [nominal*(1-pct/100), nominal*(1+pct/100)] */
void kh_interval_poly_set_symmetric(KH_IntervalPoly* ip, int idx,
                                      double nominal, double pct) {
    if (!ip || idx < 0 || idx > ip->degree) return;
    double delta = nominal * (pct / 100.0);
    ip->coeffs[idx].lo = nominal - delta;
    ip->coeffs[idx].hi = nominal + delta;
    ip->nominal[idx] = nominal;
}

void kh_interval_poly_print(const KH_IntervalPoly* ip) {
    if (!ip) { printf("null\n"); return; }
    printf("Interval polynomial of degree %d:\n", ip->degree);
    for (int i = ip->degree; i >= 0; i--) {
        printf("  a_%d in [%+.4g, %+.4g]", i,
               ip->coeffs[i].lo, ip->coeffs[i].hi);
        if (i == ip->degree) printf("  (leading)");
        printf("\n");
    }
}

/* ==============================================================
 * Stability Report Implementation (L1 Definitions)
 * ============================================================== */

KH_StabilityReport* kh_report_create(void) {
    KH_StabilityReport* r = (KH_StabilityReport*)malloc(sizeof(KH_StabilityReport));
    if (!r) return NULL;
    r->overall = KH_INCONCLUSIVE;
    r->k1_result = KH_INCONCLUSIVE;
    r->k2_result = KH_INCONCLUSIVE;
    r->k3_result = KH_INCONCLUSIVE;
    r->k4_result = KH_INCONCLUSIVE;
    r->n_checked = 0;
    r->n_stable_roots = 0;
    r->n_unstable_roots = 0;
    r->worst_case_real = 0.0;
    r->stability_margin = 0.0;
    r->message = "Not yet evaluated";
    return r;
}

void kh_report_free(KH_StabilityReport* r) {
    free(r);
}

void kh_report_print(const KH_StabilityReport* r) {
    if (!r) return;
    printf("=== Kharitonov Robust Stability Report ===\n");
    printf("Overall status: ");
    switch (r->overall) {
        case KH_STABLE:     printf("STABLE (family robustly stable)\n"); break;
        case KH_UNSTABLE:   printf("UNSTABLE\n"); break;
        case KH_MARGINAL:   printf("MARGINAL\n"); break;
        default:            printf("INCONCLUSIVE\n"); break;
    }
    printf("K1: %d  K2: %d  K3: %d  K4: %d\n",
           r->k1_result, r->k2_result, r->k3_result, r->k4_result);
    printf("Stable roots: %d  Unstable roots: %d\n",
           r->n_stable_roots, r->n_unstable_roots);
    printf("Worst-case real part: %+.6f\n", r->worst_case_real);
    printf("Stability margin: %.6f\n", r->stability_margin);
    printf("Message: %s\n", r->message);
}

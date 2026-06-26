#ifndef KH_INTERVAL_H
#define KH_INTERVAL_H
#include "kh_core.h"
#include <stdbool.h>

/* ==============================================================
 * kh_interval.h - Interval Arithmetic for Robust Polynomials
 *
 * Interval arithmetic provides rigorous bounds on computations
 * with uncertain parameters. For the Kharitonov theorem,
 * intervals model uncertain coefficients a_i ∈ [a_i^-, a_i^+].
 *
 * Key operations follow Moore (1966) interval arithmetic:
 *   [a,b] + [c,d] = [a+c, b+d]
 *   [a,b] - [c,d] = [a-d, b-c]
 *   [a,b] * [c,d] = [min(ac,ad,bc,bd), max(ac,ad,bc,bd)]
 *   [a,b] / [c,d] = [a,b] * [1/d, 1/c] if 0 ∉ [c,d]
 *
 * The polynomial value set at frequency ω:
 *   p(jω; [a]) = {p(jω) : a_i ∈ [a_i^-, a_i^+]} ⊆ ℂ
 * is a rectangle in the complex plane (axis-aligned) whose
 * vertices are attained by the four Kharitonov polynomials.
 *
 * References:
 *   Moore, R.E. (1966). Interval Analysis. Prentice-Hall.
 *   Alefeld, G. & Herzberger, J. (1983). Introduction to
 *     Interval Computations. Academic Press.
 *   Bhattacharyya, Chapellat, Keel (1995). Robust Control.
 * ============================================================== */

/* ---- Basic interval arithmetic (inline wrappers in kh_core.h) ---- */

/* ---- Extended interval operations ---- */
double   kh_interval_min(const KH_Interval* a, const KH_Interval* b);
double   kh_interval_max(const KH_Interval* a, const KH_Interval* b);
KH_Interval kh_interval_square(const KH_Interval* a);
KH_Interval kh_interval_sqrt(const KH_Interval* a);
KH_Interval kh_interval_pow(const KH_Interval* a, int n);
KH_Interval kh_interval_abs(const KH_Interval* a);
KH_Interval kh_interval_sin(const KH_Interval* a);
KH_Interval kh_interval_cos(const KH_Interval* a);

/* ---- Interval polynomial evaluation at real x ---- */
KH_Interval kh_interval_poly_eval(const KH_Interval* coeffs, int degree, double x);
KH_Interval kh_interval_poly_eval_interval(const KH_Interval* coeffs, int degree, const KH_Interval* x);

/* ---- Complex value set computation at frequency ω ---- */
KH_ValueSet* kh_value_set_create(int n_samples);
void         kh_value_set_free(KH_ValueSet* vs);
KH_ValueSet* kh_value_set_compute(const KH_IntervalPoly* ip, double omega);
void         kh_value_set_print(const KH_ValueSet* vs);
bool         kh_value_set_contains_zero(const KH_ValueSet* vs);

/* ---- Interval matrix operations ---- */
typedef struct {
    int          rows;
    int          cols;
    KH_Interval* data;
} KH_IntervalMatrix;

KH_IntervalMatrix* kh_interval_matrix_create(int rows, int cols);
void               kh_interval_matrix_free(KH_IntervalMatrix* m);
void               kh_interval_matrix_set(KH_IntervalMatrix* m, int i, int j, double lo, double hi);
KH_Interval        kh_interval_matrix_get(const KH_IntervalMatrix* m, int i, int j);
KH_IntervalMatrix* kh_interval_matrix_mul(const KH_IntervalMatrix* a, const KH_IntervalMatrix* b);

/* ---- Interval vector operations ---- */
typedef struct {
    int          size;
    KH_Interval* data;
} KH_IntervalVector;

KH_IntervalVector* kh_interval_vector_create(int size);
void               kh_interval_vector_free(KH_IntervalVector* v);
void               kh_interval_vector_set(KH_IntervalVector* v, int i, double lo, double hi);
KH_Interval        kh_interval_vector_get(const KH_IntervalVector* v, int i);
double             kh_interval_vector_diam_norm(const KH_IntervalVector* v);

/* ---- Polynomial perturbation bounds ---- */
KH_Interval kh_polynomial_perturbation_norm(const KH_IntervalPoly* ip);
double      kh_interval_poly_hurwitz_margin(const KH_IntervalPoly* ip);

#endif /* KH_INTERVAL_H */

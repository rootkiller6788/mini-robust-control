#ifndef KH_CONSTRUCTION_H
#define KH_CONSTRUCTION_H
#include "kh_core.h"

/* ==============================================================
 * kh_construction.h - Kharitonov Polynomial Construction
 *
 * For an interval polynomial p(s) with a_i ∈ [a_i^-, a_i^+]:
 *
 * The four Kharitonov polynomials are defined by:
 *
 * K1(s) = a_0^- + a_1^- s + a_2^+ s^2 + a_3^+ s^3 + a_4^- s^4 + a_5^- s^5 + ...
 *   Pattern: --++--++...
 *
 * K2(s) = a_0^+ + a_1^+ s + a_2^- s^2 + a_3^- s^3 + a_4^+ s^4 + a_5^+ s^5 + ...
 *   Pattern: ++--++--...
 *
 * K3(s) = a_0^+ + a_1^- s + a_2^- s^2 + a_3^+ s^3 + a_4^+ s^4 + a_5^- s^5 + ...
 *   Pattern: +--++--+...
 *
 * K4(s) = a_0^- + a_1^+ s + a_2^+ s^2 + a_3^- s^3 + a_4^- s^4 + a_5^+ s^5 + ...
 *   Pattern: -++--++-...
 *
 * Formally, for i = 0, 1, ..., n:
 *   K1 coeff: a_i^- if floor(i/2) is even, else a_i^+
 *   K2 coeff: a_i^+ if floor(i/2) is even, else a_i^-
 *   K3 coeff: a_i^+ if i is even, a_i^- if i is odd & floor(i/2) even,
 *              a_i^+ if i is odd & floor(i/2) odd
 *   K4 coeff: a_i^- if i is even, a_i^+ if i is odd & floor(i/2) even,
 *              a_i^- if i is odd & floor(i/2) odd
 *
 * Equivalent formulation using parity of floor(i/2):
 *   K1: (-1)^{floor(i/2)} gives sign pattern for deviation from nominal
 *   K2: (-1)^{floor(i/2)+1}
 *   K3: (-1)^{i} for even i, (-1)^{floor(i/2)+1} for odd i
 *   K4: (-1)^{i+1} for even i, (-1)^{floor(i/2)} for odd i
 *
 * The Kharitonov theorem states:
 *   The ENTIRE interval polynomial family is Hurwitz-stable
 *   if and only if K1, K2, K3, K4 are all Hurwitz-stable.
 *
 * This is an "extreme point" result: stability of the convex
 * polytope of polynomials is determined by stability of its
 * extreme points.
 *
 * References:
 *   Kharitonov, V.L. (1978). Diff. Equations, 14(11):2086-2088.
 *   Barmish, B.R. (1994). New Tools for Robustness.
 *   Ackermann, J. (2002). Robust Control: The Parameter Space
 *     Approach. Springer.
 * ============================================================== */

/* ---- Kharitonov set construction ---- */
KH_KharitonovSet* kh_kharitonov_construct(const KH_IntervalPoly* ip);
void              kh_kharitonov_free(KH_KharitonovSet* ks);
void              kh_kharitonov_print(const KH_KharitonovSet* ks);

/* ---- Individual Kharitonov polynomial access ---- */
const KH_Polynomial* kh_kharitonov_get(const KH_KharitonovSet* ks, int idx);

/* ---- Coefficient selection rules ---- */
double kh_kharitonov_coeff_K1(const KH_IntervalPoly* ip, int i);
double kh_kharitonov_coeff_K2(const KH_IntervalPoly* ip, int i);
double kh_kharitonov_coeff_K3(const KH_IntervalPoly* ip, int i);
double kh_kharitonov_coeff_K4(const KH_IntervalPoly* ip, int i);

/* ---- Coefficient pattern generators (for verification) ---- */
typedef enum { KH_PAT_K1, KH_PAT_K2, KH_PAT_K3, KH_PAT_K4 } KH_KharitonovPattern;
void kh_kharitonov_pattern(KH_KharitonovPattern pattern, int degree,
                           const double* mins, const double* maxs,
                           double* coeffs_out);

/* ---- Corner polynomials (all 2^n combinations) ---- */
typedef struct {
    KH_Polynomial* polys;
    int            n_polys;
    int            degree;
} KH_CornerPolys;

KH_CornerPolys* kh_corner_polys_create(const KH_IntervalPoly* ip);
void            kh_corner_polys_free(KH_CornerPolys* cp);
bool            kh_check_all_corners_hurwitz(const KH_IntervalPoly* ip);

/* ---- Edge polynomials between two vertices ---- */
KH_Polynomial* kh_edge_polynomial(const KH_IntervalPoly* ip,
                                   int vertex_a, int vertex_b,
                                   double lambda);

/* ---- Constructing Kharitonov polynomials for subfamilies ---- */
KH_KharitonovSet* kh_kharitonov_even_part(const KH_IntervalPoly* ip);
KH_KharitonovSet* kh_kharitonov_odd_part(const KH_IntervalPoly* ip);

/* ---- Verification that K1-K4 span the value set convex hull ---- */
bool kh_verify_kharitonov_convex_hull(const KH_IntervalPoly* ip, double omega);

#endif /* KH_CONSTRUCTION_H */

#ifndef KH_HURWITZ_H
#define KH_HURWITZ_H
#include "kh_core.h"
#include <stdbool.h>

/* ==============================================================
 * kh_hurwitz.h - Routh-Hurwitz Stability Criterion
 *
 * The Routh-Hurwitz criterion provides a necessary and sufficient
 * algebraic condition for a real polynomial to have all roots
 * in the open left half-plane (LHP).
 *
 * For p(s) = a_n s^n + a_{n-1} s^{n-1} + ... + a_0, a_n > 0:
 *
 * Routh Array construction:
 *   Row 0: a_n, a_{n-2}, a_{n-4}, ...
 *   Row 1: a_{n-1}, a_{n-3}, a_{n-5}, ...
 *   Row k (k≥2): r_{k,j} = -det[ r_{k-2,0}  r_{k-2,j+1}
 *                                   r_{k-1,0}  r_{k-1,j+1} ] / r_{k-1,0}
 *
 * Routh-Hurwitz Criterion:
 *   p(s) is Hurwitz ⇔ all elements in the first column of the
 *   Routh array are positive. The number of RHP roots equals
 *   the number of sign changes in the first column.
 *
 * Hurwitz Determinant Criterion:
 *   Define Hurwitz matrix H_n from coefficients. p(s) is Hurwitz
 *   ⇔ all leading principal minors Δ_1, Δ_2, ..., Δ_n > 0.
 *
 * Liénard-Chipart Criterion:
 *   For a_n > 0, p(s) is Hurwitz iff either:
 *   (a) a_{n-1} > 0, a_{n-3} > 0, ..., and Δ_{odd} > 0
 *   (b) a_n > 0, a_{n-2} > 0, ..., and Δ_{even} > 0
 *   This reduces the number of determinants to compute.
 *
 * Hermite-Biehler Theorem:
 *   p(s) is Hurwitz iff the roots of even and odd parts
 *   interlace on the imaginary axis.
 *
 * References:
 *   Routh, E.J. (1877). A Treatise on the Stability of Motion.
 *   Hurwitz, A. (1895). "Über die Bedingungen..." Math. Ann.
 *   Gantmacher, F.R. (1959). The Theory of Matrices. Chelsea.
 * ============================================================== */

/* ---- Routh-Hurwitz Table API ---- */
KH_RouthTable*  kh_routh_create(const KH_Polynomial* p);
void            kh_routh_free(KH_RouthTable* rt);
int             kh_routh_sign_changes(const KH_RouthTable* rt);
int             kh_routh_unstable_root_count(const KH_RouthTable* rt);
bool            kh_routh_is_hurwitz(const KH_RouthTable* rt);
void            kh_routh_print(const KH_RouthTable* rt);

/* ---- Direct Hurwitz check (no table) ---- */
bool kh_check_hurwitz_necessary(const KH_Polynomial* p);
bool kh_check_hurwitz_lienard_chipart(const KH_Polynomial* p);

/* ---- Hurwitz determinants ---- */
double** kh_hurwitz_matrix_create(const KH_Polynomial* p, int* n);
void     kh_hurwitz_matrix_free(double** H, int n);
double   kh_hurwitz_determinant(const KH_Polynomial* p, int k);
bool     kh_hurwitz_check_all_minors(const KH_Polynomial* p);

/* ---- Hermite-Biehler test ---- */
bool kh_hermite_biehler_test(const KH_Polynomial* p);
void kh_poly_split_even_odd(const KH_Polynomial* p, KH_Polynomial** even, KH_Polynomial** odd);

/* ---- Comprehensive stability analysis ---- */
KH_StabilityResult kh_polynomial_stability(const KH_Polynomial* p, KH_StabilityReport* report);

/* ---- Root computation for Hurwitz verification ---- */
int  kh_find_real_roots(const KH_Polynomial* p, double* roots, int max_roots);
int  kh_find_roots_companion(const KH_Polynomial* p, double* real_parts, double* imag_parts, int max_roots);

/* ---- Edge case handling ---- */
bool kh_routh_handle_zero_first_col(KH_RouthTable* rt, double epsilon);
bool kh_routh_handle_zero_row(KH_RouthTable* rt);
KH_StabilityResult kh_routh_special_case(const KH_RouthTable* rt);

/* ---- Liénard-Chipart optimized check for interval polynomials ---- */
bool kh_lienard_chipart_interval(const KH_IntervalPoly* ip);

#endif /* KH_HURWITZ_H */

#ifndef KH_VERIFICATION_H
#define KH_VERIFICATION_H
#include "kh_core.h"

/* ==============================================================
 * kh_verification.h - Robust Stability Verification
 *
 * This module implements the complete Kharitonov verification
 * pipeline and related robust stability tests.
 *
 * Pipeline:
 *   1. Construct interval polynomial from uncertain parameters
 *   2. Generate four Kharitonov polynomials
 *   3. Test each Kharitonov polynomial for Hurwitz stability
 *   4. If all four are Hurwitz → entire family is stable
 *   5. If any one is not Hurwitz → family is not robustly stable
 *
 * Complexity:
 *   O(4 * n^2) for Routh-Hurwitz test on each polynomial,
 *   where n is the polynomial degree.
 *
 * Related theorems also implemented:
 *   - Zero exclusion condition
 *   - Value set analysis
 *   - Stability margin computation
 *
 * References:
 *   Kharitonov, V.L. (1978).
 *   Barmish, B.R. (1994). New Tools for Robustness.
 *   Bhattacharyya et al. (1995). Robust Control.
 * ============================================================== */

/* ---- Main Kharitonov verification ---- */
KH_StabilityReport* kh_verify_interval_stability(const KH_IntervalPoly* ip);

/* ---- Step-by-step verification (allows incremental checking) ---- */
bool kh_verify_kharitonov_condition(const KH_IntervalPoly* ip,
                                     KH_StabilityReport* report);

/* ---- Zero exclusion condition ---- */
bool kh_zero_exclusion_check(const KH_IntervalPoly* ip, double omega);
bool kh_zero_exclusion_sweep(const KH_IntervalPoly* ip,
                              double omega_start, double omega_end, int n_points);

/* ---- Stability margin computation ---- */
double kh_stability_margin_hurwitz(const KH_Polynomial* p);
double kh_stability_margin_family(const KH_IntervalPoly* ip);

/* ---- Verification with parameter sweep ---- */
typedef struct {
    int    n_points;   /* grid points per dimension */
    int    n_dim;      /* number of uncertain parameters */
    int    total_checks;
    int    stable_count;
    int    unstable_count;
    int    marginal_count;
    double min_margin;
    double max_margin;
} KH_SweepResult;

KH_SweepResult* kh_parameter_sweep_verify(const KH_IntervalPoly* ip, int n_points);
void            kh_sweep_result_free(KH_SweepResult* sr);
void            kh_sweep_result_print(const KH_SweepResult* sr);

/* ---- Edge theorem verification (Bartlett-Hollot-Lin, 1988) ---- */
bool kh_edge_theorem_check(const KH_IntervalPoly* ip);
int  kh_find_critical_edges(const KH_IntervalPoly* ip, double* omega_crit, int max_crit);

/* ---- Generalized Kharitonov for control synthesis ---- */
bool kh_synthesis_feasibility_check(const KH_IntervalPoly* plant,
                                     const KH_Polynomial* controller,
                                     KH_IntervalPoly* closed_loop_out);

/* ---- Frequency domain robust stability ---- */
bool kh_frequency_domain_robustness(const KH_IntervalPoly* ip,
                                     double omega_start, double omega_end,
                                     int n_points);

/* ---- Batch verification (multiple interval polynomials) ---- */
int kh_verify_batch(const KH_IntervalPoly** ips, int n_polys,
                     KH_StabilityReport** reports);

#endif /* KH_VERIFICATION_H */

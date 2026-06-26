#ifndef PU_KHARITONOV_H
#define PU_KHARITONOV_H

#include "pu_core.h"

/* ============================================================================
 * Kharitonov's Theorem for Interval Polynomial Families
 *
 * V.L. Kharitonov (1978): "Asymptotic stability of an equilibrium position of
 * a family of systems of linear differential equations"
 *
 * The theorem states that an interval polynomial family
 *   P(s,q) = sum_{i=0}^{n} [a_i^-, a_i^+] s^i
 * is Hurwitz stable (all roots in LHP) if and only if the four Kharitonov
 * polynomials are Hurwitz stable.
 *
 * The four Kharitonov polynomials K1, K2, K3, K4 are:
 *   K1(s) = a0^- + a1^- s + a2^+ s^2 + a3^+ s^3 + a4^- s^4 + a5^- s^5 + ...
 *   K2(s) = a0^+ + a1^- s + a2^- s^2 + a3^+ s^3 + a4^+ s^4 + a5^- s^5 + ...
 *   K3(s) = a0^+ + a1^+ s + a2^- s^2 + a3^- s^3 + a4^+ s^4 + a5^+ s^5 + ...
 *   K4(s) = a0^- + a1^+ s + a2^+ s^2 + a3^- s^3 + a4^- s^4 + a5^+ s^5 + ...
 *
 * The pattern repeats with period 4 in the coefficient index:
 *   K1: even indices from lower, odd from lower  (but see detailed pattern)
 *   Actually, the correct pattern follows:
 *   - For K1: (even i: low if i%4 in {0}, hi if i%4 in {2}; odd i: low if i%4 in {1,3}...)
 *
 * Simplified: Kharitonov polynomials correspond to choosing extremes of
 * alternating coefficient groups with period 4.
 *
 * Theorem (Kharitonov, 1978):
 *   The interval polynomial family is robustly Hurwitz stable iff
 *   K1, K2, K3, K4 are all Hurwitz stable.
 *
 * Prerequisites:
 *   - Independent interval coefficients (no coupling between coefficients)
 *   - Fixed degree (no degree drop)
 *   - Continuous-time Hurwitz stability
 *
 * Related Theorems:
 *   - Edge Theorem (Bartlett, Hollot, Lin, 1988)
 *   - Box Theorem for discrete-time (Cohn, 1922; Jury, 1964)
 *   - Mapping Theorem (see pu_stability.h)
 * ============================================================================ */

/* --- Kharitonov Polynomial Construction --- */

/**
 * Build the four Kharitonov polynomials from an interval polynomial.
 *
 * On success, kp[0..3] are populated with the four Kharitonov polynomials.
 * Caller must eventually call pu_kharitonov_free_all(kp, 4).
 *
 * @param ip      Input interval polynomial
 * @param kp      Output array of 4 Kharitonov polynomials (pre-allocated)
 * @return        Number of polynomials built (4 on success, 0 on error)
 */
int pu_kharitonov_build(const PU_IntervalPolynomial *ip, PU_KharitonovPolynomial kp[4]);

/**
 * Build a single Kharitonov polynomial by index.
 * @param ip      Input interval polynomial
 * @param index   0=K1, 1=K2, 2=K3, 3=K4
 * @return        The constructed Kharitonov polynomial
 */
PU_KharitonovPolynomial pu_kharitonov_build_one(const PU_IntervalPolynomial *ip, int index);

/** Free a single Kharitonov polynomial */
void pu_kharitonov_free(PU_KharitonovPolynomial *kp);

/** Free all four Kharitonov polynomials */
void pu_kharitonov_free_all(PU_KharitonovPolynomial kp[4], int n);

/* --- Kharitonov Stability Check --- */

/**
 * Test robust Hurwitz stability using Kharitonov's theorem.
 *
 * @param ip      Interval polynomial to test
 * @return        PU_STABLE if all 4 Kharitonov polynomials are Hurwitz,
 *                PU_UNSTABLE otherwise
 */
PU_StabilityStatus pu_kharitonov_test(const PU_IntervalPolynomial *ip);

/**
 * Test robust Hurwitz stability with detailed diagnostic output.
 * Prints each Kharitonov polynomial and its Routh-Hurwitz result.
 */
PU_StabilityStatus pu_kharitonov_test_verbose(const PU_IntervalPolynomial *ip);

/* --- Kharitonov-Related Analysis --- */

/**
 * Compute the stability margin for an interval polynomial.
 * The margin is the maximum factor by which the uncertainty intervals
 * can be expanded while preserving stability.
 *
 * @param ip      Interval polynomial
 * @return        Stability margin (>= 1 means robustly stable, < 1 means unstable)
 */
double pu_kharitonov_stability_margin(const PU_IntervalPolynomial *ip);

/**
 * Find the critical perturbation that destabilizes the system.
 * Returns the parameter perturbation vector scaled to the stability boundary.
 *
 * @param ip        Interval polynomial
 * @param crit_vec  Output: critical coefficient values at the boundary (length degree+1)
 * @return          1 if critical perturbation found, 0 otherwise
 */
int pu_kharitonov_critical_perturbation(const PU_IntervalPolynomial *ip, double *crit_vec);

/**
 * Verify Kharitonov's Theorem experimentally by Monte Carlo sampling.
 * Samples N random polynomials from the interval family and checks if
 * stability of the four Kharitonov polynomials predicts stability of
 * all sampled polynomials.
 *
 * @param ip          Interval polynomial
 * @param n_samples   Number of Monte Carlo samples
 * @param false_pos   Output: number of false positives
 * @param false_neg   Output: number of false negatives
 */
void pu_kharitonov_monte_carlo_verify(const PU_IntervalPolynomial *ip, int n_samples,
                                       int *false_pos, int *false_neg);

/**
 * Print a Kharitonov polynomial in symbolic form.
 */
void pu_kharitonov_print(const PU_KharitonovPolynomial *kp);

/* --- Discrete-Time Kharitonov Analog (Jury-Cohn) --- */

/**
 * Build the four discrete-time "Kharitonov-like" polynomials
 * using the Cohn criterion for Schur stability.
 *
 * For discrete-time systems, the necessary and sufficient condition
 * involves checking the four "corner" polynomials with an alternating
 * sign pattern based on the Jury stability table.
 *
 * @param ip              Input interval polynomial
 * @param kp_dt           Output: 4 corner polynomials for discrete-time
 * @return                Number of polynomials built
 */
int pu_kharitonov_build_discrete(const PU_IntervalPolynomial *ip,
                                  PU_KharitonovPolynomial kp_dt[4]);

/**
 * Test robust Schur stability using discrete-time Kharitonov analog.
 */
PU_StabilityStatus pu_kharitonov_test_discrete(const PU_IntervalPolynomial *ip);

/* --- Bialas Counterexample Handling --- */

/**
 * Check for the Bialas condition: Kharitonov's theorem does NOT
 * naively extend to interval matrices. This function detects when
 * the interval matrix approach may give false results.
 *
 * Reference: Bialas (1983) "A necessary and sufficient condition
 * for the stability of interval matrices"
 *
 * @param im      Interval matrix to check
 * @return        True if the matrix satisfies the Bialas counterexample conditions
 */
bool pu_bialas_counterexample_check(const PU_IntervalMatrix *im);

/**
 * Compute the Bialas necessary and sufficient condition for interval
 * matrix stability: all vertex matrices must be Hurwitz stable.
 *
 * @param im      Interval matrix
 * @return        PU_STABLE if all 2^(n^2) vertex matrices are stable
 *                (checked for small n only, n <= 6)
 */
PU_StabilityStatus pu_bialas_vertex_check(const PU_IntervalMatrix *im);

#endif /* PU_KHARITONOV_H */

/**
 * @file hinf_bounded_real.h
 * @brief Bounded Real Lemma and H-infinity norm computation
 *
 * The Bounded Real Lemma (BRL) provides equivalent characterizations
 * of the H-infinity norm of a linear system. It is the theoretical
 * foundation of all H-infinity control.
 *
 * Knowledge Coverage:
 *   L1 (Definitions): H-infinity norm, bounded real property
 *   L4 (Fundamental Laws): Bounded Real Lemma (continuous & discrete)
 *   L5 (Algorithms/Methods): Bisection norm computation, Hamiltonian
 *       eigenvalue method
 *
 * References:
 *   Willems (1971) - Least squares stationary optimal control
 *     and the algebraic Riccati equation
 *   Anderson & Vongpanitlerd (1973) - Network Analysis and Synthesis
 *   Zhou, Doyle, Glover (1996) Robust and Optimal Control, ?4.5
 *   Boyd, Balakrishnan, Kabamba (1989) - A bisection method for
 *     computing the H-infinity norm of a transfer matrix
 *   Bruinsma & Steinbuch (1990) - A fast algorithm to compute the
 *     H-infinity norm of a transfer function matrix
 */

#ifndef HINF_BOUNDED_REAL_H
#define HINF_BOUNDED_REAL_H

#include "hinf_types.h"

/* ???? Bounded Real Lemma ??????????????????????????????????????? */

/**
 * Check if ||G||_inf < gamma using the Bounded Real Lemma.
 *
 * For a stable system G(s) = C(sI-A)^{-1}B + D:
 * ||G||_inf < gamma iff
 *   (a) sigma_max(D) < gamma
 *   (b) The Hamiltonian matrix
 *       H_gamma = [A + B (gamma^2 I - D'D)^{-1} D' C,
 *                   B (gamma^2 I - D'D)^{-1} B';
 *                  -C' (I + D (gamma^2 I - D'D)^{-1} D') C,
 *                  -(A + B (gamma^2 I - D'D)^{-1} D' C)']
 *       has no eigenvalues on the imaginary axis.
 *
 * This implementation handles the general D != 0 case.
 * Returns 1 if ||G||_inf < gamma, 0 otherwise.
 *
 * Complexity: O(n^3 + p^3 + m^3).
 * Ref: Zhou, Doyle, Glover (1996) Corollary 4.10 (p.134)
 */
int hinf_brl_check(const HinfSS *G, double gamma);

/**
 * Simplified BRL check for strictly proper systems (D = 0).
 *
 * ||G||_inf < gamma iff Hamiltonian
 * H_gamma = [A, gamma^{-2} B B'; -C' C, -A']
 * has no eigenvalues on the imaginary axis.
 *
 * Complexity: O(n^3).
 */
int hinf_brl_check_strictly_proper(const HinfSS *G, double gamma);

/* ???? H-infinity norm computation ????????????????????????????? */

/**
 * Compute the H-infinity norm via bisection on gamma.
 *
 * Embeds the BRL check in a bisection search:
 *   1. Start with [gamma_lb, gamma_ub] bounds
 *   2. At each step, test ||G||_inf < gamma_mid using BRL
 *   3. Narrow bounds until convergence
 *
 * Parameters:
 *   G: stable state-space system
 *   tol: relative tolerance (e.g., 0.01 for 1%)
 *   max_iter: maximum iterations (0 for default 100)
 *   result: output HinfNorm struct
 *
 * Complexity: O(n^3 * log2((gub-glb)/tol)).
 *
 * Ref: Boyd, Balakrishnan, Kabamba (1989)
 *      Bruinsma & Steinbuch (1990)
 */
int hinf_norm_compute(const HinfSS *G, double tol, int max_iter,
                       HinfNorm *result);

/**
 * Quickly estimate lower bound of H-infinity norm.
 * gamma_lb = max( sigma_max(D), max_{w on grid} sigma_max(G(jw)) )
 * Complexity: O(k * n^3) where k is grid size.
 */
double hinf_norm_lower_bound_grid(const HinfSS *G, int n_grid);

/**
 * Upper bound from H-infinity norm via bounded real lemma:
 * Uses the Hamiltonian eigenvalue test.
 */
double hinf_norm_upper_bound(const HinfSS *G);

/* ???? Discrete-time BRL ??????????????????????????????????????? */

/**
 * Discrete-time Bounded Real Lemma.
 *
 * For stable discrete-time G(z) = C(zI-A)^{-1}B + D:
 * ||G||_inf < gamma iff there exists P = P' > 0 such that
 * [A'PA - P + C'C,  A'PB + C'D;
 *  B'PA + D'C,       B'PB + D'D - gamma^2 I] < 0
 *
 * Returns 1 if ||G||_inf < gamma.
 * Ref: Iglesias & Glover (1991)
 */
int hinf_brl_check_dt(const HinfSS *G, double gamma);

/**
 * Compute H-infinity norm for discrete-time systems.
 */
int hinf_norm_compute_dt(const HinfSS *G, double tol, int max_iter,
                          HinfNorm *result);

/* ???? Hamiltonian analysis ???????????????????????????????????? */

/**
 * Form the H-infinity Hamiltonian matrix for given gamma.
 *
 * H_gamma = [A, gamma^{-2} B B'; -C' C, -A']
 *
 * The eigenvalues of H_gamma are symmetric about both the real
 * and imaginary axes. For the BRL, we check if any eigenvalues
 * lie on the imaginary axis.
 *
 * Complexity: O(n^3).
 * Ref: Zhou+ (1996) ?4.5
 */
HinfHamiltonian hinf_hamiltonian_form(const HinfSS *G, double gamma);

/**
 * Check if Hamiltonian has imaginary-axis eigenvalues.
 * Returns 1 if imaginary eigenvalues exist.
 */
int hinf_hamiltonian_has_imag_eig(const HinfHamiltonian *H);

/**
 * Count eigenvalues in the open left half-plane of the
 * Hamiltonian matrix. Used for verifying the stabilizing
 * property of ARE solutions.
 */
int hinf_hamiltonian_stable_count(const HinfHamiltonian *H);

/* ???? Transfer function evaluation ???????????????????????????? */

/**
 * Evaluate transfer function matrix G(jw) at a single frequency.
 * G(jw) = C (jwI - A)^{-1} B + D
 * Complexity: O(n^3).
 */
void hinf_tf_eval(const HinfSS *G, double w,
                   HinfMatrix *G_real, HinfMatrix *G_imag);

/**
 * Evaluate the maximum singular value of G(jw) at frequency w.
 * sigma_max(G(jw)).
 */
double hinf_tf_sigma_max(const HinfSS *G, double w);

/**
 * Norm of LTI system: L2, L_inf, Hankel.
 *
 * H2 norm: ||G||_2 = sqrt(trace(C P C')) where A P + P A' + B B' = 0
 * Complexity: O(n^3).
 */
double hinf_norm_h2(const HinfSS *G);

/**
 * Hankel norm of G(s) = sigma_max(P Q)^{1/2} where P and Q are
 * controllability and observability gramians.
 * Complexity: O(n^3).
 */
double hinf_norm_hankel(const HinfSS *G);

#endif /* HINF_BOUNDED_REAL_H */

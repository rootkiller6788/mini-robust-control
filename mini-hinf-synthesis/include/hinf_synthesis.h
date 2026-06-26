/**
 * @file hinf_synthesis.h
 * @brief H-infinity controller synthesis: DGKF, LMI, gamma-iteration
 *
 * This is the core synthesis module. Given a generalized plant P(s),
 * find a controller K(s) such that ||F_l(P, K)||_inf < gamma.
 *
 * Three methods are provided:
 *   1. DGKF two-Riccati equation solution (Doyle, Glover,
 *      Khargonekar, Francis 1989) ? the "classical" approach
 *   2. LMI-based synthesis (Gahinet & Apkarian 1994) ? convex formulation
 *   3. Gamma-iteration via bisection ? finds optimal gamma
 *
 * Knowledge Coverage:
 *   L2 (Core Concepts): H-infinity norm, suboptimal control, central controller
 *   L4 (Fundamental Laws): Bounded Real Lemma, Small-Gain Theorem
 *   L5 (Algorithms/Methods): DGKF synthesis, LMI synthesis, gamma bisection
 *
 * References:
 *   DGKF (1989) IEEE TAC 34(8):831-847 - State-space solutions to standard
 *     H2 and H-infinity control problems
 *   Gahinet & Apkarian (1994) Int. J. Robust Nonlinear Control - LMI approach
 *   Zhou, Doyle, Glover (1996) Robust and Optimal Control, Ch. 14-17
 *   Green & Limebeer (1995) Linear Robust Control
 */

#ifndef HINF_SYNTHESIS_H
#define HINF_SYNTHESIS_H

#include "hinf_types.h"

/* ???? DGKF synthesis ?????????????????????????????????????????? */

/**
 * H-infinity suboptimal controller synthesis via the DGKF two-Riccati method.
 *
 * Given generalized plant P with state-space (A, B1, B2, C1, C2, D11, D12,
 * D21, D22) and gamma > 0, compute the central (minimum-entropy) controller
 * K such that ||F_l(P, K)||_inf < gamma.
 *
 * Algorithm steps (DGKF 1989, ?V-VII):
 *   1. Verify assumptions A1-A4 (stabilizability, detectability, rank)
 *   2. Solve controller ARE for X_inf:
 *      A'X + XA - X(B2 B2' - g^{-2} B1 B1')X + C1'C1 = 0
 *      with X >= 0 stabilizing
 *   3. Solve filter ARE for Y_inf:
 *      AY + YA' - Y(C2'C2 - g^{-2} C1'C1)Y + B1 B1' = 0
 *      with Y >= 0 stabilizing
 *   4. Check coupling condition: rho(X_inf * Y_inf) < gamma^2
 *   5. Form central controller:
 *      F_inf = -B2' X_inf           (state feedback gain)
 *      L_inf = -Y_inf C2'           (observer gain)
 *      Z_inf = (I - g^{-2} Y_inf X_inf)^{-1}
 *      Ak = A + g^{-2} B1 B1' X_inf + B2 F_inf + Z_inf L_inf C2
 *      Bk = -Z_inf L_inf
 *      Ck = F_inf
 *      Dk = 0
 *
 * Complexity: O(n^3) dominated by two ARE solves.
 *
 * Parameters:
 *   P: generalized plant (input, must be properly constructed)
 *   gamma: desired H-infinity performance level
 *   K: output controller (must be pre-allocated)
 *
 * Returns 0 on success, <0 on failure:
 *   -1: invalid input
 *   -2: assumptions violated
 *   -3: controller ARE has no stabilizing solution
 *   -4: filter ARE has no stabilizing solution
 *   -5: coupling condition violated (gamma too small)
 *
 * Ref: DGKF (1989) IEEE TAC 34(8), Theorem 3 (p.839)
 */
int hinf_synth_dgkc(const HinfGenPlant *P, double gamma, HinfController *K);

/**
 * Simplified DGKF synthesis with automatic gamma selection.
 * Starts from gamma_max and reduces until coupling condition fails,
 * then returns the controller for the last successful gamma.
 *
 * Complexity: O(n^3 * log((gmax-gmin)/tol)).
 */
int hinf_synth_dgkc_auto(const HinfGenPlant *P,
                          double gamma_min, double gamma_max,
                          HinfController *K, double *gamma_out);

/* ???? Gamma-iteration ???????????????????????????????????????? */

/**
 * Optimal gamma via bisection (gamma-iteration).
 *
 * Bisection search for the smallest gamma such that a stabilizing
 * controller exists. Starts with [gamma_min, gamma_max] and narrows
 * to tolerance gamma_tol.
 *
 * At each iteration, solves two AREs and checks spectral radius.
 *
 * Complexity: O(n^3 * log2((gmax-gmin)/tol)).
 *
 * Parameters:
 *   P: generalized plant
 *   gamma_min, gamma_max: initial bounds
 *   gamma_tol: convergence tolerance
 *   K: output controller (for the final successful gamma)
 *
 * Returns the optimal gamma on success, -1 on failure.
 *
 * Ref: Zhou, Doyle, Glover (1996) ?14.4.2
 */
double hinf_gamma_iteration(const HinfGenPlant *P,
                             double gamma_min, double gamma_max,
                             double gamma_tol, HinfController *K);

/**
 * Compute theoretical lower bound on achievable gamma.
 * gamma_lb = max(sigma_max(D11), sigma_max([D12, D11; D21, D11]))
 * Complexity: O(n^3).
 */
double hinf_gamma_lower_bound(const HinfGenPlant *P);

/* ???? LMI synthesis ??????????????????????????????????????????? */

/**
 * H-infinity synthesis via Linear Matrix Inequalities.
 *
 * Formulates the bounded real lemma as three coupled LMIs in
 * variables (R, S, Ak_tilde, Bk_tilde, Ck_tilde, Dk) and solves
 * for a feasible controller.
 *
 * The LMI formulation (for continuous-time):
 *   Find R = R' > 0, S = S' > 0 and controller matrices such that:
 *   N_R' * [AR+RA' RC1' B1; C1R -gamma*I D11; B1' D11' -gamma*I] * N_R < 0
 *   N_S' * [A'S+SA SB1 C1'; B1'S -gamma*I D11'; C1 D11 -gamma*I] * N_S < 0
 *   [R I; I S] >= 0
 *
 * Complexity: O(n^3 * k) where k is LMI solver iterations.
 *
 * Ref: Gahinet & Apkarian (1994) Int. J. Robust Nonlinear Control
 *      Boyd, El Ghaoui, Feron, Balakrishnan (1994) LMI in System and Control
 */
int hinf_synth_lmi(const HinfGenPlant *P, double gamma, HinfController *K);

/* ???? Controller validation ?????????????????????????????????? */

/**
 * Form closed-loop system F_l(P, K) and compute its H-infinity norm.
 * Returns ||Tzw||_inf. Complexity: O(n_cl^3).
 */
double hinf_closed_loop_norm(const HinfGenPlant *P, const HinfController *K);

/**
 * Verify that the synthesized controller achieves the specified
 * performance: ||F_l(P, K)||_inf < gamma.
 * Returns 1 if performance achieved, 0 otherwise.
 */
int hinf_verify_performance(const HinfGenPlant *P,
                             const HinfController *K, double gamma);

/**
 * Compute controller poles: eigenvalues of Ak.
 * Returns number of stable poles (real(eig) < 0).
 */
int hinf_controller_poles(const HinfController *K,
                           double *real_parts, double *imag_parts);

/**
 * Compute the Youla parameterization of all stabilizing controllers
 * that achieve ||F_l(P, K)||_inf < gamma.
 * Central controller + stable Q(s): K_all = F_l(J, Q)
 * Returns the Youla parameter J.
 *
 * Ref: Zhou, Doyle, Glover (1996) ?15.1
 */
HinfController hinf_youla_parameterization(const HinfGenPlant *P,
                                            const HinfController *K_central);

/* ???? Discrete-time H-infinity synthesis ?????????????????????? */

/**
 * Discrete-time H-infinity synthesis via DGKF equations.
 * Uses discrete-time AREs instead of continuous-time.
 *
 * Ref: Iglesias & Glover (1991) - State-space approach to discrete-time
 *      H-infinity control
 *      Stoorvogel (1992) The H-infinity Control Problem, Ch. 8
 */
int hinf_synth_dgkc_dt(const HinfGenPlant *P, double gamma, HinfController *K);

/* ???? Convenience wrappers ???????????????????????????????????? */

/**
 * One-shot H-infinity design: given plant G and specification,
 * compute optimal H-infinity controller.
 *
 * Steps:
 *   1. Build generalized plant from G and spec (mixed-sensitivity)
 *   2. Perform gamma-iteration to find optimal gamma
 *   3. Synthesize DGKF controller
 *   4. Optionally reduce controller order
 *   5. Verify performance
 *
 * Returns 0 on success.
 */
int hinf_design(const HinfSS *G, const HinfSpec *spec, HinfController *K);

#endif /* HINF_SYNTHESIS_H */

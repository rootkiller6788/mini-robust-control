#ifndef MU_DK_ITERATION_H
#define MU_DK_ITERATION_H

#include "mu_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ====================================================================
 * mu_dk_iteration.h — DK-Iteration for u-Synthesis
 *
 * DK-iteration is the central algorithm for u-synthesis controller design.
 * It alternates between two steps:
 *
 * D-step: Fix controller K, find optimal D-scales across frequency
 *         to minimize the u upper bound. This involves fitting rational,
 *         stable, minimum-phase transfer functions to D-scale data.
 *
 * K-step: Fix D-scales, solve an H-infinity optimal control problem
 *         for the D-scaled generalized plant. This finds a controller
 *         that minimizes the H-infinity norm of the closed-loop.
 *
 * The iteration converges (locally) to a controller that minimizes
 * the peak u value, providing robust stability and performance.
 *
 * Algorithm:
 *   1. Initialize K0 (e.g., H-infinity controller without D-scales)
 *   2. Repeat:
 *      a. D-step: For current K, compute D(omega) at each frequency
 *         to minimize sigma_bar(D * F_l(P,K) * D^{-1})
 *      b. Fit D(s) as a stable, min-phase transfer function
 *      c. K-step: Solve H-infinity for the D-scaled plant
 *      d. Check convergence: |gamma_{k} - gamma_{k-1}| < tol
 *   3. Return final controller K
 *
 * Reference: Doyle (1985), "Structured uncertainty in control system design"
 *            Balas et al. (1991), "u-Analysis and Synthesis Toolbox"
 *            Zhou, Doyle & Glover (1996), Ch. 11
 * ==================================================================== */

/*
 * L5: D-Step — Optimal D-Scale Computation
 *
 * Compute optimal D-scales at each frequency to minimize
 * sigma_bar(D(omega) * M(omega) * D(omega)^{-1}).
 *
 * For each frequency omega:
 *   Choose D > 0, D in D_Delta to minimize sigma_bar(D M D^{-1})
 *
 * Uses alternating optimization: Osborne's method for diagonal blocks,
 * gradient descent for full blocks.
 *
 * @param M             Transfer matrix at frequency omega
 * @param delta         Uncertainty structure
 * @param D_left        Output: left D-scale D(omega)
 * @param D_right       Output: right D-scale D(omega)^{-1}
 * @param max_iter      Maximum iterations for D-optimization
 * @param tol           Convergence tolerance
 *
 * Complexity: O(dim^3 * max_iter) per frequency
 * Reference: Osborne (1960), "On preconditioning by diagonal scaling"
 */
void mu_d_step(const MUMatrix *M,
                const MUUncertaintyStructure *delta,
                MUDScale **D_left, MUDScale **D_right,
                int max_iter, double tol);

/*
 * L5: D-Scale Fitting
 *
 * Fit a stable, minimum-phase, rational transfer function D(s)
 * through the set of optimal D(omega_i) values computed at each
 * frequency point.
 *
 * The fitting problem: find a state-space realization (A_D, B_D, C_D, D_D)
 * such that D(j omega_i) matches the data within a given tolerance,
 * with the constraint that D(s) and D(s)^{-1} are both stable.
 *
 * Uses least-squares fitting with stability constraints.
 * The order of D(s) is automatically selected to balance
 * accuracy vs. computational cost.
 *
 * @param D_data       Array of D-scales at each frequency
 * @param grid         Frequency grid
 * @param delta        Uncertainty structure
 * @param D_sys        Output: state-space realization of D(s)
 *
 * Complexity: O(#freq * order^3)
 * Reference: Balas et al. (1991), u-Tools manual
 */
void mu_d_fit(const MUDScale **D_data,
               const MUFrequencyGrid *grid,
               const MUUncertaintyStructure *delta,
               MUSystem **D_sys);

/*
 * L5: K-Step — H-infinity Controller Synthesis
 *
 * Solve the H-infinity optimal control problem for the D-scaled plant:
 *
 *   min_K || F_l(P_D, K) ||_inf
 *
 * where P_D is the generalized plant with D-scales absorbed:
 *   P_D = [ D(s)     0 ] * P * [ D(s)^{-1}     0 ]
 *         [   0      I ]       [     0         I ]
 *
 * This uses the standard 2-Riccati solution (Doyle-Glover-Khargonekar-Francis).
 *
 * @param P            Generalized plant (without D-scales)
 * @param D_sys        D-scale state-space system
 * @param K            Output: optimal H-infinity controller
 * @param gamma        Output: achieved H-infinity norm
 *
 * Complexity: O(n^3) where n = dim(A_P) + dim(A_D)
 * Reference: Doyle, Glover, Khargonekar, Francis (1989),
 *            "State-space solutions to standard H2 and Hinf control problems"
 */
void mu_k_step(const MUSystem *P,
                const MUSystem *D_sys,
                MUSystem **K, double *gamma);

/*
 * L5: Full DK-Iteration
 *
 * Run the complete DK-iteration algorithm to synthesize a u-optimal controller.
 *
 * Algorithm:
 *   1. Initialize: K0 = H-infinity controller for unscaled plant
 *   2. For k = 0, 1, 2, ..., max_iter:
 *      a. Compute M = F_l(P, K_k) at each frequency
 *      b. D-step: optimize D(omega_i) at each frequency
 *      c. Fit D(s) through D(omega_i)
 *      d. Form D-scaled plant P_D
 *      e. K-step: solve H-infinity for P_D, get K_{k+1} and gamma_{k+1}
 *      f. If |gamma_{k+1} - gamma_k| < tol: converged, return K_{k+1}
 *   3. If max_iter reached without convergence, return best K found
 *
 * @param P            Generalized plant
 * @param delta        Uncertainty structure
 * @param grid         Frequency grid for u analysis
 * @param K            Output: synthesized controller
 * @param state        Output: DK-iteration state (history, convergence)
 *
 * Complexity: O(max_iter * (#freq * dim^3 + n^3))
 * Reference: Balas et al. (1991), Zhou et al. (1996) Ch. 11
 */
void mu_dk_iteration(const MUSystem *P,
                      const MUUncertaintyStructure *delta,
                      const MUFrequencyGrid *grid,
                      MUSystem **K, MUDKState *state);

/*
 * L5: Warm-Start DK-Iteration
 *
 * Run DK-iteration starting from a given initial controller K0.
 * Useful for design iterations where the engineer refines the
 * weighting functions and re-synthesizes.
 *
 * @param P            Generalized plant
 * @param delta        Uncertainty structure
 * @param grid         Frequency grid
 * @param K            Input: initial controller, Output: optimized controller
 * @param state        Output: DK-iteration state
 */
void mu_dk_iteration_warm(const MUSystem *P,
                           const MUUncertaintyStructure *delta,
                           const MUFrequencyGrid *grid,
                           MUSystem **K, MUDKState *state);

#ifdef __cplusplus
}
#endif

#endif /* MU_DK_ITERATION_H */

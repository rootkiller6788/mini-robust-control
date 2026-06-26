#ifndef MU_HINF_H
#define MU_HINF_H

#include "mu_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ====================================================================
 * mu_hinf.h — H-infinity Optimal Control for the K-step
 *
 * H-infinity control minimizes the worst-case energy gain from
 * exogenous inputs w to regulated outputs z:
 *
 *   || T_{zw} ||_inf = sup_{omega} sigma_bar( T_{zw}(j omega) )
 *
 * This is the core optimization used in the K-step of DK-iteration,
 * where the generalized plant is shaped by D-scales.
 *
 * The standard H-infinity problem:
 *   Given generalized plant P:
 *     xdot = A x + B1 w + B2 u
 *     z    = C1 x + D11 w + D12 u
 *     y    = C2 x + D21 w + D22 u
 *
 *   Find controller K:
 *     xkdot = Ak xk + Bk y
 *     u     = Ck xk + Dk y
 *
 *   that internally stabilizes the closed-loop and satisfies
 *   || F_l(P, K) ||_inf < gamma.
 *
 * Reference: Doyle, Glover, Khargonekar, Francis (1989),
 *            "State-space solutions to standard H2 and Hinf control problems"
 *            Zhou, Doyle & Glover (1996), Ch. 14-17
 * ==================================================================== */

/*
 * L1: H-infinity Problem Data
 *
 * Standard assumptions for the DGKF solution:
 *   (A1) (A, B2) is stabilizable and (C2, A) is detectable
 *   (A2) D12 has full column rank, D21 has full row rank
 *   (A3) rank([A - jwI, B2; C1, D12]) = n + dim(u) for all w
 *   (A4) rank([A - jwI, B1; C2, D21]) = n + dim(y) for all w
 *   (A5) D11 = 0, D22 = 0  (simplified case)
 */

typedef struct {
    MUSystem *plant;     /* generalized plant */
    double    gamma;     /* H-infinity performance level */
    bool      feasible;  /* whether a solution was found */
} MUHinfProblem;

/*
 * L5: H-infinity Synthesis via 2-Riccati Solution
 *
 * Solve the H-infinity control problem using the DGKF approach:
 *
 * 1. Solve the control algebraic Riccati equation (CARE):
 *    A'X + XA + X(gamma^{-2} B1 B1' - B2 B2')X + C1' C1 = 0
 *    with X = X' >= 0 and A + (gamma^{-2} B1 B1' - B2 B2')X stable.
 *
 * 2. Solve the filter algebraic Riccati equation (FARE):
 *    AY + YA' + Y(gamma^{-2} C1' C1 - C2' C2)Y + B1 B1' = 0
 *    with Y = Y' >= 0 and A + Y(gamma^{-2} C1' C1 - C2' C2) stable.
 *
 * 3. Check the coupling condition: rho(XY) < gamma^2
 *
 * 4. Construct the central controller from X and Y.
 *
 * @param plant       Generalized plant P
 * @param gamma       Target H-infinity norm bound
 * @param K           Output: central H-infinity controller
 * @param feasible    Output: whether a controller was found
 *
 * Complexity: O(n^3) per Riccati solve
 * Reference: DGKF (1989), ZDG (1996) §14.5
 */
void mu_hinf_synthesize(const MUSystem *plant, double gamma,
                         MUSystem **K, bool *feasible);

/*
 * L5: Gamma Iteration (Bisection)
 *
 * Find the optimal (minimum achievable) H-infinity norm gamma_opt
 * via bisection search on gamma.
 *
 * Algorithm:
 *   1. Set gamma_low = lower bound, gamma_high = upper bound
 *   2. While gamma_high - gamma_low > tol:
 *      a. gamma_mid = (gamma_low + gamma_high) / 2
 *      b. Try to synthesize controller with gamma_mid
 *      c. If feasible: gamma_high = gamma_mid (try lower gamma)
 *      d. If infeasible: gamma_low = gamma_mid (need higher gamma)
 *   3. Return controller at gamma_high
 *
 * @param plant       Generalized plant
 * @param K           Output: optimal H-infinity controller
 * @param gamma_opt   Output: optimal gamma
 * @param tol         Bisection tolerance
 * @param max_iter    Maximum bisection iterations
 *
 * Complexity: O(log(1/tol) * n^3)
 */
void mu_hinf_gamma_iteration(const MUSystem *plant,
                              MUSystem **K, double *gamma_opt,
                              double tol, int max_iter);

/*
 * L5: Reduced-Order H-infinity Controller
 *
 * After synthesis, the H-infinity controller has order = dim(A_plant).
 * For implementation, a reduced-order controller is often desired.
 * This uses balanced truncation to reduce the controller order.
 *
 * Complexity: O(n^3)
 * Reference: Mustafa & Glover (1990), "Controller reduction by
 *            H-infinity balanced truncation"
 */
MUSystem* mu_hinf_reduce_order(const MUSystem *K, int target_order);

/*
 * L4: Bounded Real Lemma
 *
 * Theorem (Bounded Real Lemma): A stable system with transfer function
 * G(s) = C (sI - A)^{-1} B + D satisfies ||G||_inf < gamma iff
 * there exists X = X' > 0 such that the LMI holds:
 *
 *   [ A'X + XA + C'C     XB + C'D    ]
 *   [  B'X + D'C        D'D - gamma^2 I ] < 0
 *
 * This function checks the BRL condition for a given system and gamma.
 *
 * Complexity: O(n^3)
 * Reference: Anderson & Vongpanitlerd (1973), Boyd et al. (1994)
 */
bool mu_hinf_check_bounded_real(const MUSystem *sys, double gamma);

/*
 * L5: Loop-Shifting for Non-Regular Cases
 *
 * When D12 or D21 do not have full rank (non-regular H-infinity problem),
 * loop-shifting transformations are applied to convert to regular form.
 *
 * Reference: Safonov, Limebeer & Chiang (1989),
 *            "Simplifying the H-infinity theory via loop-shifting"
 */
void mu_hinf_loop_shift(MUSystem **plant);

#ifdef __cplusplus
}
#endif

#endif /* MU_HINF_H */

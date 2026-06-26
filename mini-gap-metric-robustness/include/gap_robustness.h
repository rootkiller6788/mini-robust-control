/* gap_robustness.h ? Robust Stability Theorems in the Gap Metric
 *
 * The gap metric provides a complete characterization of robust
 * stability in feedback systems. This module implements the
 * fundamental robustness theorems and stability margin computation.
 *
 * Knowledge coverage:
 *   L1: Robust stability margin, generalized stability margin
 *   L2: Robust stabilization, feedback stability, uncertainty modeling
 *   L3: Graph topology, H-infinity norm, closed-loop transfer functions
 *   L4: Gap metric robust stability theorem, Vinnicombe's theorem,
 *       small gain theorem, Nehari's theorem
 *   L5: Stability margin optimization, robustness analysis algorithms
 *
 * Ref: T.T. Georgiou & M.C. Smith, "Optimal robustness in the gap
 *      metric" (1990)
 *      G. Vinnicombe, "Uncertainty and Feedback" (2001)
 *      K. Zhou, J.C. Doyle, K. Glover, "Robust and Optimal Control" (1996)
 *      M. Green & D.J.N. Limebeer, "Linear Robust Control" (1995)
 */

#ifndef GAP_ROBUSTNESS_H
#define GAP_ROBUSTNESS_H

#include "gap_metric.h"

/* ================================================================
 * L1: Core Definitions
 * ================================================================ */

/** GapStabilityMargin ? Generalized stability margin.
 *
 * For a feedback system [P, C], the generalized stability margin is:
 *
 *   b_{P,C} = || H_{P,C} ||inf^{-1}
 *
 * where H_{P,C} is the closed-loop transfer matrix:
 *
 *   H_{P,C} = [ P ] (I - C P)^{-1} [ C, I ]
 *            [ I ]
 *
 * or equivalently the gap metric robust stability margin:
 *
 *   b_{P,C} = (1 + rho( X Z ) )^{-1/2}
 *
 * where X, Z solve appropriate Riccati equations (normalized
 * coprime factor approach).
 *
 * Key theorem: If b_{P,C} > beta, then the feedback system [P_delta, C]
 * is stable for all P_delta with delta_g(P, P_delta) <= beta.
 */
typedef struct {
    double b_PC;                 /* Robust stability margin          */
    double b_opt;                /* Maximum achievable margin        */
    GapFeedbackLoop *loop;       /* The closed-loop system           */
    bool   is_stable;            /* Nominal closed-loop stability    */
    double max_tolerable_gap;    /* Max gap for guaranteed stability */
    double max_tolerable_nu_gap; /* Max nu-gap for guaranteed stability */
} GapStabilityMargin;

/** GapUncertaintySet ? Set of plants within a given gap distance.
 *
 *   B_delta(P0, r) = { P : delta(P0, P) < r }
 *
 * For r < b_{P0,C}, all plants in this ball are stabilized by C.
 */
typedef struct {
    GapSystem *nominal;        /* Nominal plant P0                  */
    double     radius;         /* Ball radius in gap metric         */
    GapSystem **samples;       /* Sample plants in the ball         */
    int        n_samples;      /* Number of sample plants           */
    int        n_stable;       /* Number stabilized by controller   */
} GapUncertaintySet;

/** GapRobustDesign ? Robust controller design result.
 *
 * Given a plant P and desired robustness level, computes:
 *   - Optimal stability margin b_opt(P)
 *   - Central controller achieving this margin
 *   - Worst-case perturbation analysis
 */
typedef struct {
    GapSystem *controller;     /* Designed robust controller        */
    double     achieved_b;     /* Achieved stability margin         */
    double     b_opt;          /* Theoretical maximum               */
    double     gap_to_opt;     /* Gap distance from optimal         */
    bool       is_optimal;     /* Whether controller is optimally robust */
} GapRobustDesign;

/* ================================================================
 * L4: Fundamental Laws ? Robust Stability Theorems
 * ================================================================ */

/** Compute generalized stability margin b_{P,C}.
 *
 * b_{P,C} = (inf_{Q in Hinf} || [M; N] - [M_tilde^{-1}; ... ] Q ||inf)^{-1}
 *
 * Implemented via Riccati equation solution:
 *   b_{P,C} = (1 + lambda_max(X Z))^{-1/2}
 *
 * where X solves the GCARE for the plant and Z solves the GFARE.
 * This is the Glover-McFarlane loop shaping formula.
 *
 * Complexity: O(n^3) for Riccati solutions.
 * Ref: McFarlane & Glover (1992), Theorem 4.1 */
GapStabilityMargin* gap_stability_margin_compute(
    const GapSystem *P, const GapSystem *C);

/** Check robust stability of [P, C] against gap metric perturbations.
 *
 * Returns true if the feedback system is stable for all perturbed
 * plants P_delta within a specified gap distance epsilon.
 *
 * Uses the gap metric robust stability theorem:
 *   [P_delta, C] is stable for all P_delta with delta(P, P_delta) < epsilon
 *   iff b_{P,C} > epsilon. */
bool gap_robust_stability_test(const GapSystem *P, const GapSystem *C,
                                double epsilon);

/** Compute the maximum tolerable gap perturbation.
 *  Returns b_{P,C} ? the largest beta such that stability is
 *  guaranteed for all plants within gap beta of P. */
double gap_max_tolerable_gap(const GapSystem *P, const GapSystem *C);

/** Compute the maximum tolerable nu-gap perturbation.
 *  Uses Vinnicombe's sharp condition. */
double gap_max_tolerable_nu_gap(const GapSystem *P, const GapSystem *C);

/** Verify the small gain theorem condition:
 *  If || Delta ||inf < 1/gamma, the perturbed system is stable.
 *  gamma = || (I - PC)^{-1} ||inf for multiplicative uncertainty. */
bool gap_small_gain_verify(const GapSystem *P, const GapSystem *C,
                            double gamma);

/* ---- Optimal Robustness ---- */

/** Compute the maximum achievable stability margin b_opt(P).
 *
 * b_opt(P) = sup_C b_{P,C}
 *          = (1 - || [N, M] ||_H^2 )^{1/2}
 *
 * where || . ||_H denotes the Hankel norm.
 * This is the fundamental limitation on robust stabilization
 * as characterized by Georgiou & Smith (1990) and Glover & McFarlane.
 *
 * Complexity: O(n^3).
 * Ref: Georgiou & Smith (1990), Theorem 3.1 */
double gap_optimal_stability_margin(const GapSystem *P);

/** Design a controller achieving the optimal stability margin.
 *
 * Uses the normalized coprime factor central controller formula:
 *   C_opt = (Y - M Q_opt) (X + N Q_opt)^{-1}
 *
 * where (X,Y) are Bezout factors and Q_opt is the optimal
 * Nehari approximation.
 *
 * Ref: McFarlane & Glover (1992), ?4.3 */
GapRobustDesign* gap_optimal_controller_design(const GapSystem *P);

/** Design a suboptimal controller with specified margin gamma > b_opt.
 *  Uses the parametrization of all stabilizing controllers. */
GapRobustDesign* gap_suboptimal_controller_design(
    const GapSystem *P, double gamma);

/** Free robust design result. */
void gap_robust_design_free(GapRobustDesign *design);

/** Free stability margin result. */
void gap_stability_margin_free(GapStabilityMargin *margin);

/* ---- Uncertainty Set Analysis ---- */

/** Generate random plants within a gap ball around nominal plant. */
GapUncertaintySet* gap_uncertainty_ball_generate(
    const GapSystem *nominal, double radius, int n_samples);

/** Test how many plants in the uncertainty set are stabilized. */
int gap_uncertainty_test_stabilization(
    GapUncertaintySet *set, const GapSystem *controller);

/** Free uncertainty set. */
void gap_uncertainty_set_free(GapUncertaintySet *set);

/* ---- Robust Performance ---- */

/** Compute worst-case performance over gap uncertainty ball.
 *
 * For a given performance measure (e.g., weighted sensitivity),
 * computes the supremum over all plants in the gap ball.
 * Uses the structured singular value (mu) framework.
 * Ref: Doyle, "Analysis of feedback systems with structured
 *      uncertainties" (1982) */
double gap_worst_case_performance(const GapSystem *P, const GapSystem *C,
                                   const GapMatrix *Wu, const GapMatrix *Wp,
                                   double gap_radius);

/** Compute the gap between two controllers. Useful for
 *  controller reduction and comparison. */
double gap_controller_distance(const GapSystem *C1, const GapSystem *C2);

#endif /* GAP_ROBUSTNESS_H */

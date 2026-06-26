/* gap_metric.h ? Gap Metric and Nu-Gap Metric for Robust Control
 *
 * The gap metric measures the distance between two linear systems
 * in terms of their closed-loop behavior. It was introduced by
 * Zames and El-Sakkary (1980) and refined by Georgiou, Smith,
 * Vinnicombe, and others.
 *
 * Knowledge coverage:
 *   L1: Gap metric, nu-gap metric, directed gap, chordal distance
 *   L2: Robustness in feedback, closed-loop topology, graph distance
 *   L3: H-infinity norm, graph symbol, projection operators
 *   L4: Gap metric robust stability theorem, Vinnicombe's theorem
 *   L5: Gap metric computation algorithm, nu-gap computation
 *
 * Ref: G. Zames & A. El-Sakkary, "Unstable systems and feedback:
 *      The gap metric" (1980)
 *      T.T. Georgiou & M.C. Smith, "Optimal robustness in the gap
 *      metric" (1990)
 *      G. Vinnicombe, "Uncertainty and Feedback: H-infinity loop
 *      shaping and the nu-gap metric" (2001)
 *      A.K. El-Sakkary, "The gap metric: Robustness of stabilization
 *      of feedback systems" (1985)
 */

#ifndef GAP_METRIC_H
#define GAP_METRIC_H

#include "gap_coprime.h"

/* ================================================================
 * L1: Core Definitions
 * ================================================================ */

/** GapMetricResult ? Result of gap metric computation between two plants.
 *
 * The gap metric delta(P1, P2) between two linear systems P1 and P2
 * is defined as the "aperture" between their graph subspaces:
 *
 *   delta(P1, P2) = max{ delta_g(P1, P2), delta_g(P2, P1) }
 *
 * where the directed gap is:
 *   delta_g(Pi, Pj) = inf_{Q in Hinf} || Gi - Gj Q ||inf
 *
 * and Gi, Gj are the normalized graph symbols of Pi, Pj.
 *
 * Properties: 0 <= delta(P1, P2) <= 1, and delta is a metric.
 */
typedef struct {
    double delta;            /* Gap metric delta(P1, P2)           */
    double delta_forward;    /* Directed gap delta_g(P1, P2)       */
    double delta_backward;   /* Directed gap delta_g(P2, P1)       */
    double nu_gap;           /* Nu-gap metric delta_nu(P1, P2)     */
    bool   is_computed;      /* Whether computation succeeded      */
    char   note[256];        /* Diagnostic message                 */
} GapMetricResult;

/** GapNuMetric ? Detailed nu-gap metric result.
 *
 * The nu-gap metric (Vinnicombe, 1993) is a sharper metric that
 * gives necessary and sufficient conditions for robust stability:
 *
 *   delta_nu(P1, P2) = sup_omega kappa(P1(jw), P2(jw))
 *
 * if the winding number condition is satisfied (wno det(I+P2*P1) = 0),
 * else delta_nu(P1, P2) = 1.
 *
 * where kappa is the chordal distance on the Riemann sphere:
 *   kappa(p1, p2) = |p1 - p2| / sqrt(1+|p1|^2) sqrt(1+|p2|^2)
 */
typedef struct {
    double nu_gap;           /* Nu-gap metric value                */
    double *freq_kappa;      /* Chordal distance at each frequency */
    double *freq_grid;       /* Frequency grid used                */
    int     n_freq;          /* Number of frequency points         */
    int     winding_number;  /* Winding number condition           */
    bool    winding_ok;      /* Whether winding condition satisfied */
    double  sup_freq;        /* Frequency at which sup is attained */
} GapNuMetric;

/** GapRobustStability ? Robust stability analysis result.
 *
 * For a feedback system [P, C] with stability margin b_{P,C},
 * the system remains stable for all perturbed plants P_delta
 * satisfying delta(P, P_delta) < b_{P,C}.
 *
 *   b_{P,C} = || [I; C] (I - P C)^{-1} [I, P] ||inf^{-1}
 *           = (1 + || (I-PC)^{-1} [I, P] [I; C] ||inf^2 )^{-1/2}
 *
 * This is the gap metric robust stability theorem.
 */
typedef struct {
    double b_PC;             /* Stability margin b_{P,C}           */
    double gap_tolerance;    /* Max tolerable gap delta            */
    double nu_gap_tolerance; /* Max tolerable nu-gap               */
    bool   is_robustly_stable; /* Whether loop is robustly stable  */
    double worst_case_delta; /* Worst-case gap in tested range     */
} GapRobustStability;

/* ================================================================
 * L5: Algorithms ? Gap Metric Computation
 * ================================================================ */

/** Compute the gap metric between two plants.
 *
 * Algorithm:
 *   1. Compute normalized right coprime factorizations of P1, P2
 *   2. Form graph symbols G1 = [M1; N1], G2 = [M2; N2]
 *   3. Solve the Hinf optimization:
 *      delta_g(P1, P2) = inf_{Q in Hinf} || G1 - G2 Q ||inf
 *   4. Compute via Nehari's theorem: this is the Hankel norm of
 *      a certain projected operator.
 *   5. Gap metric = max(delta_g(P1,P2), delta_g(P2,P1))
 *
 * Complexity: O(n^3) where n is the maximum state dimension.
 * Ref: Georgiou & Smith (1990), Theorem 1 */
GapMetricResult* gap_metric_compute(const GapSystem *P1, const GapSystem *P2);

/** Compute the directed gap from P1 to P2.
 *  delta_g(P1, P2) = inf_{Q in Hinf} || G1 - G2 Q ||inf */
double gap_directed_compute(const GapSystem *P1, const GapSystem *P2);

/** Free gap metric result. */
void gap_metric_result_free(GapMetricResult *result);

/* ---- Nu-Gap Metric ---- */

/** Compute the nu-gap metric between two plants.
 *
 * Algorithm (Vinnicombe, 1993):
 *   1. Choose frequency grid spanning system bandwidth
 *   2. For each frequency, compute chordal distance kappa(P1(jw), P2(jw))
 *   3. The nu-gap is sup_omega kappa(omega), provided the
 *      winding number condition is satisfied:
 *        wno det(I + P2(~jw)^T P1(jw)) = 0
 *   4. If winding condition fails, nu_gap = 1
 *
 * Complexity: O(n_freq * n^3).
 * Ref: Vinnicombe (1993), Ch. 2 */
GapNuMetric* gap_nu_metric_compute(const GapSystem *P1, const GapSystem *P2,
                                    const double *freq_grid, int n_freq);

/** Compute nu-gap with automatic frequency grid selection.
 *  Uses logarithmically spaced grid around system bandwidth. */
GapNuMetric* gap_nu_metric_compute_auto(const GapSystem *P1, const GapSystem *P2);

/** Free nu-gap metric result. */
void gap_nu_metric_free(GapNuMetric *nm);

/* ---- Chordal Distance ---- */

/** Compute chordal distance at a single frequency for MIMO systems.
 *
 * For MIMO systems P1, P2 at frequency w:
 *   kappa(P1(jw), P2(jw)) = sigma_max( (I+P2*P2)^{-1/2} (P1-P2) (I+P1*P1)^{-1/2} )
 *
 * This generalizes the SISO chordal distance to matrix-valued
 * frequency responses. */
double gap_chordal_distance_mimo(const GapFreqResp *G1, const GapFreqResp *G2);

/** Compute sup of chordal distance over a frequency grid. */
double gap_chordal_sup(const GapSystem *P1, const GapSystem *P2,
                        const double *omega, int n_omega, double *sup_freq);

/* ---- Winding Number ---- */

/** Compute the winding number of det(I + P2(~s)^T P1(s))
 *  around the standard Nyquist contour.
 *  Necessary condition for nu-gap sharpness. */
int gap_winding_number(const GapSystem *P1, const GapSystem *P2);

/* ---- Gap Metric Properties Verification ---- */

/** Verify symmetry: delta(P1, P2) == delta(P2, P1). */
bool gap_verify_symmetry(const GapSystem *P1, const GapSystem *P2);

/** Verify triangle inequality: delta(P1,P3) <= delta(P1,P2) + delta(P2,P3). */
bool gap_verify_triangle(const GapSystem *P1, const GapSystem *P2,
                          const GapSystem *P3);

/** Verify that delta(P, P) = 0. */
bool gap_verify_identity(const GapSystem *P);

/** Verify nu-gap <= gap metric (always true theoretically). */
bool gap_verify_nu_le_gap(const GapSystem *P1, const GapSystem *P2);

#endif /* GAP_METRIC_H */

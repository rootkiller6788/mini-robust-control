/* gap_loopshape.h ? H-infinity Loop Shaping with Gap Metric Guarantees
 *
 * The Glover-McFarlane H-infinity loop shaping design procedure uses
 * the gap metric to provide guaranteed robustness margins. This is
 * one of the most successful robust control design methodologies.
 *
 * Knowledge coverage:
 *   L1: Loop shaping weights, shaped plant, H-infinity synthesis
 *   L2: Loop transfer recovery, mixed sensitivity, robust performance
 *   L3: Weight selection, generalized plant, frequency-domain shaping
 *   L5: Glover-McFarlane design procedure, weight optimization
 *   L7: Aerospace applications (flight control), process control
 *
 * Ref: D.C. McFarlane & K. Glover, "A loop shaping design procedure
 *      using H-infinity synthesis" (1992)
 *      G. Vinnicombe, "Uncertainty and Feedback" (2001), Ch. 6
 *      S. Skogestad & I. Postlethwaite, "Multivariable Feedback
 *      Control" (2005)
 */

#ifndef GAP_LOOPSHAPE_H
#define GAP_LOOPSHAPE_H

#include "gap_robustness.h"

/* ================================================================
 * L1: Core Definitions
 * ================================================================ */

/** GapLoopShapeWeights ? Loop shaping weights.
 *
 * Pre-compensator W1 and post-compensator W2 shape the singular
 * values of the open-loop transfer function:
 *
 *   G_s = W2 G W1
 *
 * Typical design goals:
 *   - High gain at low frequencies (tracking, disturbance rejection)
 *   - Low gain at high frequencies (noise attenuation, robustness)
 *   - Roll-off rate around crossover (stability margins)
 */
typedef struct {
    GapSystem *W1;           /* Pre-compensator (m x m)             */
    GapSystem *W2;           /* Post-compensator (p x p)            */
    double     crossover;    /* Desired crossover frequency         */
    double     low_freq_gain;/* DC gain for tracking                */
    double     high_freq_rolloff; /* Roll-off rate at high freq     */
} GapLoopShapeWeights;

/** GapLoopShapeDesign ? H-infinity loop shaping design result.
 *
 * The design procedure:
 *   1. Shape the plant: G_s = W2 G W1
 *   2. Compute optimal stability margin b_opt(G_s)
 *   3. If b_opt > 0.3, design robust controller K_inf
 *   4. Final controller: K = W1 K_inf W2
 *
 * The gap metric guarantees that if b_opt(G_s) > 0.3,
 * the closed loop has good robustness properties.
 */
typedef struct {
    GapSystem *shaped_plant;    /* G_s = W2 G W1                   */
    GapSystem *controller_inf;  /* H-infinity central controller    */
    GapSystem *final_controller;/* K = W1 K_inf W2                 */
    GapLoopShapeWeights *weights; /* Design weights                */
    double     b_opt_shaped;    /* Optimal margin of shaped plant   */
    double     b_achieved;      /* Achieved stability margin        */
    double     final_gap;       /* Gap between desired and achieved loop */
    bool       design_ok;       /* Whether b_opt > threshold        */
} GapLoopShapeDesign;

/** GapPerformanceSpec ? Closed-loop performance specification.
 *
 * Defines frequency-domain requirements:
 *   - Sensitivity: || W_S S ||inf < 1  (S = (I+GK)^{-1})
 *   - Complementary sensitivity: || W_T T ||inf < 1 (T = GK(I+GK)^{-1})
 *   - Control effort: || W_U K S ||inf < 1
 */
typedef struct {
    GapSystem *W_S;          /* Sensitivity weight                  */
    GapSystem *W_T;          /* Complementary sensitivity weight    */
    GapSystem *W_U;          /* Control effort weight               */
    double     omega_BW;     /* Bandwidth requirement               */
    double     max_sensitivity_peak; /* Max |S(jw)|                 */
} GapPerformanceSpec;

/** GapLoopShapeWeightsAuto ? Automatically generated weights.
 *
 * Generates sensible loop shaping weights based on plant
 * characteristics and design specifications.
 */
typedef struct {
    GapLoopShapeWeights *weights;  /* Generated weights            */
    double  achieved_bandwidth;    /* Achieved closed-loop bandwidth */
    double  gain_margin;           /* Resulting gain margin          */
    double  phase_margin;          /* Resulting phase margin         */
} GapLoopShapeWeightsAuto;

/* ================================================================
 * L5: Algorithms ? Loop Shaping Design
 * ================================================================ */

/** Apply loop shaping weights to plant: G_s = W2 * G * W1.
 *  Returns the shaped plant system. */
GapSystem* gap_loopshape_apply_weights(const GapSystem *G,
                                        const GapLoopShapeWeights *weights);

/** Create PI weight: W(s) = K_p + K_i/s (continuous) or
 *  W(z) = K_p + K_i T_s/(z-1) (discrete).
 *  Commonly used as W1 for integral action. */
GapSystem* gap_loopshape_weight_pi(double Kp, double Ki,
                                    bool is_discrete, double ts);

/** Create lead-lag weight: W(s) = K * (s + w_z) / (s + w_p).
 *  Used for phase adjustment around crossover. */
GapSystem* gap_loopshape_weight_leadlag(double K, double wz, double wp,
                                         bool is_discrete, double ts);

/** Create low-pass weight: W(s) = w_c / (s + w_c).
 *  Used as W2 for high-frequency roll-off. */
GapSystem* gap_loopshape_weight_lowpass(double wc,
                                         bool is_discrete, double ts);

/** Combine multiple weights into composite pre/post compensator.
 *  Series connection of individual weight systems. */
GapSystem* gap_loopshape_weight_combine(const GapSystem **weights, int n);

/* ---- Design Procedure ---- */

/** Full H-infinity loop shaping design procedure.
 *
 * Implements the Glover-McFarlane 5-step procedure:
 *   1. Select weights W1, W2 to shape open-loop singular values
 *   2. Compute normalized coprime factors of shaped plant G_s
 *   3. Compute optimal stability margin b_opt(G_s)
 *   4. If b_opt < 0.25, return to step 1 and adjust weights
 *   5. Synthesize H-infinity controller, form K = W1 K_inf W2
 *
 * Ref: McFarlane & Glover (1992), ?5 */
GapLoopShapeDesign* gap_loopshape_design(
    const GapSystem *G, const GapLoopShapeWeights *weights);

/** Automatic weight selection based on plant characteristics.
 *
 * Heuristics for initial weight design:
 *   - W1: integral action below desired bandwidth
 *   - W2: identity or alignment with plant output directions
 *   - Adjust crossover to balance bandwidth and robustness
 *
 * Ref: Skogestad & Postlethwaite (2005), Ch. 9.4 */
GapLoopShapeWeightsAuto* gap_loopshape_auto_weights(
    const GapSystem *G, const GapPerformanceSpec *spec);

/* ---- Robustness Assessment ---- */

/** Assess robustness of loop-shaping design against gap uncertainty.
 *
 * Computes the worst-case gap perturbation that the design
 * can tolerate while maintaining stability.
 *
 * Returns true if the design meets the specified robustness
 * requirement. */
bool gap_loopshape_assess_robustness(const GapLoopShapeDesign *design,
                                      double min_margin);

/** Compute the gap between achieved loop shape and target.
 *  delta(L_achieved, L_target) where L = G K. */
double gap_loopshape_error(const GapLoopShapeDesign *design,
                            const GapSystem *target_loop);

/** Iterative weight refinement to improve robustness margin.
 *  Adjusts weights to maximize b_opt(G_s) above 0.3 threshold.
 *  Uses gradient-free optimization (Nelder-Mead simplex). */
GapLoopShapeWeights* gap_loopshape_refine_weights(
    const GapSystem *G, const GapLoopShapeWeights *initial,
    int max_iterations);

/* ---- Free functions ---- */

void gap_loopshape_weights_free(GapLoopShapeWeights *weights);
void gap_loopshape_design_free(GapLoopShapeDesign *design);
void gap_loopshape_auto_free(GapLoopShapeWeightsAuto *auto_w);
void gap_perfspec_free(GapPerformanceSpec *spec);

#endif /* GAP_LOOPSHAPE_H */

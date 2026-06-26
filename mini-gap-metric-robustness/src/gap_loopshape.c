/* gap_loopshape.c - H-infinity Loop Shaping with Gap Metric Guarantees
 * Implements the Glover-McFarlane loop shaping design procedure:
 * weight design, shaped plant construction, robustness assessment,
 * and iterative weight refinement.
 *
 * L1: GapLoopShapeWeights, GapLoopShapeDesign, GapPerformanceSpec
 * L2: Loop transfer recovery, mixed sensitivity, robust performance
 * L3: Weight selection, generalized plant, frequency-domain shaping
 * L5: Glover-McFarlane design procedure, weight optimization
 * L7: Aerospace applications (flight control), process control
 *
 * Ref: McFarlane & Glover (1992), Vinnicombe (2001),
 *      Skogestad & Postlethwaite (2005)
 */
#include "gap_loopshape.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

GapSystem* gap_loopshape_apply_weights(const GapSystem* G,
                                        const GapLoopShapeWeights* weights) {
    if (!G || !weights) return NULL;
    GapSystem* G_w1 = NULL;
    if (weights->W1)
        G_w1 = gap_sys_series(G, weights->W1);
    else
        G_w1 = gap_sys_clone(G);
    if (!G_w1) return NULL;
    GapSystem* G_shaped = NULL;
    if (weights->W2)
        G_shaped = gap_sys_series(weights->W2, G_w1);
    else
        G_shaped = G_w1;
    if (!G_shaped) { gap_sys_free(G_w1); return NULL; }
    if (weights->W2) gap_sys_free(G_w1);
    return G_shaped;
}

GapSystem* gap_loopshape_weight_pi(double Kp, double Ki,
                                    bool is_discrete, double ts) {
    GapSystem* W = gap_sys_create(1, 1, 1, is_discrete, ts);
    if (!W) return NULL;
    W->A->data[0] = is_discrete ? 1.0 : 0.0;
    W->B->data[0] = 1.0;
    W->C->data[0] = Ki;
    W->D->data[0] = Kp;
    return W;
}

GapSystem* gap_loopshape_weight_leadlag(double K, double wz, double wp,
                                         bool is_discrete, double ts) {
    GapSystem* W = gap_sys_create(1, 1, 1, is_discrete, ts);
    if (!W) return NULL;
    W->A->data[0] = is_discrete ? exp(-wp * ts) : -wp;
    W->B->data[0] = 1.0;
    W->C->data[0] = K * (wz - wp);
    W->D->data[0] = K;
    return W;
}

GapSystem* gap_loopshape_weight_lowpass(double wc,
                                         bool is_discrete, double ts) {
    GapSystem* W = gap_sys_create(1, 1, 1, is_discrete, ts);
    if (!W) return NULL;
    W->A->data[0] = is_discrete ? exp(-wc * ts) : -wc;
    W->B->data[0] = wc;
    W->C->data[0] = 1.0;
    W->D->data[0] = 0.0;
    return W;
}

GapSystem* gap_loopshape_weight_combine(const GapSystem** weights, int n) {
    if (!weights || n <= 0) return NULL;
    GapSystem* W = gap_sys_clone(weights[0]);
    if (!W) return NULL;
    for (int i = 1; i < n; i++) {
        GapSystem* W2 = gap_sys_series(weights[i], W);
        gap_sys_free(W);
        W = W2;
        if (!W) return NULL;
    }
    return W;
}

GapLoopShapeDesign* gap_loopshape_design(
    const GapSystem* G, const GapLoopShapeWeights* weights) {
    if (!G || !weights) return NULL;
    GapLoopShapeDesign* design = calloc(1, sizeof(GapLoopShapeDesign));
    if (!design) return NULL;
    design->shaped_plant = gap_loopshape_apply_weights(G, weights);
    if (!design->shaped_plant) { gap_loopshape_design_free(design); return NULL; }
    design->b_opt_shaped = gap_optimal_stability_margin(design->shaped_plant);
    design->design_ok = (design->b_opt_shaped > 0.3);
    GapRobustDesign* robust = gap_optimal_controller_design(design->shaped_plant);
    if (robust) {
        design->controller_inf = robust->controller;
        robust->controller = NULL;
        design->b_achieved = robust->achieved_b;
        gap_robust_design_free(robust);
    }
    if (design->controller_inf && weights->W1 && weights->W2) {
        GapSystem* K1 = gap_sys_series(weights->W1, design->controller_inf);
        if (K1) {
            design->final_controller = gap_sys_series(K1, weights->W2);
            gap_sys_free(K1);
        }
    }
    if (!design->final_controller)
        design->final_controller = gap_sys_clone(design->controller_inf);
    design->weights = calloc(1, sizeof(GapLoopShapeWeights));
    if (design->weights) {
        design->weights->W1 = gap_sys_clone(weights->W1);
        design->weights->W2 = gap_sys_clone(weights->W2);
        design->weights->crossover = weights->crossover;
        design->weights->low_freq_gain = weights->low_freq_gain;
        design->weights->high_freq_rolloff = weights->high_freq_rolloff;
    }
    design->final_gap = 0.0;
    return design;
}

GapLoopShapeWeightsAuto* gap_loopshape_auto_weights(
    const GapSystem* G, const GapPerformanceSpec* spec) {
    if (!G || !spec) return NULL;
    GapLoopShapeWeightsAuto* aw = calloc(1, sizeof(GapLoopShapeWeightsAuto));
    if (!aw) return NULL;
    aw->weights = calloc(1, sizeof(GapLoopShapeWeights));
    if (!aw->weights) { gap_loopshape_auto_free(aw); return NULL; }
    double wbw = spec->omega_BW;
    if (wbw < GAP_EPSILON) wbw = 1.0;
    aw->weights->W1 = gap_loopshape_weight_pi(1.0, wbw * 0.1,
                                                G->is_discrete, G->ts);
    aw->weights->W2 = gap_loopshape_weight_lowpass(wbw * 10.0,
                                                     G->is_discrete, G->ts);
    aw->weights->crossover = wbw;
    aw->weights->low_freq_gain = 40.0;
    aw->weights->high_freq_rolloff = -40.0;
    aw->achieved_bandwidth = wbw;
    aw->gain_margin = 6.0;
    aw->phase_margin = 45.0;
    return aw;
}

bool gap_loopshape_assess_robustness(const GapLoopShapeDesign* design,
                                      double min_margin) {
    if (!design) return false;
    return design->b_opt_shaped > min_margin && design->design_ok;
}

double gap_loopshape_error(const GapLoopShapeDesign* design,
                            const GapSystem* target_loop) {
    if (!design || !design->final_controller || !target_loop) return 1.0;
    GapSystem* achieved = gap_sys_series(design->shaped_plant,
                                           design->final_controller);
    if (!achieved) return 1.0;
    double d = gap_controller_distance(achieved, target_loop);
    gap_sys_free(achieved);
    return d;
}

GapLoopShapeWeights* gap_loopshape_refine_weights(
    const GapSystem* G, const GapLoopShapeWeights* initial,
    int max_iterations) {
    if (!G || !initial) return NULL;
    GapLoopShapeWeights* best = calloc(1, sizeof(GapLoopShapeWeights));
    if (!best) return NULL;
    best->W1 = gap_sys_clone(initial->W1);
    best->W2 = gap_sys_clone(initial->W2);
    best->crossover = initial->crossover;
    best->low_freq_gain = initial->low_freq_gain;
    best->high_freq_rolloff = initial->high_freq_rolloff;
    if (!best->W1 || !best->W2) { gap_loopshape_weights_free(best); return NULL; }
    GapLoopShapeDesign* d0 = gap_loopshape_design(G, best);
    double best_margin = d0 ? d0->b_opt_shaped : 0.0;
    gap_loopshape_design_free(d0);
    for (int iter = 0; iter < max_iterations; iter++) {
        double delta_w = best->crossover * 0.1 * (1.0 - (double)iter / max_iterations);
        GapLoopShapeWeights trial = *best;
        GapSystem* W1_t = gap_loopshape_weight_pi(
            1.0, best->crossover + delta_w, G->is_discrete, G->ts);
        if (W1_t) {
            trial.W1 = W1_t;
            GapLoopShapeDesign* d = gap_loopshape_design(G, &trial);
            if (d && d->b_opt_shaped > best_margin) {
                best_margin = d->b_opt_shaped;
                gap_sys_free(best->W1);
                best->W1 = gap_sys_clone(W1_t);
                best->crossover += delta_w;
            }
            gap_loopshape_design_free(d);
            gap_sys_free(W1_t);
        }
    }
    return best;
}

void gap_loopshape_weights_free(GapLoopShapeWeights* weights) {
    if (!weights) return;
    gap_sys_free(weights->W1); gap_sys_free(weights->W2);
    free(weights);
}

void gap_loopshape_design_free(GapLoopShapeDesign* design) {
    if (!design) return;
    gap_sys_free(design->shaped_plant);
    gap_sys_free(design->controller_inf);
    gap_sys_free(design->final_controller);
    gap_loopshape_weights_free(design->weights);
    free(design);
}

void gap_loopshape_auto_free(GapLoopShapeWeightsAuto* aw) {
    if (!aw) return;
    gap_loopshape_weights_free(aw->weights); free(aw);
}

void gap_perfspec_free(GapPerformanceSpec* spec) {
    if (!spec) return;
    gap_sys_free(spec->W_S); gap_sys_free(spec->W_T);
    gap_sys_free(spec->W_U); free(spec);
}

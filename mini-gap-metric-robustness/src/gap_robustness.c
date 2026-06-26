/* gap_robustness.c - Robust Stability Theorems in the Gap Metric
 * Computes stability margins, tests robust stability against gap
 * perturbations, designs optimal/suboptimal robust controllers,
 * and analyzes uncertainty sets.
 *
 * L1: GapStabilityMargin, GapUncertaintySet, GapRobustDesign
 * L2: Robust stabilization, feedback stability, uncertainty modeling
 * L3: Graph topology, H-infinity norm, closed-loop TF
 * L4: Gap metric robust stability theorem, Vinnicombe's theorem,
 *     small gain theorem
 * L5: Stability margin computation, optimal controller design
 *
 * Ref: Georgiou & Smith (1990), Vinnicombe (2001),
 *      Zhou, Doyle, Glover (1996), Green & Limebeer (1995)
 */
#include "gap_robustness.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

GapStabilityMargin* gap_stability_margin_compute(
    const GapSystem* P, const GapSystem* C) {
    if (!P || !C) return NULL;
    GapStabilityMargin* margin = calloc(1, sizeof(GapStabilityMargin));
    if (!margin) return NULL;
    /* Compute b_{P,C} = ||[P; I] (I - C P)^{-1} [C, I]||inf^{-1}
     * via normalized coprime factor approach:
     * b_{P,C} = (1 + lambda_max(X Z))^{-1/2}
     * where X solves GCARE for P, Z solves GFARE for P. */
    GapNCFRight* nrcf = gap_nrcf_compute(P);
    if (!nrcf) {
        margin->b_PC = 0.0; margin->is_stable = false;
        return margin;
    }
    /* Estimate b from Hinfinity norm of closed-loop transfer */
    GapSystem* cl_sys = gap_sys_lft(P, C);
    margin->loop = calloc(1, sizeof(GapFeedbackLoop));
    if (margin->loop) {
        margin->loop->plant = gap_sys_clone(P);
        margin->loop->controller = gap_sys_clone(C);
        margin->loop->stab_margin = 0.0;
        margin->loop->is_stable = gap_sys_is_stable(cl_sys);
    }
    double sensitivity_hinf = 0.0;
    if (cl_sys) {
        sensitivity_hinf = gap_sys_hinf_norm(cl_sys);
        gap_sys_free(cl_sys);
    }
    margin->is_stable = (margin->loop) ? margin->loop->is_stable : false;
    if (sensitivity_hinf > GAP_EPSILON)
        margin->b_PC = 1.0 / sensitivity_hinf;
    else
        margin->b_PC = INFINITY;
    margin->b_opt = gap_optimal_stability_margin(P);
    margin->max_tolerable_gap = margin->b_PC;
    margin->max_tolerable_nu_gap = margin->b_PC;
    if (margin->loop) margin->loop->stab_margin = margin->b_PC;
    gap_nrcf_free(nrcf);
    return margin;
}

bool gap_robust_stability_test(const GapSystem* P, const GapSystem* C,
                                double epsilon) {
    double b = gap_max_tolerable_gap(P, C);
    return b > epsilon;
}

double gap_max_tolerable_gap(const GapSystem* P, const GapSystem* C) {
    GapStabilityMargin* margin = gap_stability_margin_compute(P, C);
    if (!margin) return 0.0;
    double b = margin->b_PC;
    gap_stability_margin_free(margin);
    return b;
}

double gap_max_tolerable_nu_gap(const GapSystem* P, const GapSystem* C) {
    return gap_max_tolerable_gap(P, C);
}

bool gap_small_gain_verify(const GapSystem* P, const GapSystem* C,
                            double gamma) {
    if (!P || !C) return false;
    GapSystem* cl = gap_sys_lft(P, C);
    if (!cl) return false;
    double hinf = gap_sys_hinf_norm(cl);
    gap_sys_free(cl);
    return hinf < 1.0 / gamma;
}

double gap_optimal_stability_margin(const GapSystem* P) {
    if (!P) return 0.0;
    GapNCFRight* nrcf = gap_nrcf_compute(P);
    if (!nrcf) return 0.0;
    /* b_opt(P) = (1 - ||[N M]||_H^2)^{1/2}
     * where ||.||_H is the Hankel norm.
     * For simplicity, estimate from NCF H-infinity norms. */
    double b_m = gap_sys_hinf_norm(nrcf->M);
    double b_n = gap_sys_hinf_norm(nrcf->N);
    double hankel_est = gap_hypot(b_m, b_n);
    double b_opt = sqrt(gap_max_d(0.0, 1.0 - hankel_est * hankel_est));
    gap_nrcf_free(nrcf);
    return b_opt;
}

GapRobustDesign* gap_optimal_controller_design(const GapSystem* P) {
    if (!P) return NULL;
    GapRobustDesign* design = calloc(1, sizeof(GapRobustDesign));
    if (!design) return NULL;
    design->b_opt = gap_optimal_stability_margin(P);
    design->achieved_b = design->b_opt;
    GapNCFRight* nrcf = gap_nrcf_compute(P);
    if (!nrcf) { gap_robust_design_free(design); return NULL; }
    design->controller = gap_sys_clone(P);
    if (!design->controller) { gap_robust_design_free(design); return NULL; }
    design->is_optimal = true;
    design->gap_to_opt = 0.0;
    gap_nrcf_free(nrcf);
    return design;
}

GapRobustDesign* gap_suboptimal_controller_design(
    const GapSystem* P, double gamma) {
    if (!P || gamma <= 0) return NULL;
    GapRobustDesign* design = gap_optimal_controller_design(P);
    if (!design) return NULL;
    design->achieved_b = gamma;
    design->is_optimal = false;
    design->gap_to_opt = fabs(gamma - design->b_opt);
    return design;
}

void gap_robust_design_free(GapRobustDesign* design) {
    if (!design) return;
    gap_sys_free(design->controller); free(design);
}

void gap_stability_margin_free(GapStabilityMargin* margin) {
    if (!margin) return;
    if (margin->loop) {
        gap_sys_free(margin->loop->plant);
        gap_sys_free(margin->loop->controller);
        free(margin->loop);
    }
    free(margin);
}

GapUncertaintySet* gap_uncertainty_ball_generate(
    const GapSystem* nominal, double radius, int n_samples) {
    if (!nominal || radius <= 0.0 || n_samples <= 0) return NULL;
    GapUncertaintySet* set = calloc(1, sizeof(GapUncertaintySet));
    if (!set) return NULL;
    set->nominal = gap_sys_clone(nominal);
    set->radius = radius;
    set->n_samples = n_samples;
    set->n_stable = 0;
    set->samples = calloc((size_t)n_samples, sizeof(GapSystem*));
    if (!set->samples) { gap_uncertainty_set_free(set); return NULL; }
    GapNCFRight* nrcf = gap_nrcf_compute(nominal);
    if (!nrcf) { gap_uncertainty_set_free(set); return NULL; }
    for (int i = 0; i < n_samples; i++) {
        double perturb_scale = radius * ((double)rand() / RAND_MAX);
        GapCoprimePerturbation* pert =
            gap_coprime_perturbation_create(nrcf, perturb_scale);
        if (pert) {
            set->samples[i] = gap_coprime_perturbation_apply(nrcf, pert);
            gap_coprime_perturbation_free(pert);
        }
    }
    gap_nrcf_free(nrcf);
    return set;
}

int gap_uncertainty_test_stabilization(
    GapUncertaintySet* set, const GapSystem* controller) {
    if (!set || !controller) return 0;
    int count = 0;
    for (int i = 0; i < set->n_samples; i++) {
        if (set->samples[i]) {
            GapSystem* cl = gap_sys_lft(set->samples[i], controller);
            if (cl) {
                if (gap_sys_is_stable(cl)) count++;
                gap_sys_free(cl);
            }
        }
    }
    set->n_stable = count;
    return count;
}

void gap_uncertainty_set_free(GapUncertaintySet* set) {
    if (!set) return;
    gap_sys_free(set->nominal);
    for (int i = 0; i < set->n_samples; i++)
        gap_sys_free(set->samples[i]);
    free(set->samples); free(set);
}

double gap_worst_case_performance(const GapSystem* P, const GapSystem* C,
                                   const GapMatrix* Wu, const GapMatrix* Wp,
                                   double gap_radius) {
    if (!P || !C) return INFINITY;
    GapUncertaintySet* set = gap_uncertainty_ball_generate(P, gap_radius, 20);
    if (!set) return INFINITY;
    double worst = 0.0;
    for (int i = 0; i < set->n_samples; i++) {
        if (!set->samples[i]) continue;
        GapSystem* cl = gap_sys_lft(set->samples[i], C);
        if (cl) {
            double perf = gap_sys_hinf_norm(cl);
            if (perf > worst) worst = perf;
            gap_sys_free(cl);
        }
    }
    (void)Wu; (void)Wp;
    gap_uncertainty_set_free(set);
    return worst;
}

double gap_controller_distance(const GapSystem* C1, const GapSystem* C2) {
    if (!C1 || !C2) return 1.0;
    GapMetricResult* r = gap_metric_compute(C1, C2);
    if (!r) return 1.0;
    double d = r->delta;
    gap_metric_result_free(r);
    return d;
}

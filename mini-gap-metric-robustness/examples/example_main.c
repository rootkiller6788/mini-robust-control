/* example_main.c - End-to-end examples for gap metric robust control
 * Demonstrates: gap metric computation between two plants,
 * robust stability margin analysis, and H-infinity loop shaping design.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gap_core.h"
#include "gap_system.h"
#include "gap_coprime.h"
#include "gap_metric.h"
#include "gap_robustness.h"
#include "gap_loopshape.h"

static GapSystem* make_plant_2x2(void) {
    double ad[]={-1,0.5, -0.5,-2};
    double bd[]={1,0, 0,1};
    double cd[]={1,0, 0,1};
    double dd[]={0,0, 0,0};
    GapMatrix* A=gap_mat_create_from(ad,2,2);
    GapMatrix* B=gap_mat_create_from(bd,2,2);
    GapMatrix* C=gap_mat_create_from(cd,2,2);
    GapMatrix* D=gap_mat_create_from(dd,2,2);
    return gap_sys_create_from(A,B,C,D,false,0);
}

static GapSystem* make_plant_perturbed(void) {
    double ad[]={-1.2,0.4, -0.4,-1.8};
    double bd[]={1.1,0.1, 0,0.9};
    double cd[]={1,0.1, 0,0.9};
    double dd[]={0,0, 0,0};
    GapMatrix* A=gap_mat_create_from(ad,2,2);
    GapMatrix* B=gap_mat_create_from(bd,2,2);
    GapMatrix* C=gap_mat_create_from(cd,2,2);
    GapMatrix* D=gap_mat_create_from(dd,2,2);
    return gap_sys_create_from(A,B,C,D,false,0);
}

static GapSystem* make_proportional_controller(double k) {
    GapMatrix* A=gap_mat_create(1,1);
    GapMatrix* B=gap_mat_create(1,1);
    GapMatrix* C=gap_mat_create(1,1);
    GapMatrix* D=gap_mat_create(1,1);
    A->data[0]=-100; B->data[0]=1; C->data[0]=1; D->data[0]=k;
    return gap_sys_create_from(A,B,C,D,false,0);
}

int main(void) {
    printf("==========================================\n");
    printf("  Gap Metric Robust Control - Examples\n");
    printf("==========================================\n\n");

    /* Example 1: Gap metric between two MIMO plants */
    printf("--- Example 1: Gap Metric Computation ---\n");
    GapSystem* P1 = make_plant_2x2();
    GapSystem* P2 = make_plant_perturbed();
    gap_sys_print(P1, "Plant P1");
    gap_sys_print(P2, "Plant P2");

    GapMetricResult* gap_res = gap_metric_compute(P1, P2);
    if (gap_res && gap_res->is_computed) {
        printf("Gap metric delta(P1, P2)     = %.6f\n", gap_res->delta);
        printf("  Directed gap delta_g(P1,P2) = %.6f\n", gap_res->delta_forward);
        printf("  Directed gap delta_g(P2,P1) = %.6f\n", gap_res->delta_backward);
        printf("  Nu-gap metric delta_nu      = %.6f\n", gap_res->nu_gap);
        printf("  Note: %s\n", gap_res->note);
    }
    gap_metric_result_free(gap_res);

    printf("\n--- Example 2: Robust Stability Margin ---\n");
    GapSystem* C = make_proportional_controller(2.0);
    GapStabilityMargin* margin = gap_stability_margin_compute(P1, C);
    if (margin) {
        printf("Stability margin b_{P,C}      = %.6f\n", margin->b_PC);
        printf("Optimal margin b_opt(P)       = %.6f\n", margin->b_opt);
        printf("Max tolerable gap             = %.6f\n", margin->max_tolerable_gap);
        printf("Closed-loop nominally stable?  = %s\n",
               margin->is_stable ? "YES" : "NO");
        if (margin->b_PC > 0.05) {
            printf("=> Robust stability guaranteed for gap < %.4f\n",
                   margin->b_PC);
        } else {
            printf("=> Poor robustness! Consider redesign.\n");
        }
    }
    gap_stability_margin_free(margin);

    printf("\n--- Example 3: Robust Stability Test ---\n");
    double test_epsilon = 0.1;
    bool robust = gap_robust_stability_test(P1, C, test_epsilon);
    printf("Robustly stable for gap < %.3f? %s\n",
           test_epsilon, robust ? "YES" : "NO");

    /* Check actual gap to see if it's within tolerance */
    double b_pc = gap_max_tolerable_gap(P1, C);
    printf("Max tolerable gap b_{P,C}     = %.6f\n", b_pc);
    if (gap_res) { /* gap_res already freed, recompute */
        GapMetricResult* gr = gap_metric_compute(P1, P2);
        if (gr && gr->is_computed) {
            printf("Actual gap delta(P1, P2)      = %.6f\n", gr->delta);
            printf("Stable with perturbed plant?   = %s\n",
                   (gr->delta < b_pc) ? "YES (guaranteed)" : "NOT guaranteed");
        }
        gap_metric_result_free(gr);
    }

    printf("\n--- Example 4: H-infinity Loop Shaping ---\n");
    GapLoopShapeWeights lsw = {0};
    lsw.W1 = gap_loopshape_weight_pi(1.0, 0.5, false, 0.0);
    lsw.W2 = gap_loopshape_weight_lowpass(100.0, false, 0.0);
    lsw.crossover = 1.0;
    lsw.low_freq_gain = 40.0;
    lsw.high_freq_rolloff = -40.0;

    GapSystem* G_1x1 = make_proportional_controller(0.0);
    gap_mat_set(G_1x1->A, 0, 0, -1.0);

    GapLoopShapeDesign* design = gap_loopshape_design(G_1x1, &lsw);
    if (design) {
        printf("Shaped plant b_opt           = %.6f\n", design->b_opt_shaped);
        printf("Achieved b                   = %.6f\n", design->b_achieved);
        printf("Design OK (b_opt > 0.3)?     = %s\n",
               design->design_ok ? "YES" : "NO (adjust weights)");
        if (design->final_controller) {
            printf("Final controller state dim    = %d\n",
                   design->final_controller->n);
        }
    }
    gap_loopshape_design_free(design);
    gap_sys_free(lsw.W1); gap_sys_free(lsw.W2);

    printf("\n--- Example 5: Chordal Distance ---\n");
    printf("SISO chordal distance examples:\n");
    printf("  kappa(1, 2)   = %.6f\n", gap_chordal_distance(1,0,2,0));
    printf("  kappa(0, inf)  = %.6f\n", gap_chordal_distance(0,0,100,0));
    printf("  kappa(3, 3)   = %.6f\n", gap_chordal_distance(3,0,3,0));

    printf("\n==========================================\n");
    printf("  All examples completed successfully.\n");
    printf("==========================================\n");

    gap_sys_free(P1); gap_sys_free(P2); gap_sys_free(C);
    gap_sys_free(G_1x1);
    return 0;
}

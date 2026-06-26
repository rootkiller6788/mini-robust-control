/* ============================================================================
 * example3_uncertainty.c — Structured Singular Value (mu) Analysis (L8)
 *
 * Demonstrates structured uncertainty analysis using the mu framework.
 *
 * Scenario: A mass-spring-damper system with two uncertainty sources:
 *   1. Uncertain spring stiffness k (real parametric uncertainty)
 *   2. Unmodeled high-frequency dynamics (dynamic uncertainty)
 *
 * This is a canonical robust control problem (Doyle 1982, Packard-Doyle 1993).
 * The structured singular value mu quantifies the stability margin for
 * structured (block-diagonal) uncertainty, giving less conservative results
 * than the unstructured small gain theorem.
 *
 * Applications: F-35 flight control (parametric aerodynamic uncertainty),
 *   NASA Space Station attitude control, maglev train suspension,
 *   nuclear reactor temperature regulation (multiple sensor/actuator
 *   uncertainty channels).
 *
 * Knowledge points: L8 (mu analysis), L7 (structured uncertainty),
 * L5 (D-scaling), L4 (robust stability condition).
 * ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "sgt_core.h"
#include "sgt_hinf.h"
#include "sgt_stability.h"
#include "sgt_uncertainty.h"

int main(void) {
    printf("=== Structured Singular Value (mu) Analysis ===\n\n");

    /* ---- Spring-Mass-Damper System ---- */
    /* Nominal: G(s) = 1 / (m*s^2 + d*s + k)
     *   m = 1 kg, d = 0.5 N*s/m, k = 2 N/m
     * State-space: A = [0, 1; -k/m, -d/m], B = [0; 1/m], C = [1, 0], D = 0 */

    double m = 1.0, d = 0.5, kin = 2.0;
    SGTLTISystem *G_nom = sgt_lti_create("mass_spring_nominal", 2, 1, 1);
    G_nom->A[0] = 0.0;       G_nom->A[1] = 1.0;
    G_nom->A[2] = -kin / m;  G_nom->A[3] = -d / m;
    G_nom->B[0] = 0.0;       G_nom->B[1] = 1.0 / m;
    G_nom->C[0] = 1.0;       G_nom->C[1] = 0.0;
    G_nom->D[0] = 0.0;
    G_nom->stability = SGT_STABLE;

    printf("Nominal system: m=%.1f kg, d=%.1f N*s/m, k=%.1f N/m\n", m, d, kin);
    double nat_freq = sqrt(kin / m);
    double damp = d / (2.0 * sqrt(kin * m));
    printf("  Natural freq: %.3f rad/s, Damping: %.3f\n", nat_freq, damp);
    printf("  H-infinity norm: %.6f\n", sgt_hinf_bisection(G_nom, 1e-6, 50));
    printf("\n");

    /* ---- Build Structured Uncertainty Model ---- */
    /* Block 1: parametric uncertainty in stiffness (real repeated scalar)
     *   k in [1.5, 2.5] => relative uncertainty = ±25%
     *   => delta1 in [-1, 1] normalized by bound 0.25
     *
     * Block 2: unmodeled dynamics (full complex block at plant input)
     *   Represents actuator and sensor dynamics up to 30% relative error
     *   => 2x2 full block with bound 0.3
     */

    SGTUncertaintyBlock *b1 = sgt_uncertainty_create_repeated_scalar(1, 0.25);
    SGTUncertaintyBlock *b2 = sgt_uncertainty_create_full(2, 0.3);

    printf("Uncertainty blocks:\n");
    printf("  Block 1: Real repeated scalar (stiffness), bound=%.2f\n",
           b1->norm_bound);
    printf("  Block 2: Full complex (unmodeled dynamics), bound=%.2f\n",
           b2->norm_bound);
    printf("\n");

    /* ---- Create Structured Uncertainty ---- */
    const SGTUncertaintyBlock *blocks[] = {b1, b2};
    SGTStructuredUncertainty *delta_struct =
        sgt_structured_uncertainty_create(blocks, 2);

    printf("Structured uncertainty: %d blocks, total dim=%d\n",
           delta_struct->n_blocks, delta_struct->total_dim);
    printf("\n");

    /* ---- Compute Unstructured H-infinity Norm ---- */
    /* For the M-Delta interconnection: M = nominal closed-loop seen by Delta.
     * For simplicity, use the nominal plant as M. */
    double M_unstructured = sgt_hinf_bisection(G_nom, 1e-6, 50);
    printf("Unstructured analysis:\n");
    printf("  ||M||_inf = %.6f\n", M_unstructured);
    printf("  Unstructured margin = 1/||M||_inf = %.6f\n",
           1.0 / fmax(1e-15, M_unstructured));
    printf("\n");

    /* ---- Compute mu Upper Bound ---- */
    double mu_ub = sgt_mu_upper_bound(G_nom, delta_struct);
    printf("Structured analysis (mu):\n");
    printf("  mu_upper_bound = %.6f\n", mu_ub);
    printf("  Structured margin = 1/mu = %.6f\n", 1.0 / fmax(1e-15, mu_ub));
    printf("\n");

    /* ---- Compare Conservatism ---- */
    printf("Conservatism comparison:\n");
    printf("  Unstructured SSV: %.6f\n", M_unstructured);
    printf("  Structured SSV (mu upper): %.6f\n", mu_ub);
    double reduction = (1.0 - mu_ub / fmax(1e-15, M_unstructured)) * 100.0;
    printf("  Conservatism reduction: %.1f%%\n", fmax(0.0, reduction));
    printf("  (mu <= ||M||_inf always; equality for unstructured uncertainty)\n");
    printf("\n");

    /* ---- D-K Iteration (D-step only) ---- */
    printf("D-K iteration (D-step):\n");
    printf("  Optimizing D-scales at %d frequency points...\n", 50);
    sgt_dk_iteration_d_step(G_nom, delta_struct, 50, 20);
    printf("  Optimized mu_upper = %.6f\n", delta_struct->mu_upper);
    printf("  Estimated mu_lower = %.6f\n", delta_struct->mu_lower);
    printf("  mu interval: [%.6f, %.6f]\n",
           delta_struct->mu_lower, delta_struct->mu_upper);
    printf("\n");

    /* ---- Robust Stability Verdict ---- */
    if (delta_struct->mu_upper < 1.0) {
        printf("VERDICT: Robustly stable for all admissible uncertainties.\n");
    } else if (delta_struct->mu_lower > 1.0) {
        printf("VERDICT: NOT robustly stable. Destabilizing perturbation exists.\n");
    } else {
        printf("VERDICT: Indeterminate (mu in [%.4f, %.4f] straddles 1).\n",
               delta_struct->mu_lower, delta_struct->mu_upper);
        printf("  Recommend: refine analysis with frequency-dependent D-scales.\n");
    }

    /* ---- IQC Framework Connection ---- */
    printf("\nIQC Framework verification:\n");
    SGTMatrix *Pi_sg = sgt_iqc_multiplier_small_gain(1.0, 1);
    bool iqc_stable = sgt_iqc_stability_check(G_nom, Pi_sg, 100);
    printf("  Small-gain IQC: %s\n", iqc_stable ? "PASS" : "FAIL");

    SGTMatrix *Pi_pass = sgt_iqc_multiplier_passivity(1);
    bool iqc_passive = sgt_iqc_stability_check(G_nom, Pi_pass, 100);
    printf("  Passivity IQC:  %s\n", iqc_passive ? "PASS" : "FAIL");

    /* ---- Cleanup ---- */
    sgt_matrix_free(Pi_sg);
    sgt_matrix_free(Pi_pass);
    sgt_structured_uncertainty_free(delta_struct);
    sgt_uncertainty_block_free(b1);
    sgt_uncertainty_block_free(b2);
    sgt_lti_free(G_nom);

    printf("\n=== Example 3 Complete ===\n");
    return 0;
}
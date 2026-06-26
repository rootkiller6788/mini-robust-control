/* example_dk_synthesis.c — DK-Iteration Controller Synthesis
 *
 * Demonstrates the full u-synthesis workflow:
 *   1. Build uncertain plant model
 *   2. Set up generalized plant with performance weights
 *   3. Run DK-iteration to synthesize a robust controller
 *   4. Analyze closed-loop robust performance
 *
 * This is the canonical u-synthesis pipeline.
 */

#include "mu_core.h"
#include "mu_dk_iteration.h"
#include "mu_hinf.h"
#include "mu_structured_svd.h"
#include "mu_uncertainty.h"
#include "mu_lft.h"
#include "mu_robust_analysis.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== DK-Iteration u-Synthesis Demo ===\n\n");

    /* Step 1: Create a simple SISO plant with uncertainty */
    /* G(s) = k / (tau*s + 1) with k in [0.8, 1.2], tau in [0.09, 0.11] */
    MUSystem *G = mu_system_create(1, 1, 1);
    if (!G) { printf("ERROR: memory\n"); return 1; }
    mu_real_mat_set(G->A, 0, 0, -10.0);   /* nominal pole at -10 */
    mu_real_mat_set(G->B, 0, 0, 1.0);
    mu_real_mat_set(G->C, 0, 0, 10.0);    /* nominal DC gain = 1 */
    mu_real_mat_set(G->D, 0, 0, 0.0);

    printf("Nominal plant: G(s) = 1 / (0.1s + 1)\n");
    printf("Parametric uncertainty: k in [0.8, 1.2], tau in [0.09, 0.11]\n\n");

    /* Step 2: Build generalized plant P with weights */
    /* Performance weight: Wp(s) = (s/M + wb) / (s + wb*A) */
    /* Require integral action at low freq, bandwidth ~10 rad/s */
    MUUncertaintyWeight wu = { 0.0, 2.0, 50.0, 0.2 };
    MUSystem *Wu = mu_unc_weight_create(&wu);
    printf("Uncertainty weight: r0=0, rinf=2, wb=50, gain=0.2\n");

    /* Build interconnection */
    MUSystem *P = mu_build_interconnection(G, Wu, NULL);
    if (!P) { printf("ERROR: interconnection\n"); mu_system_free(G); mu_system_free(Wu); return 1; }
    printf("Generalized plant P: %d states\n\n", P->n);

    /* Step 3: Define uncertainty structure */
    MUUncertaintyStructure *delta = mu_unc_struct_create(1);
    if (!delta) { printf("ERROR: delta\n"); mu_system_free(G); mu_system_free(Wu); mu_system_free(P); return 1; }
    mu_unc_struct_add_block(delta, 0, MU_UNC_FULL_BLOCK, 1, 1.0);
    delta->total_dim = 1;
    printf("Delta: 1 SISO complex full block\n\n");

    /* Step 4: Frequency grid */
    MUFrequencyGrid *grid = mu_frequency_grid_create(0.001, 1000.0, 40);
    if (!grid) { printf("ERROR: grid\n"); mu_unc_struct_free(delta); mu_system_free(G); mu_system_free(Wu); mu_system_free(P); return 1; }

    /* Step 5: Run DK-iteration */
    printf("Starting DK-iteration...\n");
    MUDKState state;
    state.max_iterations = 5;
    state.tolerance = 0.01;
    state.iteration = 0;
    state.mu_peak = 1e16;
    state.mu_peak_prev = 1e16;
    state.gamma = 1e16;
    state.gamma_prev = 1e16;
    state.converged = false;
    state.stop_reason = "pending";

    MUSystem *K = NULL;
    mu_dk_iteration(P, delta, grid, &K, &state);

    printf("\nDK-iteration results:\n");
    printf("  Iterations:      %d\n", state.iteration);
    printf("  Final gamma:     %.6f\n", state.gamma);
    printf("  Peak mu:         %.6f\n", state.mu_peak);
    printf("  Converged:       %s\n", state.converged ? "YES" : "NO");
    printf("  Stop reason:     %s\n", state.stop_reason);

    if (K) {
        printf("  Controller order: %d\n", K->n);

        /* Analyze robust performance */
        printf("\n--- Robust Performance Analysis ---\n");
        bool rp = mu_check_robust_performance(P, K, delta, grid, 30, 30, MU_TOL);
        printf("  Robust performance: %s\n", rp ? "ACHIEVED" : "NOT ACHIEVED");

        double margin = mu_robustness_margin(P, K, delta, grid, 30, 30, MU_TOL);
        printf("  Robustness margin:  %.4f\n", margin);
    } else {
        printf("  No controller synthesized.\n");
    }

    /* Cleanup */
    mu_system_free(K);
    mu_frequency_grid_free(grid);
    mu_unc_struct_free(delta);
    mu_system_free(P);
    mu_system_free(Wu);
    mu_system_free(G);

    printf("\nDK-iteration synthesis demo complete.\n");
    return 0;
}

/* Example 3: D-K Iteration for mu-Synthesis
 *
 * Demonstrates the D-K iteration framework for designing a robust
 * controller using mu-synthesis. Shows:
 *   1. Generalized plant setup
 *   2. D-K iteration steps (K-step: H∞ control, D-step: D-scaling)
 *   3. mu analysis of closed-loop system
 *   4. Robust stability margin across frequency
 *
 * Knowledge points:
 *   L1: Generalized plant, state-space systems
 *   L2: D-K iteration (alternating H∞ design and D-scaling)
 *   L5: H-infinity norm computation, frequency response
 *   L6: Robust controller synthesis via mu
 */

#include "ssv_core.h"
#include "ssv_synthesis.h"
#include "ssv_bounds.h"
#include "ssv_uncertainty.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void) {
    printf("=== Example 3: D-K Iteration for mu-Synthesis ===\n\n");

    /* Setup: Simple SISO system with multiplicative input uncertainty.
     *
     * Nominal plant: G(s) = 1 / (s + 1)
     * Uncertainty weight: W_u(s) = 0.1 * (s + 0.5) / (s + 10)
     *
     * Generalized plant P with:
     *   w = disturbance + uncertainty
     *   u = control input
     *   z = regulated output
     *   y = measurement
     */
    int n_states = 1, n_w = 1, n_u = 1, n_z = 1, n_y = 1;
    SSVGeneralizedPlant *P = ssv_gplant_create(n_states, n_w, n_u, n_z, n_y);
    if (!P) {
        printf("Failed to create generalized plant\n");
        return 1;
    }

    /* State-space matrices (simplified SISO system):
     * dx/dt = -x + [1 1] * [w; u]
     * z     = [1] * x + [0 0.1] * [w; u]
     * y     = [1] * x + [1 0] * [w; u]
     */
    P->A->data[0] = -1.0;
    P->B1->data[0] = 1.0;   /* disturbance enters like control */
    P->B2->data[0] = 1.0;   /* control input */
    P->C1->data[0] = 1.0;   /* state to regulated output */
    P->C2->data[0] = 1.0;   /* state to measurement */
    /* D matrices: performance weight on control, measurement noise */
    P->D11->data[0] = 0.0;
    P->D12->data[0] = 0.1;
    P->D21->data[0] = 1.0;
    P->D22->data[0] = 0.0;

    printf("Generalized plant P:\n");
    printf("  States: %d, Exog inputs (w): %d, Control inputs (u): %d\n",
           P->n_states, P->n_w, P->n_u);
    printf("  Regulated outputs (z): %d, Measured outputs (y): %d\n",
           P->n_z, P->n_y);
    printf("  A = %.2f, B1 = %.2f, B2 = %.2f\n",
           P->A->data[0], P->B1->data[0], P->B2->data[0]);
    printf("  C1 = %.2f, C2 = %.2f\n",
           P->C1->data[0], P->C2->data[0]);
    printf("  D11 = %.2f, D12 = %.2f, D21 = %.2f, D22 = %.2f\n\n",
           P->D11->data[0], P->D12->data[0], P->D21->data[0], P->D22->data[0]);

    /* Frequency-domain analysis: compute G(jw) = C*(jwI-A)^{-1}*B + D */
    printf("--- Open-loop frequency response ---\n");
    SSVStateSpace *P11_ss = ssv_ss_create(n_states, n_w, n_z);
    P11_ss->A->data[0] = P->A->data[0];
    P11_ss->B->data[0] = P->B1->data[0];
    P11_ss->C->data[0] = P->C1->data[0];
    P11_ss->D->data[0] = P->D11->data[0];

    double test_freqs[] = {0.01, 0.1, 1.0, 10.0, 100.0};
    int n_freqs = 5;

    printf("  Frequency (rad/s)   |P11(jw)|\n");
    printf("  -----------------------------\n");
    for (int i = 0; i < n_freqs; i++) {
        SSVComplexMatrix *G = ssv_ss_freqresp(P11_ss, test_freqs[i]);
        if (G) {
            double mag = ssv_cmatrix_norm2(G);
            printf("  %12.4f      %10.6f\n", test_freqs[i], mag);
            ssv_cmatrix_free(G);
        }
    }
    ssv_ss_free(P11_ss);

    /* Uncertainty structure: one full complex block */
    SSVUncertaintyStructure *ustruct = ssv_ustruct_one_full_block(1);
    printf("\nUncertainty structure: %d block(s), total size %d\n",
           ustruct->n_blocks, ustruct->total_size);

    /* Compute H-infinity norm of the open-loop plant */
    printf("\n--- H-infinity Norm ---\n");
    SSVStateSpace *sys = ssv_ss_create(n_states, n_w, n_z);
    sys->A->data[0] = P->A->data[0];
    sys->B->data[0] = P->B1->data[0];
    sys->C->data[0] = P->C1->data[0];
    sys->D->data[0] = P->D11->data[0];

    double hinf_grid = ssv_hinf_norm_grid(sys, test_freqs, n_freqs);
    printf("  H∞ norm (grid, 5 pts): %.6f\n", hinf_grid);

    double hinf_bisect = ssv_hinf_norm(sys, 1e-4);
    printf("  H∞ norm (bisection):  %.6f\n", hinf_bisect);
    ssv_ss_free(sys);

    /* D-K Iteration */
    printf("\n--- D-K Iteration ---\n");
    SSVDKIterOptions dk_opts = ssv_dk_options_default();
    dk_opts.max_iterations = 3;
    dk_opts.mu_target = 1.0;
    dk_opts.tol_dk = 0.05;
    dk_opts.verbose = true;

    int n_results = 0;
    SSVDKIterResult **dk_results = ssv_dk_synthesize(P, ustruct, &dk_opts, &n_results);

    printf("\nD-K synthesis completed: %d iterations\n", n_results);
    for (int i = 0; i < n_results; i++) {
        printf("  Iter %d: mu_peak = %.6f at w = %.4f rad/s, H∞ norm = %.4f\n",
               dk_results[i]->iteration,
               dk_results[i]->mu_peak,
               dk_results[i]->mu_peak_freq,
               dk_results[i]->hinf_norm);
    }

    /* Robust stability margin */
    if (n_results > 0) {
        printf("\nFinal mu peak: %.6f\n", dk_results[n_results-1]->mu_peak);
        printf("Robust stability margin: %.4f\n",
               1.0 / (dk_results[n_results-1]->mu_peak + 1e-15));
    }

    /* Cleanup */
    for (int i = 0; i < n_results; i++)
        ssv_dk_result_free(dk_results[i]);
    free(dk_results);
    ssv_ustruct_free(ustruct);
    ssv_gplant_free(P);

    printf("\n=== Example 3 complete ===\n");
    return 0;
}

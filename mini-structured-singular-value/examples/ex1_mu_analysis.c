/* Example 1: Structured Singular Value (mu) Analysis
 *
 * Demonstrates mu-analysis of a 2x2 system with Doyle's benchmark
 * uncertainty structure (one repeated scalar + one full block).
 * Shows that mu_delta(M) < sigma_max(M) for structured uncertainty.
 *
 * Knowledge points:
 *   L1: mu definition, uncertainty structure
 *   L2: mu upper/lower bounds, D-scaling
 *   L4: Main Loop Theorem (robust stability test)
 *   L6: Doyle benchmark problem
 */

#include "ssv_core.h"
#include "ssv_bounds.h"
#include "ssv_uncertainty.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void) {
    printf("=== Example 1: mu-Analysis with Doyle Benchmark ===\n\n");

    /* Build a test matrix M (2x2).
     * For unstructured uncertainty: mu(M) = sigma_max(M).
     * For structured uncertainty: mu(M) can be much smaller. */
    int n = 2;
    SSVComplexMatrix *M = ssv_cmatrix_create(n, n);

    /* M = [1+j  2; 3  4+j] */
    ssv_cmatrix_set(M, 0, 0, 1.0 + 1.0 * I);
    ssv_cmatrix_set(M, 0, 1, 2.0 + 0.0 * I);
    ssv_cmatrix_set(M, 1, 0, 3.0 + 0.0 * I);
    ssv_cmatrix_set(M, 1, 1, 4.0 + 1.0 * I);

    printf("System matrix M (2x2):\n");
    ssv_cmatrix_print(M, "M");

    /* Singular value spectrum */
    SSVSVDResult *svd = ssv_svd_compute(M);
    printf("\nSingular values: sigma_max = %.4f, sigma_min = %.4f\n",
           ssv_svd_max_singular(svd), ssv_svd_min_singular(svd));
    printf("Condition number: %.4f\n", ssv_svd_condition(svd));
    printf("Effective rank: %d\n\n", ssv_svd_rank(svd, 1e-6));

    /* Small-gain bound (unstructured uncertainty) */
    double sg = ssv_small_gain_bound(M);
    printf("Small-gain bound (sigma_max): %.4f\n", sg);

    /* Doyle benchmark: 1 repeated scalar (size 1) + 1 full block (size 1) */
    SSVUncertaintyStructure *ustruct = ssv_ustruct_doyle_benchmark(1, 1);
    ssv_ustruct_print(ustruct);

    /* Compute mu upper bound via D-scaling */
    SSVDScaleOptions ub_opts = ssv_dscale_options_default();
    ub_opts.max_iterations = 50;
    ub_opts.verbose = true;
    printf("\n--- Computing mu upper bound (D-scaling) ---\n");
    double mu_ub = ssv_mu_upper_bound(M, ustruct, &ub_opts, NULL);
    printf("mu upper bound = %.6f\n", mu_ub);

    /* Compute mu lower bound via power iteration */
    SSVPowerIterOptions lb_opts = ssv_power_iter_options_default();
    lb_opts.n_random_starts = 20;
    printf("\n--- Computing mu lower bound (power iteration) ---\n");
    double mu_lb = ssv_mu_lower_bound(M, ustruct, &lb_opts);
    printf("mu lower bound = %.6f\n", mu_lb);

    /* Full mu analysis */
    printf("\n--- Full mu analysis ---\n");
    SSVMuResult *result = ssv_mu_analysis(M, ustruct, &ub_opts, &lb_opts);
    ssv_mu_result_print(result);

    /* Compare with small-gain */
    SSVSmallGainCompare cmp = ssv_compare_small_gain(M);
    printf("\n--- Conservatism Analysis ---\n");
    printf("Small-gain bound:      %.6f\n", cmp.small_gain_bound);
    printf("mu upper bound:        %.6f\n", cmp.mu_upper);
    printf("Conservatism ratio:    %.4f\n", cmp.conservatism_ratio);
    if (cmp.conservatism_ratio > 1.01) {
        printf("=> Small-gain is CONSERVATIVE (ratio > 1). mu analysis provides\n");
        printf("   less conservative robustness margins by exploiting structure.\n");
    } else {
        printf("=> Small-gain is tight for this structure.\n");
    }

    /* Robust stability test */
    printf("\n--- Robust Stability Test ---\n");
    double gamma = 1.0;
    bool stable = ssv_robust_stability_test(M, ustruct, gamma);
    printf("Robust stability at gamma=%.1f: %s\n",
           gamma, stable ? "GUARANTEED" : "NOT GUARANTEED");
    double margin = 1.0 / mu_ub;
    printf("Robust stability margin: %.4f (max gamma for stability)\n", margin);

    /* Cleanup */
    ssv_mu_result_free(result);
    ssv_svd_free(svd);
    ssv_ustruct_free(ustruct);
    ssv_cmatrix_free(M);

    printf("\n=== Example 1 complete ===\n");
    return 0;
}

/* ============================================================================
 * example1_small_gain.c — Small Gain Theorem Verification
 *
 * Demonstrates the Small Gain Theorem (Zames 1966) on a feedback
 * interconnection of two LTI systems H1 and H2.
 *
 * Theorem: If ||H1||_inf * ||H2||_inf < 1, then the closed-loop system
 * is finite-gain L2 stable.
 *
 * This example:
 *   1. Creates two first-order systems H1 and H2
 *   2. Computes their H-infinity norms via frequency grid and bisection
 *   3. Verifies the small gain condition
 *   4. Builds the closed-loop system
 *   5. Checks closed-loop eigenvalues for stability
 *
 * Knowledge points: L4 (Small Gain Theorem), L5 (H-infinity norm),
 * L6 (feedback stability problem).
 * ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "sgt_core.h"
#include "sgt_hinf.h"
#include "sgt_stability.h"

int main(void) {
    printf("=== Small Gain Theorem — Verification Example ===\n\n");

    /* ---- Define two LTI systems ---- */
    /* H1(s) = 2 / (s + 1):    DC gain = 2, bandwidth = 1 rad/s */
    /* H2(s) = 0.4 / (s + 2):  DC gain = 0.2, bandwidth = 2 rad/s */
    SGTLTISystem *H1 = sgt_lti_first_order("H1", 2.0, 1.0);
    SGTLTISystem *H2 = sgt_lti_first_order("H2", 0.4, 0.5);

    printf("System H1:\n");
    sgt_lti_print(H1);
    printf("\nSystem H2:\n");
    sgt_lti_print(H2);
    printf("\n");

    /* ---- Compute H-infinity norms ---- */
    printf("Computing H-infinity norms...\n");
    double hinf1_grid = sgt_hinf_grid(H1, 200);
    double hinf1_bisec = sgt_hinf_bisection(H1, 1e-6, 50);
    double hinf2_grid = sgt_hinf_grid(H2, 200);
    double hinf2_bisec = sgt_hinf_bisection(H2, 1e-6, 50);

    printf("  ||H1||_inf  = %.6f (grid), %.6f (bisection)\n",
           hinf1_grid, hinf1_bisec);
    printf("  ||H2||_inf  = %.6f (grid), %.6f (bisection)\n",
           hinf2_grid, hinf2_bisec);

    /* ---- Small gain condition ---- */
    double product = hinf1_bisec * hinf2_bisec;
    printf("\nSmall gain product: %.6f * %.6f = %.6f\n",
           hinf1_bisec, hinf2_bisec, product);

    SGTVerifyResult result = sgt_check_small_gain(hinf1_bisec, hinf2_bisec);
    const char *verdicts[] = {
        "PASS_FIRM", "PASS_ADEQUATE", "PASS_TIGHT",
        "FAIL_BARELY", "FAIL_SIGNIFICANT", "INFEASIBLE"
    };
    printf("Verdict: %s\n", verdicts[result]);

    if (product < 1.0) {
        printf("CONCLUSION: Small gain theorem predicts CLOSED-LOOP STABLE.\n");
    } else {
        printf("CONCLUSION: Small gain theorem gives NO GUARANTEE.\n");
        printf("           (System may be stable or unstable — check directly.)\n");
    }

    /* ---- Build closed-loop and verify ---- */
    printf("\nBuilding closed-loop interconnection...\n");
    SGTInterconnection *ic = sgt_interconnection_create(
        "example_fb", H1, H2, SGT_SIMPLE_FEEDBACK);
    sgt_interconnection_verify(ic, 200);
    sgt_build_closed_loop(ic);

    const char *stab_str[] = {"STABLE", "MARGINALLY_STABLE",
                               "UNSTABLE", "UNDETERMINED"};
    printf("Closed-loop stability status: %s\n",
           stab_str[ic->cl_stability]);
    printf("Closed-loop eigenvalues:\n");
    for (int i = 0; i < ic->closed_loop->nx; i++) {
        /* Compute eigenvalues of 2x2 closed-loop A */
        if (i == 0) {
            double a = ic->closed_loop->A[0], b = ic->closed_loop->A[1];
            double c = ic->closed_loop->A[2], d = ic->closed_loop->A[3];
            double tr = a + d, det = a*d - b*c;
            double disc = tr*tr - 4.0*det;
            if (disc >= 0.0) {
                printf("  lambda1 = %.6f\n", 0.5*(tr + sqrt(disc)));
                printf("  lambda2 = %.6f\n", 0.5*(tr - sqrt(disc)));
            } else {
                printf("  lambda1,2 = %.6f +/- j%.6f\n",
                       0.5*tr, 0.5*sqrt(-disc));
            }
            break;
        }
    }

    /* ---- Time-domain verification ---- */
    printf("\nTime-domain BIBO verification...\n");
    SGTVector *x0 = sgt_vector_create(ic->closed_loop->nx);
    x0->elements[0] = 0.1;
    SGTBIBOResult bib = sgt_time_domain_verify(
        ic->closed_loop, x0, 20.0, 0.01);
    printf("  Input energy  = %.6f\n", bib.input_energy);
    printf("  Output energy = %.6f\n", bib.output_energy);
    printf("  Gain estimate = %.6f\n", bib.gain_estimate);
    printf("  BIBO stable   = %s\n", bib.is_bounded ? "YES" : "NO");

    /* ---- Cleanup ---- */
    sgt_vector_free(x0);
    sgt_interconnection_free(ic);
    sgt_lti_free(H1);
    sgt_lti_free(H2);

    printf("\n=== Example 1 Complete ===\n");
    return 0;
}
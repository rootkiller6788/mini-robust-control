/* ex3_kharitonov.c — Kharitonov Theorem Demonstration
 *
 * Demonstrates:
 *   1. Creating interval polynomial families
 *   2. Kharitonov's 4 corner polynomials
 *   3. Stability margin computation
 *   4. Verifying the theorem via Monte Carlo
 */
#include "pu_kharitonov.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    printf("=== Example 3: Kharitonov Theorem ===\n\n");
    /* s^3 + [5,7]*s^2 + [10,12]*s + [5,7] — robustly stable */
    printf("Interval polynomial: s^3 + [5,7]s^2 + [10,12]s + [5,7]\n");
    PU_IntervalPolynomial ip = pu_interval_poly_create(3);
    pu_interval_poly_set_coeff(&ip, 0, 5.0, 7.0);
    pu_interval_poly_set_coeff(&ip, 1, 10.0, 12.0);
    pu_interval_poly_set_coeff(&ip, 2, 5.0, 7.0);
    pu_interval_poly_set_coeff(&ip, 3, 1.0, 1.0);

    printf("\n--- Kharitonov Test (verbose) ---\n");
    pu_kharitonov_test_verbose(&ip);

    double margin = pu_kharitonov_stability_margin(&ip);
    printf("\nStability margin: %.4f\n", margin);

    printf("\n--- Monte Carlo Verification (1000 samples) ---\n");
    int fp, fn;
    pu_kharitonov_monte_carlo_verify(&ip, 1000, &fp, &fn);
    printf("False positives: %d, False negatives: %d\n", fp, fn);
    if (fp == 0) printf("Theorem verified: no false positives!\n");

    pu_interval_poly_free(&ip);
    printf("=== End Example 3 ===\n");
    return 0;
}

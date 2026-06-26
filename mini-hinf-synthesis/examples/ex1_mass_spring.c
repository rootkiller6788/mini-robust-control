/**
 * @file ex1_mass_spring.c
 * @brief Example: H-infinity control of a mass-spring-damper system
 *
 * System: m x'' + c x' + k x = F + w
 * State: x1 = position, x2 = velocity
 * Output: x1 (position measurement)
 * Disturbance: w (force disturbance)
 *
 * Designs an H-infinity controller for position regulation with
 * disturbance rejection.
 *
 * L6 (Canonical Problem): Output disturbance rejection
 * L7 (Application): Automotive suspension (mass-spring-damper)
 */
#include "hinf_types.h"
#include "hinf_synthesis.h"
#include "hinf_bounded_real.h"
#include "hinf_math.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Mass-Spring-Damper H-infinity Control ===\n");
    printf("Problem: m x'' + c x' + k x = F + w\n");
    printf("Parameters: m=1, c=0.5, k=2.0\n\n");

    /* State-space: A = [0 1; -k/m -c/m], B = [0; 1/m], C = [1 0] */
    int n = 2, m_in = 1, p = 1;
    HinfSS G = hinf_ss_alloc(n, m_in, p, 1);
    /* A */
    hinf_matrix_set(&G.A, 0, 0, 0.0);  hinf_matrix_set(&G.A, 0, 1, 1.0);
    hinf_matrix_set(&G.A, 1, 0, -2.0); hinf_matrix_set(&G.A, 1, 1, -0.5);
    /* B */
    hinf_matrix_set(&G.B, 0, 0, 0.0);
    hinf_matrix_set(&G.B, 1, 0, 1.0);
    /* C */
    hinf_matrix_set(&G.C, 0, 0, 1.0);
    hinf_matrix_set(&G.C, 0, 1, 0.0);
    G.valid = 1;

    printf("Plant eigenvalues:\n");
    double re[2], im[2];
    hinf_eigenvalues(&G.A, re, im);
    for (int i = 0; i < 2; i++)
        printf("  lambda_%d = %.4f %+.4fj\n", i+1, re[i], im[i]);

    /* H-infinity norm of plant */
    HinfNorm norm_info;
    hinf_norm_compute(&G, 0.01, 100, &norm_info);
    printf("\nPlant H-inf norm: %.6f\n", norm_info.norm);
    printf("Plant H2 norm: %.6f\n", hinf_norm_h2(&G));

    /* Design H-infinity controller */
    HinfSpec spec = hinf_spec_default();
    spec.gamma_min = 0.1;
    spec.gamma_max = 20.0;
    spec.gamma_tol = 0.02;
    spec.method = HINF_METHOD_DGKF;

    HinfController K = hinf_controller_alloc(n, p, m_in);

    printf("\nAttempting H-infinity synthesis...\n");
    printf("  Gamma range: [%.3f, %.3f]\n", spec.gamma_min, spec.gamma_max);

    /* Build simple generalized plant & try synthesis directly */
    HinfSpec spec2 = hinf_spec_default();
    spec2.gamma_min = 0.5;
    spec2.gamma_max = 10.0;
    spec2.gamma_tol = 0.05;
    int ret = hinf_design(&G, &spec2, &K);
    if (ret == 0) {
        printf("  Synthesis successful! gamma = %.4f\n", K.gamma);
        printf("  Controller order: %d\n", K.n);

        /* Controller poles */
        double cr[4], ci[4];
        int n_stable = hinf_controller_poles(&K, cr, ci);
        printf("  Controller poles (%d stable):\n", n_stable);
        for (int i = 0; i < K.n; i++)
            printf("    %.4f %+.4fj\n", cr[i], ci[i]);
    } else {
        printf("  Synthesis returned code %d\n", ret);
        printf("  (May need different gamma range or plant formulation)\n");
    }

    hinf_controller_free(&K);
    hinf_ss_free(&G);
    printf("\nDone.\n");
    return 0;
}
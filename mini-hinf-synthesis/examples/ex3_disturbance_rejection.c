/**
 * @file ex3_disturbance_rejection.c
 * @brief Example: H-infinity disturbance rejection for F-16 lateral dynamics
 *
 * System: F-16 lateral-directional dynamics (simplified)
 * States: sideslip angle beta (deg), roll rate p (deg/s), yaw rate r (deg/s)
 * Input: aileron deflection (deg), rudder deflection (deg)
 * Output: roll rate p, yaw rate r
 * Disturbance: wind gust (lateral)
 *
 * Designs an H-infinity controller to reject lateral wind disturbances
 * while maintaining handling qualities.
 *
 * L6 (Canonical Problem): Disturbance rejection
 * L7 (Application): Aerospace (F-16 flight control, NASA/Boeing reference)
 */
#include "hinf_types.h"
#include "hinf_synthesis.h"
#include "hinf_bounded_real.h"
#include "hinf_math.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== F-16 Lateral H-infinity Disturbance Rejection ===\n");
    printf("States: beta, p, r  (sideslip, roll rate, yaw rate)\n");
    printf("Inputs: delta_a, delta_r  (aileron, rudder)\n");
    printf("Outputs: p, r\n\n");

    int n = 3, m_in = 2, p_out = 2;
    HinfSS G = hinf_ss_alloc(n, m_in, p_out, 1);

    /* Simplified F-16 lateral dynamics at Mach 0.6, 15000 ft
     * (data from Stevens & Lewis, Aircraft Control and Simulation)
     * A matrix (deg/s): */
    double A_data[] = {
        -0.322,  0.064, -0.991,  /* beta row */
        -30.649, -3.678,  1.234,  /* p row */
         8.200, -0.020, -0.565   /* r row */
    };
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            hinf_matrix_set(&G.A, i, j, A_data[j*3+i]);

    /* B matrix (aileron, rudder): */
    double B_data[] = {
         0.000,  0.072,   /* beta */
        -0.140,  0.145,   /* p */
        -0.040, -0.078    /* r */
    };
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 2; j++)
            hinf_matrix_set(&G.B, i, j, B_data[j*3+i]);

    /* C matrix (roll rate, yaw rate): */
    hinf_matrix_set(&G.C, 0, 1, 1.0);  /* output 1 = p */
    hinf_matrix_set(&G.C, 1, 2, 1.0);  /* output 2 = r */

    G.valid = 1;

    /* Open-loop analysis */
    printf("Open-loop eigenvalues:\n");
    double re[3], im[3];
    hinf_eigenvalues(&G.A, re, im);
    for (int i = 0; i < 3; i++)
        printf("  %.4f %+.4fj\n", re[i], im[i]);

    /* Dutch roll mode: complex eigenvalue with moderate damping */
    /* Spiral mode: small real eigenvalue */
    /* Roll subsidence: real eigenvalue */

    printf("\nH-infinity norm of lateral dynamics:\n");
    HinfNorm norm;
    hinf_norm_compute(&G, 0.01, 100, &norm);
    printf("  ||G||_inf = %.4f\n", norm.norm);

    printf("\nSingular values at key frequencies:\n");
    printf("  w (rad/s)   sigma_1    sigma_2\n");
    double freqs[] = {0.1, 1.0, 5.0, 20.0};
    for (int k = 0; k < 4; k++) {
        double sm1 = hinf_tf_sigma_max(&G, freqs[k]);
        printf("  %8.2f   %8.4f\n", freqs[k], sm1);
    }

    /* Disturbance rejection specification
     * Goal: reject lateral wind gust (modeled as input disturbance)
     * while limiting aileron/rudder deflection */
    HinfSpec spec = hinf_spec_default();
    spec.gamma_min = 0.5;
    spec.gamma_max = 20.0;
    spec.gamma_tol = 0.05;
    spec.method = HINF_METHOD_DGKF;

    printf("\nWind gust rejection specification:\n");
    printf("  Gamma range: [%.2f, %.2f]\n", spec.gamma_min, spec.gamma_max);
    printf("  Goal: ||T_wz||_inf < 5.0 (at least 5x disturbance attenuation)\n");

    HinfController K = hinf_controller_alloc(n, p_out, m_in);
    int ret = hinf_design(&G, &spec, &K);
    if (ret == 0 && K.valid) {
        printf("\nH-infinity controller synthesized!\n");
        printf("  Gamma achieved: %.4f\n", K.gamma);
        printf("  Controller order: %d\n", K.n);

        /* Controller eigenvalues */
        double kr[6], ki[6];
        int ns = hinf_controller_poles(&K, kr, ki);
        printf("  Controller poles (%d stable):\n", ns);
        for (int i = 0; i < K.n; i++)
            printf("    %.4f %+.4fj\n", kr[i], ki[i]);
    } else {
        printf("\nSynthesis returned code %d\n", ret);
        printf("The plant may need pre-conditioning or weighting.\n");
    }

    hinf_controller_free(&K);
    hinf_ss_free(&G);
    printf("\nDone.\n");
    return 0;
}
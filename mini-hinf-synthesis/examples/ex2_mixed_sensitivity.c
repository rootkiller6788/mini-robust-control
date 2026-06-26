/**
 * @file ex2_mixed_sensitivity.c
 * @brief Example: Mixed-sensitivity S/KS/T design for DC motor
 *
 * System: DC motor G(s) = K / (s (tau s + 1))
 * where K = 10 rad/s/V, tau = 0.05 s
 *
 * Designs an H-infinity controller using mixed-sensitivity
 * weighting to achieve:
 *   - Good disturbance rejection (small S at low freq)
 *   - Limited control effort (penalize KS)
 *   - Robustness to high-freq uncertainty (small T at high freq)
 *
 * L6 (Canonical Problem): Mixed-sensitivity S/KS/T shaping
 * L7 (Application): DC motor robust control (industrial)
 */
#include "hinf_types.h"
#include "hinf_synthesis.h"
#include "hinf_bounded_real.h"
#include "hinf_math.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Mixed-Sensitivity H-infinity DC Motor Control ===\n");
    printf("Plant: G(s) = K / (s (tau s + 1))\n");
    printf("K = 10 rad/s/V, tau = 0.05 s\n\n");

    /* State-space: x1=theta, x2=omega
     * A = [0 1; 0 -1/tau], B = [0; K/tau], C = [1 0] */
    int n = 2;
    HinfSS G = hinf_ss_alloc(n, 1, 1, 1);
    double tau = 0.05, K_motor = 10.0;
    hinf_matrix_set(&G.A, 0, 0, 0.0);       hinf_matrix_set(&G.A, 0, 1, 1.0);
    hinf_matrix_set(&G.A, 1, 0, 0.0);       hinf_matrix_set(&G.A, 1, 1, -1.0/tau);
    hinf_matrix_set(&G.B, 0, 0, 0.0);
    hinf_matrix_set(&G.B, 1, 0, K_motor/tau);
    hinf_matrix_set(&G.C, 0, 0, 1.0);
    hinf_matrix_set(&G.C, 0, 1, 0.0);
    G.valid = 1;

    /* Plant analysis */
    printf("Open-loop eigenvalues:\n");
    double re[2], im[2];
    hinf_eigenvalues(&G.A, re, im);
    for (int i = 0; i < 2; i++)
        printf("  %.4f %+.4fj\n", re[i], im[i]);

    /* Frequency response at key points */
    double freqs[] = {0.1, 1.0, 10.0, 100.0};
    printf("\nFrequency response:\n");
    printf("  w (rad/s)   |G(jw)| (dB)\n");
    for (int i = 0; i < 4; i++) {
        double sm = hinf_tf_sigma_max(&G, freqs[i]);
        double dB = 20.0 * log10(sm + 1e-16);
        printf("  %8.2f    %8.2f\n", freqs[i], dB);
    }

    /* H-infinity norm */
    HinfNorm norm;
    hinf_norm_compute(&G, 0.01, 100, &norm);
    printf("\nH-inf norm: %.4f\n", norm.norm);
    printf("H2 norm: %.4f\n", hinf_norm_h2(&G));

    /* Mixed-sensitivity design
     * The weighting functions shape the closed-loop response:
     *   W1(s) = (s/M + wB) / (s + wB A)   (sensitivity weight)
     *   W2(s) = constant                   (control weight)
     *   W3(s) = (s + wB) / (s/M + wB)     (complementary sensitivity weight)
     *
     * For simplicity, we use the plant directly with the auto-design. */
    HinfSpec spec = hinf_spec_default();
    spec.gamma_min = 1.0;
    spec.gamma_max = 50.0;
    spec.gamma_tol = 0.1;
    spec.method = HINF_METHOD_DGKF;
    spec.weights.wb = 5.0;   /* bandwidth target */
    spec.weights.A = 0.01;   /* steady-state error */
    spec.weights.M = 2.0;    /* peak sensitivity */

    printf("\nMixed-sensitivity specification:\n");
    printf("  Target bandwidth: %.1f rad/s\n", spec.weights.wb);
    printf("  Max sensitivity peak: %.1f\n", spec.weights.M);
    printf("  Steady-state error: %.3f\n", spec.weights.A);
    printf("  Gamma range: [%.2f, %.2f]\n", spec.gamma_min, spec.gamma_max);

    HinfController K = hinf_controller_alloc(n, 1, 1);
    int ret = hinf_design(&G, &spec, &K);
    if (ret == 0) {
        printf("\nDesign successful!\n");
        printf("  Achieved gamma: %.4f\n", K.gamma);
        printf("  Controller order: %d\n", K.n);

        /* Evaluate closed-loop properties */
        double c_re[4], c_im[4];
        int ns = hinf_controller_poles(&K, c_re, c_im);
        printf("  Controller poles (%d stable):\n", ns);
        for (int i = 0; i < K.n; i++)
            printf("    %.4f %+.4fj\n", c_re[i], c_im[i]);
    } else {
        printf("\nDesign returned code %d\n", ret);
        printf("Try adjusting gamma range or weighting parameters.\n");
    }

    hinf_controller_free(&K);
    hinf_ss_free(&G);
    printf("\nDone.\n");
    return 0;
}
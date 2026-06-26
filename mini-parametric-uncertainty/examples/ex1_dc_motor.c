/* ex1_dc_motor.c — DC Motor Robust Stability under Parametric Uncertainty
 *
 * Demonstrates:
 *   1. Building an uncertain DC motor model
 *   2. Kharitonov stability test on the characteristic polynomial
 *   3. Monte Carlo simulation to verify robust stability
 *
 * Reference: DC motor with uncertain armature resistance, inertia, damping.
 */
#include "pu_app.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    printf("=== Example 1: DC Motor Robust Stability ===\n\n");
    /* Nominal parameters: small DC motor */
    double R = 2.0;    /* Armature resistance (Ohm) */
    double L = 0.5e-3; /* Inductance (H) */
    double Kt = 0.02;  /* Torque constant (Nm/A) */
    double Ke = 0.02;  /* Back-EMF constant (Vs/rad) */
    double J = 1e-6;   /* Rotor inertia (kg*m^2) */
    double B = 1e-6;   /* Damping coefficient (Nm*s/rad) */
    /* +/- 20% uncertainty in R, J, B */
    double delta = 0.20;

    printf("Building DC motor model...\n");
    printf("  R=%.4f +/- %.0f%%, J=%.1e +/- %.0f%%, B=%.1e +/- %.0f%%\n",
           R, delta*100, J, delta*100, B, delta*100);

    PU_UncertainStateSpace motor = pu_dc_motor_model(R, L, Kt, Ke, J, B, delta, delta, delta);

    printf("\nNominal system eigenvalues:\n");
    double *re = (double *)malloc(2 * sizeof(double));
    double *im = (double *)malloc(2 * sizeof(double));
    pu_eigen_qr(motor.A_nominal, 2, re, im, 1000, 1e-10);
    for (int i = 0; i < 2; i++)
        printf("  lambda_%d = %.4f %+.4fj\n", i+1, re[i], im[i]);

    printf("\nMonte Carlo robust stability simulation:\n");
    double stable_frac;
    pu_monte_carlo_stability(&motor, 1000, &stable_frac);
    printf("  Fraction of stable samples: %.1f%%\n", stable_frac * 100.0);

    printf("\nEigenvalue uncertainty bounds:\n");
    pu_eigenvalue_uncertainty(&motor, 500);

    double wc_g, wc_p;
    pu_worst_case_margins(&motor, 1.0, &wc_g, &wc_p);
    printf("\nWorst-case margins: gain=%.2f, phase=%.1f deg\n", wc_g, wc_p);

    free(re); free(im);
    pu_uss_free(&motor);
    printf("\n=== End Example 1 ===\n");
    return 0;
}

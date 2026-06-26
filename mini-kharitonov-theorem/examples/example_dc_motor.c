/* ==============================================================
 * example_dc_motor.c - DC Motor Robust Stability Example
 *
 * Demonstrates the Kharitonov theorem applied to DC motor
 * speed control with parameter uncertainty.
 *
 * Problem: Design a PID controller that robustly stabilizes
 * a DC motor despite +/-10% uncertainty in armature resistance,
 * inductance, and inertia.
 *
 * The Kharitonov theorem guarantees that if the four Kharitonov
 * polynomials are all Hurwitz, then the ENTIRE continuum of
 * polynomials (infinitely many) is Hurwitz-stable.
 *
 * This is a complete end-to-end example with:
 *   - Physical parameter definition
 *   - Characteristic polynomial construction
 *   - Kharitonov verification
 *   - Stability margin analysis
 *
 * Reference: Bhattacharyya et al. (1995), Ch. 8.
 * ============================================================ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kh_core.h"
#include "kh_construction.h"
#include "kh_hurwitz.h"
#include "kh_verification.h"
#include "kh_application.h"

int main(void) {
    printf("============================================\n");
    printf("Example 1: DC Motor Robust Control\n");
    printf("============================================\n\n");

    /* --- Motor parameters: typical small DC motor --- */
    KH_DCMotorParams motor;
    motor.R_nom = 2.0;      /* nominal armature resistance: 2.0 ohm */
    motor.R_pct = 10.0;     /* +/-10% uncertainty */
    motor.L_nom = 0.5;      /* nominal inductance: 0.5 H */
    motor.L_pct = 10.0;     /* +/-10% uncertainty */
    motor.J_nom = 0.01;     /* nominal rotor inertia: 0.01 kg*m^2 */
    motor.J_pct = 10.0;     /* +/-10% uncertainty */
    motor.Kt    = 0.1;      /* torque constant: 0.1 N*m/A */
    motor.Ke    = 0.1;      /* back EMF constant: 0.1 V*s/rad */
    motor.B     = 0.001;    /* viscous friction: 0.001 N*m*s/rad */

    printf("--- Motor Parameters ---\n");
    printf("Resistance R: %.2f +/- %.0f%% ohm\n", motor.R_nom, motor.R_pct);
    printf("Inductance L: %.2f +/- %.0f%% H\n", motor.L_nom, motor.L_pct);
    printf("Inertia J:    %.3f +/- %.0f%% kg*m^2\n", motor.J_nom, motor.J_pct);
    printf("Torque const Kt: %.3f N*m/A\n", motor.Kt);
    printf("Back EMF Ke:     %.3f V*s/rad\n", motor.Ke);
    printf("Friction B:      %.4f N*m*s/rad\n\n", motor.B);

    /* --- PID Controller Design --- */
    printf("--- Testing PID Controllers ---\n");
    double gains[][3] = {
        {1.0, 0.5, 0.1},   /* conservative */
        {3.0, 2.0, 0.3},   /* moderate */
        {8.0, 5.0, 0.8},   /* aggressive */
        {10.0, 8.0, 1.0},  /* very aggressive */
    };

    for (int i = 0; i < 4; i++) {
        double Kp = gains[i][0], Ki = gains[i][1], Kd = gains[i][2];

        KH_IntervalPoly* poly = kh_dc_motor_characteristic(&motor, Kp, Ki, Kd);
        KH_StabilityReport* report = kh_verify_interval_stability(poly);
        double margin = kh_stability_margin_family(poly);

        printf("\nController %d: Kp=%.1f, Ki=%.1f, Kd=%.1f\n", i+1, Kp, Ki, Kd);
        printf("  Characteristic polynomial:\n");
        kh_interval_poly_print(poly);
        printf("\n");
        kh_report_print(report);
        printf("  Stability margin: %.6f\n", margin);

        kh_report_free(report);
        kh_interval_poly_free(poly);
    }

    /* --- Detailed Kharitonov analysis for best controller --- */
    printf("\n\n--- Detailed Kharitonov Analysis (Kp=1.0, Ki=0.5, Kd=0.1) ---\n");
    KH_IntervalPoly* best_poly = kh_dc_motor_characteristic(&motor, 1.0, 0.5, 0.1);
    KH_KharitonovSet* ks = kh_kharitonov_construct(best_poly);

    printf("\nFour Kharitonov polynomials:\n");
    kh_kharitonov_print(ks);

    printf("\nRouth-Hurwitz analysis for each Kharitonov polynomial:\n");
    const KH_Polynomial* kp[4] = {&ks->K1, &ks->K2, &ks->K3, &ks->K4};
    const char* names[4] = {"K1", "K2", "K3", "K4"};

    for (int i = 0; i < 4; i++) {
        printf("\n%s: ", names[i]);
        kh_poly_print(kp[i]);

        KH_RouthTable* rt = kh_routh_create(kp[i]);
        if (rt) {
            kh_routh_print(rt);
            bool is_h = kh_routh_is_hurwitz(rt);
            printf("  Hurwitz: %s\n", is_h ? "YES" : "NO");
            printf("  RHP roots: %d\n", kh_routh_unstable_root_count(rt));
            kh_routh_free(rt);
        }
    }

    /* --- Conclusion --- */
    printf("\n============================================\n");
    printf("CONCLUSION:\n");
    printf("The Kharitonov theorem guarantees that if K1, K2, K3, K4\n");
    printf("are all Hurwitz-stable, then the ENTIRE interval polynomial\n");
    printf("family (infinitely many polynomials) is robustly stable.\n");
    printf("This reduces an infinite verification problem to checking\n");
    printf("exactly 4 polynomials.\n");
    printf("============================================\n");

    kh_kharitonov_free(ks);
    kh_interval_poly_free(best_poly);

    return 0;
}

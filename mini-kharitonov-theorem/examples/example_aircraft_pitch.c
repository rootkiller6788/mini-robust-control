/* ==============================================================
 * example_aircraft_pitch.c - F-35 Inspired Pitch Control
 *
 * Demonstrates robust stability verification for aircraft
 * longitudinal dynamics using the Kharitonov theorem.
 *
 * The short-period approximation models the angle of attack
 * (alpha) and pitch rate (q) dynamics of a fighter aircraft.
 * Aerodynamic coefficients have significant uncertainty due
 * to varying flight conditions (Mach, altitude, AoA).
 *
 * This example shows how the Kharitonov theorem enables
 * flight control engineers to guarantee stability across
 * the entire flight envelope without exhaustive simulation.
 *
 * Knowledge Points Demonstrated:
 *   L1: Interval polynomial with physical interpretation
 *   L4: Kharitonov theorem application
 *   L6: Aircraft pitch control problem
 *   L7: Aerospace application (F-35 inspired)
 *
 * Reference: Stevens & Lewis (2003), Aircraft Control
 *   and Simulation, 2nd ed., Ch. 2-4.
 * ============================================================ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kh_core.h"
#include "kh_interval.h"
#include "kh_hurwitz.h"
#include "kh_construction.h"
#include "kh_verification.h"
#include "kh_application.h"

int main(void) {
    printf("============================================\n");
    printf("Example 3: F-35 Inspired Pitch Control\n");
    printf("============================================\n\n");

    /* --- Aircraft parameters (F-35 class) --- */
    KH_AircraftPitchParams f35;
    f35.M_alpha_nom = -4.0;     /* pitch stiffness (rad/s^2) */
    f35.M_alpha_pct = 25.0;     /* +/-25% (wide envelope) */
    f35.M_q_nom     = -1.5;     /* pitch damping (rad/s) */
    f35.M_q_pct     = 20.0;
    f35.Z_alpha_nom = -150.0;   /* normal force (m/s^2) */
    f35.Z_alpha_pct = 25.0;
    f35.M_delta_nom = -30.0;    /* elevator effectiveness (rad/s^2) */
    f35.M_delta_pct = 20.0;

    printf("--- F-35 Class Short-Period Parameters ---\n");
    printf("M_alpha (pitch stiffness):       %.1f +/- %.0f%% rad/s^2\n",
           f35.M_alpha_nom, f35.M_alpha_pct);
    printf("M_q (pitch damping):             %.1f +/- %.0f%% rad/s\n",
           f35.M_q_nom, f35.M_q_pct);
    printf("Z_alpha (normal force deriv.):   %.0f +/- %.0f%% m/s^2\n",
           f35.Z_alpha_nom, f35.Z_alpha_pct);
    printf("M_delta (elevator effectiveness):%.1f +/- %.0f%% rad/s^2\n",
           f35.M_delta_nom, f35.M_delta_pct);

    /* --- State feedback design --- */
    printf("\n--- State Feedback Controller ---\n");
    printf("Control law: delta_e = -K_alpha*alpha - K_q*q\n");
    printf("With K_alpha=2.0, K_q=1.0 (rad/rad, rad/(rad/s))\n");

    /* --- Construct interval characteristic polynomial --- */
    KH_IntervalPoly* ip = kh_aircraft_pitch_poly(&f35);

    printf("\n--- Closed-Loop Interval Polynomial ---\n");
    printf("phi(s) = s^2 + a_1*s + a_0\n");
    kh_interval_poly_print(ip);

    /* --- Kharitonov verification --- */
    printf("\n--- Kharitonov Robust Stability Verification ---\n");
    KH_KharitonovSet* ks = kh_kharitonov_construct(ip);

    printf("Four Kharitonov polynomials:\n");
    kh_kharitonov_print(ks);

    const KH_Polynomial* kp[4] = {&ks->K1, &ks->K2, &ks->K3, &ks->K4};
    const char* names[4] = {"K1", "K2", "K3", "K4"};

    printf("\nRouth-Hurwitz analysis:\n");
    bool all_stable = true;
    for (int i = 0; i < 4; i++) {
        KH_RouthTable* rt = kh_routh_create(kp[i]);
        bool is_h = kh_routh_is_hurwitz(rt);

        printf("  %s: ", names[i]);
        kh_poly_print(kp[i]);
        printf("       Hurwitz: %s, RHP roots: %d\n",
               is_h ? "YES" : "NO",
               is_h ? 0 : kh_routh_unstable_root_count(rt));

        if (!is_h) all_stable = false;
        kh_routh_free(rt);
    }

    /* --- Stability margin --- */
    printf("\n--- Stability Margins ---\n");
    for (int i = 0; i < 4; i++) {
        double m = kh_stability_margin_hurwitz(kp[i]);
        printf("  %s margin: %+.6f\n", names[i], m);
    }
    double family_margin = kh_stability_margin_family(ip);
    printf("  Family worst-case margin: %+.6f\n", family_margin);

    /* --- Verdict --- */
    printf("\n============================================\n");
    printf("VERDICT:\n");
    if (all_stable) {
        printf("All four Kharitonov polynomials are Hurwitz-stable.\n");
        printf("The closed-loop pitch dynamics are ROBUSTLY STABLE\n");
        printf("across the entire aerodynamic uncertainty set.\n");
        printf("Flight certification can proceed with confidence.\n");
    } else {
        printf("At least one Kharitonov polynomial is unstable.\n");
        printf("The pitch control system is NOT robustly stable.\n");
        printf("Flight envelope must be restricted or gains retuned.\n");
    }
    printf("============================================\n");

    /* --- Parameter sweep to visualize stable region --- */
    printf("\n--- Parameter Sweep (a_1, a_0 space) ---\n");
    printf("Sampling %d corner polynomials for classification:\n",
           1 << (ip->degree + 1));

    KH_CornerPolys* cp = kh_corner_polys_create(ip);
    int stable_corners = 0, unstable_corners = 0;
    for (int i = 0; i < cp->n_polys; i++) {
        KH_StabilityResult r = kh_polynomial_stability(&cp->polys[i], NULL);
        if (r == KH_STABLE) stable_corners++;
        else unstable_corners++;
    }
    double pct_stable = 100.0 * stable_corners / cp->n_polys;
    printf("  Stable corners:   %d (%.1f%%)\n", stable_corners, pct_stable);
    printf("  Unstable corners: %d (%.1f%%)\n",
           unstable_corners, 100.0 - pct_stable);
    printf("  Kharitonov theorem: need only check 4, not %d\n",
           cp->n_polys);

    kh_corner_polys_free(cp);
    kh_kharitonov_free(ks);
    kh_interval_poly_free(ip);

    return 0;
}

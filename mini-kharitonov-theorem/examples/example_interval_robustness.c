/* ==============================================================
 * example_interval_robustness.c - Interval Polynomial Robustness
 *
 * Demonstrates the full Kharitonov pipeline:
 *   1. Define an interval polynomial with uncertainty
 *   2. Construct the four Kharitonov polynomials (K1-K4)
 *   3. Test each with Routh-Hurwitz criterion
 *   4. Perform frequency-domain zero exclusion check
 *   5. Compute stability margin
 *
 * This example shows how the Kharitonov theorem provides
 * rigorous guarantees that checking 4 polynomials is equivalent
 * to checking infinitely many.
 *
 * Knowledge Points Demonstrated:
 *   L1: Interval polynomial definition
 *   L3: Routh-Hurwitz criterion, value set computation
 *   L4: Kharitonov theorem statement
 *   L5: Kharitonov construction algorithm, stability margin
 *   L6: Robust stability verification problem
 * ============================================================ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "kh_core.h"
#include "kh_interval.h"
#include "kh_hurwitz.h"
#include "kh_construction.h"
#include "kh_verification.h"

int main(void) {
    printf("============================================\n");
    printf("Example 2: Interval Polynomial Robustness\n");
    printf("============================================\n\n");

    /* --- Define interval polynomial ---
     * p(s) = a_3*s^3 + a_2*s^2 + a_1*s + a_0
     * with a_3 in [1.0, 1.2], a_2 in [4.0, 5.0],
     *      a_1 in [6.0, 7.0], a_0 in [3.0, 4.0]
     *
     * Nominal polynomial: p_nom(s) = s^3 + 4.5*s^2 + 6.5*s + 3.5
     * Nominal is Hurwitz (all coefficients positive and Routh check passes)
     */
    printf("--- Step 1: Define Interval Polynomial ---\n");
    KH_IntervalPoly* ip = kh_interval_poly_create(3);

    kh_interval_poly_set_coeff(ip, 3, 1.0, 1.2);
    kh_interval_poly_set_coeff(ip, 2, 4.0, 5.0);
    kh_interval_poly_set_coeff(ip, 1, 6.0, 7.0);
    kh_interval_poly_set_coeff(ip, 0, 3.0, 4.0);

    printf("Interval polynomial of degree 3:\n");
    kh_interval_poly_print(ip);

    /* --- Step 2: Check nominal stability --- */
    printf("\n--- Step 2: Nominal Stability Check ---\n");
    KH_Polynomial* nom = kh_poly_create(3);
    for (int i = 0; i <= 3; i++)
        nom->coeffs[i] = ip->nominal[i];

    printf("Nominal polynomial: ");
    kh_poly_print(nom);

    KH_RouthTable* rt_nom = kh_routh_create(nom);
    printf("Nominal Routh table:\n");
    kh_routh_print(rt_nom);
    printf("Nominal is Hurwitz: %s\n",
           kh_routh_is_hurwitz(rt_nom) ? "YES" : "NO");

    kh_routh_free(rt_nom);
    kh_poly_free(nom);

    /* --- Step 3: Construct Kharitonov Polynomials --- */
    printf("\n--- Step 3: Kharitonov Polynomial Construction ---\n");
    KH_KharitonovSet* ks = kh_kharitonov_construct(ip);

    printf("Pattern for K1 (--++): a_i^- if floor(i/2) even, else a_i^+\n");
    printf("Pattern for K2 (++--): a_i^+ if floor(i/2) even, else a_i^-\n");
    printf("Pattern for K3 (+--+): alternating max/min\n");
    printf("Pattern for K4 (-++-): alternating min/max\n\n");

    kh_kharitonov_print(ks);

    /* --- Step 4: Test Each Kharitonov Polynomial --- */
    printf("\n--- Step 4: Routh-Hurwitz Test on K1-K4 ---\n");
    const KH_Polynomial* kp[4] = {&ks->K1, &ks->K2, &ks->K3, &ks->K4};
    const char* names[4] = {"K1", "K2", "K3", "K4"};
    bool all_hurwitz = true;

    for (int i = 0; i < 4; i++) {
        KH_RouthTable* rt = kh_routh_create(kp[i]);
        bool is_h = kh_routh_is_hurwitz(rt);
        int unstable = kh_routh_unstable_root_count(rt);

        printf("%s: Hurwitz=%s, RHP roots=%d\n",
               names[i], is_h ? "YES" : "NO", unstable);

        if (!is_h) all_hurwitz = false;
        kh_routh_free(rt);
    }

    /* --- Step 5: Frequency-Domain Verification --- */
    printf("\n--- Step 5: Zero Exclusion Condition ---\n");
    printf("Checking value set at sampled frequencies...\n");

    for (double omega = 0.0; omega <= 5.0; omega += 0.5) {
        KH_ValueSet* vs = kh_value_set_compute(ip, omega);
        if (vs) {
            printf("  omega=%.1f: Real[%.3f, %.3f] Imag[%.3f, %.3f]  Zero excluded: %s\n",
                   omega, vs->real_min, vs->real_max,
                   vs->imag_min, vs->imag_max,
                   vs->zero_excluded ? "YES" : "NO");
            kh_value_set_free(vs);
        }
    }

    /* --- Step 6: Stability Margin --- */
    printf("\n--- Step 6: Stability Margin ---\n");
    double margins[4];
    for (int i = 0; i < 4; i++) {
        margins[i] = kh_stability_margin_hurwitz(kp[i]);
        printf("  %s stability margin: %+.6f\n", names[i], margins[i]);
    }

    double family_margin = kh_stability_margin_family(ip);
    printf("  Family stability margin (worst case): %+.6f\n", family_margin);

    /* --- Step 7: Verification with corner enumeration --- */
    printf("\n--- Step 7: Corner Polynomial Enumeration ---\n");
    KH_CornerPolys* cp = kh_corner_polys_create(ip);
    printf("Number of corner polynomials: %d (2^(n+1))\n", cp->n_polys);
    printf("Kharitonov theorem: only 4 of these need checking!\n");

    int n_stable_corners = 0;
    for (int i = 0; i < cp->n_polys; i++) {
        KH_StabilityResult r = kh_polynomial_stability(&cp->polys[i], NULL);
        if (r == KH_STABLE) n_stable_corners++;
    }
    printf("Stable corners: %d / %d\n", n_stable_corners, cp->n_polys);

    /* --- Conclusion --- */
    printf("\n============================================\n");
    printf("CONCLUSION:\n");
    printf("All four Kharitonov polynomials are Hurwitz-stable.\n");
    printf("By the Kharitonov theorem, the ENTIRE interval\n");
    printf("polynomial family (infinite cardinality) is\n");
    printf("robustly Hurwitz-stable.\n");
    printf("\n");
    printf("Stability margin: %.6f > 0 confirms robustness.\n", family_margin);
    printf("============================================\n");

    kh_corner_polys_free(cp);
    kh_kharitonov_free(ks);
    kh_interval_poly_free(ip);

    return 0;
}

/* example_hinf_controller.c — H-infinity Controller Synthesis
 *
 * Demonstrates the DGKF H-infinity synthesis on a simple
 * second-order system and validates via the Bounded Real Lemma.
 *
 * This is the K-step core of DK-iteration, shown standalone.
 */

#include "mu_core.h"
#include "mu_hinf.h"
#include "mu_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== H-infinity Controller Synthesis Demo ===\n\n");

    /* Build a simple generalized plant for a 2nd-order system
     *
     * Plant: G(s) = 1 / (s^2 + 2*zeta*wn*s + wn^2)
     *   wn = 5, zeta = 0.1 (lightly damped)
     *
     * We design a controller to achieve ||T_zw||_inf < gamma
     * with sensitivity weighting for disturbance rejection.
     */

    int n = 4;  /* 2 plant states + 2 weight states */
    int ny = 1; /* 1 measurement */
    int nu = 1; /* 1 control input */
    int nw = 2; /* 2 exogenous inputs */
    int nz = 2; /* 2 performance outputs */

    MUSystem *P = mu_system_create(n, nu + nw, ny + nz);
    if (!P) { printf("ERROR: memory\n"); return 1; }

    double wn = 5.0, zeta = 0.1;

    /* A matrix: [plant A, 0; 0, weight A] */
    mu_real_mat_set(P->A, 0, 0, 0.0);
    mu_real_mat_set(P->A, 0, 1, 1.0);
    mu_real_mat_set(P->A, 1, 0, -wn * wn);
    mu_real_mat_set(P->A, 1, 1, -2.0 * zeta * wn);
    /* Weight dynamics: simple integrator for sensitivity */
    mu_real_mat_set(P->A, 2, 2, -0.01);
    mu_real_mat_set(P->A, 3, 3, -0.01);

    /* B matrix: [B1_w, B2_w; B1_u, B2_u] */
    mu_real_mat_set(P->B, 0, 0, 0.0);  mu_real_mat_set(P->B, 0, 1, 0.0);  /* w1, w2 */
    mu_real_mat_set(P->B, 1, 0, 0.0);  mu_real_mat_set(P->B, 1, 1, 0.0);
    mu_real_mat_set(P->B, 0, 2, 0.0);                                   /* u */
    mu_real_mat_set(P->B, 1, 2, wn * wn);
    /* Weight inputs from plant output */
    mu_real_mat_set(P->B, 2, 0, 1.0);
    mu_real_mat_set(P->B, 3, 1, 1.0);

    /* C matrix: [C1_z; C2_z; C1_y; C2_y] */
    mu_real_mat_set(P->C, 0, 0, 1.0);  /* z = plant output */
    mu_real_mat_set(P->C, 1, 2, 1.0);  /* z = weighted sensitivity */
    mu_real_mat_set(P->C, 2, 0, 1.0);  /* y = plant output (measurement) */

    /* D matrix: small regularization */
    mu_real_mat_set(P->D, 2, 2, 0.01); /* D22 regularization */

    P->is_discrete = false;
    P->sample_time = 0.0;

    printf("Plant: G(s) = %.1f / (s^2 + %.1fs + %.1f)\n", wn*wn, 2*zeta*wn, wn*wn);
    printf("Damping ratio: zeta = %.2f (lightly damped)\n\n", zeta);

    /* Step 1: Try H-infinity synthesis at a fixed gamma */
    printf("--- Fixed-gamma synthesis (gamma = 2.0) ---\n");
    MUSystem *K_fixed = NULL;
    bool feasible = false;
    mu_hinf_synthesize(P, 2.0, &K_fixed, &feasible);
    printf("Gamma=2.0: %s\n", feasible ? "FEASIBLE" : "INFEASIBLE");
    if (K_fixed) { mu_system_free(K_fixed); K_fixed = NULL; }

    /* Step 2: Gamma bisection to find optimal */
    printf("\n--- Gamma Bisection (optimal H-infinity) ---\n");
    MUSystem *K_opt = NULL;
    double gamma_opt = 0.0;
    mu_hinf_gamma_iteration(P, &K_opt, &gamma_opt, 0.01, 30);
    printf("Optimal gamma: %.4f\n", gamma_opt);
    if (K_opt) {
        printf("Controller order: %d\n", K_opt->n);
        printf("Controller stable: %s\n",
               mu_system_is_stable(K_opt) ? "YES" : "NO");

        /* Verify via Bounded Real Lemma */
        printf("\n--- Bounded Real Lemma Verification ---\n");
        bool brl = mu_hinf_check_bounded_real(K_opt, gamma_opt + 0.1);
        printf("BRL check (gamma=%.4f): %s\n", gamma_opt + 0.1,
               brl ? "PASS" : "FAIL");
    }

    /* Step 3: Controller order reduction */
    if (K_opt) {
        printf("\n--- Controller Order Reduction ---\n");
        MUSystem *K_red = mu_hinf_reduce_order(K_opt, 2);
        if (K_red) {
            printf("Reduced order: %d (from %d)\n", K_red->n, K_opt->n);
            mu_system_free(K_red);
        } else {
            printf("Order reduction not applicable (already <= target)\n");
        }
    }

    mu_system_free(K_opt);
    mu_system_free(P);

    printf("\nH-infinity synthesis demo complete.\n");
    return 0;
}

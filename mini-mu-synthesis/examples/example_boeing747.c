/* example_boeing747.c — Boeing 747 Lateral u-Analysis
 *
 * Demonstrates u-synthesis workflow on a real aircraft model:
 *   1. Build Boeing 747 lateral-directional model
 *   2. Set up uncertainty structure (real scalars + unmodeled dynamics)
 *   3. Compute u across frequency grid
 *   4. Assess robust stability margin
 *
 * Reference: Boeing Commercial Airplane Group lateral-directional model
 *            NASA TP-1990
 */

#include "mu_core.h"
#include "mu_uncertainty.h"
#include "mu_structured_svd.h"
#include "mu_lft.h"
#include "mu_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Boeing 747 Lateral-Directional u-Analysis ===\n\n");

    /* Step 1: Build nominal model */
    MUSystem *G = mu_unc_boeing747_model();
    if (!G) { printf("ERROR: Failed to create Boeing 747 model\n"); return 1; }
    printf("Boeing 747 model: %d states, %d inputs, %d outputs\n",
           G->n, G->m, G->p);
    printf("States: sideslip(beta), roll-rate(p), yaw-rate(r), roll-angle(phi)\n");
    printf("Inputs: aileron(delta_a), rudder(delta_r)\n\n");

    /* Step 2: Define uncertainty structure */
    MUUncertaintyStructure *delta = mu_unc_boeing747_structure();
    if (!delta) { printf("ERROR: Failed to create uncertainty structure\n"); mu_system_free(G); return 1; }
    printf("Uncertainty structure: %d blocks, total dim %d\n",
           delta->num_blocks, delta->total_dim);
    printf("  Block 0: Real scalar — roll-rate damping L_p (20%% variation)\n");
    printf("  Block 1: Real scalar — yaw-rate damping N_r (25%% variation)\n");
    printf("  Block 2: Complex full 2x2 — unmodeled actuator dynamics\n\n");

    /* Step 3: Check nominal stability */
    printf("Nominal stability check... ");
    bool stable = mu_system_is_stable(G);
    printf("%s\n", stable ? "STABLE" : "UNSTABLE");

    /* Step 4: Frequency grid for analysis */
    MUFrequencyGrid *grid = mu_frequency_grid_create(0.01, 100.0, 30);
    if (!grid) { printf("ERROR: Failed to create frequency grid\n"); mu_unc_struct_free(delta); mu_system_free(G); return 1; }
    printf("Frequency grid: %.2f to %.2f rad/s, %d points\n",
           grid->w_min, grid->w_max, grid->num_points);

    /* Step 5: Compute u across frequency (peak-finding) */
    printf("\nu-analysis across frequency:\n");
    printf("  omega (rad/s)    mu_upper    mu_lower     gap\n");
    printf("  ------------    --------    --------    ------\n");

    MUMuResult *results = (MUMuResult*)calloc((size_t)grid->num_points,
        sizeof(MUMuResult));
    double peak_mu = 0.0;
    double peak_freq = 0.0;

    for (int i = 0; i < grid->num_points; i++) {
        double w = grid->omega[i];
        MUMatrix *Mjw = mu_system_freqresp(G, w);
        if (!Mjw) continue;

        /* Build M for u analysis */
        MUMatrix *M = mu_system_to_m(G, w, delta->total_dim, delta->total_dim);
        if (!M) { mu_mat_free(Mjw); continue; }

        results[i] = mu_compute(M, delta, w, 30, 30, MU_TOL);
        double gap = results[i].mu_upper - results[i].mu_lower;

        printf("  %8.4f      %8.4f    %8.4f    %8.4f\n",
               w, results[i].mu_upper, results[i].mu_lower, gap);

        if (results[i].mu_upper > peak_mu) {
            peak_mu = results[i].mu_upper;
            peak_freq = w;
        }

        mu_mat_free(Mjw); mu_mat_free(M);
    }

    /* Step 6: Summary */
    printf("\n=== Analysis Summary ===\n");
    printf("Peak u value:     %.4f at omega = %.4f rad/s\n", peak_mu, peak_freq);
    printf("Robust stability threshold: u_peak < 1.0 required\n");

    if (peak_mu < 1.0) {
        printf("Robustness margin: %.4f (system is ROBUSTLY STABLE)\n",
               1.0 / (peak_mu > MU_EPSILON ? peak_mu : MU_EPSILON));
        printf("The Boeing 747 lateral dynamics tolerate %.0f%% larger\n",
               (1.0 / peak_mu - 1.0) * 100.0);
        printf("uncertainty than modeled before instability occurs.\n");
    } else {
        printf("WARNING: Peak u = %.4f >= 1.0 — robust stability NOT guaranteed.\n", peak_mu);
        printf("The critical frequency is %.4f rad/s. Consider adding\n", peak_freq);
        printf("a yaw damper or redesigning the roll-rate feedback loop.\n");
    }

    /* Step 7: Structured vs unstructured comparison */
    printf("\n--- Conservatism Analysis ---\n");
    {
        MUMatrix *M_crit = mu_system_freqresp(G, peak_freq);
        if (M_crit) {
            double sn = mu_mat_spectral_norm(M_crit);
            printf("Unstructured bound (sigma_bar): %.4f\n", sn);
            printf("Structured u (upper):          %.4f\n", peak_mu);
            printf("Conservatism reduction:         %.1f%%\n",
                   (sn > MU_EPSILON) ? (sn - peak_mu) / sn * 100.0 : 0.0);
            mu_mat_free(M_crit);
        }
    }

    free(results);
    mu_frequency_grid_free(grid);
    mu_unc_struct_free(delta);
    mu_system_free(G);

    printf("\nBoeing 747 u-analysis complete.\n");
    return 0;
}

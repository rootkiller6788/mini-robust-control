/* Example 2: Uncertainty Modeling with LFTs
 *
 * Demonstrates how different uncertainty types (additive, multiplicative
 * input/output, feedback) are modeled as Linear Fractional Transformations
 * (LFTs) for mu-analysis.
 *
 * Knowledge points:
 *   L1: Upper/Lower LFT definitions
 *   L3: LFT as interconnection formalism
 *   L5: Additive, multiplicative, feedback uncertainty LFT construction
 *   L6: Robustness analysis via LFT + mu
 */

#include "ssv_core.h"
#include "ssv_lft.h"
#include "ssv_bounds.h"
#include "ssv_uncertainty.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void) {
    printf("=== Example 2: Uncertainty Modeling with LFTs ===\n\n");

    /* Nominal plant G_nom = [2  0; 0  1] (diagonal 2x2 system) */
    int p = 2, m = 2;
    SSVComplexMatrix *G_nom = ssv_cmatrix_create(p, m);
    ssv_cmatrix_set(G_nom, 0, 0, 2.0);
    ssv_cmatrix_set(G_nom, 1, 1, 1.0);
    printf("Nominal plant G_nom:\n");
    ssv_cmatrix_print(G_nom, "G_nom");

    /* Uncertainty weights:
     * W_u: input weight, shapes how uncertainty enters the plant
     * W_y: output weight, shapes how uncertainty affects output */
    SSVComplexMatrix *W_u = ssv_cmatrix_create(p, 1);
    ssv_cmatrix_set(W_u, 0, 0, 0.2);
    ssv_cmatrix_set(W_u, 1, 0, 0.1);
    printf("\nInput uncertainty weight W_u:\n");
    ssv_cmatrix_print(W_u, "W_u");

    SSVComplexMatrix *W_y = ssv_cmatrix_eye(1);
    printf("Output uncertainty weight W_y = I_1\n");

    /* 1. Additive uncertainty: G = G_nom + W_u * Delta * W_y */
    printf("\n--- Additive Uncertainty ---\n");
    SSVComplexMatrix *M_add = ssv_lft_additive_uncertainty(G_nom, W_u, W_y);
    printf("M_add (%dx%d):\n", M_add->rows, M_add->cols);
    ssv_cmatrix_print(M_add, "M_add");

    /* 2. Input multiplicative uncertainty: G = G_nom*(I + W_u*Delta*W_y) */
    printf("\n--- Input Multiplicative Uncertainty ---\n");
    SSVComplexMatrix *M_in = ssv_lft_input_multiplicative(G_nom, W_u, W_y);
    printf("M_in (%dx%d):\n", M_in->rows, M_in->cols);
    ssv_cmatrix_print(M_in, "M_in");

    /* 3. Output multiplicative uncertainty: G = (I + W_u*Delta*W_y)*G_nom */
    printf("\n--- Output Multiplicative Uncertainty ---\n");
    SSVComplexMatrix *M_out = ssv_lft_output_multiplicative(G_nom, W_u, W_y);
    printf("M_out (%dx%d):\n", M_out->rows, M_out->cols);
    ssv_cmatrix_print(M_out, "M_out");

    /* 4. Feedback uncertainty: G = G_nom*(I + W_u*Delta*W_y)^{-1} */
    printf("\n--- Feedback Uncertainty ---\n");
    SSVComplexMatrix *M_fb = ssv_lft_feedback_uncertainty(G_nom, W_u, W_y);
    printf("M_fb (%dx%d):\n", M_fb->rows, M_fb->cols);
    ssv_cmatrix_print(M_fb, "M_fb");

    /* mu analysis for each uncertainty type */
    SSVUncertaintyStructure *u = ssv_ustruct_one_full_block(1);
    printf("\n--- mu Analysis Comparison ---\n");

    SSVDScaleOptions opts = ssv_dscale_options_default();
    opts.max_iterations = 30;

    SSVComplexMatrix *models[] = {M_add, M_in, M_out, M_fb};
    const char *names[] = {"Additive", "Input Mult", "Output Mult", "Feedback"};

    for (int i = 0; i < 4; i++) {
        /* Extract M11 block (connects to Delta) */
        int blk_size = models[i]->cols - m;
        if (blk_size > 0) {
            /* The upper-left block is M11 for mu analysis */
            SSVComplexMatrix *M11 = ssv_cmatrix_create(blk_size, blk_size);
            for (int r = 0; r < blk_size; r++)
                for (int c = 0; c < blk_size; c++)
                    M11->data[c * blk_size + r] = models[i]->data[c * models[i]->rows + r];

            double mu_ub = ssv_mu_upper_bound(M11, u, &opts, NULL);
            printf("  %-15s: mu_ub = %.6f\n", names[i], mu_ub);
            ssv_cmatrix_free(M11);
        }
    }

    /* Demonstrate LFT well-posedness check:
     * Compute Fu(M, Delta) and verify well-posedness */
    printf("\n--- LFT Well-Posedness Demo ---\n");
    SSVLFTMatrix *lft = ssv_lft_create(M_in, 1, 1, p, m);
    SSVComplexMatrix *Delta = ssv_cmatrix_create(1, 1);
    ssv_cmatrix_set(Delta, 0, 0, 0.5);

    bool wp = ssv_lft_upper_is_wellposed(lft, Delta, 1e-8);
    printf("Upper LFT well-posed (delta=0.5): %s\n", wp ? "YES" : "NO");

    SSVComplexMatrix *Fu = ssv_lft_upper(lft, Delta);
    if (Fu) {
        printf("Fu(M, Delta) result (%dx%d):\n", Fu->rows, Fu->cols);
        ssv_cmatrix_print(Fu, "Fu");
    }

    /* Cleanup */
    ssv_cmatrix_free(Fu);
    ssv_cmatrix_free(Delta);
    ssv_lft_free(lft);
    ssv_cmatrix_free(G_nom);
    ssv_cmatrix_free(W_u);
    ssv_cmatrix_free(W_y);
    ssv_cmatrix_free(M_add);
    ssv_cmatrix_free(M_in);
    ssv_cmatrix_free(M_out);
    ssv_cmatrix_free(M_fb);
    ssv_ustruct_free(u);

    printf("\n=== Example 2 complete ===\n");
    return 0;
}

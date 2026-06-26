/* ex2_msd.c — Mass-Spring-Damper Robust Controller Design
 *
 * Demonstrates:
 *   1. Building uncertain MSD model
 *   2. Robust pole placement controller design
 *   3. Closed-loop eigenvalue analysis under uncertainty
 */
#include "pu_app.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    printf("=== Example 2: MSD Robust Controller Design ===\n\n");
    double m=1.0, c=0.5, k=10.0;
    double dm=0.25, dc=0.25, dk=0.25; /* 25% uncertainty */
    printf("MSD: m=%.2f +/- %.0f%%, c=%.2f +/- %.0f%%, k=%.2f +/- %.0f%%\n",
           m, dm*100, c, dc*100, k, dk*100);

    PU_UncertainStateSpace msd = pu_msd_model(m, c, k, dm, dc, dk);
    printf("Open-loop eigenvalues:\n");
    double *re=(double*)malloc(2*sizeof(double)),*im=(double*)malloc(2*sizeof(double));
    pu_eigen_qr(msd.A_nominal, 2, re, im, 1000, 1e-10);
    for (int i = 0; i < 2; i++) printf("  %.4f %+.4fj\n", re[i], im[i]);

    double wn=5.0, xi=0.7, K[2];
    printf("\nDesign: wn=%.1f, xi=%.1f\n", wn, xi);
    PU_StabilityStatus s = pu_msd_robust_controller(&msd, wn, xi, K);
    printf("Feedback gains: K=[%.4f, %.4f]\n", K[0], K[1]);
    printf("Robust stability: %s\n", s==PU_STABLE?"STABLE":"NOT STABLE");

    double sf;
    pu_monte_carlo_stability(&msd, 500, &sf);
    printf("MonteCarlo stable fraction: %.1f%%\n", sf*100.0);

    free(re); free(im); pu_uss_free(&msd);
    printf("=== End Example 2 ===\n");
    return 0;
}

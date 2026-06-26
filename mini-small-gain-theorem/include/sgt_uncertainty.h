#ifndef SGT_UNCERTAINTY_H
#define SGT_UNCERTAINTY_H

#include "sgt_core.h"

/* ============================================================================
 * Small Gain Theorem — Uncertainty Modeling and Robustness Analysis
 *
 * Robust control theory models the discrepancy between a nominal plant
 * model and the true physical system as "uncertainty". The small gain
 * theorem provides the fundamental stability condition for feedback
 * systems with norm-bounded uncertainty.
 *
 * Key concepts:
 *   1. Uncertainty modeling (additive, multiplicative, coprime factor)
 *   2. Robust stability condition: ||M * Delta||_inf < 1
 *   3. Structured singular value (mu) for structured uncertainty
 *   4. D-K iteration for mu-synthesis
 *
 * Reference:
 *   K. Zhou, J.C. Doyle, K. Glover, "Robust and Optimal Control," 1996.
 *   A. Packard, J. Doyle, "The complex structured singular value,"
 *     Automatica, 1993.
 * ============================================================================ */

/* --- Uncertainty Block Construction --- */

SGTUncertaintyBlock* sgt_uncertainty_create_full(int dim, double bound);
SGTUncertaintyBlock* sgt_uncertainty_create_repeated_scalar(int dim, double bound);
SGTUncertaintyBlock* sgt_uncertainty_create_dynamic(const SGTLTISystem *weight, double bound);
void sgt_uncertainty_block_free(SGTUncertaintyBlock *block);

/* --- Structured Uncertainty --- */

SGTStructuredUncertainty* sgt_structured_uncertainty_create(const SGTUncertaintyBlock **blocks, int n_blocks);
void sgt_structured_uncertainty_free(SGTStructuredUncertainty *delta);

/* --- Robust Stability Analysis --- */

double sgt_robust_stability_unstructured(const SGTLTISystem *M, double gamma_uncertainty);
double sgt_mu_lower_bound(const SGTLTISystem *M, const SGTStructuredUncertainty *delta, int max_iter, double tol);
double sgt_mu_upper_bound(const SGTLTISystem *M, const SGTStructuredUncertainty *delta);
double sgt_mu_upper_bound_freq(const SGTMatrix *M_real, const SGTMatrix *M_imag, const SGTStructuredUncertainty *delta, int max_iter, double tol);

/* --- D-K Iteration --- */

void sgt_dk_iteration_d_step(const SGTLTISystem *M, SGTStructuredUncertainty *delta, int n_omega, int max_iter);
SGTLTISystem* sgt_fit_dscale_tf(const double *omega, const double *d_mag, int n_omega, int order);

/* --- Uncertainty Modeling --- */

SGTLTISystem* sgt_multiplicative_weight(const SGTLTISystem **plants, int n_plants, int nominal_idx, int n_omega);
SGTLTISystem* sgt_additive_weight(const SGTLTISystem **plants, int n_plants, int nominal_idx, int n_omega);

/* --- IQC Framework Connection --- */

SGTMatrix* sgt_iqc_multiplier_small_gain(double gamma, int signal_dim);
SGTMatrix* sgt_iqc_multiplier_passivity(int signal_dim);
bool sgt_iqc_stability_check(const SGTLTISystem *G, const SGTMatrix *Pi, int n_omega);

#endif /* SGT_UNCERTAINTY_H */
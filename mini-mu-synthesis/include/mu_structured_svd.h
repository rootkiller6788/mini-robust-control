#ifndef MU_STRUCTURED_SVD_H
#define MU_STRUCTURED_SVD_H

#include "mu_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ====================================================================
 * mu_structured_svd.h — Structured Singular Value (u) Computation
 *
 * The structured singular value u_Delta(M) measures the "size" of the
 * smallest structured perturbation Delta that makes I - M*Delta singular.
 *
 *   u_Delta(M) = 1 / min{ sigma_bar(Delta) : Delta in Delta, det(I - M*Delta) = 0 }
 *
 * Key property: u_Delta(alpha * M) = |alpha| * u_Delta(M)
 *
 * This module implements:
 *   1. Upper bound via D-scale optimization (convex, LMI-based)
 *   2. Lower bound via power iteration (non-convex, local convergence)
 *   3. Osborne balancing for numerical stability
 *
 * Reference: Packard & Doyle (1993), "The complex structured singular value"
 *            Young & Doyle (1990), "Computation of u with real and complex uncertainties"
 * ==================================================================== */

/*
 * L4: Small-u Theorem (Robust Stability)
 *
 * Theorem: Let M(s) be a stable transfer matrix and Delta be a structured
 * uncertainty set with ||Delta||_inf <= 1. The closed-loop system
 * F_u(M, Delta) is robustly stable iff
 *
 *   sup_{omega in R} u_Delta(M(j omega)) < 1
 *
 * This is the generalization of the small-gain theorem to structured uncertainty.
 * Proof: Zhou, Doyle & Glover (1996), Theorem 11.7
 */
#define MU_ROBUST_STABILITY_THRESHOLD  1.0

/*
 * L5: u Upper Bound Computation (D-Scale Method)
 *
 * u_Delta(M) <= inf_{D in D_Delta} sigma_bar(D M D^{-1})
 *
 * This is a convex optimization problem in D. At each frequency,
 * we minimize sigma_bar(D M D^{-1}) over D > 0, D in D_Delta.
 *
 * The optimization alternates between:
 *   1. Fix D, compute sigma_bar(D M D^{-1}) via SVD
 *   2. Update D to reduce the spectral norm
 *
 * Convergence is guaranteed but potentially slow.
 *
 * @param M         Complex matrix to analyze
 * @param delta     Uncertainty structure defining D_Delta
 * @param max_iter  Maximum iterations for D-optimization
 * @param tol       Convergence tolerance
 * @return          Upper bound on u_Delta(M)
 *
 * Complexity: O(dim^3 * max_iter)
 */
double mu_upper_bound(const MUMatrix *M,
                       const MUUncertaintyStructure *delta,
                       int max_iter, double tol);

/*
 * L5: u Lower Bound (Power Iteration)
 *
 * u_Delta(M) >= max_{Q in Q_Delta} rho(M Q)
 *
 * where Q_Delta = { Q in Delta : sigma_bar(Q) <= 1 }
 * and rho(.) is the spectral radius.
 *
 * Power iteration: iteratively find the worst-case Q that maximizes
 * the spectral radius of M Q, providing a lower bound on u.
 *
 * This is a non-convex optimization; local maxima are common.
 * Multiple random restarts improve the chance of finding the global max.
 *
 * @param M         Complex matrix to analyze
 * @param delta     Uncertainty structure
 * @param max_iter  Maximum power iterations
 * @param tol       Convergence tolerance
 * @return          Lower bound on u_Delta(M) (0 if computation fails)
 *
 * Complexity: O(dim^3 * max_iter)
 * Reference: Packard & Doyle (1993), §4.3
 *            Young & Doyle (1990)
 */
double mu_lower_bound(const MUMatrix *M,
                       const MUUncertaintyStructure *delta,
                       int max_iter, double tol);

/*
 * L5: u Computation with Both Bounds
 *
 * Compute both upper and lower bounds for u at a given frequency.
 * The true u lies between these bounds. A small gap indicates
 * an accurate result (the gap is zero for certain structures).
 *
 * Complexity: O(dim^3 * (upper_iter + lower_iter))
 */
MUMuResult mu_compute(const MUMatrix *M,
                       const MUUncertaintyStructure *delta,
                       double frequency,
                       int upper_iter, int lower_iter, double tol);

/*
 * L5: u Analysis Across Frequency Grid
 *
 * Compute u across a grid of frequencies. This is the core analysis
 * step: evaluate u_Delta(M(j omega)) at each omega and find the peak.
 *
 * The peak u value determines the robust stability/performance margin.
 *
 * @param M_sys     State-space system to analyze
 * @param delta     Uncertainty structure
 * @param grid      Frequency grid
 * @param results   Output array of MUMuResult (pre-allocated, grid->num_points)
 * @return          Peak u value across all frequencies
 *
 * Complexity: O(#freq * dim^3 * max_iter)
 */
double mu_analysis_grid(const MUSystem *M_sys,
                         const MUUncertaintyStructure *delta,
                         const MUFrequencyGrid *grid,
                         MUMuResult *results);

/*
 * L4: Robust Stability Test
 *
 * Check if the closed-loop system is robustly stable by verifying
 * that the peak u value is less than 1 across all frequencies.
 *
 * Theorem (Robust Stability via u): The interconnection F_u(M, Delta)
 * is robustly stable for all ||Delta||_inf <= 1 iff:
 *   sup_{omega} u_Delta(M(j omega)) < 1
 */
bool mu_is_robustly_stable(const MUSystem *M_sys,
                            const MUUncertaintyStructure *delta,
                            const MUFrequencyGrid *grid,
                            int upper_iter, int lower_iter, double tol);

/*
 * L4: Robust Performance Test
 *
 * Theorem (Main Loop Theorem): The interconnection achieves robust
 * performance (RP) iff the augmented u (with a performance block)
 * is less than 1 across all frequencies.
 *
 * RP <=> u_{Delta_a}(M(j omega)) < 1 for all omega
 *
 * where Delta_a = diag(Delta, Delta_p) with an additional
 * "full" complex performance block Delta_p.
 */
bool mu_is_robust_performance(const MUSystem *M_sys,
                               const MUUncertaintyStructure *delta,
                               const MUFrequencyGrid *grid,
                               int upper_iter, int lower_iter, double tol);

/*
 * L8: Mixed-u Upper Bound
 *
 * For systems with mixed real/complex uncertainties, the pure complex-u
 * upper bound can be conservative. Mixed-u uses a combination of D-scales
 * (for complex blocks) and G-scales (for real blocks):
 *
 *   u_Delta(M) <= inf_{D,G} { alpha : M* D M - j(G M - M* G) <= alpha^2 D }
 *
 * This LMI condition provides a tighter upper bound for real uncertainties.
 *
 * Complexity: O(dim^6) (LMI solution via interior-point)
 * Reference: Fan, Tits & Doyle (1991), "Robustness in the presence of
 *            mixed parametric uncertainty and unmodeled dynamics"
 */
double mu_mixed_upper_bound(const MUMatrix *M,
                             const MUUncertaintyStructure *delta,
                             int max_iter, double tol);

#ifdef __cplusplus
}
#endif

#endif /* MU_STRUCTURED_SVD_H */

#ifndef MU_ROBUST_ANALYSIS_H
#define MU_ROBUST_ANALYSIS_H

#include "mu_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ====================================================================
 * mu_robust_analysis.h — Robust Stability and Performance Analysis
 *
 * Four key properties for a closed-loop system with uncertainty:
 *
 *   NS: Nominal Stability — system stable with Delta = 0
 *   NP: Nominal Performance — performance specs met with Delta = 0
 *   RS: Robust Stability — system stable for all ||Delta||_inf <= 1
 *   RP: Robust Performance — specs met for all ||Delta||_inf <= 1
 *
 * Relationships (Skogestad & Postlethwaite 2005, Table 8.1):
 *   RP => RS, NP
 *   RS and NP does NOT imply RP
 *
 *   NS <=> M internally stable
 *   NP <=> sigma_bar(M22) < 1 for all omega
 *   RS <=> u_Delta(M11) < 1 for all omega
 *   RP <=> u_{Delta_a}(M) < 1 for all omega (Main Loop Theorem)
 *
 * Reference: Skogestad & Postlethwaite (2005), Ch. 8
 *            Zhou, Doyle & Glover (1996), Ch. 11
 * ==================================================================== */

/*
 * L2: Nominal Stability Analysis
 *
 * Check if the closed-loop system M = F_l(P, K) is internally stable.
 * This requires:
 *   a) All eigenvalues of A_cl have Re(s) < 0 (continuous)
 *   b) No pole-zero cancellations in the RHP
 *
 * Complexity: O(n^3)
 */
bool mu_check_nominal_stability(const MUSystem *P, const MUSystem *K);

/*
 * L2: Nominal Performance Analysis
 *
 * Check if the nominal closed-loop satisfies performance specs:
 *   sigma_bar(M22(j omega)) < 1 for all omega
 *
 * This is an H-infinity norm check on the nominal performance channel.
 *
 * Complexity: O(#freq * n^3)
 */
bool mu_check_nominal_performance(const MUSystem *P, const MUSystem *K,
                                   const MUFrequencyGrid *grid);

/*
 * L4: Robust Stability via u-Analysis
 *
 * Theorem: The system F_u(M, Delta) is robustly stable for all
 * ||Delta||_inf <= 1 iff u_Delta(M11(j omega)) < 1 for all omega.
 *
 * This is the generalization of Nyquist to structured uncertainty.
 */
bool mu_check_robust_stability(const MUSystem *P, const MUSystem *K,
                                const MUUncertaintyStructure *delta,
                                const MUFrequencyGrid *grid,
                                int upper_iter, int lower_iter, double tol);

/*
 * L4: Robust Performance via Augmented u (Main Loop Theorem)
 *
 * Theorem (Main Loop Theorem): The system achieves robust performance
 * for all ||Delta||_inf <= 1 iff:
 *
 *   mu_{Delta_a}(M(j omega)) < 1 for all omega
 *
 * where Delta_a = [ Delta     0    ]
 *                  [   0    Delta_p ]
 *
 * with Delta_p being a full complex block of appropriate dimensions
 * (representing the performance specification).
 *
 * Reference: Packard & Doyle (1993), Main Loop Theorem
 *            Zhou et al. (1996), Theorem 11.9
 */
bool mu_check_robust_performance(const MUSystem *P, const MUSystem *K,
                                  const MUUncertaintyStructure *delta,
                                  const MUFrequencyGrid *grid,
                                  int upper_iter, int lower_iter, double tol);

/*
 * L4: Robustness Margin
 *
 * Compute the structured robustness margin: the largest k such that
 * the system remains stable for all ||Delta||_inf <= k.
 *
 *   k_max = 1 / sup_{omega} u_Delta(M11(j omega))
 *
 * A margin > 1 means the system tolerates larger uncertainty than designed.
 * A margin < 1 means the system fails to meet the robustness spec.
 *
 * Complexity: O(#freq * dim^3 * max_iter)
 */
double mu_robustness_margin(const MUSystem *P, const MUSystem *K,
                             const MUUncertaintyStructure *delta,
                             const MUFrequencyGrid *grid,
                             int upper_iter, int lower_iter, double tol);

/*
 * L4: Worst-Case Uncertainty
 *
 * Find the worst-case perturbation Delta_wc (within the structure)
 * that causes instability at a given frequency.
 *
 * This is found by extracting the principal directions from the
 * power iteration that produces the u lower bound.
 *
 * @param M           Transfer matrix at the critical frequency
 * @param delta       Uncertainty structure
 * @param Delta_wc    Output: worst-case perturbation
 *
 * Complexity: O(dim^3 * max_iter)
 * Reference: Young & Doyle (1990)
 */
void mu_worst_case_perturbation(const MUMatrix *M,
                                 const MUUncertaintyStructure *delta,
                                 MUMatrix **Delta_wc);

/*
 * L2: Structured vs Unstructured Comparison
 *
 * Compare the structured singular value u with the unstructured
 * bound (spectral norm). The ratio u / sigma_bar(M) shows how much
 * conservatism is removed by accounting for uncertainty structure.
 *
 * sigma_bar(M) >= u(M) >= rho(M)   (where rho is spectral radius)
 *
 * For block-diagonal uncertainty with S repeated scalars and F full blocks:
 *   u(M) <= sigma_bar(M), with equality when Delta is a single full block.
 *
 * Complexity: O(dim^3)
 */
typedef struct {
    double spectral_norm;     /* sigma_bar(M) -- unstructured bound */
    double spectral_radius;   /* rho(M) -- theoretical lower bound */
    double mu_upper;          /* u upper bound */
    double mu_lower;          /* u lower bound */
    double conservatism;      /* (sigma_bar - mu_lower) / sigma_bar */
} MUComparisonResult;

MUComparisonResult mu_compare_structured_unstructured(
    const MUMatrix *M, const MUUncertaintyStructure *delta,
    int upper_iter, int lower_iter, double tol);

#ifdef __cplusplus
}
#endif

#endif /* MU_ROBUST_ANALYSIS_H */

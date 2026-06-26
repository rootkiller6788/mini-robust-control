#ifndef SSV_BOUNDS_H
#define SSV_BOUNDS_H

#include "ssv_core.h"
#include "ssv_uncertainty.h"
#include <complex.h>
#include <stdbool.h>

/* ============================================================================
 * Structured Singular Value — Upper and Lower Bounds
 *
 * Computing mu exactly is NP-hard for general uncertainty structures.
 * In practice, upper and lower bounds are computed:
 *
 * Upper bound (convex):
 *   mu_Delta(M) <= inf_{D in D} sigma_max(D * M * D^{-1})
 *   This is a convex optimization (LMI) and provides a guaranteed bound.
 *
 * Lower bound (non-convex):
 *   mu_Delta(M) = max_{Q in Q} rho(Q * M)
 *   where rho is the spectral radius and Q is the set of unitary block-diagonal
 *   matrices commuting with Delta. This is a non-convex optimization.
 *
 *   For purely complex uncertainty with 3 or fewer blocks, the upper bound
 *   is exact: mu = inf_D sigma_max(D M D^{-1}).
 *
 * References:
 *   Packard & Doyle, "The complex structured singular value" (Automatica, 1993)
 *   Fan & Tits, "A measure of worst case H∞ performance" (1986)
 *   Young & Doyle, "A lower bound for the mixed mu problem" (1990)
 * ============================================================================ */

/* --- Upper Bound via D-Scaling --- */

/** Options for the D-scaling upper bound optimization. */
typedef struct {
    int max_iterations;        /* maximum number of iterations (default 100) */
    double tolerance;          /* convergence tolerance (default 1e-6) */
    bool use_osborne_init;     /* use Osborne balancing for initial D (default true) */
    double step_size;          /* step size for D update (0 < step <= 1, default 0.5) */
    bool verbose;              /* print iteration progress */
} SSVDScaleOptions;

/** Default D-scaling options. */
SSVDScaleOptions ssv_dscale_options_default(void);

/** Compute the mu upper bound via D-scaling optimization.
 *
 *  Algorithm:
 *    1. Initialize D = I
 *    2. Compute M_scaled = D * M * D^{-1}
 *    3. Compute SVD: M_scaled = U * Sigma * V^H
 *    4. Update D to reduce sigma_max(D M D^{-1})
 *    5. Iterate until convergence
 *
 *  This implements the Osborne-like iteration for the optimal D-scaling,
 *  which finds a local minimum of the convex upper bound.
 *
 *  @param M complex matrix to analyze
 *  @param ustruct uncertainty structure (defines D set)
 *  @param options optimization options
 *  @param out_dscale (optional, may be NULL) output optimal D-scaling
 *  @return mu upper bound
 */
double ssv_mu_upper_bound(const SSVComplexMatrix *M,
                           const SSVUncertaintyStructure *ustruct,
                           const SSVDScaleOptions *options,
                           SSVDScaleMatrices **out_dscale);

/** Compute the mu upper bound at multiple frequency points.
 *  This is the standard mu-analysis: robust stability margin across frequency.
 *
 *  @param M_array array of complex matrices, one per frequency
 *  @param frequencies array of frequency values in rad/s
 *  @param n_freq number of frequency points
 *  @param ustruct uncertainty structure
 *  @param options optimization options
 *  @return array of SSVMuResult (length n_freq), caller must free each
 */
SSVMuResult** ssv_mu_upper_bound_frequency(
    const SSVComplexMatrix **M_array,
    const double *frequencies,
    int n_freq,
    const SSVUncertaintyStructure *ustruct,
    const SSVDScaleOptions *options);

/* --- Lower Bound via Power Iteration --- */

/** Options for the power iteration lower bound. */
typedef struct {
    int max_iterations;        /* maximum power iterations (default 500) */
    double tolerance;          /* convergence tolerance (default 1e-8) */
    int n_random_starts;       /* number of random restarts (default 10) */
    bool use_complex_perturb;  /* use complex perturbation for real problems */
    bool verbose;              /* print iteration progress */
} SSVPowerIterOptions;

/** Default power iteration options. */
SSVPowerIterOptions ssv_power_iter_options_default(void);

/** Compute the mu lower bound via power iteration.
 *
 *  The power iteration (Packard & Doyle, 1993) computes a local maximum of
 *  rho(Q * M) over block-diagonal unitary Q. At each fixed point,
 *  |Q^H M a|_i provides a lower bound.
 *
 *  Algorithm:
 *    For each random start:
 *      1. Initialize b(0) randomly, a(0) = 0, theta(0) = 0
 *      2. Iterate:
 *         a. a(k+1) = M * (theta(k) * b(k)), normalize per-block
 *         b. b(k+1) = M^H * a(k+1), normalize per-block
 *         c. theta(k+1) = (b(k+1))^H * M^H * a(k+1) / |...|
 *         d. mu_lb_k = |theta(k+1)|
 *      3. Keep best lower bound across starts
 *
 *  @param M complex matrix to analyze
 *  @param ustruct uncertainty structure
 *  @param options power iteration options
 *  @return mu lower bound (0 <= mu_lb <= mu_true)
 */
double ssv_mu_lower_bound(const SSVComplexMatrix *M,
                           const SSVUncertaintyStructure *ustruct,
                           const SSVPowerIterOptions *options);

/** Compute both upper and lower bounds for mu, returning the full result.
 *  This is the main entry point for mu-analysis of a single frequency.
 */
SSVMuResult* ssv_mu_analysis(const SSVComplexMatrix *M,
                              const SSVUncertaintyStructure *ustruct,
                              const SSVDScaleOptions *ub_options,
                              const SSVPowerIterOptions *lb_options);

/** Compute mu analysis across a frequency grid.
 *  @return array of SSVMuResult (length n_freq), caller must free each
 */
SSVMuResult** ssv_mu_analysis_frequency(
    const SSVComplexMatrix **M_array,
    const double *frequencies,
    int n_freq,
    const SSVUncertaintyStructure *ustruct,
    const SSVDScaleOptions *ub_options,
    const SSVPowerIterOptions *lb_options);

/* --- Robust Stability and Performance Tests --- */

/** Test robust stability: mu(M) < 1/gamma?
 *
 *  Robust stability condition (Main Loop Theorem):
 *    The perturbed system is robustly stable for all ||Delta||_inf <= gamma
 *    iff mu_Delta(M(jw)) < 1/gamma for all frequencies w.
 *
 *  @param M nominal closed-loop transfer matrix at a frequency
 *  @param ustruct uncertainty structure
 *  @param gamma uncertainty level (typically 1.0)
 *  @return true if robustly stable at this frequency
 */
bool ssv_robust_stability_test(const SSVComplexMatrix *M,
                                const SSVUncertaintyStructure *ustruct,
                                double gamma);

/** Test robust performance: Is the worst-case performance level <= beta?
 *
 *  Robust performance is tested by augmenting M with a performance block
 *  Delta_p (full complex) and testing mu_augmented(M_aug) < 1.
 *
 *  Main Loop Theorem for performance:
 *    The system achieves robust performance level <= 1/gamma
 *    iff mu_{Delta_aug}(M_aug(jw)) < 1/gamma for all w,
 *    where Delta_aug = diag(Delta, Delta_p).
 *
 *  @param M nominal closed-loop (with performance channels)
 *  @param ustruct uncertainty structure (without performance block)
 *  @param n_perf number of performance channels (size of Delta_p)
 *  @return true if robust performance condition holds
 */
bool ssv_robust_performance_test(const SSVComplexMatrix *M,
                                  const SSVUncertaintyStructure *ustruct,
                                  int n_perf);

/** Compute the robust stability margin: maximum gamma such that
 *  stability is guaranteed for all ||Delta||_inf <= gamma.
 *
 *  margin = 1 / max_w mu_Delta(M(jw))
 */
double ssv_robust_stability_margin(const SSVMuResult **mu_results, int n_freq);

/** Find the worst-case frequency where mu peaks.
 *  @return index of the peak frequency in mu_results
 */
int ssv_mu_peak_frequency(const SSVMuResult **mu_results, int n_freq);

/* --- Mixed mu (Real + Complex) --- */

/** Compute upper bound for mixed real/complex mu using D,G-scaling.
 *
 *  For real uncertainties, the upper bound requires both D-scaling and
 *  G-scaling (skew-symmetric). This is more conservative but necessary
 *  because real mu is discontinuous.
 *
 *  mu_{mixed}(M) <= inf_{D in D, G in G} alpha such that:
 *    M^H D M + j(G M - M^H G) <= alpha^2 D
 *
 *  @param M complex matrix
 *  @param ustruct uncertainty structure (must include real blocks)
 *  @param options D-scale options
 *  @return mixed mu upper bound
 */
double ssv_mixed_mu_upper_bound(const SSVComplexMatrix *M,
                                 const SSVUncertaintyStructure *ustruct,
                                 const SSVDScaleOptions *options);

#endif /* SSV_BOUNDS_H */

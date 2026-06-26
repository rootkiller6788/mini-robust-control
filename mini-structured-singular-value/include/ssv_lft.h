#ifndef SSV_LFT_H
#define SSV_LFT_H

#include "ssv_core.h"
#include <complex.h>
#include <stdbool.h>

/* ============================================================================
 * Linear Fractional Transformations (LFTs)
 *
 * LFTs are the fundamental building block for uncertainty modeling and
 * controller interconnection in robust control.
 *
 * Given a partitioned matrix M:
 *
 *   M = [ M11  M12 ]
 *       [ M21  M22 ]
 *
 * Upper LFT:  Fu(M, Delta_u) = M22 + M21 * Delta_u * (I - M11*Delta_u)^{-1} * M12
 *   Used to pull out uncertainty from the upper channels.
 *
 * Lower LFT:  Fl(M, Delta_l) = M11 + M12 * Delta_l * (I - M22*Delta_l)^{-1} * M21
 *   Used to close the loop with a controller on the lower channels.
 *
 * Star Product (Redheffer): S(M, P) = generalized interconnection.
 *
 * References:
 *   Zhou, Doyle & Glover — Robust and Optimal Control, Ch. 9 (1996)
 *   Redheffer, "On a certain LFT..." (1960)
 *   Packard, "Gain scheduling via LFTs" (1994)
 * ============================================================================ */

/* --- LFT Data Structure --- */

/** Partitioned matrix for LFT operations.
 *
 *  Stored as a single complex matrix with partition information:
 *    matrix = [ M11  M12 ]
 *             [ M21  M22 ]
 *
 *  where M11 is n1 x m1, M12 is n1 x m2,
 *        M21 is n2 x m1, M22 is n2 x m2.
 */
typedef struct {
    SSVComplexMatrix *M;       /* full complex matrix */
    int n1;                    /* rows/cols of M11 block */
    int m1;                    /* cols of M11 (inputs from / outputs to) */
    int n2;                    /* rows of M22 block */
    int m2;                    /* cols of M22 */
} SSVLFTMatrix;

/** Create an LFT matrix from a partitioned complex matrix.
 *  @param M full matrix (caller retains ownership, internal copy made)
 *  @param n1 number of rows/cols in M11
 *  @param m1 number of cols in M11 (= rows in Delta_u / cols in Delta_l)
 *  @param n2 number of rows in M22
 *  @param m2 number of cols in M22
 */
SSVLFTMatrix* ssv_lft_create(const SSVComplexMatrix *M, int n1, int m1, int n2, int m2);

/** Free LFT matrix. */
void ssv_lft_free(SSVLFTMatrix *lft);

/** Create an LFT from a state-space system (frequency-domain).
 *  Given G(s) = C*(sI-A)^{-1}*B + D, at frequency w, return the complex
 *  matrix that can be used in LFT operations.
 */
SSVLFTMatrix* ssv_lft_from_state_space(const SSVRealMatrix *A,
                                         const SSVRealMatrix *B,
                                         const SSVRealMatrix *C,
                                         const SSVRealMatrix *D,
                                         double w);

/* --- Upper LFT: Fu(M, Delta_u) --- */

/** Compute the upper LFT: Fu(M, Delta_u) = M22 + M21 * Delta_u * (I - M11*Delta_u)^{-1} * M12
 *
 *  This represents the system seen from below when uncertainty Delta_u
 *  is pulled out from the upper channels.
 *
 *  Well-posedness condition: I - M11*Delta_u must be invertible.
 *
 *  @param lft partitioned matrix M
 *  @param Delta_u upper uncertainty/parameter block (m1 x n1)
 *  @return Fu(M, Delta_u) as a complex matrix (n2 x m2), NULL if not well-posed
 */
SSVComplexMatrix* ssv_lft_upper(const SSVLFTMatrix *lft,
                                  const SSVComplexMatrix *Delta_u);

/** Check well-posedness of upper LFT. */
bool ssv_lft_upper_is_wellposed(const SSVLFTMatrix *lft,
                                  const SSVComplexMatrix *Delta_u,
                                  double tol);

/* --- Lower LFT: Fl(M, Delta_l) --- */

/** Compute the lower LFT: Fl(M, Delta_l) = M11 + M12 * Delta_l * (I - M22*Delta_l)^{-1} * M21
 *
 *  This represents the closed-loop system when controller Delta_l (= K)
 *  is connected in feedback from the lower channels.
 *
 *  Well-posedness condition: I - M22*Delta_l must be invertible.
 *
 *  @param lft partitioned matrix M
 *  @param Delta_l lower feedback block (m2 x n2), typically a controller
 *  @return Fl(M, Delta_l) as a complex matrix (n1 x m1), NULL if not well-posed
 */
SSVComplexMatrix* ssv_lft_lower(const SSVLFTMatrix *lft,
                                  const SSVComplexMatrix *Delta_l);

/** Check well-posedness of lower LFT. */
bool ssv_lft_lower_is_wellposed(const SSVLFTMatrix *lft,
                                  const SSVComplexMatrix *Delta_l,
                                  double tol);

/* --- Redheffer Star Product --- */

/** Compute the Redheffer star product S = M * P.
 *
 *  Given two partitioned systems, the star product interconnects them:
 *
 *    S = [ S11  S12 ]
 *        [ S21  S22 ]
 *
 *  where:
 *    S11 = M11 + M12 * P11 * (I - M22*P11)^{-1} * M21
 *    S12 = M12 * (I - P11*M22)^{-1} * P12
 *    S21 = P21 * (I - M22*P11)^{-1} * M21
 *    S22 = P22 + P21 * M22 * (I - P11*M22)^{-1} * P12
 *
 *  This corresponds to the series interconnection of two LFT systems.
 *  The Redheffer star product is the generalization of transfer function
 *  multiplication to multi-port LFT representations.
 *
 *  @param M first LFT system
 *  @param P second LFT system (must have compatible dimensions)
 *  @return star product LFT, NULL if dimensions incompatible
 */
SSVLFTMatrix* ssv_lft_star_product(const SSVLFTMatrix *M, const SSVLFTMatrix *P);

/* --- LFT Chain (Multiple Uncertainty Blocks) --- */

/** Chain multiple LFT systems in series: M1 * M2 * ... * Mn
 *  Uses the Redheffer star product repeatedly.
 *
 *  @param lfts array of LFT matrices
 *  @param n_lfts number of LFT matrices
 *  @return chained LFT, NULL on error
 */
SSVLFTMatrix* ssv_lft_chain(const SSVLFTMatrix **lfts, int n_lfts);

/* --- LFT for Uncertainty Modeling --- */

/** Create the M-Delta interconnection from a nominal model with uncertainty.
 *
 *  Given a nominal state-space model with additive uncertainty:
 *    G = G_nom + W_u * Delta * W_y
 *
 *  The LFT form is:
 *    [ z ]   [ 0      W_y    ] [ w ]
 *    [ y ] = [ W_u    G_nom  ] [ u ]
 *
 *  where w = Delta * z (uncertainty feedback).
 *
 *  @param G_nom nominal transfer matrix at frequency w
 *  @param W_u input uncertainty weight at frequency w
 *  @param W_y output uncertainty weight at frequency w
 *  @return M matrix for mu-analysis [M11 M12; M21 M22]
 */
SSVComplexMatrix* ssv_lft_additive_uncertainty(
    const SSVComplexMatrix *G_nom,
    const SSVComplexMatrix *W_u,
    const SSVComplexMatrix *W_y);

/** Create M-Delta for multiplicative input uncertainty:
 *    G = G_nom * (I + W_u * Delta * W_y)
 *
 *  M = [ -W_y*G_nom*W_u   -W_y*G_nom ]
 *      [ W_u               G_nom      ]
 */
SSVComplexMatrix* ssv_lft_input_multiplicative(
    const SSVComplexMatrix *G_nom,
    const SSVComplexMatrix *W_u,
    const SSVComplexMatrix *W_y);

/** Create M-Delta for multiplicative output uncertainty:
 *    G = (I + W_u * Delta * W_y) * G_nom
 *
 *  M = [ -W_y*G_nom*W_u   -W_y*G_nom ]
 *      [ G_nom*W_u          G_nom     ]
 */
SSVComplexMatrix* ssv_lft_output_multiplicative(
    const SSVComplexMatrix *G_nom,
    const SSVComplexMatrix *W_u,
    const SSVComplexMatrix *W_y);

/** Create M-Delta for feedback uncertainty (coprime factor type):
 *    G = N * (M + Delta_M)^{-1}  (left coprime factor uncertainty)
 *
 *  Simplified: G = G_nom * (I + Delta_feeback)^{-1}
 */
SSVComplexMatrix* ssv_lft_feedback_uncertainty(
    const SSVComplexMatrix *G_nom,
    const SSVComplexMatrix *W_u,
    const SSVComplexMatrix *W_y);

/* --- Augmented M for Robust Performance --- */

/** Augment M with a performance block for robust performance analysis.
 *
 *  Standard N-Delta structure with performance:
 *
 *    [ z_Delta ]   [ M11   M12 ] [ w_Delta ]
 *    [ z_perf  ] = [ M21   M22 ] [ w_perf  ]
 *
 *  For robust performance test via mu, we treat performance as an additional
 *  full complex uncertainty block Delta_p with ||Delta_p||_inf <= 1.
 *
 *  The augmented structure is:
 *    Delta_aug = diag(Delta, Delta_p)
 *
 *  and we test mu_{Delta_aug}(M) < 1.
 *
 *  @param M original M matrix [n x n]
 *  @param n_delta number of uncertainty channels
 *  @param n_perf number of performance channels
 *  @return augmented M (n+n_perf) x (n+n_perf)
 */
SSVComplexMatrix* ssv_lft_augment_performance(
    const SSVComplexMatrix *M, int n_delta, int n_perf);

#endif /* SSV_LFT_H */

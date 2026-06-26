/* gap_coprime.h — Normalized Coprime Factorization for Gap Metric Theory
 *
 * Central to the gap metric is the normalized coprime factorization
 * (NCF) of a plant G. This module computes left and right NCFs and
 * the graph symbol G_sym = [M; N].
 *
 * Knowledge coverage:
 *   L1: Left/right coprime factorization, normalized coprime factors
 *   L2: Bezout identity, graph symbol, perturbed coprime factors
 *   L3: Riccati equation for NCF, stabilizing solution
 *   L4: Small gain theorem connection, robust stability condition
 *   L5: CARE solver for NCF computation
 *
 * Ref: D.C. McFarlane & K. Glover, "Robust Controller Design Using
 *      Normalized Coprime Factor Plant Descriptions" (1990)
 *      M. Vidyasagar, "Control System Synthesis: A Factorization Approach" (1985)
 *      C.N. Nett, C.A. Jacobson, M.J. Balas, "A connection between state-space
 *      and doubly coprime fractional representations" (1984)
 */

#ifndef GAP_COPRIME_H
#define GAP_COPRIME_H

#include "gap_system.h"

/* ================================================================
 * L1: Core Definitions — Coprime Factorizations
 * ================================================================ */

/** GapNCF — Normalized Right Coprime Factorization (NRCF).
 *
 * For a stabilizable and detectable plant G(s) with state-space (A,B,C,D),
 * the normalized right coprime factorization G = N M^{-1} is obtained by:
 *
 *   [M(s)]  =  [ A + B F,   B R^{-1/2} ]
 *   [N(s)]     [ C + D F,   D R^{-1/2} ]
 *             [     F,      R^{-1/2}    ]
 *
 * where F = -R^{-1}(B^T X + D^T C), R = I + D^T D,
 * and X is the stabilizing solution to the Generalized Control
 * Algebraic Riccati Equation (GCARE):
 *   (A - B R^{-1} D^T C)^T X + X (A - B R^{-1} D^T C)
 *   - X B R^{-1} B^T X + C^T (I - D R^{-1} D^T) C = 0
 */
typedef struct {
    GapSystem *M;          /* Denominator factor (p × p system)     */
    GapSystem *N;          /* Numerator factor (p × m system)       */
    GapMatrix *X;          /* Stabilizing Riccati solution          */
    GapMatrix *F;          /* State feedback gain                   */
    bool       is_normalized; /* Whether N*N + M*M = I on jω-axis  */
} GapNCFRight;

/** GapNCFLeft — Normalized Left Coprime Factorization (NLCF).
 *
 * G = M̃^{-1} Ñ where M̃, Ñ are left coprime and satisfy:
 *   M̃ M̃* + Ñ Ñ* = I  on the jω-axis.
 *
 * Dual construction using the observer Riccati equation (GFARE).
 */
typedef struct {
    GapSystem *Mtilde;     /* Denominator factor (m × m system)     */
    GapSystem *Ntilde;     /* Numerator factor (m × p system)       */
    GapMatrix *Y;          /* Stabilizing Riccati solution          */
    GapMatrix *L;          /* Observer gain                        */
    bool       is_normalized; /* Whether M̃ M̃* + Ñ Ñ* = I on jω-axis */
} GapNCFLeft;

/** GapCoprimePerturbation — Additive perturbation of coprime factors.
 *
 * A perturbed plant G_Δ is represented as:
 *   G_Δ = (N + Δ_N) (M + Δ_M)^{-1}
 *
 * where [Δ_N; Δ_M] ∈ RH∞ is a stable perturbation.
 * The gap metric essentially measures the minimum H∞ norm
 * of such a perturbation that connects two plants.
 */
typedef struct {
    GapSystem *Delta_N;    /* Numerator perturbation               */
    GapSystem *Delta_M;    /* Denominator perturbation             */
    double     norm_bound; /* ||[Δ_N; Δ_M]||∞ bound               */
} GapCoprimePerturbation;

/* ================================================================
 * L5: Algorithms — NCF Computation
 * ================================================================ */

/** Compute normalized right coprime factorization (NRCF) of plant G.
 *
 * Algorithm:
 *   1. Form Hamiltonian matrix for GCARE
 *   2. Solve Riccati equation for X (stabilizing solution)
 *   3. Compute F = -R^{-1}(B^T X + D^T C)
 *   4. Build state-space realizations of M and N
 *
 * Returns NULL if G is not stabilizable or detectable.
 * Complexity: O(n^3) where n = state dimension.
 * Ref: McFarlane & Glover (1990), §3.2 */
GapNCFRight* gap_nrcf_compute(const GapSystem *G);

/** Compute normalized left coprime factorization (NLCF) of plant G.
 *  Dual construction using the observer Riccati equation.
 *  Ref: McFarlane & Glover (1990), §3.3 */
GapNCFLeft* gap_nlcf_compute(const GapSystem *G);

/** Free right coprime factorization result. */
void gap_nrcf_free(GapNCFRight *nrcf);

/** Free left coprime factorization result. */
void gap_nlcf_free(GapNCFLeft *nlcf);

/* ---- Graph Symbol ---- */

/** Build the graph symbol G_sym = [M; N] from NRCF.
 *
 * The graph symbol is the column concatenation of M and N,
 * forming a (p+m)×m system that encodes the entire plant geometry:
 *
 *   G_sym(s) = [ M(s) ]
 *              [ N(s) ]
 *
 * For normalized factors, G_sym is inner: G_sym* G_sym = I. */
GapGraphSymbol* gap_graph_symbol_from_nrcf(const GapNCFRight *nrcf);

/** Build graph symbol from left coprime factorization. */
GapGraphSymbol* gap_graph_symbol_from_nlcf(const GapNCFLeft *nlcf);

/** Free graph symbol. */
void gap_graph_symbol_free(GapGraphSymbol *gs);

/* ---- Bezout Identity Verification ---- */

/** Verify the Bezout identity for right coprime factors:
 *  There exist X_r, Y_r ∈ RH∞ such that X_r M + Y_r N = I.
 *  Returns the H∞ norm of the residual. */
double gap_verify_bezout_right(const GapNCFRight *nrcf);

/** Verify the Bezout identity for left coprime factors:
 *  M̃ X_l + Ñ Y_l = I for some X_l, Y_l ∈ RH∞. */
double gap_verify_bezout_left(const GapNCFLeft *nlcf);

/* ---- Perturbation Construction ---- */

/** Construct a coprime factor perturbation with specified norm bound.
 *  Useful for testing robust stability margins against
 *  gap metric perturbations. */
GapCoprimePerturbation* gap_coprime_perturbation_create(
    const GapNCFRight *nrcf, double epsilon);

/** Apply perturbation to plant: G_Δ = (N+Δ_N)(M+Δ_M)^{-1}.
 *  Returns the perturbed system. */
GapSystem* gap_coprime_perturbation_apply(
    const GapNCFRight *nrcf, const GapCoprimePerturbation *pert);

/** Free coprime perturbation. */
void gap_coprime_perturbation_free(GapCoprimePerturbation *pert);

/* ---- Geometric Distance ---- */

/** Compute the chordal distance between two points on the Riemann sphere.
 *
 * κ(P1(jω), P2(jω)) = |P1 - P2| / sqrt(1 + |P1|^2) / sqrt(1 + |P2|^2)
 *
 * This is the pointwise gap metric for SISO systems at frequency ω.
 * For MIMO, the chordal distance generalizes to the subspace gap.
 * Ref: G. Vinnicombe, "Measuring the Robustness of Feedback Systems"
 *      PhD Thesis, Cambridge (1993), Ch. 1 */
double gap_chordal_distance(double P1_re, double P1_im,
                             double P2_re, double P2_im);

/** Compute chordal distance between two transfer function evaluations. */
double gap_chordal_distance_matrix(const GapFreqResp *G1,
                                    const GapFreqResp *G2, double omega);

#endif /* GAP_COPRIME_H */

/**
 * @file hinf_riccati.h
 * @brief Algebraic Riccati Equation (ARE) solver for H-infinity control
 *
 * The ARE is the computational core of H-infinity synthesis.
 * For the DGKF algorithm, we solve two AREs:
 *
 *   Controller ARE: A^T X + X A - X (B2 B2^T - g^{-2} B1 B1^T) X + C1^T C1 = 0
 *   Filter ARE:     A Y + Y A^T - Y (C2^T C2 - g^{-2} C1^T C1) Y + B1 B1^T = 0
 *
 * Both must have stabilizing solutions (A-GX stable) and
 * satisfy coupling condition rho(X Y) < gamma^2.
 *
 * Two solution methods:
 *   1. Schur vector method (Laub 1979): build Hamiltonian, compute stable
 *      invariant subspace via ordered real Schur decomposition.
 *   2. Matrix sign function method (Roberts 1980): iterative sign function
 *      of the Hamiltonian yields the stabilizing solution.
 *
 * Knowledge Coverage:
 *   L3 (Math Structures): Hamiltonian matrix, stable invariant subspace
 *   L4 (Fundamental Laws): ARE solvability conditions, coupling condition
 *   L5 (Algorithms): Schur method, sign function method, Newton iteration
 *
 * References:
 *   Laub (1979) A Schur Method for Solving ARE, IEEE TAC 24(6):913-921
 *   Roberts (1980) Linear model reduction and solution of the ARE,
 *     Int. J. Control 32(4):677-687
 *   Lancaster & Rodman (1995) Algebraic Riccati Equations, Oxford
 *   Bini, Iannazzo, Meini (2012) Numerical Solution of ARE, SIAM
 *   DGKF (1989) IEEE TAC 34(8):831-847
 */

#ifndef HINF_RICCATI_H
#define HINF_RICCATI_H

#include "hinf_types.h"

/* ========== ARE Solution =============================================== */

/**
 * Solve the continuous-time algebraic Riccati equation (CARE):
 *   A^T X + X A - X R X + Q = 0
 *
 * with the stabilizing property: Re(λ(A - R X)) < 0.
 *
 * Uses the Schur vector method (Laub 1979):
 *   1. Form the 2n x 2n Hamiltonian matrix:
 *      H = [A,   -R;  -Q,  -A^T]
 *   2. Compute ordered real Schur decomposition: U^T H U = [T11 T12; 0 T22]
 *      where T11 has all stable eigenvalues and T22 all unstable.
 *   3. Partition U = [U11; U21] and solve X = U21 U11^{-1}.
 *
 * Parameters:
 *   A: n x n system matrix
 *   R: n x n symmetric (may be indefinite for H-inf)
 *   Q: n x n symmetric (may be indefinite for H-inf)
 *   X: output n x n solution matrix (must be pre-allocated)
 *
 * Returns 0 on success, <0 on failure:
 *   -1: no stabilizing solution exists
 *   -2: U11 singular (invariant subspace not graph subspace)
 *   -3: invalid input
 *
 * Complexity: O(n^3) dominated by Schur decomposition.
 * Ref: Laub (1979) Theorem 2, Algorithm on p.917
 */
int hinf_care_solve(const HinfMatrix *A, const HinfMatrix *R,
                     const HinfMatrix *Q, HinfMatrix *X,
                     HinfCARE *info);

/**
 * Solve the ARE using the matrix sign function method.
 *
 *   H = [A, -R; -Q, -A^T]
 *   sign(H) = H * H^{-1/2}
 *
 * The stabilizing solution X is obtained from:
 *   sign(H) * [I; X] = -[I; X]   =>   (I + sign(H)) [I; X] = 0
 *
 * This method is quadratically convergent but requires H to have
 * no eigenvalues on the imaginary axis.
 *
 * Returns 0 on success.
 * Ref: Roberts (1980), Higham (2008) Ch. 5
 */
int hinf_care_sign(const HinfMatrix *A, const HinfMatrix *R,
                    const HinfMatrix *Q, HinfMatrix *X,
                    HinfCARE *info);

/* ========== H-infinity Specific AREs =================================== */

/**
 * Build the Hamiltonian matrix for the controller ARE in DGKF.
 *
 * H_X = [ A,
 *         -(B2 B2^T - gamma^{-2} B1 B1^T);
 *         -C1^T C1,
 *         -A^T ]
 *
 * Complexity: O(n^3) for matrix products.
 * Ref: DGKF (1989) eq. (23)
 */
int hinf_hamiltonian_controller(const HinfGenPlant *P, double gamma,
                                  HinfHamiltonian *H);

/**
 * Build the Hamiltonian matrix for the filter ARE in DGKF.
 *
 * H_Y = [ A^T,
 *         -(C2^T C2 - gamma^{-2} C1^T C1);
 *         -B1 B1^T,
 *         -A ]^T
 *
 * (Hamiltonian for the dual system.)
 *
 * Complexity: O(n^3) for matrix products.
 * Ref: DGKF (1989) eq. (24)
 */
int hinf_hamiltonian_filter(const HinfGenPlant *P, double gamma,
                              HinfHamiltonian *H);

/**
 * Solve the controller ARE for DGKF:
 *   A^T X + X A - X (B2 B2^T - g^{-2} B1 B1^T) X + C1^T C1 = 0
 *
 * Returns the stabilizing solution X >= 0.
 * Ref: DGKF (1989) Theorem 3, eq. (25)
 */
int hinf_care_controller(const HinfGenPlant *P, double gamma,
                          HinfMatrix *X, HinfCARE *info);

/**
 * Solve the filter ARE for DGKF:
 *   A Y + Y A^T - Y (C2^T C2 - g^{-2} C1^T C1) Y + B1 B1^T = 0
 *
 * Returns the stabilizing solution Y >= 0.
 * Ref: DGKF (1989) Theorem 3, eq. (27)
 */
int hinf_care_filter(const HinfGenPlant *P, double gamma,
                      HinfMatrix *Y, HinfCARE *info);

/* ========== DGKF Assumption Checking =================================== */

/**
 * Verify the standard DGKF assumptions required for the two-Riccati
 * H-infinity solution:
 *
 * (A1) (A, B2) is stabilizable and (C2, A) is detectable.
 * (A2) D12 has full column rank: [D12, D_perp] is unitary where
 *      D_perp completes the orthogonal basis.
 * (A3) D21 has full row rank: [D21; D_perp] is unitary.
 * (A4) rank([A-jwI, B2; C1, D12]) = n + n_u for all w.
 *
 * Returns 0 if all satisfied, bitmask of failed assumptions otherwise.
 * Complexity: O(n^3).
 * Ref: DGKF (1989) Section IV, Assumptions A1-A4
 */
int hinf_dgkf_check_assumptions(const HinfGenPlant *P, double gamma);

/**
 * Check the DGKF coupling condition:
 *   rho(X_inf * Y_inf) < gamma^2
 *
 * This is necessary and sufficient for the existence of
 * a gamma-suboptimal controller.
 *
 * Returns 1 if satisfied, 0 otherwise.
 * Ref: DGKF (1989) Lemma 4, Theorem 3
 */
int hinf_dgkf_coupling_check(const HinfMatrix *X, const HinfMatrix *Y,
                               double gamma);

/* ========== Discrete-time ARE ========================================== */

/**
 * Solve the discrete-time ARE:
 *   X = A^T X A - A^T X B (R + B^T X B)^{-1} B^T X A + Q
 *
 * with stabilizing property. Used in discrete-time H-inf.
 *
 * Uses the symplectic pencil approach.
 *
 * Complexity: O(n^3).
 * Ref: Pappas, Laub, Sandell (1980) IEEE TAC 25(4):631-641
 *      Ionescu & Weiss (1993) IMA J. Math. Control & Info 10
 */
int hinf_dare_solve(const HinfMatrix *A, const HinfMatrix *R,
                     const HinfMatrix *Q, HinfMatrix *X,
                     HinfCARE *info);

#endif /* HINF_RICCATI_H */
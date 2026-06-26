#ifndef MU_MATRIX_H
#define MU_MATRIX_H

#include "mu_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ====================================================================
 * mu_matrix.h — Extended Matrix Operations for u-Synthesis
 *
 * Implements linear algebra operations needed for:
 *   - Riccati equation solvers (Hinf synthesis)
 *   - Matrix sign function and spectral factorization
 *   - QR decomposition for minimum realization
 *   - Lyapunov equation solver
 *
 * Reference: Golub & Van Loan (2013), "Matrix Computations", 4th ed.
 *            Laub (1979), "A Schur method for solving algebraic Riccati equations"
 * ==================================================================== */

/*
 * L3: QR Decomposition
 *
 * A = Q * R where Q is unitary (Q* Q = I) and R is upper triangular.
 * Uses Householder reflections.
 *
 * Complexity: O(n^3)
 * Reference: Golub & Van Loan (2013), §5.2
 */
void mu_mat_qr(const MUMatrix *A, MUMatrix **Q, MUMatrix **R);

/*
 * L3: Schur Decomposition
 *
 * A = U * T * U* where U is unitary and T is upper triangular
 * (or quasi-upper triangular for real matrices).
 *
 * The eigenvalues of A appear on the diagonal of T.
 *
 * Complexity: O(n^3)
 * Reference: Golub & Van Loan (2013), §7.4
 */
void mu_mat_schur(const MUMatrix *A, MUMatrix **U, MUMatrix **T);

/*
 * L5: Algebraic Riccati Equation (ARE)
 *
 * Solve the continuous-time algebraic Riccati equation:
 *
 *   A' X + X A + X R X + Q = 0
 *
 * where A, Q, R are n×n real matrices, Q = Q', R = R'.
 * The stabilizing solution (A + R X has eigenvalues in open LHP).
 *
 * Uses the Schur vector method (Laub 1979):
 *   1. Form Hamiltonian matrix H = [A, R; -Q, -A']
 *   2. Compute Schur decomposition of H
 *   3. Reorder to have stable eigenvalues first
 *   4. Extract X = U21 * U11^{-1}
 *
 * Complexity: O(n^3)
 * Reference: Laub (1979), "A Schur method for solving algebraic Riccati equations"
 *            Arnold & Laub (1984), "Generalized eigenproblem algorithms"
 */
typedef struct {
    MUMatrix *X;        /* Solution matrix (may be NULL if no stabilizing sol.) */
    bool      found;    /* true if a stabilizing solution exists */
    int       iterations;
} MURiccatiResult;

MURiccatiResult mu_solve_riccati(const MUMatrix *A, const MUMatrix *R,
                                  const MUMatrix *Q);

/*
 * L5: Lyapunov Equation
 *
 * Solve A X + X A' + Q = 0 (continuous-time Lyapunov equation)
 *
 * Uses Bartels-Stewart algorithm with Schur decomposition.
 *
 * Complexity: O(n^3)
 * Reference: Bartels & Stewart (1972),
 *            "Solution of the matrix equation AX + XB = C"
 */
MUMatrix* mu_solve_lyapunov(const MUMatrix *A, const MUMatrix *Q);

/*
 * L3: Cholesky Decomposition
 *
 * A = L * L' where A = A' > 0 and L is lower triangular.
 * Verifies positive definiteness of D-scales and Riccati solutions.
 *
 * Complexity: O(n^3)
 */
bool mu_mat_cholesky(const MUMatrix *A, MUMatrix **L);

/*
 * L3: Matrix Exponential
 *
 * expm(A) = sum_{k=0}^{inf} A^k / k!
 *
 * Uses scaling-and-squaring with Pade approximation.
 * Needed for discrete-time conversions and time-domain simulation.
 *
 * Complexity: O(n^3)
 * Reference: Higham (2005), "The scaling and squaring method for
 *            the matrix exponential revisited"
 */
MUMatrix* mu_mat_exp(const MUMatrix *A);

/*
 * L3: Controllability Gramian
 *
 * W_c = integral_0^inf exp(A t) B B' exp(A' t) dt
 *
 * Solved via Lyapunov: A W_c + W_c A' + B B' = 0
 *
 * Complexity: O(n^3)
 */
MUMatrix* mu_gramian_controllability(const MUMatrix *A, const MUMatrix *B);

/*
 * L3: Observability Gramian
 *
 * W_o = integral_0^inf exp(A' t) C' C exp(A t) dt
 *
 * Solved via Lyapunov: A' W_o + W_o A + C' C = 0
 *
 * Complexity: O(n^3)
 */
MUMatrix* mu_gramian_observability(const MUMatrix *A, const MUMatrix *C);

/*
 * L3: Hankel Singular Values
 *
 * sigma_i = sqrt(lambda_i(W_c * W_o))
 *
 * Measure the joint controllability and observability of each state.
 * Used for balanced realization and model reduction.
 *
 * Complexity: O(n^3)
 * Reference: Moore (1981), "Principal component analysis in linear systems"
 */
int mu_hankel_singular_values(const MUMatrix *A, const MUMatrix *B,
                               const MUMatrix *C, double *hsv);

/*
 * L5: Minimum Realization
 *
 * Find a minimal (controllable and observable) realization of a system.
 * Removes states that are uncontrollable or unobservable.
 *
 * Complexity: O(n^3)
 */
MUSystem* mu_minreal(const MUSystem *sys);

/*
 * L8: Balanced Realization
 *
 * Transform system to balanced coordinates where controllability
 * and observability gramians are equal and diagonal.
 *
 * This reveals which states are most important for input-output behavior
 * and enables safe truncation for model reduction.
 *
 * Complexity: O(n^3)
 * Reference: Moore (1981), "Principal component analysis in linear systems"
 */
MUSystem* mu_balance(const MUSystem *sys, double *hsv);

/*
 * L8: Balanced Truncation (Model Reduction)
 *
 * Reduce system order by truncating states with small Hankel singular values.
 *
 * Error bound: ||G(s) - G_r(s)||_inf <= 2 * sum_{i=r+1}^{n} sigma_i
 *
 * Complexity: O(n^3)
 * Reference: Enns (1984), Glover (1984)
 */
MUSystem* mu_balanced_truncation(const MUSystem *sys, int target_order);

/*
 * L3: Condition Number
 *
 * kappa(A) = sigma_max(A) / sigma_min(A)
 *
 * Measures sensitivity of linear system solutions to perturbations.
 * Large condition numbers indicate numerical difficulty in mu computation.
 *
 * Complexity: O(n^3)
 */
double mu_mat_condition(const MUMatrix *A);

/*
 * L3: Osborne Balancing for D-Scales
 *
 * Find diagonal similarity transform D such that D A D^{-1} has
 * minimal Frobenius norm (or nearly equal row/column norms).
 *
 * This is the core iteration in the D-step for scalar uncertainties.
 *
 * Algorithm: iteratively scale rows/columns to balance norms.
 *
 * Complexity: O(n^2 * max_iter)
 * Reference: Osborne (1960), "On preconditioning by diagonal scaling"
 *            Parlett & Reinsch (1969), "Balancing a matrix for
 *            calculation of eigenvalues and eigenvectors"
 */
double* mu_osborne_balancing(const MUMatrix *A, int n, int max_iter, double tol);

#ifdef __cplusplus
}
#endif

#endif /* MU_MATRIX_H */

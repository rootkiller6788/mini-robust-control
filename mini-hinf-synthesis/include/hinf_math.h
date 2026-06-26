/**
 * @file hinf_math.h
 * @brief Advanced matrix mathematics: eigenvalues, QR, Schur, spectral radius
 *
 * Provides the linear algebra backbone for Algebraic Riccati Equation
 * solvers, the Bounded Real Lemma eigenvalue test, and H-infinity synthesis.
 *
 * Knowledge Coverage:
 *   L3 (Math Structures): Schur form, real block-diagonal decomposition
 *   L5 (Algorithms): QR algorithm with Wilkinson shifts, eigenvalue ordering
 *
 * References:
 *   Golub & Van Loan (2013) Matrix Computations, 4th ed., Ch. 7-8
 *   Laub (1979) A Schur Method for Solving ARE, IEEE TAC 24(6)
 *   Francis (1961) The QR Transformation, Parts I-II, Comp. J. 4
 */

#ifndef HINF_MATH_H
#define HINF_MATH_H

#include "hinf_types.h"

/* ========== Eigenvalue Computation ================================= */

/** Compute all eigenvalues of a real square matrix via double-shift QR.
 *  Output: real_parts[n], imag_parts[n]. Complexity: O(n^3).
 *  Ref: Golub & Van Loan (2013) Algorithm 7.5.2, Francis (1961) */
int hinf_eigenvalues(const HinfMatrix *A, double *real_parts, double *imag_parts);

/** Symmetric eigenvalue problem via tridiagonalization + implicit QL.
 *  All eigenvalues real. Complexity: O(n^3).
 *  Ref: Golub & Van Loan (2013) Alg 8.3.3 */
int hinf_eigenvalues_symm(const HinfMatrix *A, double *eigvals);

/** Dominant eigenvalue (largest |lambda|) via power iteration.
 *  Complexity: O(k * n^2). Ref: Wilkinson (1965) Ch. 9 */
double hinf_eigenvalue_dominant(const HinfMatrix *A, double *eigvec, int max_iter);

/* ========== Real Schur Decomposition =============================== */

/** Real Schur: A = Q T Q^T. T quasi-upper-triangular, Q orthogonal.
 *  Output: T overwrites A; Q optional. Complexity: O(n^3).
 *  Ref: Golub & Van Loan (2013) Alg 7.5.2, Laub (1979) IEEE TAC */
int hinf_schur_decomp(HinfMatrix *A, HinfMatrix *Q);

/** Reorder eigenvalue block at position select to top-left.
 *  Ref: Bai & Demmel (1993) Lin. Alg. Appl. 186 */
int hinf_schur_ord_select(HinfMatrix *T, HinfMatrix *Q, int select);

/** Order stable eigenvalues (Re < 0) first. Ref: Laub (1979) */
int hinf_schur_ord_stable(HinfMatrix *T, HinfMatrix *Q);

/** Count eigenvalues with Re < 0 in Schur form. */
int hinf_schur_count_stable(const HinfMatrix *T);

/** Count eigenvalues on imaginary axis (|Re| < tol). */
int hinf_schur_count_imag(const HinfMatrix *T, double tol);

/* ========== Spectral Radius ======================================== */

/** rho(A) = max_i |lambda_i(A)|. Used in coupling condition.
 *  Complexity: O(k * n^2) via power iteration. */
double hinf_spectral_radius(const HinfMatrix *A);

/** rho(A * B) without forming product explicitly. */
double hinf_spectral_radius_product(const HinfMatrix *A, const HinfMatrix *B);

/* ========== Matrix Sign Function =================================== */

/** Matrix sign via Newton-Schulz: X_{k+1} = 0.5 X_k (3I - X_k^2).
 *  Quadratic convergence. Used for ARE. Complexity: O(k*n^3).
 *  Ref: Roberts (1980) Int. J. Control, Higham (2008) */
int hinf_matrix_sign(HinfMatrix *A, HinfMatrix *signA, int max_iter, double tol);

/** Matrix sqrt A^{1/2} for SPD via Denman-Beavers iteration. */
int hinf_matrix_sqrt(HinfMatrix *A, HinfMatrix *sqrtA, int max_iter);

/* ========== QR Factorization ====================================== */

/** QR via Householder reflections: A = Q R. Complexity: O(m*n^2).
 *  Ref: Golub & Van Loan (2013) Alg 5.2.1 */
int hinf_qr_decomp(HinfMatrix *A, HinfMatrix *Q);

/** Least-squares solve min||Ax-b||_2 via QR. Returns residual norm. */
double hinf_qr_ls_solve(const HinfMatrix *A, const double *b, double *x);

/* ========== Singular Value Decomposition =========================== */

/** Max singular value via power iteration on A^T A.
 *  Ref: Golub & Van Loan (2013) Sec 8.6 */
double hinf_svd_sigma_max(const HinfMatrix *A, double *u, double *v, int max_iter);

/** Full SVD via Golub-Kahan bidiagonalization + QR.
 *  Ref: Golub & Kahan (1965) J. SIAM B, Golub & Van Loan Alg 8.6.2 */
int hinf_svd(const HinfMatrix *A, double *sigma, HinfMatrix *U, HinfMatrix *VT);

/* ========== Matrix Balancing ====================================== */

/** Balance via diagonal similarity for eigenvalue accuracy.
 *  Ref: Osborne (1960), Parlett & Reinsch (1969) */
int hinf_balance(HinfMatrix *A, double *scale);

/* ========== Symplectic Structure =================================== */

/** Check Hamiltonian: H = [A R; Q -A^T] with R=RT, Q=QT.
 *  Ref: Paige & Van Loan (1981) Lin. Alg. Appl. 41 */
int hinf_check_hamiltonian(const HinfMatrix *H, double tol);

/** Check symplectic: U^T J U = J. */
int hinf_check_symplectic(const HinfMatrix *U, double tol);

/** Form J = [0 I_n; -I_n 0] (canonical symplectic matrix). */
void hinf_symplectic_J(HinfMatrix *J, int n);

/* ========== Lyapunov & Sylvester ================================== */

/** CT Lyapunov: A X + X A^T + Q = 0 (A stable, Q=QT).
 *  Bartels-Stewart via Schur. Complexity: O(n^3). */
int hinf_lyapunov_ct(const HinfMatrix *A, const HinfMatrix *Q, HinfMatrix *X);

/** Sylvester: A X + X B = C. Ref: Bartels & Stewart (1972) */
int hinf_sylvester(const HinfMatrix *A, const HinfMatrix *B,
                    const HinfMatrix *C, HinfMatrix *X);

/* ========== Matrix Exponential ==================================== */

/** exp(A*t) via scaling-and-squaring + Pade. Ref: Higham (2005) SIREV. */
int hinf_matrix_exp(const HinfMatrix *A, double t, HinfMatrix *expAt);

#endif /* HINF_MATH_H */
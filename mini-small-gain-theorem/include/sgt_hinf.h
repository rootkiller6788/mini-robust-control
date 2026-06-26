#ifndef SGT_HINF_H
#define SGT_HINF_H

#include "sgt_core.h"

/* ============================================================================
 * Small Gain Theorem — H-infinity Norm Computation
 *
 * The H-infinity norm is the induced L2 gain of an LTI system:
 *   ||G||_inf = sup_{omega in R} sigma_max(G(j omega))
 *
 * For the small gain theorem, we need to compute ||H1||_inf and ||H2||_inf.
 * This header provides algorithms for computing H-infinity norms,
 * solving algebraic Riccati equations, and evaluating frequency responses.
 *
 * Key algorithms:
 *   1. Frequency grid method (simple but approximate)
 *   2. Bisection method (exact, via Hamiltonian matrix eigenvalues)
 *   3. Riccati-based method (via the Bounded Real Lemma)
 * ============================================================================ */

/* --- H-infinity Norm Computation --- */

/**
 * Compute the H-infinity norm of an LTI system using a frequency grid.
 * Evaluates sigma_max(G(j omega)) at n_omega logarithmically spaced
 * frequency points and returns the peak.
 *
 * This method is computationally simple but may miss sharp peaks between
 * grid points. For exact computation, use sgt_hinf_bisection().
 *
 * Complexity: O(n_omega * (nx^2 * ny * nu)) per frequency point.
 * Reference: Zhou, Doyle, Glover (1996), Section 4.3. */
double sgt_hinf_grid(const SGTLTISystem *sys, int n_omega);

/**
 * Compute the H-infinity norm exactly using the bisection method.
 *
 * Algorithm (Boyd, Balakrishnan, Kabamba 1989):
 *   The H-infinity norm of G(s) = (A,B,C,D) is < gamma iff:
 *     (i)  gamma > sigma_max(D)
 *     (ii) The Hamiltonian matrix
 *            H_gamma = [A + B*R^{-1}*D^T*C,   B*R^{-1}*B^T;
 *                       -C^T*(I+D*R^{-1}*D^T)*C, -(A + B*R^{-1}*D^T*C)^T]
 *          where R = gamma^2*I - D^T*D, has no eigenvalues on the
 *          imaginary axis.
 *
 * Bisection refines a lower bound gamma_lo (initially sigma_max(D)) and
 * an upper bound gamma_hi (from norm or frequency grid) until
 * gamma_hi - gamma_lo < tol.
 *
 * Complexity: O(log2(range/tol) * nx^3).
 * Reference: S. Boyd, V. Balakrishnan, P. Kabamba, "A bisection method
 *   for computing the H-infinity norm of a transfer matrix and related
 *   problems," Math. Control Signals Systems, 1989. */
double sgt_hinf_bisection(const SGTLTISystem *sys, double tol, int max_iter);

/**
 * Compute sigma_max(G(j omega)) at a single frequency.
 * This is the maximum singular value of the complex matrix
 * G(jw) = C*(jwI - A)^{-1}*B + D.
 *
 * Uses LAPACK-style approach: builds the real 2x2 block representation
 * and computes singular values via eigenvalue decomposition.
 * Complexity: O(nx^3) due to the (jwI-A)^{-1} solve. */
double sgt_freq_response_max_sv(const SGTLTISystem *sys, double omega);

/* --- Hamiltonian Matrix --- */

/**
 * Build the Hamiltonian matrix associated with the H-infinity problem:
 *   H_gamma = [A + B*Z1*C,     B*Z2*B^T;
 *              -C^T*Z3*C,      -(A + B*Z1*C)^T]
 * where Z1, Z2, Z3 are functions of gamma and D.
 *
 * The matrix is stored in a 2nx x 2nx array.
 * Returns NULL if gamma <= sigma_max(D) (infeasible).
 * Complexity: O(nx^3) for matrix multiplications. */
SGTMatrix* sgt_build_hamiltonian(const SGTLTISystem *sys, double gamma);

/**
 * Check if the Hamiltonian matrix has eigenvalues on the imaginary axis.
 * This is the key test in the bisection method: if H_gamma has imaginary
 * axis eigenvalues, then ||G||_inf >= gamma.
 *
 * Uses a symplectic QR-like algorithm specialized for Hamiltonian
 * eigenvalue structure.
 * Complexity: O(nx^3). */
bool sgt_hamiltonian_has_imag_eigs(const SGTMatrix *H, double tol);

/* --- Bounded Real Lemma --- */

/**
 * Check if ||G||_inf < gamma using the Bounded Real Lemma.
 *
 * Bounded Real Lemma (BRL):
 *   ||G||_inf < gamma  iff  there exists X = X^T > 0 solving the ARE:
 *     A^T*X + X*A + C^T*C + X*B*(gamma^2*I - D^T*D)^{-1}*B^T*X = 0
 *   and (A + B*(gamma^2*I - D^T*D)^{-1}*(B^T*X + D^T*C)) is stable.
 *
 * This implementation solves the ARE via the associated Hamiltonian
 * matrix and the stable invariant subspace method.
 *
 * Complexity: O(nx^3).
 * Reference: P. Gahinet, P. Apkarian, "A linear matrix inequality approach
 *   to H-infinity control," Int. J. Robust Nonlinear Control, 1994. */
bool sgt_bounded_real_lemma(const SGTLTISystem *sys, double gamma);

/* --- Algebraic Riccati Equation --- */

/**
 * Solve the continuous-time algebraic Riccati equation (CARE):
 *   A^T*X + X*A - X*B*R^{-1}*B^T*X + Q = 0
 *
 * Solved via the Hamiltonian matrix method:
 *   H = [A, -B*R^{-1}*B^T; -Q, -A^T]
 * The stabilizing solution X is obtained from the Schur vectors
 * corresponding to the stable eigenvalues of H.
 *
 * Returns the solution X as a newly allocated matrix, or NULL if no
 * stabilizing solution exists.
 * Complexity: O(nx^3).
 * Reference: A.J. Laub, "A Schur method for solving algebraic Riccati
 *   equations," IEEE Trans. AC, 1979. */
SGTMatrix* sgt_solve_care(const SGTMatrix *A, const SGTMatrix *B,
                           const SGTMatrix *R, const SGTMatrix *Q);

/**
 * Solve a small (2x2 or 3x3) Riccati equation analytically.
 * For nx <= 3, uses closed-form solutions derived from the
 * Hamiltonian eigenvector method applied directly.
 * Returns the stabilizing solution X, or NULL if none exists. */
SGTMatrix* sgt_solve_care_small(const SGTMatrix *A, const SGTMatrix *B,
                                 const SGTMatrix *R, const SGTMatrix *Q);

/* --- H2 Norm --- */

/**
 * Compute the H2 norm of an LTI system (assuming D=0 for strict properness).
 * ||G||_2 = sqrt( (1/2*pi) * integral_{-inf}^{inf} trace(G(jw)*G(jw)^H) dw )
 *
 * Computed via the controllability Gramian Wc:
 *   ||G||_2^2 = trace(C * Wc * C^T)
 * where A*Wc + Wc*A^T + B*B^T = 0.
 *
 * Complexity: O(nx^3) for the Lyapunov equation. */
double sgt_h2_norm(const SGTLTISystem *sys);

/**
 * Compute the Hankel norm of an LTI system.
 * The Hankel singular values are the square roots of the eigenvalues
 * of the product of controllability and observability Gramians.
 * The Hankel norm is the largest Hankel singular value.
 * Complexity: O(nx^3). */
double sgt_hankel_norm(const SGTLTISystem *sys);

/* --- Structured Norm Computation --- */

/**
 * Compute the H-infinity norm with D-scaling.
 * For structured uncertainty with block-diagonal Delta, the structured
 * singular value mu is bounded by:
 *   inf_{D in D} sigma_max(D * M(jw) * D^{-1})
 * where D is a set of scaling matrices that commute with Delta.
 *
 * This function computes the upper bound via frequency-dependent
 * D-scaling optimization at each frequency point.
 *
 * Complexity: O(n_omega * n_d_iter * n_total_dim^3). */
double sgt_hinf_with_dscale(const SGTLTISystem *sys,
                             const SGTStructuredUncertainty *delta_struct,
                             int n_omega, int n_d_iter);

#endif /* SGT_HINF_H */

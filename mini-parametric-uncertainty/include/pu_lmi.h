#ifndef PU_LMI_H
#define PU_LMI_H

#include "pu_core.h"

/* ============================================================================
 * Linear Matrix Inequality (LMI) Methods for Robust Control
 *
 * LMI formulations for parametric uncertainty problems:
 *
 * 1. Quadratic Stability LMI:
 *    Find P = P' > 0 such that A_i'*P + P*A_i < 0 for all i
 *
 * 2. Bounded Real Lemma (H-infinity):
 *    ||G(s)||_inf < gamma iff exists P > 0 such that
 *    [ A'*P + P*A   P*B   C'  ]
 *    [    B'*P    -gamma*I D'  ] < 0
 *    [     C         D   -gamma*I ]
 *
 * 3. H2 Performance LMI:
 *    ||G(s)||_2 < nu iff exists P > 0, W such that
 *    A'*P + P*A + C'*C < 0, trace(B'*P*B) < nu^2
 *
 * 4. Generalized Eigenvalue Problem (GEVP):
 *    Minimize lambda subject to A(x) < lambda*B(x), B(x) > 0, C(x) > 0
 *
 * Reference: Boyd, El Ghaoui, Feron, Balakrishnan (1994)
 *   "Linear Matrix Inequalities in System and Control Theory", SIAM.
 * ============================================================================ */

/* --- LMI Data Structures --- */

/** An LMI constraint: F(x) = F0 + sum x_i * F_i < 0 */
typedef struct {
    int      n;                /* Matrix dimension */
    int      n_vars;           /* Number of decision variables */
    double **F0;               /* Constant term (n x n, symmetric) */
    double ***F;               /* Coefficient matrices: n_vars x n x n, each symmetric */
} PU_LMI_Constraint;

/** A set of LMI constraints */
typedef struct {
    int               n_constraints;
    PU_LMI_Constraint *constraints;
} PU_LMI_Problem;

/** Solution of an LMI problem */
typedef struct {
    int      n_vars;
    double  *x;                /* Decision variables */
    double   objective;        /* Objective value */
    bool     feasible;         /* True if feasible solution found */
    int      iterations;       /* Iterations used by solver */
    double   duality_gap;      /* Final duality gap */
} PU_LMI_Solution;

/* --- LMI Construction --- */

/** Create an LMI constraint of dimension n with n_vars variables */
PU_LMI_Constraint pu_lmi_constraint_create(int n, int n_vars);

/** Set the constant term F0 */
void pu_lmi_constraint_set_F0(PU_LMI_Constraint *c, double **F0);

/** Set the coefficient matrix for variable k */
void pu_lmi_constraint_set_Fk(PU_LMI_Constraint *c, int k, double **Fk);

/** Free LMI constraint */
void pu_lmi_constraint_free(PU_LMI_Constraint *c);

/** Create an LMI problem */
PU_LMI_Problem pu_lmi_problem_create(int n_constraints);

/** Add a constraint to the problem */
void pu_lmi_problem_add(PU_LMI_Problem *prob, const PU_LMI_Constraint *c);

/** Free LMI problem */
void pu_lmi_problem_free(PU_LMI_Problem *prob);

/* --- LMI Solvers --- */

/**
 * Solve LMI feasibility using the alternating projection method
 * (also known as the "method of alternating projections" or MAP).
 *
 * This is a simple iterative method suitable for small problems.
 * For larger problems, interior-point methods (e.g., SeDuMi, SDPT3)
 * should be used.
 *
 * @param prob       LMI problem
 * @param max_iter   Maximum iterations
 * @param tol        Convergence tolerance
 * @return           LMI solution
 */
PU_LMI_Solution pu_lmi_solve_ap(const PU_LMI_Problem *prob, int max_iter, double tol);

/**
 * Solve LMI feasibility using a simplified primal-dual interior point method.
 *
 * Minimizes t subject to F_i(x) <= t*I.
 *
 * @param prob       LMI problem
 * @param max_iter   Maximum iterations
 * @param tol        Tolerance
 * @return           LMI solution
 */
PU_LMI_Solution pu_lmi_solve_ip(const PU_LMI_Problem *prob, int max_iter, double tol);

/** Free LMI solution */
void pu_lmi_solution_free(PU_LMI_Solution *sol);

/** Print LMI solution */
void pu_lmi_solution_print(const PU_LMI_Solution *sol);

/* --- Specialized Control LMIs --- */

/**
 * Build and solve the Quadratic Stability LMI:
 *   Find P = P' > 0 such that A_i'*P + P*A_i < 0 for all i = 1..N
 *
 * @param A_vertices  Array of N vertex matrices (each n x n)
 * @param N           Number of vertices
 * @param n           System dimension
 * @param P_out       Output Lyapunov matrix (n x n), caller-allocated
 * @return            True if feasible
 */
bool pu_lmi_quadratic_stability(double ***A_vertices, int N, int n, double **P_out);

/**
 * Build and solve the Bounded Real Lemma LMI:
 *   Find P = P' > 0 and gamma > 0 minimizing gamma such that
 *   [ A'*P + P*A   P*B   C'  ]
 *   [    B'*P    -gamma*I D'  ] < 0
 *   [     C         D   -gamma*I ]
 *
 * @param A          System A matrix (n x n)
 * @param B          System B matrix (n x m)
 * @param C          System C matrix (p x n)
 * @param D          System D matrix (p x m)
 * @param n, m, p    Dimensions
 * @param gamma_opt  Output: optimal H-infinity norm bound
 * @return           True if LMI solved
 */
bool pu_lmi_hinf_norm(double **A, double **B, double **C, double **D,
                       int n, int m, int p, double *gamma_opt);

/**
 * Solve robust H2 performance via LMI.
 *
 * @param A_vertices  Vertex A matrices for polytopic uncertainty
 * @param B_vertices  Vertex B matrices
 * @param C_vertices  Vertex C matrices
 * @param N           Number of vertices
 * @param n, m, p     Dimensions
 * @param nu_opt      Output: optimal H2 norm bound
 * @return            True if feasible
 */
bool pu_lmi_h2_robust(double ***A_vertices, double ***B_vertices,
                       double ***C_vertices, int N, int n, int m, int p,
                       double *nu_opt);

/**
 * State-feedback robust stabilization via LMI.
 * Find K such that (A_i + B_i*K) is stable for all i.
 *
 * Uses the change of variables Y = K*Q where Q = P^{-1}.
 * LMI: A_i*Q + Q*A_i' + B_i*Y + Y'*B_i' < 0 for all i, Q > 0.
 *
 * @param A_vertices  Vertex A matrices
 * @param B_vertices  Vertex B matrices
 * @param N           Number of vertices
 * @param n, m        Dimensions
 * @param K_out       Output feedback gain (m x n), caller-allocated
 * @return            True if feasible
 */
bool pu_lmi_state_feedback_robust(double ***A_vertices, double ***B_vertices,
                                    int N, int n, int m, double **K_out);

/**
 * Build a parameter-dependent Lyapunov function LMI.
 * Find P(q) = P0 + sum q_i * P_i > 0 such that
 * A(q)'*P(q) + P(q)*A(q) < 0 for all q in the polytope.
 *
 * This is less conservative than quadratic stability.
 */
bool pu_lmi_parameter_dependent_lyapunov(double ***A_vertices,
                                           const double **param_vertices,
                                           int N, int n, int n_params,
                                           double **P0_out,
                                           double ***Pk_out);

#endif /* PU_LMI_H */

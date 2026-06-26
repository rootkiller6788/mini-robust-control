#ifndef PU_STABILITY_H
#define PU_STABILITY_H

#include "pu_core.h"

/* ============================================================================
 * Robust Stability Analysis under Parametric Uncertainty
 *
 * Key theoretical foundations:
 *
 * 1. Zero Exclusion Principle:
 *    An interval polynomial family is robustly stable iff 0 is not in
 *    the value set P(j*omega, Q) for any omega >= 0, and one member
 *    of the family is stable.
 *
 * 2. Mapping Theorem (Barmish, 1994):
 *    For an interval polynomial, the value set at s = j*omega is a
 *    rectangle in the complex plane whose corners are given by the
 *    four Kharitonov polynomials evaluated at s = j*omega.
 *
 * 3. Edge Theorem (Bartlett-Hollot-Lin, 1988):
 *    For a polytope of polynomials, stability of all edges implies
 *    stability of the entire polytope.
 *
 * 4. Lyapunov Quadratic Stability:
 *    If there exists P = P' > 0 such that A_i'*P + P*A_i < 0 for all
 *    vertex systems A_i, then the polytopic system is quadratically stable.
 *
 * 5. Small-Gain Theorem:
 *    For systems with norm-bounded uncertainty Delta:
 *    Robust stability holds if ||M||_infinity * ||Delta||_infinity < 1.
 * ============================================================================ */

/* --- Value Set Computation --- */

/**
 * Compute the value set of an interval polynomial at s = j*omega.
 *
 * Uses the Mapping Theorem: the value set is a rectangle in the
 * complex plane, axis-aligned only in special cases.
 *
 * For general interval polynomials, we compute the convex hull
 * by evaluating the 4 Kharitonov polynomials at s = j*omega.
 *
 * @param ip         Interval polynomial
 * @param omega      Frequency (rad/s)
 * @return           Value set structure (caller must call pu_value_set_free)
 */
PU_ValueSet pu_value_set_compute(const PU_IntervalPolynomial *ip, double omega);

/**
 * Compute value set by dense gridding of the parameter space.
 * More accurate than Mapping Theorem alone for non-interval structures.
 *
 * @param ip            Interval polynomial
 * @param omega         Frequency
 * @param grid_points   Number of grid points per dimension
 * @return              Value set from gridding
 */
PU_ValueSet pu_value_set_grid(const PU_IntervalPolynomial *ip, double omega, int grid_points);

/** Free a value set structure */
void pu_value_set_free(PU_ValueSet *vs);

/** Print a value set summary */
void pu_value_set_print(const PU_ValueSet *vs);

/* --- Zero Exclusion Principle --- */

/**
 * Test robust stability using the Zero Exclusion Principle.
 *
 * Sweeps over a frequency grid and checks if 0 is in the value set
 * at any frequency. If 0 is never in the value set and the nominal
 * system is stable, the family is robustly stable.
 *
 * @param ip          Interval polynomial
 * @param omega_min   Minimum frequency to check
 * @param omega_max   Maximum frequency to check
 * @param n_points    Number of frequency grid points
 * @return            PU_STABLE or PU_UNSTABLE
 */
PU_StabilityStatus pu_zero_exclusion_test(const PU_IntervalPolynomial *ip,
                                           double omega_min, double omega_max,
                                           int n_points);

/**
 * Find the critical frequency where the value set first touches the origin.
 *
 * @param ip            Interval polynomial
 * @param omega_crit    Output: critical frequency
 * @param omega_min     Search start
 * @param omega_max     Search end
 * @param n_points      Grid points
 * @return              Distance to origin at critical frequency (> 0 => stable)
 */
double pu_zero_exclusion_margin(const PU_IntervalPolynomial *ip,
                                 double *omega_crit,
                                 double omega_min, double omega_max,
                                 int n_points);

/* --- Polytopic Stability (Edge Theorem) --- */

/**
 * Create a polytopic polynomial family from vertex polynomials.
 *
 * @param degree       Degree of all vertex polynomials
 * @param n_vertices   Number of vertex polynomials
 * @param vertices     Array of coefficient vectors (n_vertices x (degree+1))
 * @return             Polynomial family (polytopic)
 */
PU_PolytopicPolynomial pu_polyfamily_create(int degree, int n_vertices, double **vertices);

/** Free polytopic polynomial */
void pu_polyfamily_free(PU_PolytopicPolynomial *pp);

/**
 * Test robust stability using the Edge Theorem.
 *
 * The edge from vertex i to vertex j is parameterized as:
 *   P(s, lambda) = (1-lambda)*P_i(s) + lambda*P_j(s),  lambda in [0,1]
 *
 * The edge is stable iff P(s, lambda) is Hurwitz for all lambda in [0,1].
 * The polytope is stable iff all edges are stable and one vertex is stable.
 *
 * @param pp          Polynomial family (polytopic)
 * @param n_divs      Number of divisions along each edge for gridding
 * @return            PU_STABLE or PU_UNSTABLE
 */
PU_StabilityStatus pu_edge_theorem_test(const PU_PolytopicPolynomial *pp, int n_divs);

/**
 * Check if a single edge of the polytope is Hurwitz stable.
 *
 * @param pp          Polynomial family (polytopic)
 * @param i           Index of first vertex
 * @param j           Index of second vertex
 * @param n_divs      Number of divisions along the edge
 * @return            True if edge is stable
 */
bool pu_edge_is_stable(const PU_PolytopicPolynomial *pp, int i, int j, int n_divs);

/* --- Lyapunov-Based Robust Stability --- */

/**
 * Construct the uncertain state-space system.
 *
 * @param n_states      Number of states
 * @param n_inputs      Number of inputs
 * @param n_outputs     Number of outputs
 * @param A_nom         Nominal A matrix (n_states x n_states)
 * @param B_nom         Nominal B matrix (n_states x n_inputs)
 * @param C_nom         Nominal C matrix (n_outputs x n_states)
 * @param D_nom         Nominal D matrix (n_outputs x n_inputs)
 * @param unc_type      Uncertainty type
 * @return              Uncertain state-space system
 */
PU_UncertainStateSpace pu_uss_create(int n_states, int n_inputs, int n_outputs,
                                      double **A_nom, double **B_nom,
                                      double **C_nom, double **D_nom,
                                      PU_UncertaintyType unc_type);

/** Free uncertain state-space system */
void pu_uss_free(PU_UncertainStateSpace *uss);

/**
 * Add a perturbation direction matrix for parameter q_i.
 * The actual A matrix is: A(q) = A_nominal + sum_i q_i * A_pert_i
 *
 * @param uss         Uncertain state-space
 * @param A_pert_i    Perturbation direction for A (n_states x n_states)
 * @param B_pert_i    Perturbation direction for B (n_states x n_inputs)
 * @return            0 on success, -1 on error
 */
int pu_uss_add_perturbation(PU_UncertainStateSpace *uss,
                             double **A_pert_i, double **B_pert_i);

/**
 * Evaluate state-space matrices at a specific parameter value.
 *
 * @param uss         Uncertain state-space
 * @param q           Parameter vector (length n_unc_params)
 * @param A_out       Output A matrix (n_states x n_states)
 * @param B_out       Output B matrix (n_states x n_inputs)
 */
void pu_uss_eval(const PU_UncertainStateSpace *uss, const double *q,
                 double **A_out, double **B_out);

/**
 * Test quadratic stability: find P > 0 such that A_i'*P + P*A_i < 0
 * for all vertex systems A_i.
 *
 * Uses a simple alternating projection method to solve the LMI feasibility.
 *
 * @param uss         Uncertain state-space (must be polytopic type)
 * @param result      Output: LMI feasibility result
 * @return            PU_STABLE if quadratically stable, PU_UNSTABLE or INDETERMINATE otherwise
 */
PU_StabilityStatus pu_quadratic_stability_test(const PU_UncertainStateSpace *uss,
                                                PU_QuadraticStabilityLMI *result);

/**
 * Compute the quadratic stability margin: the largest gamma such that
 * A_i'*P + P*A_i < 0 for all A_i and some P > 0, with A_i scaled by gamma.
 */
double pu_quadratic_stability_margin(const PU_UncertainStateSpace *uss);

/**
 * Solve the continuous-time Lyapunov equation and return results.
 *
 * For a stable nominal A, find P such that A'*P + P*A + Q = 0.
 *
 * @param A         Hurwitz stable matrix (n x n)
 * @param Q         Positive definite matrix (n x n)
 * @param n         Matrix dimension
 * @return          Lyapunov data with P matrix
 */
PU_LyapunovData pu_lyapunov_solve_ct(const double **A, double **Q, int n);

/** Free Lyapunov data */
void pu_lyapunov_data_free(PU_LyapunovData *ld);

/* --- Robust Stability Margin --- */

/**
 * Compute the robust stability radius: maximum uncertainty under which
 * the system remains stable.
 *
 * Uses frequency gridding and eigenvalue analysis.
 *
 * @param uss             Uncertain state-space system
 * @param omega_min       Minimum frequency
 * @param omega_max       Maximum frequency
 * @param n_freq_points   Number of frequency points
 * @return                Robust stability result
 */
PU_RobustStabilityResult pu_robust_stability_radius(const PU_UncertainStateSpace *uss,
                                                     double omega_min, double omega_max,
                                                     int n_freq_points);

/** Free robust stability result */
void pu_robust_stability_result_free(PU_RobustStabilityResult *result);

/** Print robust stability result */
void pu_robust_stability_result_print(const PU_RobustStabilityResult *result);

/* --- Interval Matrix Stability --- */

/**
 * Create interval matrix from lower and upper bounds.
 *
 * @param rows        Number of rows
 * @param cols        Number of cols
 * @param lower       Lower bound matrix
 * @param upper       Upper bound matrix
 * @return            Interval matrix
 */
PU_IntervalMatrix pu_interval_matrix_create(int rows, int cols,
                                              double **lower, double **upper);

/** Free interval matrix */
void pu_interval_matrix_free(PU_IntervalMatrix *im);

/**
 * Check stability of interval matrix.
 * Uses necessary condition (vertex check) and sufficient condition
 * (Gershgorin circles / norm bounds).
 */
PU_StabilityStatus pu_interval_matrix_stability(const PU_IntervalMatrix *im);

/**
 * Sample a matrix uniformly from the interval matrix.
 */
void pu_interval_matrix_sample(const PU_IntervalMatrix *im, double **out);

/**
 * Print interval matrix.
 */
void pu_interval_matrix_print(const PU_IntervalMatrix *im);

#endif /* PU_STABILITY_H */

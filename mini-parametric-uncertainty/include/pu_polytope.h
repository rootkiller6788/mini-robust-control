#ifndef PU_POLYTOPE_H
#define PU_POLYTOPE_H

#include "pu_core.h"

/* ============================================================================
 * Polytopic Uncertainty Models and Analysis
 *
 * Polytopic uncertainty: the uncertain parameters q belong to a convex
 * polytope Q in R^p. The system matrices are affine in q:
 *
 *   A(q) = A0 + sum_{i=1}^{p} q_i * A_i
 *   B(q) = B0 + sum_{i=1}^{p} q_i * B_i
 *
 * where q in conv{q^1, q^2, ..., q^N} and each q^j is a vertex.
 *
 * The image of the polytope under an affine map is also a polytope,
 * which greatly simplifies robust analysis.
 *
 * Vertex Enumeration Problem: Given an affine system, the stability
 * of all vertex systems implies stability of the entire polytope
 * ONLY for certain properties (e.g., quadratic stability, H-infinity norm).
 *
 * Key References:
 *   Boyd, El Ghaoui, Feron, Balakrishnan (1994) — Linear Matrix Inequalities
 *   Barmish (1994) — New Tools for Robustness
 * ============================================================================ */

/* --- Polytope in Parameter Space --- */

typedef struct {
    int      dim;              /* Dimension of parameter space (p) */
    int      n_vertices;       /* Number of polytope vertices (N) */
    double **vertices;         /* Vertex coordinates (N x p) */
    double  *center;           /* Center of polytope (dim entries) */
    double **A_vertices;       /* System matrices at vertices, size N * (n*n) packed */
    double **B_vertices;       /* Input matrices at vertices */
    int      n_states;         /* System order */
    int      n_inputs;         /* Number of inputs */
} PU_Polytope;

/** Create a polytope from vertices */
PU_Polytope pu_polytope_create(int dim, int n_vertices, double **vertices);

/** Free polytope */
void pu_polytope_free(PU_Polytope *p);

/** Compute the center of the polytope */
void pu_polytope_compute_center(PU_Polytope *p);

/** Check if a point is inside the polytope (via barycentric coordinates) */
bool pu_polytope_contains(const PU_Polytope *p, const double *point);

/** Sample a point uniformly from the polytope (Dirichlet distribution) */
void pu_polytope_sample(const PU_Polytope *p, double *point);

/** Compute the diameter of the polytope */
double pu_polytope_diameter(const PU_Polytope *p);

/** Print polytope information */
void pu_polytope_print(const PU_Polytope *p);

/* --- Affine Parameter-Dependent State-Space --- */

typedef struct {
    int       n_states;
    int       n_inputs;
    int       n_outputs;
    int       n_params;            /* Number of uncertain parameters */
    double  **A0;                  /* Nominal A */
    double  **B0;                  /* Nominal B */
    double  **C0;                  /* Nominal C */
    double  **D0;                  /* Nominal D */
    double ***A_param;             /* A_i matrices: n_params x n x n */
    double ***B_param;             /* B_i matrices: n_params x n x m */
    double  **param_bounds;        /* Parameter bounds: n_params x 2 (min,max) */
} PU_AffineStateSpace;

/** Create affine parameter-dependent state-space system */
PU_AffineStateSpace pu_affine_ss_create(int n_states, int n_inputs, int n_outputs,
                                          int n_params,
                                          double **A0, double **B0,
                                          double **C0, double **D0);

/** Set the parameter matrices A_i and B_i for parameter i */
void pu_affine_ss_set_param(PU_AffineStateSpace *ass, int param_idx,
                              double **Ai, double **Bi);

/** Set parameter bounds */
void pu_affine_ss_set_bounds(PU_AffineStateSpace *ass, int param_idx,
                               double lower, double upper);

/** Evaluate A(q) and B(q) at parameter q */
void pu_affine_ss_eval(const PU_AffineStateSpace *ass, const double *q,
                        double **A_out, double **B_out);

/** Free affine state-space */
void pu_affine_ss_free(PU_AffineStateSpace *ass);

/** Print affine state-space system */
void pu_affine_ss_print(const PU_AffineStateSpace *ass);

/* --- Polytope Operations --- */

/** Compute convex hull vertices from a set of points (Gift Wrapping for dim <= 3) */
int pu_convex_hull_2d(const double *points, int n_points,
                       int *hull_indices, int max_hull);

/** Compute bounding box of a polytope */
void pu_polytope_bounding_box(const PU_Polytope *p,
                                double *lower, double *upper);

/** Check if two polytopes intersect */
bool pu_polytope_intersect(const PU_Polytope *p1, const PU_Polytope *p2);

/** Compute the Chebyshev center (maximum inscribed ball) of a polytope
 *  described by A*x <= b. Uses linear programming via simplified simplex. */
int pu_chebyshev_center(int m, int n, double **A, double *b,
                         double *center, double *radius);

/* --- Branch and Bound for Robust Stability --- */

/**
 * Use branch-and-bound to determine if a polytopic system is robustly stable.
 * Recursively subdivides the parameter box and checks stability at each
 * sub-box using a sufficient condition.
 *
 * @param ass          Affine state-space system
 * @param max_depth    Maximum subdivision depth
 * @param tol          Tolerance for subdivision
 * @return             PU_STABLE if all sub-boxes are stable
 */
PU_StabilityStatus pu_branch_and_bound_stability(const PU_AffineStateSpace *ass,
                                                   int max_depth, double tol);

#endif /* PU_POLYTOPE_H */

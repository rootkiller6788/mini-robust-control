#ifndef PU_CORE_H
#define PU_CORE_H

#include <stdbool.h>
#include <stddef.h>

/* ============================================================================
 * Parametric Uncertainty Core Types and Utilities
 *
 * Foundations of robust control under parametric uncertainty.
 * Key contributions from:
 *   Kharitonov (1978) — Interval polynomial stability
 *   Bartlett, Hollot, Lin (1988) — Edge Theorem
 *   Barmish (1994) — New Tools for Robustness of Linear Systems
 *   Doyle, Francis, Tannenbaum (1992) — Feedback Control Theory
 *   Zhou, Doyle, Glover (1996) — Robust and Optimal Control
 *   Ackermann (1993) — Robust Control: The Parameter Space Approach
 *   Bhattacharyya, Chapellat, Keel (1995) — Robust Control: The Parametric Approach
 * ============================================================================ */

/* --- Fundamental Enumerations --- */

/** Type of parametric uncertainty structure */
typedef enum {
    PU_UNCERTAINTY_INTERVAL = 0,    /* Each parameter in [q_min, q_max] independently */
    PU_UNCERTAINTY_POLYTOPE = 1,    /* Convex hull of vertex systems */
    PU_UNCERTAINTY_NORM_BOUNDED = 2,/* ||Δ|| ≤ γ for some norm */
    PU_UNCERTAINTY_AFFINE = 3,      /* Affine dependence A(q) = A0 + Σ q_i A_i */
    PU_UNCERTAINTY_MULTILINEAR = 4, /* Multilinear in parameters */
    PU_UNCERTAINTY_LFT = 5          /* Linear Fractional Transformation representation */
} PU_UncertaintyType;

/** Stability status for uncertain systems */
typedef enum {
    PU_STABLE = 0,                  /* Robustly stable for all uncertainty */
    PU_UNSTABLE = 1,                /* Unstable for some admissible uncertainty */
    PU_MARGINALLY_STABLE = 2,       /* On the stability boundary */
    PU_INDETERMINATE = 3            /* Cannot determine with current test */
} PU_StabilityStatus;

/** Polynomial family type */
typedef enum {
    PU_POLY_INTERVAL = 0,           /* Independent interval coefficients */
    PU_POLY_POLYTOPE = 1,           /* Convex combination of vertex polynomials */
    PU_POLY_AFFINE = 2,             /* Coefficients affine in parameters */
    PU_POLY_MULTILINEAR = 3         /* Coefficients multilinear in parameters */
} PU_PolynomialFamilyType;

/** Robust performance measure type */
typedef enum {
    PU_PERF_H2 = 0,                 /* H2 norm performance */
    PU_PERF_HINF = 1,               /* H-infinity norm performance */
    PU_PERF_WORST_CASE = 2,         /* Worst-case over uncertainty set */
    PU_PERF_GUARANTEED_COST = 3     /* Guaranteed cost bound */
} PU_PerformanceMeasure;

/** Parameter distribution type */
typedef enum {
    PU_PARAM_UNIFORM = 0,           /* Uniform distribution */
    PU_PARAM_NORMAL = 1,            /* Normal/gaussian distribution */
    PU_PARAM_BETA = 2,              /* Beta distribution (bounded) */
    PU_PARAM_DISCRETE = 3           /* Discrete set of values */
} PU_ParameterDistribution;

/* --- Core Data Structures --- */

/** Single uncertain parameter */
typedef struct {
    char   *name;                   /* Parameter identifier */
    double  nominal;                /* Nominal value */
    double  lower;                  /* Lower bound (interval) or min */
    double  upper;                  /* Upper bound (interval) or max */
    double  relative_uncertainty;   /* Relative uncertainty: (upper-lower)/|nominal| */
    int     is_fixed;               /* 1 if parameter is known exactly */
    PU_ParameterDistribution dist;  /* Distribution type for stochastic analysis */
    double  dist_param1;            /* Distribution parameter 1 */
    double  dist_param2;            /* Distribution parameter 2 */
} PU_Parameter;

/** Vector of uncertain parameters */
typedef struct {
    int          n_params;          /* Number of uncertain parameters */
    PU_Parameter *params;           /* Array of parameters */
    char         *description;      /* Human-readable description */
} PU_ParameterVector;

/** Interval polynomial: P(s) = sum_{i=0}^{n} [a_i^-, a_i^+] s^i */
typedef struct {
    int      degree;                /* Polynomial degree n */
    double  *coeff_lower;           /* Lower bound coefficients a_i^- (length n+1) */
    double  *coeff_upper;           /* Upper bound coefficients a_i^+ (length n+1) */
    double  *coeff_nominal;         /* Nominal coefficients */
    bool     is_monic;              /* True if leading coefficient = 1 */
} PU_IntervalPolynomial;

/** Kharitonov polynomial (one of 4 corner polynomials) */
typedef struct {
    int      degree;                /* Polynomial degree */
    double  *coeff;                 /* Fixed coefficients (length n+1) */
    int      kharitonov_index;      /* 0=K1(--), 1=K2(-+), 2=K3(+-), 3=K4(++) */
    bool     is_hurwitz;            /* Result of Routh-Hurwitz test */
} PU_KharitonovPolynomial;

/** Vertex polynomial for polytopic family */
typedef struct {
    int      degree;
    double  *coeff;
    int      vertex_index;
    double  *barycentric;           /* Barycentric coordinates */
} PU_VertexPolynomial;

/** Polytopic polynomial family: P(s,lambda) = Sum lambda_i P_i(s), lambda_i >= 0, Sum lambda_i = 1 */
typedef struct {
    int                    degree;
    int                    n_vertices;    /* Number of vertex polynomials */
    PU_VertexPolynomial   *vertices;      /* Vertex polynomials */
    int                    n_edges;       /* Number of edges (combinatorial) */
    int                   *edge_pairs;    /* Pairs (i,j) of connected vertices, length 2*n_edges */
    double                *lambda_center; /* Center of the polytope */
} PU_PolytopicPolynomial;

/** Continuous-time uncertain state-space system: dx = A(q)x + B(q)u, y = C(q)x + D(q)u */
typedef struct {
    int      n_states;                /* Number of states */
    int      n_inputs;                /* Number of inputs */
    int      n_outputs;               /* Number of outputs */
    double **A_nominal;               /* Nominal A matrix (n x n) */
    double **B_nominal;               /* Nominal B matrix (n x m) */
    double **C_nominal;               /* Nominal C matrix (p x n) */
    double **D_nominal;               /* Nominal D matrix (p x m) */
    PU_UncertaintyType unc_type;      /* Type of uncertainty structure */
    int      n_unc_params;            /* Number of uncertain parameters */
    double ***A_pert;                 /* Perturbation direction matrices: n_params x n x n */
    double ***B_pert;                 /* Perturbation direction matrices: n_params x n x m */
    PU_ParameterVector params;        /* The uncertain parameters */
    bool     is_continuous;           /* True for continuous-time, false for discrete */
} PU_UncertainStateSpace;

/** Interval matrix: A in [A_lower, A_upper] elementwise */
typedef struct {
    int      rows;
    int      cols;
    double **lower;                   /* Lower bound matrix */
    double **upper;                   /* Upper bound matrix */
} PU_IntervalMatrix;

/** Robust stability margin result */
typedef struct {
    double   stability_radius;        /* Max gamma such that stable for ||Delta|| <= gamma */
    double   worst_case_real_part;    /* Max real part of any eigenvalue in uncertainty set */
    double   h_infinity_norm;         /* H-infinity norm of the worst-case transfer function */
    double  *mu_bounds;               /* Upper/lower bounds for structured singular value */
    int      n_frequency_points;      /* Number of frequency grid points used */
    double  *frequency_grid;          /* Frequency grid for analysis */
    PU_StabilityStatus status;        /* Final stability verdict */
    char    *method_used;             /* Description of method */
} PU_RobustStabilityResult;

/** Value set of a polynomial family at a specific frequency s=j*omega */
typedef struct {
    double   omega;                   /* Frequency omega */
    int      n_points;                /* Number of points in the value set */
    double  *real_part;               /* Real parts of P(j*omega) */
    double  *imag_part;               /* Imaginary parts of P(j*omega) */
    double   center_real;             /* Center of value set (real) */
    double   center_imag;             /* Center of value set (imag) */
    double   radius;                  /* Approximate radius */
    bool     contains_origin;         /* True if 0+j0 in value set */
    double   min_distance_to_origin;  /* Minimum distance from origin to value set */
} PU_ValueSet;

/** Routh array for stability testing */
typedef struct {
    int      degree;                  /* Polynomial degree */
    int      n_rows;                  /* Number of rows in Routh array */
    double **array;                   /* Routh array (n_rows x variable width) */
    int      sign_changes;            /* Number of sign changes in first column */
    bool     is_stable;               /* True if all first column elements > 0 */
} PU_RouthArray;

/** Lyapunov-based robust stability data */
typedef struct {
    int      n_states;
    double **P;                       /* Lyapunov matrix P > 0 */
    double **Q;                       /* Q = -(A'*P + P*A) > 0 */
    double   min_eigenvalue_P;        /* lambda_min(P) for stability margin */
    double   condition_number;        /* kappa(P) = lambda_max/lambda_min */
} PU_LyapunovData;

/** Quadratic stability LMI data: find P > 0 s.t. A_i'*P + P*A_i < 0 for all i */
typedef struct {
    int      n_states;
    int      n_vertices;              /* Number of vertex systems */
    double **P;                       /* Common Lyapunov matrix */
    double   t_min;                   /* Optimal value of auxiliary variable */
    bool     is_feasible;             /* True if LMI is feasible */
    int      iterations;              /* Number of iterations for solver */
} PU_QuadraticStabilityLMI;

/* --- Core API --- */

/** Create a new uncertain parameter */
PU_Parameter pu_param_create(const char *name, double nominal,
                              double lower, double upper);

/** Free a parameter (free internal string) */
void pu_param_free(PU_Parameter *p);

/** Create parameter vector from array of parameters */
PU_ParameterVector pu_param_vector_create(int n, const PU_Parameter *params);

/** Free parameter vector */
void pu_param_vector_free(PU_ParameterVector *pv);

/** Get effective range of a parameter */
void pu_param_effective_range(const PU_Parameter *p, double *eff_lower, double *eff_upper);

/** Sample parameter uniformly from its range */
double pu_param_sample(const PU_Parameter *p);

/** Gaussian sample with given mean and std */
double pu_gaussian_sample(double mean, double std);

/** Print parameter info */
void pu_param_print(const PU_Parameter *p);

/* --- Interval Polynomial API --- */

/** Create an interval polynomial of given degree */
PU_IntervalPolynomial pu_interval_poly_create(int degree);

/** Set coefficient bounds for a specific index */
void pu_interval_poly_set_coeff(PU_IntervalPolynomial *poly, int index,
                                 double lower, double upper);

/** Set monic (leading coefficient fixed to 1) */
void pu_interval_poly_set_monic(PU_IntervalPolynomial *poly, bool monic);

/** Free interval polynomial */
void pu_interval_poly_free(PU_IntervalPolynomial *poly);

/** Evaluate interval polynomial at s = a + jb, returns complex value bounds */
void pu_interval_poly_eval(const PU_IntervalPolynomial *poly, double a, double b,
                            double *real_lower, double *real_upper,
                            double *imag_lower, double *imag_upper);

/** Sample a specific polynomial from the interval family */
void pu_interval_poly_sample(const PU_IntervalPolynomial *poly, double *coeff_out);

/** Print interval polynomial (symbolic form) */
void pu_interval_poly_print(const PU_IntervalPolynomial *poly);

/* --- Matrix Utility API --- */

/** Allocate an n x m double matrix, zero-initialized */
double** pu_matrix_alloc(int rows, int cols);

/** Free an n x m matrix */
void pu_matrix_free(double **mat, int rows);

/** Allocate identity matrix of size n */
double** pu_matrix_eye(int n);

/** Copy matrix src to dst (both n x m) */
void pu_matrix_copy(double **dst, double **src, int rows, int cols);

/** Matrix addition C = A + B */
void pu_matrix_add(double **A, double **B, double **C, int rows, int cols);

/** Matrix subtraction C = A - B */
void pu_matrix_sub(double **A, double **B, double **C, int rows, int cols);

/** Matrix multiplication C = A * B (A: m x k, B: k x n, C: m x n) */
void pu_matrix_mul(double **A, double **B, double **C, int m, int k, int n);

/** Matrix scaling B = alpha * A */
void pu_matrix_scale(double **A, double **B, double alpha, int rows, int cols);

/** Transpose B = A' (A: rows x cols, B: cols x rows) */
void pu_matrix_transpose(double **A, double **B, int rows, int cols);

/** Print matrix to stdout */
void pu_matrix_print(double **mat, int rows, int cols, const char *name);

/** Compute eigenvalues of a 2x2 real matrix (real parts in re, imag in im) */
void pu_eigen_2x2(double **A, double *re1, double *im1, double *re2, double *im2);

/** Compute eigenvalues using QR algorithm (simple, for small matrices) */
int  pu_eigen_qr(double **A, int n, double *real_part, double *imag_part,
                 int max_iter, double tol);

/** Compute spectral radius rho(A) = max |lambda_i| */
double pu_spectral_radius(double **A, int n);

/** Compute matrix norm (Frobenius) */
double pu_matrix_norm_fro(double **A, int rows, int cols);

/** Solve Lyapunov equation A'*P + P*A + Q = 0 (Bartels-Stewart simplified for n <= 10) */
int  pu_lyapunov_solve(double **A, double **Q, double **P, int n);

/** Check if matrix P is positive definite (via Cholesky attempt) */
bool pu_is_positive_definite(double **P, int n);

/** Routh-Hurwitz stability test for a fixed-coefficient polynomial */
PU_RouthArray pu_routh_hurwitz(const double *coeff, int degree);

/** Print Routh array */
void pu_routh_print(const PU_RouthArray *ra);

/** Free Routh array */
void pu_routh_free(PU_RouthArray *ra);

/** Check if all roots have negative real part (continuous-time Hurwitz stable) */
bool pu_is_hurwitz_stable(const double *coeff, int degree);

/** Check if all roots are inside unit circle (discrete-time Schur stable) */
bool pu_is_schur_stable(const double *coeff, int degree);

/* --- Constants --- */
#define PU_PI         3.14159265358979323846
#define PU_EPS        1e-12
#define PU_INF        1e308
#define PU_MAX_DEGREE 128
#define PU_MAX_STATES 128
#define PU_MAX_ITER_QR 1000
#define PU_QR_TOL     1e-10

#endif /* PU_CORE_H */

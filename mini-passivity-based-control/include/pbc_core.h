#ifndef PBC_CORE_H
#define PBC_CORE_H

#include <stdbool.h>
#include <stddef.h>

/* ==============================================================
 * pbc_core.h — Passivity-Based Control Core Definitions
 *
 * Passivity-Based Control (PBC) exploits the energy properties
 * of physical systems to design stabilizing controllers.
 *
 * Core concept (Willems 1972, Hill-Moylan 1976):
 *   A system is passive if there exists a non-negative storage
 *   function S(x) such that:
 *     S(x(t)) - S(x(0)) ≤ ∫₀ᵗ uᵀ(s)y(s) ds
 *
 * Equivalent differential form:
 *     dS/dt ≤ uᵀ y         (dissipation inequality)
 *
 * For strictly passive systems:
 *     dS/dt ≤ uᵀ y - φ(x)  with φ(x) > 0 for x ≠ 0
 *
 * Supply rate: w(t) = uᵀ(t) y(t) (input-output power)
 *
 * Important special cases:
 *   - Electrical: u = voltage, y = current, S = stored energy
 *   - Mechanical: u = force/torque, y = velocity
 *   - Robotics: u = joint torque, y = joint velocity
 *
 * References:
 *   Willems (1972) "Dissipative Dynamical Systems"
 *   Hill & Moylan (1976) "Stability of Nonlinear Dissipative Systems"
 *   Ortega et al. (1998) "Passivity-Based Control of Euler-Lagrange Systems"
 *   van der Schaft (2000) "L2-Gain and Passivity Techniques"
 *   Khalil (2002) "Nonlinear Systems", Chapter 6
 *   Brogliato et al. (2007) "Dissipative Systems Analysis and Control"
 * ============================================================== */

/* --- Vector in R^n --- */
typedef struct {
    double* data;
    int    size;
} PBC_Vector;

/* --- Dense matrix (row-major, flat array) --- */
typedef struct {
    double* data;
    int    rows;
    int    cols;
} PBC_Matrix;

/* --- Passivity type classification --- */
typedef enum {
    PBC_PASSIVE_NONE      = 0,  /* Not passive */
    PBC_PASSIVE           = 1,  /* Passive (lossless or dissipative) */
    PBC_PASSIVE_STRICT_INPUT   = 2,  /* Input strictly passive */
    PBC_PASSIVE_STRICT_OUTPUT  = 3,  /* Output strictly passive */
    PBC_PASSIVE_STRICT         = 4,  /* Strictly passive (both) */
    PBC_PASSIVE_LOSSLESS       = 5,  /* Lossless (dS/dt = u^T y) */
    PBC_PASSIVE_FINITE_GAIN    = 6   /* Finite L2-gain */
} PBC_PassivityType;

/* --- Storage function types --- */
typedef enum {
    PBC_STORAGE_QUADRATIC      = 0,  /* S(x) = x^T P x, P > 0 */
    PBC_STORAGE_KINETIC        = 1,  /* S(q̇) = ½ q̇ᵀ M(q) q̇ */
    PBC_STORAGE_TOTAL_ENERGY   = 2,  /* S = KE + PE */
    PBC_STORAGE_HAMILTONIAN    = 3,  /* H(q,p) = ½ pᵀ M⁻¹(q) p + V(q) */
    PBC_STORAGE_CUSTOM         = 4   /* User-defined storage function */
} PBC_StorageType;

/* --- Storage function descriptor --- */
typedef struct {
    PBC_StorageType type;
    PBC_Matrix*     P;        /* Quadratic weight matrix (n×n) */
    double          gamma;    /* Damping coefficient */
    double          min_value; /* Energy at equilibrium */
    bool            is_positive_definite;
    bool            is_radially_unbounded;
    int             n_states;
} PBC_StorageFunction;

/* --- Supply rate --- */
typedef struct {
    double* u;    /* Input vector */
    double* y;    /* Output vector */
    int     dim;  /* Dimension (inputs = outputs for square systems) */
    double  value; /* Computed supply rate w = uᵀ y */
} PBC_SupplyRate;

/* --- Dissipation inequality result --- */
typedef struct {
    double S_current;   /* Current storage value */
    double S_initial;   /* Initial storage value */
    double supply_integral;  /* ∫₀ᵗ w(s) ds */
    double dissipation; /* S(t) - S(0) - ∫ w ≤ 0 for passivity */
    bool   satisfied;   /* Whether dissipation inequality holds */
} PBC_Dissipation;

/* --- System state descriptor --- */
typedef struct {
    double* x;       /* State vector */
    double* u;       /* Current input */
    double* y;       /* Current output */
    int     n;       /* Number of states */
    int     m;       /* Number of inputs */
    int     p;       /* Number of outputs */
    double  t;       /* Current time */
    double  dt;      /* Time step */
} PBC_SystemState;

/* --- Generic passive system structure --- */
typedef struct PBC_System PBC_System;

typedef void (*PBC_StateFunc)(const PBC_System* sys, const double* x,
                               const double* u, double* xdot);
typedef void (*PBC_OutputFunc)(const PBC_System* sys, const double* x,
                                double* y);
typedef double (*PBC_StorageEval)(const PBC_System* sys, const double* x);

struct PBC_System {
    char*            name;
    int              n;        /* State dimension */
    int              m;        /* Input dimension */
    int              p;        /* Output dimension */
    PBC_PassivityType passivity_type;
    PBC_StorageFunction storage;
    PBC_StateFunc    f;        /* xdot = f(x, u) */
    PBC_OutputFunc   h;        /* y = h(x) */
    PBC_StorageEval  S;        /* S = S(x) */
    void*            data;     /* System-specific data */
    double*          x;        /* Current state */
    double           t;        /* Current time */
    bool             is_initialized;
};

/* --- Passivity test parameters (KYP-like) --- */
typedef struct {
    double epsilon;     /* Strict passivity margin */
    double delta;       /* Output strictness margin */
    double nu;          /* Input strictness margin */
    double gamma;       /* L2-gain bound */
    int    max_iter;    /* Max iterations for numerical checks */
    double tol;         /* Convergence tolerance */
} PBC_PassivityParams;

/* ==============================================================
 * Matrix / Vector API
 * ============================================================== */

/** Create a zero-initialized matrix. Returns NULL on allocation failure. */
PBC_Matrix* pbc_matrix_create(int rows, int cols);

/** Free a matrix allocated by pbc_matrix_create. */
void pbc_matrix_free(PBC_Matrix* m);

/** Create a deep copy of a matrix. */
PBC_Matrix* pbc_matrix_copy(const PBC_Matrix* m);

/** Create n×n identity matrix. */
PBC_Matrix* pbc_matrix_identity(int n);

/** Create zero-filled matrix. */
PBC_Matrix* pbc_matrix_zeros(int rows, int cols);

/** Set element at (row, col). Bounds-checked. */
void pbc_matrix_set(PBC_Matrix* m, int row, int col, double val);

/** Get element at (row, col). Returns 0.0 if out of bounds. */
double pbc_matrix_get(const PBC_Matrix* m, int row, int col);

/** Fill all elements with given value. */
void pbc_matrix_fill(PBC_Matrix* m, double val);

/** Matrix addition: C = A + B. Returns NULL if dimensions mismatch. */
PBC_Matrix* pbc_matrix_add(const PBC_Matrix* a, const PBC_Matrix* b);

/** Matrix subtraction: C = A - B. */
PBC_Matrix* pbc_matrix_sub(const PBC_Matrix* a, const PBC_Matrix* b);

/** Matrix multiplication: C = A * B. Dimensions must match. */
PBC_Matrix* pbc_matrix_mul(const PBC_Matrix* a, const PBC_Matrix* b);

/** Scalar multiplication: B = alpha * A. */
PBC_Matrix* pbc_matrix_scale(const PBC_Matrix* a, double alpha);

/** Matrix transpose: B = Aᵀ. */
PBC_Matrix* pbc_matrix_transpose(const PBC_Matrix* a);

/** Frobenius norm ||A||_F = sqrt(∑ aᵢⱼ²). */
double pbc_matrix_norm_frobenius(const PBC_Matrix* a);

/** Trace of a square matrix. Returns NAN if not square. */
double pbc_matrix_trace(const PBC_Matrix* a);

/** Determinant of 1×1, 2×2, or 3×3 matrix. Returns NAN otherwise. */
double pbc_matrix_det(const PBC_Matrix* a);

/** Check if matrix is symmetric (within tolerance). */
bool pbc_matrix_is_symmetric(const PBC_Matrix* a, double tol);

/** Check if matrix is positive definite (via Cholesky attempt). */
bool pbc_matrix_is_positive_definite(const PBC_Matrix* a);

/** Check if matrix is positive semi-definite. */
bool pbc_matrix_is_positive_semidefinite(const PBC_Matrix* a);

/** Check if matrix is skew-symmetric: Aᵀ = -A. */
bool pbc_matrix_is_skew_symmetric(const PBC_Matrix* a, double tol);

/** Print matrix with label to stdout. */
void pbc_matrix_print(const PBC_Matrix* m, const char* label);

/* --- Vector operations --- */

/** Create a zero-initialized vector. */
PBC_Vector* pbc_vector_create(int size);

/** Free a vector. */
void pbc_vector_free(PBC_Vector* v);

/** Set element. Bounds-checked. */
void pbc_vector_set(PBC_Vector* v, int idx, double val);

/** Get element. Returns 0.0 if out of bounds. */
double pbc_vector_get(const PBC_Vector* v, int idx);

/** Euclidean norm ||v||₂. */
double pbc_vector_norm(const PBC_Vector* v);

/** Dot product a·b. Returns 0.0 if sizes differ. */
double pbc_vector_dot(const PBC_Vector* a, const PBC_Vector* b);

/** Print vector with label. */
void pbc_vector_print(const PBC_Vector* v, const char* label);

/* ==============================================================
 * Storage Function API
 * ============================================================== */

/**
 * Initialize a quadratic storage function S(x) = xᵀ P x.
 * @param n  State dimension
 * @param P  n×n positive definite matrix (copied)
 * @return Initialized storage function descriptor
 */
PBC_StorageFunction* pbc_storage_quadratic_init(int n, const PBC_Matrix* P);

/**
 * Initialize a total-energy storage function.
 *   S(q, q̇) = ½ q̇ᵀ M q̇ + ½ qᵀ K q   (M = inertia, K = stiffness)
 * @param n  Number of generalized coordinates
 * @param M  n×n inertia matrix (positive definite)
 * @param K  n×n stiffness matrix (positive semi-definite)
 */
PBC_StorageFunction* pbc_storage_mechanical_init(int n,
    const PBC_Matrix* M, const PBC_Matrix* K);

/** Free a storage function descriptor. */
void pbc_storage_free(PBC_StorageFunction* sf);

/** Evaluate S(x) for given state vector. */
double pbc_storage_eval(const PBC_StorageFunction* sf, const double* x);

/** Compute gradient ∇S(x). Result is an n-vector (allocated by caller). */
void pbc_storage_gradient(const PBC_StorageFunction* sf,
                           const double* x, double* grad);

/** Compute Hessian ∇²S(x). Result is n×n (allocated by caller). */
void pbc_storage_hessian(const PBC_StorageFunction* sf,
                          const double* x, PBC_Matrix* H);

/* ==============================================================
 * System API
 * ============================================================== */

/**
 * Create a passive system descriptor.
 * @param name  Human-readable identifier
 * @param n     State dimension
 * @param m     Input dimension
 * @param p     Output dimension
 * @param f     State function: xdot = f(x, u)
 * @param h     Output function: y = h(x)
 * @param S     Storage function evaluator
 * @param data  System-specific void pointer
 */
PBC_System* pbc_system_create(const char* name, int n, int m, int p,
                               PBC_StateFunc f, PBC_OutputFunc h,
                               PBC_StorageEval S, void* data);

/** Free a system descriptor. */
void pbc_system_free(PBC_System* sys);

/** Set initial state. */
void pbc_system_set_state(PBC_System* sys, const double* x0);

/** Compute state derivative: xdot = f(x, u). xdot must be pre-allocated. */
void pbc_system_f(const PBC_System* sys, const double* x,
                   const double* u, double* xdot);

/** Compute output: y = h(x). y must be pre-allocated. */
void pbc_system_h(const PBC_System* sys, const double* x, double* y);

/** Compute supply rate: w = uᵀ y. */
double pbc_supply_rate(const double* u, const double* y, int dim);

/**
 * Numerically verify the dissipation inequality over a trajectory.
 * Uses trapezoidal integration for ∫ uᵀ y dt.
 *
 * @param sys       System descriptor
 * @param t_span     Array of time points (length n_pts)
 * @param u_traj     Input trajectory (n_pts × m, row-major)
 * @param n_pts      Number of trajectory points
 * @return Dissipation inequality result
 */
PBC_Dissipation pbc_verify_dissipation(const PBC_System* sys,
    const double* t_span, const double* u_traj, int n_pts);

/** State update with Forward Euler: x_{k+1} = x_k + dt · f(x_k, u_k). */
void pbc_step_euler(PBC_System* sys, const double* u, double dt);

/** State update with RK4. */
void pbc_step_rk4(PBC_System* sys, const double* u, double dt);

#endif /* PBC_CORE_H */

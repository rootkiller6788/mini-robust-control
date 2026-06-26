#ifndef PBC_IDA_H
#define PBC_IDA_H

#include "pbc_core.h"

/* ==============================================================
 * pbc_ida.h — Interconnection and Damping Assignment (IDA-PBC)
 *
 * IDA-PBC is a systematic methodology for stabilizing nonlinear
 * systems by assigning a desired port-Hamiltonian structure.
 *
 * Original system (Port-Hamiltonian form):
 *   dx/dt = [J(x) - R(x)] ∇H(x) + g(x) u
 *
 * Desired closed-loop system:
 *   dx/dt = [J_d(x) - R_d(x)] ∇H_d(x)
 *
 * where:
 *   J_d(x) = -J_dᵀ(x)   (desired interconnection, skew-symmetric)
 *   R_d(x) = R_dᵀ(x) ≥ 0 (desired damping, positive semi-definite)
 *   H_d(x)              (desired Hamiltonian, with minimum at x*)
 *
 * The IDA-PBC control law solves the Matching Equation:
 *   [J(x) - R(x)] ∇H(x) + g(x) u_PBC =
 *       [J_d(x) - R_d(x)] ∇H_d(x)
 *
 * Key design steps:
 *   1. Choose desired Hamiltonian H_d with minimum at equilibrium
 *   2. Assign desired interconnection J_d and damping R_d
 *   3. Solve for u_PBC from the matching equation
 *
 * Algebraic solution (Ortega et al. 2002):
 *   Let g⊥(x) be a full-rank left annihilator of g(x).
 *   Then the PDE constraint is:
 *     g⊥ [J∇H - (J_d - R_d)∇H_d] = 0
 *
 * References:
 *   Ortega, van der Schaft, Maschke, Escobar (2002)
 *     "Interconnection and Damping Assignment Passivity-Based
 *      Control of Port-Controlled Hamiltonian Systems"
 *   Ortega & Garcia-Canseco (2004)
 *     "Interconnection and Damping Assignment Passivity-Based
 *      Control: A Survey"
 *   van der Schaft & Jeltsema (2014)
 *     "Port-Hamiltonian Systems Theory: An Introductory Overview"
 *   Duindam, Macchelli, Stramigioli, Bruyninckx (2009)
 *     "Modeling and Control of Complex Physical Systems"
 * ============================================================== */

/* --- Desired interconnection and damping structure --- */
typedef struct {
    PBC_Matrix* J_d;   /* Desired interconnection matrix (skew-symmetric) */
    PBC_Matrix* R_d;   /* Desired damping matrix (positive semi-definite) */
    double*     x_star; /* Desired equilibrium point */
    int         n;     /* State dimension */
} PBC_IDA_Structure;

/* --- Desired Hamiltonian --- */
typedef struct {
    PBC_Matrix* M_d;        /* Desired inertia/mass matrix (positive definite) */
    PBC_Matrix* K_d;        /* Desired stiffness (positive definite) */
    double*     q_star;     /* Desired position */
    double*     p_star;     /* Desired momentum (typically 0) */
    int         n_dof;      /* Degrees of freedom */
    double      H_min;      /* Minimum Hamiltonian value (at equilibrium) */
} PBC_IDA_Hamiltonian;

/* --- IDA-PBC matching equation result --- */
typedef struct {
    double* u_pbc;          /* Computed control input */
    double  matching_error; /* Residual of matching equation */
    bool    solved;         /* Whether matching was successful */
    double  H_d_value;      /* Current desired Hamiltonian */
    double  H_d_dot;        /* dH_d/dt (should be ≤ 0) */
} PBC_IDA_Result;

/* --- IDA-PBC controller --- */
typedef struct {
    PBC_IDA_Structure    structure;
    PBC_IDA_Hamiltonian  hamiltonian;
    PBC_Matrix*          M;        /* Inertia matrix M(q) (model) */
    PBC_Matrix*          C;        /* Coriolis matrix C(q,q̇) (model) */
    double*              g;        /* Gravity vector g(q) (model) */
    int                  n;
    bool                 is_configured;
} PBC_IDA_Controller;

/* ==============================================================
 * IDA Structure and Hamiltonian Setup
 * ============================================================== */

/**
 * Initialize desired interconnection and damping structure.
 *
 * @param n       State dimension
 * @param J_d     Desired skew-symmetric interconnection (n×n, copied)
 * @param R_d     Desired damping (n×n, positive semi-definite, copied)
 * @param x_star  Desired equilibrium (n elements, copied)
 */
PBC_IDA_Structure* pbc_ida_structure_create(int n,
    const PBC_Matrix* J_d, const PBC_Matrix* R_d, const double* x_star);

/** Free IDA structure. */
void pbc_ida_structure_free(PBC_IDA_Structure* s);

/**
 * Initialize desired Hamiltonian for mechanical systems.
 * H_d(q, p) = ½ pᵀ M_d⁻¹ p + ½ (q - q*)ᵀ K_d (q - q*)
 *
 * @param n_dof   Degrees of freedom
 * @param M_d     Desired inertia (n×n, positive definite)
 * @param K_d     Desired stiffness (n×n, positive definite)
 * @param q_star  Desired position (n elements)
 */
PBC_IDA_Hamiltonian* pbc_ida_hamiltonian_create(int n_dof,
    const PBC_Matrix* M_d, const PBC_Matrix* K_d, const double* q_star);

/** Free IDA Hamiltonian. */
void pbc_ida_hamiltonian_free(PBC_IDA_Hamiltonian* Hd);

/* ==============================================================
 * IDA-PBC Controller
 * ============================================================== */

/**
 * Create an IDA-PBC controller for fully-actuated mechanical systems.
 *
 * Model: M(q) q̈ + C(q,q̇) q̇ + g(q) = τ
 *
 * The IDA-PBC solution for fully-actuated systems:
 *   τ = g(q) - K_d K_d_inv C(q,q̇) (q̇ - 0)
 *       - K_d (q - q*) - K_v q̇
 *
 * where K_d is the desired stiffness and K_v is velocity damping.
 *
 * @param n   Degrees of freedom
 * @param M   Inertia matrix (constant or configuration-dependent)
 * @param C   Coriolis matrix (constant for simple systems)
 * @param g   Gravity vector at current config
 * @param J_d Desired interconnection
 * @param R_d Desired damping
 * @param Hd  Desired Hamiltonian
 */
PBC_IDA_Controller* pbc_ida_controller_create(int n,
    const PBC_Matrix* M, const PBC_Matrix* C,
    const PBC_Matrix* J_d, const PBC_Matrix* R_d,
    const PBC_IDA_Hamiltonian* Hd);

/** Free IDA-PBC controller. */
void pbc_ida_controller_free(PBC_IDA_Controller* ctrl);

/**
 * Compute IDA-PBC control law.
 *
 * @param ctrl   IDA-PBC controller
 * @param q      Current position (n elements)
 * @param qdot   Current velocity (n elements)
 * @param g_vec  Gravity vector g(q) (n elements, can be NULL)
 * @param tau    Output control torque (n elements, allocated by caller)
 * @return IDA-PBC result with matching diagnostics
 */
PBC_IDA_Result pbc_ida_compute_control(PBC_IDA_Controller* ctrl,
    const double* q, const double* qdot, const double* g_vec, double* tau);

/**
 * Solve the matching equation for a simple 1-DOF system analytically.
 *
 * For 1-DOF mechanical: m q̈ + c q̇ + g(q) = u
 * Desired system: m_d q̈ + k_d (q - q*) + k_v q̇ = 0
 *
 * The control law is:
 *   u = g(q) + (m/m_d - 1)(k_v q̇ + k_d (q - q*)) - k_v q̇ - k_d (q - q*)
 *
 * @param m      Mass/inertia
 * @param c      Damping coefficient
 * @param md     Desired mass
 * @param kd     Desired stiffness
 * @param kv     Desired damping (velocity gain)
 * @param q      Current position
 * @param qdot   Current velocity
 * @param qstar  Desired position
 * @param g_val  Gravity force g(q)
 * @return Control force u
 */
double pbc_ida_1dof_solve(double m, double c, double md, double kd,
                           double kv, double q, double qdot,
                           double qstar, double g_val);

/**
 * IDA-PBC for n-DOF fully-actuated mechanical systems.
 *
 * Assumes: M(q) q̈ + C(q,q̇) q̇ + g(q) = τ
 *
 * Control law (Ortega et al. 2002, Eq. 32):
 *   τ = g(q) + [M(q) M_d⁻¹ - I] K_d (q - q*)
 *       + [M(q) M_d⁻¹ C(q,q̇) - C(q,q̇)] q̇
 *       - M(q) M_d⁻¹ K_v q̇
 *
 * @param M       Inertia matrix M(q) (n×n)
 * @param C       Coriolis matrix C(q,q̇) (n×n)
 * @param Md      Desired inertia (n×n, positive definite)
 * @param Kd      Desired stiffness (n×n, positive definite)
 * @param Kv      Velocity damping matrix (n×n, positive definite)
 * @param q       Current position
 * @param qdot    Current velocity
 * @param qstar   Desired position
 * @param g_vec   Gravity vector (n elements)
 * @param tau     Output control torque (n elements, allocated by caller)
 * @param n       Degrees of freedom
 */
void pbc_ida_mechanical_ndof(const PBC_Matrix* M, const PBC_Matrix* C,
    const PBC_Matrix* Md, const PBC_Matrix* Kd, const PBC_Matrix* Kv,
    const double* q, const double* qdot, const double* qstar,
    const double* g_vec, double* tau, int n);

/**
 * Verify the matching equation residual numerically.
 *
 * Matching condition:
 *   [J - R]∇H + g u = [J_d - R_d]∇H_d
 *
 * @param ctrl   IDA-PBC controller
 * @param q      Current position
 * @param qdot   Current velocity
 * @param u      Applied control
 * @param tol    Tolerance
 * @return Residual norm of the matching equation
 */
double pbc_ida_verify_matching(const PBC_IDA_Controller* ctrl,
    const double* q, const double* qdot, const double* u, double tol);

#endif /* PBC_IDA_H */

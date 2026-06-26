#ifndef PBC_ENERGY_SHAPING_H
#define PBC_ENERGY_SHAPING_H

#include "pbc_core.h"

/* ==============================================================
 * pbc_energy_shaping.h — Energy Shaping and Damping Injection
 *
 * The two fundamental operations of Passivity-Based Control:
 *
 * 1. ENERGY SHAPING:
 *    Modify the potential energy of the closed-loop system so
 *    the desired equilibrium becomes the global minimum.
 *
 *    For Euler-Lagrange systems:
 *      M(q) q̈ + C(q,q̇) q̇ + g(q) = τ
 *    Energy shaping control law:
 *      τ = g(q) - K_p (q - q_d) + τ_di
 *    where K_p shapes the closed-loop potential energy:
 *      V_d(q) = ½ (q - q_d)ᵀ K_p (q - q_d)
 *
 * 2. DAMPING INJECTION:
 *    Add virtual damping to ensure asymptotic stability.
 *
 *    τ_di = -K_d q̇    (linear damping injection)
 *    or more generally:
 *    τ_di = -Ψ(q̇)      (nonlinear damping, Ψ(0)=0, xᵀΨ(x) > 0)
 *
 * Combined PBC law:
 *   τ = g(q) - K_p (q - q_d) - K_d q̇
 *
 * Passivity Theorem (feedback interconnection):
 *   If H1 and H2 are passive, then their negative feedback
 *   interconnection is also passive and stable.
 *
 * References:
 *   Takegaki & Arimoto (1981) "A New Feedback Method for Dynamic
 *     Control of Manipulators"
 *   Ortega & Spong (1989) "Adaptive Motion Control of Rigid Robots:
 *     A Tutorial"
 *   Ortega et al. (1998) "Passivity-Based Control of Euler-Lagrange
 *     Systems"
 *   van der Schaft (2000) "L2-Gain and Passivity Techniques"
 *   Kelly, Santibañez, Loria (2005) "Control of Robot Manipulators
 *     in Joint Space"
 * ============================================================== */

/* --- Energy shaping controller parameters --- */
typedef struct {
    PBC_Matrix* Kp;       /* Position gain (proportional) */
    PBC_Matrix* Kd;       /* Velocity gain (derivative/damping) */
    PBC_Matrix* Ki;       /* Integral gain (optional, for PI-type) */
    double*     q_desired;  /* Desired position (n×1) */
    double*     qdot_desired; /* Desired velocity (n×1) */
    int         n_joints;   /* Number of DOF */
    bool        gravity_compensation; /* Enable g(q) compensation */
    bool        use_nonlinear_damping; /* Use Ψ(q̇) instead of Kd q̇ */
} PBC_EnergyShapingParams;

/* --- Nonlinear damping function type --- */
typedef void (*PBC_DampingFunc)(const double* qdot, int n, double* tau_damp);

/* --- Energy shaping result --- */
typedef struct {
    double* tau;         /* Computed control torque/force */
    double  kinetic_energy_cl;   /* Closed-loop kinetic energy */
    double  potential_energy_cl; /* Closed-loop potential energy */
    double  total_energy_cl;     /* Closed-loop total energy */
    double  damping_power;       /* Dissipated power = q̇ᵀ Kd q̇ */
    bool    is_stable;           /* Theoretical stability check */
} PBC_EnergyShapingResult;

/* ==============================================================
 * Basic Energy Shaping API
 * ============================================================== */

/**
 * Initialize energy shaping controller parameters.
 * Kp and Kd are n×n symmetric positive definite matrices (copied).
 *
 * @param n_joints  Number of degrees of freedom
 * @param Kp        n×n position gain matrix
 * @param Kd        n×n damping gain matrix
 * @param q_desired Desired joint positions (n elements, copied)
 */
PBC_EnergyShapingParams* pbc_es_params_create(int n_joints,
    const PBC_Matrix* Kp, const PBC_Matrix* Kd, const double* q_desired);

/** Free energy shaping parameters. */
void pbc_es_params_free(PBC_EnergyShapingParams* params);

/**
 * Compute the energy-shaping + damping-injection control law:
 *   τ = g(q) - Kp (q - q_d) - Kd q̇
 *
 * @param params    Controller parameters (with gain matrices)
 * @param q         Current position (n elements)
 * @param qdot      Current velocity (n elements)
 * @param g         Gravity vector g(q) (n elements, can be NULL for no comp)
 * @param tau       Output control torque (n elements, allocated by caller)
 */
void pbc_es_compute_control(const PBC_EnergyShapingParams* params,
    const double* q, const double* qdot, const double* g, double* tau);

/**
 * Compute control with nonlinear damping:
 *   τ = g(q) - Kp (q - q_d) - Ψ(q̇)
 *
 * @param damping   Nonlinear damping function Ψ
 */
void pbc_es_compute_control_nonlinear(const PBC_EnergyShapingParams* params,
    const double* q, const double* qdot, const double* g,
    PBC_DampingFunc damping, double* tau);

/* ==============================================================
 * Gravity Compensation
 * ============================================================== */

/**
 * Compute the gravity vector for a 2-link planar robot.
 *   g₁(q) = (m₁l_{c1} + m₂l₁) g cos(q₁) + m₂ l_{c2} g cos(q₁ + q₂)
 *   g₂(q) = m₂ l_{c2} g cos(q₁ + q₂)
 *
 * @param m1,m2   Link masses (kg)
 * @param l1,l2   Link lengths (m)
 * @param lc1,lc2 Center-of-mass distances (m)
 * @param g_acc   Gravitational acceleration (m/s²)
 * @param q       Joint angles [q1, q2] (rad)
 * @param g_out   Gravity torque [g1, g2] (allocated by caller)
 */
void pbc_gravity_2link_planar(double m1, double m2, double l1, double l2,
    double lc1, double lc2, double g_acc, const double* q, double* g_out);

/**
 * Compute gravity vector for n-link serial chain via iterative Newton-Euler
 * (simplified: each link modeled as point mass at COM).
 *
 * @param n           Number of links
 * @param masses      Link masses (n elements)
 * @param lengths     Link lengths (n elements)
 * @param com_dist    COM distances from joint (n elements)
 * @param g_acc       Gravitational acceleration
 * @param q           Joint angles (n elements)
 * @param g_out       Gravity torques (n elements, allocated by caller)
 */
void pbc_gravity_nlink(int n, const double* masses, const double* lengths,
    const double* com_dist, double g_acc, const double* q, double* g_out);

/* ==============================================================
 * Damping Functions
 * ============================================================== */

/**
 * Linear damping: Ψ(q̇) = Kd · q̇.
 * This is the standard damping injection.
 */
void pbc_damping_linear(const double* qdot, int n, const PBC_Matrix* Kd,
                         double* tau_damp);

/**
 * Hyperbolic tangent (smooth saturation) damping:
 *   Ψ_i(q̇) = Kd_ii · tanh(α_i · q̇_i)
 * Bounded damping torque, useful for actuator saturation.
 *
 * @param Kd     Diagonal damping gains (n×n, only diagonal used)
 * @param alpha  Saturation shape parameters (n elements)
 */
void pbc_damping_tanh(const double* qdot, int n, const PBC_Matrix* Kd,
                       const double* alpha, double* tau_damp);

/**
 * Power-law damping (fractional):
 *   Ψ_i(q̇) = Kd_ii · |q̇_i|^γ · sign(q̇_i), 0 < γ < 1
 * Provides finite-time convergence in some cases.
 */
void pbc_damping_power_law(const double* qdot, int n, const PBC_Matrix* Kd,
                            double gamma, double* tau_damp);

/**
 * Custom damping function wrapper for the DampingFunc type.
 * @param ctx  Pointer to damping context struct (cast in callback)
 */
typedef struct {
    int type;        /* 0=linear, 1=tanh, 2=power_law, 3=custom */
    PBC_Matrix* Kd;  /* Damping matrix */
    double* params;  /* Type-specific parameters */
} PBC_DampingContext;

void pbc_damping_custom(const double* qdot, int n, void* ctx,
                         double* tau_damp);

/* ==============================================================
 * Closed-Loop Energy Analysis
 * ============================================================== */

/**
 * Compute closed-loop energy at a given state.
 *   H_d(q,q̇) = ½ q̇ᵀ M(q) q̇ + ½ (q - q_d)ᵀ Kp (q - q_d)
 *
 * @param M       Inertia matrix M(q) (n×n)
 * @param Kp      Position gain (n×n)
 * @param q       Current position
 * @param qdot    Current velocity
 * @param qd      Desired position
 * @param n       DOF
 * @return Total closed-loop energy
 */
double pbc_closed_loop_energy(const PBC_Matrix* M, const PBC_Matrix* Kp,
    const double* q, const double* qdot, const double* qd, int n);

/**
 * Compute time derivative of closed-loop energy along trajectories.
 *   dH_d/dt = -q̇ᵀ Kd q̇ ≤ 0  (under ideal PBC)
 *
 * This verifies Lyapunov stability of the closed-loop system.
 *
 * @param qdot   Current velocity
 * @param Kd     Damping matrix
 * @param n      DOF
 * @return Energy derivative (should be ≤ 0 for stability)
 */
double pbc_energy_derivative(const double* qdot, const PBC_Matrix* Kd, int n);

/**
 * Check the LaSalle invariance principle for the closed-loop system.
 * If dH_d/dt ≤ 0 and the largest invariant set in {x | dH_d/dt = 0}
 * is the equilibrium, the system is asymptotically stable.
 *
 * @param qdot   Current velocity
 * @param Kd     Damping matrix (must be positive definite)
 * @param n      DOF
 * @return true if dH_d/dt = 0 implies q̇ = 0 (LaSalle condition satisfied)
 */
bool pbc_lasalle_condition(const double* qdot, const PBC_Matrix* Kd, int n,
                            double tol);

/**
 * Simulate closed-loop response under energy-shaping PBC.
 * Uses RK4 integration of the composite system dynamics.
 *
 * @param M       Inertia matrix function M(q) — passed as constant
 * @param C       Coriolis matrix function C(q,q̇)
 * @param g       Gravity vector function g(q)
 * @param params  Energy shaping parameters
 * @param q0      Initial position
 * @param qdot0   Initial velocity
 * @param T       Simulation horizon
 * @param dt      Time step
 * @param n_steps Output: number of steps
 * @param t_out   Output: time array (caller frees)
 * @param q_out   Output: position trajectory (caller frees)
 */
void pbc_es_simulate(const PBC_Matrix* M, const PBC_Matrix* C,
    const double* g_vec, const PBC_EnergyShapingParams* params,
    const double* q0, const double* qdot0, double T, double dt,
    int n, int* n_steps, double** t_out, double** q_out, double** qdot_out);

#endif /* PBC_ENERGY_SHAPING_H */

#ifndef PBC_EL_SYSTEMS_H
#define PBC_EL_SYSTEMS_H

#include "pbc_core.h"

/* ==============================================================
 * pbc_el_systems.h — Euler-Lagrange Mechanical Systems
 *
 * Euler-Lagrange (EL) systems are the natural domain of PBC.
 *
 * EL equation of motion:
 *   d/dt (∂L/∂q̇) - ∂L/∂q = τ
 *   where L(q, q̇) = T(q, q̇) - V(q) is the Lagrangian
 *
 * Standard form:
 *   M(q) q̈ + C(q, q̇) q̇ + g(q) = τ
 *
 * Key structural properties (exploited by PBC):
 *   1. M(q) = Mᵀ(q) > 0   (inertia is symmetric positive definite)
 *   2. Ṁ(q) - 2C(q, q̇) is skew-symmetric (passivity of mechanical map τ→q̇)
 *   3. g(q) = ∂V/∂q      (gravity is gradient of potential)
 *
 * Passivity property (fundamental for PBC):
 *   dE/dt = q̇ᵀ τ    where E = ½ q̇ᵀ M q̇ + V(q)
 *   This means the mechanical system is passive with input τ,
 *   output q̇, and storage function = total energy.
 *
 * References:
 *   Goldstein, Poole, Safko (2002) "Classical Mechanics", 3rd Ed.
 *   Spong, Hutchinson, Vidyasagar (2006) "Robot Modeling and Control"
 *   Ortega et al. (1998) "Passivity-Based Control of Euler-Lagrange Systems"
 *   Kelly, Santibañez, Loria (2005) "Control of Robot Manipulators
 *     in Joint Space"
 *   Slotine & Li (1987) "On the Adaptive Control of Robot Manipulators"
 * ============================================================== */

/* --- Euler-Lagrange system model --- */
typedef struct {
    int         n_dof;     /* Degrees of freedom */
    PBC_Matrix* M;         /* Inertia matrix M(q) (n×n) */
    PBC_Matrix* C;         /* Coriolis/centrifugal matrix C(q,q̇) (n×n) */
    double*     g;         /* Gravity vector g(q) (n elements) */
    double*     q;         /* Generalized coordinates */
    double*     qdot;      /* Generalized velocities */
    double*     qddot;     /* Generalized accelerations */
    double      T;         /* Kinetic energy T = ½ q̇ᵀ M q̇ */
    double      V;         /* Potential energy V(q) */
    double      E;         /* Total energy E = T + V */
    char*       name;
} PBC_EL_System;

/* --- Robot manipulator link parameters --- */
typedef struct {
    double mass;       /* Link mass (kg) */
    double length;     /* Link length (m) */
    double com_dist;   /* Center-of-mass distance from joint (m) */
    double inertia;    /* Moment of inertia about COM (kg·m²) */
    double I_rotor;    /* Rotor inertia (kg·m²) for motor */
    double gear_ratio; /* Gear reduction ratio */
} PBC_LinkParams;

/* --- Robot manipulator model --- */
typedef struct {
    int            n_links;
    PBC_LinkParams* links;
    PBC_EL_System  el;         /* Euler-Lagrange model */
    bool           has_gravity;
    double         gravity_acc; /* Gravitational acceleration (m/s², default 9.81) */
} PBC_RobotModel;

/* --- Trajectory data --- */
typedef struct {
    double* t;       /* Time vector */
    double* q;       /* Position trajectory (n_steps × n_dof, row-major) */
    double* qdot;    /* Velocity trajectory */
    double* tau;     /* Control torque trajectory */
    int     n_steps;
    int     n_dof;
} PBC_Trajectory;

/* ==============================================================
 * Euler-Lagrange System API
 * ============================================================== */

/**
 * Create an EL system with given DOF.
 * All matrices allocated to size n×n, vectors to size n.
 */
PBC_EL_System* pbc_el_create(int n_dof, const char* name);

/** Free an EL system. */
void pbc_el_free(PBC_EL_System* el);

/**
 * Update the EL system state: compute energies and acceleration.
 * After calling, el->qddot contains M⁻¹ (τ - C q̇ - g).
 */
void pbc_el_update(PBC_EL_System* el, const double* tau);

/**
 * Verify fundamental EL structural properties:
 * 1. M is symmetric positive definite
 * 2. Ṁ - 2C is skew-symmetric (passivity of τ → q̇ map)
 *
 * @param qdot_dot  q̈ for Ṁ computation (Ṁ = dM/dt via chain rule)
 * @return true if both properties hold within tolerance
 */
bool pbc_el_verify_properties(const PBC_EL_System* el,
                               const double* qdot_dot, double tol);

/**
 * Check passivity of the EL system: dE/dt = q̇ᵀ τ.
 * This is the fundamental passivity property of mechanical systems.
 *
 * @return Energy derivative; system is passive since dE/dt = q̇ᵀ τ
 */
double pbc_el_passivity_check(const PBC_EL_System* el, const double* tau);

/* ==============================================================
 * Robot Model API
 * ============================================================== */

/**
 * Create a robot model from link parameters.
 *
 * @param n_links  Number of links
 * @param links    Array of link parameters (copied)
 */
PBC_RobotModel* pbc_robot_create(int n_links, const PBC_LinkParams* links);

/** Free a robot model. */
void pbc_robot_free(PBC_RobotModel* robot);

/**
 * Compute the inertia matrix M(q) for a serial robot.
 * Uses the composite-rigid-body algorithm (CRBA).
 *
 * @param robot  Robot model
 * @param q      Joint positions (n_links elements)
 * @param M_out  Output inertia matrix (n×n, allocated by caller)
 */
void pbc_robot_inertia_matrix(const PBC_RobotModel* robot,
                               const double* q, PBC_Matrix* M_out);

/**
 * Compute the Coriolis/centripetal matrix C(q,q̇) using Christoffel symbols:
 *   C_{ij} = ½ ∑_k (∂M_{ij}/∂q_k + ∂M_{ik}/∂q_j - ∂M_{kj}/∂q_i) q̇_k
 *
 * @param robot  Robot model
 * @param q      Joint positions
 * @param qdot   Joint velocities
 * @param C_out  Output Coriolis matrix (n×n, allocated by caller)
 */
void pbc_robot_coriolis_matrix(const PBC_RobotModel* robot,
    const double* q, const double* qdot, PBC_Matrix* C_out);

/**
 * Compute the gravity vector g(q).
 *
 * @param robot  Robot model
 * @param q      Joint positions
 * @param g_out  Output gravity vector (n elements, allocated by caller)
 */
void pbc_robot_gravity_vector(const PBC_RobotModel* robot,
                               const double* q, double* g_out);

/**
 * Forward dynamics: given τ, compute q̈ = M⁻¹ (τ - C q̇ - g).
 */
void pbc_robot_forward_dynamics(const PBC_RobotModel* robot,
    const double* q, const double* qdot, const double* tau, double* qddot);

/**
 * Inverse dynamics: given q, q̇, q̈, compute τ = M q̈ + C q̇ + g.
 */
void pbc_robot_inverse_dynamics(const PBC_RobotModel* robot,
    const double* q, const double* qdot, const double* qddot, double* tau);

/* ==============================================================
 * Specific Robot Models
 * ============================================================== */

/**
 * Create a 1-DOF pendulum model.
 *   I q̈ + m g l_c sin(q) = τ
 */
PBC_RobotModel* pbc_robot_1dof_pendulum(double mass, double length,
                                         double com_dist, double inertia);

/**
 * Create a 2-DOF planar elbow robot (standard benchmark).
 * Link 1: base rotation (revolute)
 * Link 2: elbow rotation (revolute, same plane)
 *
 * Inertia matrix:
 *   M₁₁ = I₁ + I₂ + m₂l₁² + 2m₂l₁lc₂ cos(q₂) + m₁lc₁² + m₂lc₂²
 *   M₁₂ = M₂₁ = I₂ + m₂lc₂² + m₂l₁lc₂ cos(q₂)
 *   M₂₂ = I₂ + m₂lc₂²
 */
PBC_RobotModel* pbc_robot_2dof_planar(const PBC_LinkParams* link1,
                                       const PBC_LinkParams* link2);

/**
 * Create a 3-DOF anthropomorphic robot (shoulder, elbow, wrist).
 */
PBC_RobotModel* pbc_robot_3dof_planar(const PBC_LinkParams* link1,
                                       const PBC_LinkParams* link2,
                                       const PBC_LinkParams* link3);

/* ==============================================================
 * Trajectory Simulation
 * ============================================================== */

/**
 * Simulate robot trajectory under a given torque law using RK4.
 *
 * @param robot      Robot model
 * @param q0, qdot0  Initial state
 * @param T, dt      Simulation horizon and time step
 * @param torque_func Callback: tau = f(q, qdot, t, user_data)
 * @param user_data  User data pointer passed to torque_func
 * @return Simulated trajectory (caller must free via pbc_trajectory_free)
 */
typedef void (*PBC_TorqueFunc)(const double* q, const double* qdot,
                                double t, void* user_data, double* tau);

PBC_Trajectory* pbc_robot_simulate(PBC_RobotModel* robot,
    const double* q0, const double* qdot0, double T, double dt,
    PBC_TorqueFunc torque_func, void* user_data);

/** Free a trajectory. */
void pbc_trajectory_free(PBC_Trajectory* traj);

/** Print trajectory statistics. */
void pbc_trajectory_stats(const PBC_Trajectory* traj);

#endif /* PBC_EL_SYSTEMS_H */

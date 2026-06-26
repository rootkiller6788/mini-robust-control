#ifndef PBC_PORT_HAMILTONIAN_H
#define PBC_PORT_HAMILTONIAN_H

#include "pbc_core.h"

/* ==============================================================
 * pbc_port_hamiltonian.h — Port-Hamiltonian Systems Framework
 *
 * Port-Hamiltonian (pH) systems provide a unified framework for
 * modeling and control of multi-physics systems.
 *
 * Standard input-state-output pH form:
 *   dx/dt = [J(x) - R(x)] ∇H(x) + g(x) u
 *        y = gᵀ(x) ∇H(x)
 *
 * where:
 *   J(x) = -Jᵀ(x)   : Interconnection structure (power-conserving)
 *   R(x) = Rᵀ(x) ≥ 0 : Dissipation structure (power-dissipating)
 *   H(x) ≥ 0         : Hamiltonian / total stored energy
 *   g(x)             : Input port matrix
 *   u, y             : Power-conjugate port variables
 *
 * Power balance (passivity proof):
 *   dH/dt = ∇Hᵀ [J - R]∇H + ∇Hᵀ g u
 *         = -∇Hᵀ R ∇H + yᵀ u              (since J is skew-symmetric)
 *         ≤ yᵀ u                            (since R ≥ 0)
 *
 * This proves ALL port-Hamiltonian systems with R ≥ 0 are passive
 * with storage function H(x).
 *
 * Dirac structure generalization:
 *   A Dirac structure D ⊂ F × F* on a vector space F satisfies D = D⊥
 *   (orthogonal with respect to the pairing bilinear form).
 *
 * References:
 *   Maschke & van der Schaft (1992) "Port-Controlled Hamiltonian Systems:
 *     Modelling Origins and System-Theoretic Properties"
 *   van der Schaft (2000) "L2-Gain and Passivity Techniques in Nonlinear
 *     Control", Chapter 4
 *   van der Schaft & Jeltsema (2014) "Port-Hamiltonian Systems Theory:
 *     An Introductory Overview"
 *   Duindam et al. (2009) "Modeling and Control of Complex Physical Systems"
 *   Ortega et al. (2002) "Putting Energy Back in Control"
 * ============================================================== */

/* --- Port-Hamiltonian system descriptor --- */
typedef struct {
    int     n;         /* State dimension */
    int     m;         /* Number of input/output ports */
    PBC_Matrix* J;     /* Interconnection matrix J(x), n×n, skew-symmetric */
    PBC_Matrix* R;     /* Dissipation matrix R(x), n×n, symmetric ≥ 0 */
    PBC_Matrix* g;     /* Input matrix g(x), n×m */
    double*   x;       /* Current state */
    double*   dH;      /* Gradient of Hamiltonian ∇H(x) */
    double    H;       /* Current Hamiltonian value */
    double    H_dot;   /* dH/dt */
    bool      is_conservative; /* R = 0 (no dissipation) */
} PBC_PH_System;

/* --- Dirac structure (finite-dimensional, bond-graph compatible) --- */
typedef struct {
    PBC_Matrix* D;       /* Dirac structure matrix */
    int num_ports;       /* Number of ports */
    int dim_flow;        /* Dimension of flow space */
    int dim_effort;      /* Dimension of effort space */
} PBC_DiracStructure;

/* --- Bond graph element types --- */
typedef enum {
    PBC_BG_RESISTOR     = 0,  /* R-element: dissipation */
    PBC_BG_CAPACITOR    = 1,  /* C-element: potential energy storage */
    PBC_BG_INERTANCE    = 2,  /* I-element: kinetic energy storage */
    PBC_BG_TRANSFORMER  = 3,  /* TF-element: power-conserving transformation */
    PBC_BG_GYRATOR      = 4,  /* GY-element: power-conserving gyrator */
    PBC_BG_SOURCE_EFFORT = 5, /* Se: effort source */
    PBC_BG_SOURCE_FLOW  = 6,  /* Sf: flow source */
    PBC_BG_JUNCTION_0   = 7,  /* 0-junction: common effort */
    PBC_BG_JUNCTION_1   = 8   /* 1-junction: common flow */
} PBC_BG_ElementType;

/* --- Bond graph element (one node) --- */
typedef struct {
    int              id;
    PBC_BG_ElementType type;
    char*            name;
    double           parameter;  /* R, C, I value, or TF/GY ratio */
    int              num_ports;
    double*          effort;     /* Port efforts */
    double*          flow;       /* Port flows */
    double           power;      /* effort · flow */
} PBC_BG_Element;

/* --- Bond graph (entire system) --- */
typedef struct {
    PBC_BG_Element* elements;
    int              num_elements;
    PBC_Matrix*      junction_matrix; /* Junction structure matrix */
} PBC_BondGraph;

/* ==============================================================
 * Port-Hamiltonian System API
 * ============================================================== */

/**
 * Create a port-Hamiltonian system.
 *
 * @param n   State dimension
 * @param m   Number of ports
 * @param J   Interconnection matrix (n×n, must be skew-symmetric)
 * @param R   Dissipation matrix (n×n, must be symmetric ≥ 0)
 * @param g   Input matrix (n×m)
 */
PBC_PH_System* pbc_ph_create(int n, int m, const PBC_Matrix* J,
                              const PBC_Matrix* R, const PBC_Matrix* g);

/** Free a port-Hamiltonian system. */
void pbc_ph_free(PBC_PH_System* ph);

/**
 * Compute state derivative for a pH system:
 *   dx/dt = [J(x) - R(x)] ∇H(x) + g(x) u
 *
 * @param ph     Port-Hamiltonian system
 * @param dH     Gradient ∇H(x) (n elements)
 * @param u      Input vector (m elements)
 * @param xdot   Output derivative (n elements, allocated by caller)
 */
void pbc_ph_dynamics(const PBC_PH_System* ph, const double* dH,
                      const double* u, double* xdot);

/**
 * Compute output: y = gᵀ(x) ∇H(x)
 *
 * @param ph   Port-Hamiltonian system
 * @param dH   Gradient ∇H(x)
 * @param y    Output vector (m elements, allocated by caller)
 */
void pbc_ph_output(const PBC_PH_System* ph, const double* dH, double* y);

/**
 * Compute power balance: dH/dt = -∇Hᵀ R ∇H + yᵀ u.
 * Verifies that if R ≥ 0, then dH/dt ≤ yᵀ u (passivity).
 */
double pbc_ph_power_balance(const PBC_PH_System* ph, const double* dH,
                             const double* u);

/**
 * Check passivity of a port-Hamiltonian system.
 * Always passive if R(x) ≥ 0.
 */
bool pbc_ph_is_passive(const PBC_PH_System* ph);

/**
 * Create the standard port-Hamiltonian representation of a
 * mass-spring-damper system (1 DOF).
 *
 * State: x = [q, p]ᵀ (position, momentum)
 *   H(q, p) = p²/(2m) + ½ k q²
 *   J = [ 0   1 ]
 *       [ -1  0 ]
 *   R = [ 0   0 ]
 *       [ 0   b ]
 *   g = [ 0 ]
 *       [ 1 ]
 */
PBC_PH_System* pbc_ph_mass_spring_damper(double m, double k, double b);

/**
 * Create the pH representation of an RLC circuit.
 * State: x = [q_C, φ_L]ᵀ (capacitor charge, inductor flux)
 *   H(q, φ) = q²/(2C) + φ²/(2L)
 * For series RLC:
 *   J = [ 0   1 ]     R = [ 0   0 ]
 *       [ -1  0 ]         [ 0   R ]
 */
PBC_PH_System* pbc_ph_rlc_circuit(double R, double L, double C, bool is_series);

/**
 * Create the pH representation of a DC motor.
 *
 * State: x = [φ_f, p_ω, q_θ]ᵀ
 *   φ_f = field flux, p_ω = angular momentum, q_θ = angle
 */
PBC_PH_System* pbc_ph_dc_motor(double R_a, double L_a, double J_m,
                                double B_m, double K_t, double K_e);

/* ==============================================================
 * Interconnection of pH Systems
 * ============================================================== */

/**
 * Interconnect two pH systems in negative feedback.
 *
 * If PH1 and PH2 are passive, their interconnection (u₁ = -y₂ + e₁,
 * u₂ = y₁ + e₂) is also passive (Passivity Theorem).
 *
 * @param ph1, ph2  Port-Hamiltonian systems
 * @return Combined pH system (caller must free)
 */
PBC_PH_System* pbc_ph_feedback_interconnect(PBC_PH_System* ph1,
                                             PBC_PH_System* ph2);

/**
 * Interconnect two pH systems in parallel (common input, summed output).
 */
PBC_PH_System* pbc_ph_parallel_interconnect(PBC_PH_System* ph1,
                                              PBC_PH_System* ph2);

/**
 * Interconnect two pH systems in series (cascade).
 */
PBC_PH_System* pbc_ph_series_interconnect(PBC_PH_System* ph1,
                                           PBC_PH_System* ph2);

/* ==============================================================
 * Dirac Structure API
 * ============================================================== */

/**
 * Create a Dirac structure from a given matrix representation.
 * A Dirac structure satisfies D + D⊥ = F × F*.
 */
PBC_DiracStructure* pbc_dirac_create(const PBC_Matrix* D,
                                      int dim_flow, int dim_effort);

/** Free a Dirac structure. */
void pbc_dirac_free(PBC_DiracStructure* ds);

/**
 * Check if a matrix represents a valid Dirac structure.
 * Condition: D = -D* (skew-symmetric with respect to the pairing).
 */
bool pbc_dirac_is_valid(const PBC_DiracStructure* ds, double tol);

/* ==============================================================
 * Bond Graph API
 * ============================================================== */

/**
 * Create a bond graph element.
 */
PBC_BG_Element* pbc_bg_element_create(int id, PBC_BG_ElementType type,
                                       const char* name, double param);

/** Free a bond graph element. */
void pbc_bg_element_free(PBC_BG_Element* e);

/**
 * Compute element constitutive relations.
 *   R: e = R·f
 *   C: e = q/C = (∫f dt)/C
 *   I: f = p/I = (∫e dt)/I (or e = I·df/dt)
 *   TF: e₁ = n·e₂, f₂ = n·f₁
 *   GY: e₁ = r·f₂, e₂ = r·f₁
 */
void pbc_bg_element_constitutive(PBC_BG_Element* e);

/**
 * Extract port-Hamiltonian representation from a bond graph.
 */
PBC_PH_System* pbc_bg_to_ph(const PBC_BondGraph* bg);

/**
 * Print the bond graph structure (for debugging/visualization).
 */
void pbc_bg_print(const PBC_BondGraph* bg);

#endif /* PBC_PORT_HAMILTONIAN_H */

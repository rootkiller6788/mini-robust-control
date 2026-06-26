#ifndef PBC_PASSIVITY_H
#define PBC_PASSIVITY_H

#include "pbc_core.h"

/* ==============================================================
 * pbc_passivity.h — Passivity Analysis and KYP Lemma
 *
 * The Kalman-Yakubovich-Popov (KYP) Lemma is the bridge between
 * frequency-domain and state-space passivity conditions.
 *
 * For an LTI system:
 *   dx/dt = A x + B u
 *        y = C x + D u
 *
 * KYP Lemma (Anderson 1967, Willems 1972):
 *   The system is passive iff there exists P = Pᵀ ≥ 0 such that:
 *     [ AᵀP + PA    PB - Cᵀ ]
 *     [ BᵀP - C   -(D + Dᵀ) ] ≤ 0
 *
 * For strict passivity, the inequality is strict negative definite.
 *
 * Port-Hamiltonian form equivalence (van der Schaft 2000):
 *   Passive ⇔ exists J = -Jᵀ, R = Rᵀ ≥ 0, and H(x) ≥ 0 s.t.:
 *     dx/dt = [J(x) - R(x)] ∇H(x) + g(x)u
 *          y = gᵀ(x) ∇H(x)
 *
 * References:
 *   Anderson (1967) "A System Theory Criterion for Positive Real Matrices"
 *   Willems (1972) "Dissipative Dynamical Systems, Parts I & II"
 *   Hill & Moylan (1976) "The Stability of Nonlinear Dissipative Systems"
 *   Khalil (2002) "Nonlinear Systems", Chapters 6.5-6.7
 *   van der Schaft (2000) "L2-Gain and Passivity Techniques in Nonlinear Control"
 *   Sepulchre, Jankovic, Kokotovic (1997) "Constructive Nonlinear Control"
 * ============================================================== */

/* --- LTI system in state-space form (for KYP checks) --- */
typedef struct {
    PBC_Matrix* A;    /* n×n system matrix */
    PBC_Matrix* B;    /* n×m input matrix */
    PBC_Matrix* C;    /* p×n output matrix */
    PBC_Matrix* D;    /* p×m feedthrough matrix */
    int n;            /* State dimension */
    int m;            /* Input dimension (= output dimension for square) */
} PBC_LTI_System;

/* --- KYP matrix (the 2×2 block matrix in KYP inequality) --- */
typedef struct {
    PBC_Matrix* block11;  /* AᵀP + PA */
    PBC_Matrix* block12;  /* PB - Cᵀ */
    PBC_Matrix* block21;  /* BᵀP - C */
    PBC_Matrix* block22;  /* -(D + Dᵀ) */
    PBC_Matrix* full;     /* Full (n+m)×(n+m) block matrix */
} PBC_KYP_Matrix;

/* --- Passivity indices (input/output feedforward/passivity) --- */
typedef struct {
    double IFP;   /* Input Feedforward Passivity index ν */
    double OFP;   /* Output Feedback Passivity index ρ */
    double ISP;   /* Input Strict Passivity index ε */
    double OSP;   /* Output Strict Passivity index δ */
    double L2gain; /* L2-gain γ */
} PBC_PassivityIndices;

/* --- Numerical passivity test result --- */
typedef struct {
    bool is_passive;
    PBC_PassivityType type;
    PBC_PassivityIndices indices;
    double KYP_min_eigenvalue;  /* Smallest eigenvalue of KYP matrix */
    double max_violation;       /* Max violation of dissipation inequality */
} PBC_PassivityTest;

/* ==============================================================
 * LTI System API
 * ============================================================== */

/** Create an LTI system descriptor. All matrices are copied. */
PBC_LTI_System* pbc_lti_create(int n, int m,
    const PBC_Matrix* A, const PBC_Matrix* B,
    const PBC_Matrix* C, const PBC_Matrix* D);

/** Free an LTI system. */
void pbc_lti_free(PBC_LTI_System* sys);

/** Create a double integrator: G(s) = 1/s². */
PBC_LTI_System* pbc_lti_double_integrator(void);

/** Create a first-order system: G(s) = K/(τs + 1). */
PBC_LTI_System* pbc_lti_first_order(double K, double tau);

/** Create a mass-spring-damper: G(s) = 1/(Ms² + Bs + K). */
PBC_LTI_System* pbc_lti_mass_spring_damper(double M, double B, double K);

/** Create a DC motor transfer function (armature-controlled). */
PBC_LTI_System* pbc_lti_dc_motor(double R, double L, double J,
                                  double b, double Kt, double Ke);

/* ==============================================================
 * KYP Lemma and Passivity Verification
 * ============================================================== */

/**
 * Build the KYP block matrix for a given LTI system and P matrix.
 *   K = [ AᵀP+PA   PB-Cᵀ ]
 *       [ BᵀP-C   -(D+Dᵀ) ]
 *
 * @param sys   LTI system
 * @param P     Candidate storage matrix (n×n, symmetric)
 * @return KYP block matrix structure (caller must pbc_kyp_free)
 */
PBC_KYP_Matrix* pbc_kyp_build(const PBC_LTI_System* sys, const PBC_Matrix* P);

/** Free a KYP matrix structure. */
void pbc_kyp_free(PBC_KYP_Matrix* kyp);

/**
 * Check whether the KYP matrix is negative semi-definite.
 * This verifies passivity via the KYP lemma.
 *
 * @param kyp  KYP matrix structure
 * @param tol  Numerical tolerance
 * @return true if K ≤ 0 (passive), false otherwise
 */
bool pbc_kyp_is_passive(const PBC_KYP_Matrix* kyp, double tol);

/**
 * Find a storage matrix P satisfying the KYP inequality.
 * Uses a simple Lyapunov-equation-inspired approach for stable A.
 *
 * @param sys      LTI system (must be square: m = p)
 * @param P_out    Output storage matrix (caller allocates n×n)
 * @param max_iter Maximum iterations for iterative refinement
 * @param tol      Convergence tolerance
 * @return true if a feasible P was found
 */
bool pbc_find_storage_matrix(const PBC_LTI_System* sys, PBC_Matrix* P_out,
                              int max_iter, double tol);

/* ==============================================================
 * Passivity Indices
 * ============================================================== */

/**
 * Compute passivity indices for an LTI system.
 *
 * IFP (ν): Largest ν such that system + νI is passive
 *          (i.e., subtracting ν from the feedforward term).
 * OFP (ρ): Largest ρ such that system with output feedback ρI is passive.
 *
 * For LTI SISO systems:
 *   IFP = -inf_{ω} Re[G(jω)]
 *   OFP = +inf_{ω} 1/Re[G(jω)] (if Re[G(jω)] > 0)
 */
PBC_PassivityIndices pbc_compute_indices(const PBC_LTI_System* sys);

/**
 * Shift passivity: apply input feedforward to make system passive.
 * G'(s) = G(s) - νI   where ν = IFP(G).
 */
PBC_LTI_System* pbc_shift_passivity_input(const PBC_LTI_System* sys,
                                           double nu);

/**
 * Shift passivity: apply output feedback to make system passive.
 * G'(s) = G(s) / (1 + ρ G(s))   where ρ = OFP(G).
 */
PBC_LTI_System* pbc_shift_passivity_output(const PBC_LTI_System* sys,
                                            double rho);

/* ==============================================================
 * Relative Degree and Zero Dynamics (Passivity Implications)
 * ============================================================== */

/**
 * Compute relative degree of an LTI SISO system.
 * Relative degree = smallest k such that C A^{k-1} B ≠ 0.
 * Systems with relative degree > 1 cannot be input strictly passive.
 *
 * @return Relative degree, or -1 if error
 */
int pbc_relative_degree(const PBC_LTI_System* sys);

/**
 * Check if system is minimum-phase (all zeros in LHP).
 * For passivity: a passive system is always minimum-phase
 * (zeros must be in closed LHP for SISO).
 *
 * @return true if minimum-phase
 */
bool pbc_is_minimum_phase(const PBC_LTI_System* sys);

/**
 * Check if system is positive real.
 * A transfer function G(s) is positive real if:
 *   1. G(s) is real for real s
 *   2. Re[G(s)] ≥ 0 for Re[s] ≥ 0 (analytic in RHP)
 *
 * Positive-real ⇔ passive for LTI SISO systems.
 */
bool pbc_is_positive_real(const PBC_LTI_System* sys);

/**
 * Spectral check for passivity: compute min Re[G(jω)] over frequency grid.
 * Passive if min Re[G(jω)] ≥ 0.
 */
double pbc_min_real_part_freq(const PBC_LTI_System* sys,
                               double w_min, double w_max, int n_points);

/* ==============================================================
 * Nonlinear Passivity Test
 * ============================================================== */

/**
 * Test passivity of a general nonlinear system via numerical integration
 * and the dissipation inequality check along random trajectories.
 *
 * Generates random initial conditions and inputs, simulates,
 * and verifies the dissipation inequality.
 *
 * @param sys       System descriptor
 * @param n_trials  Number of random trials
 * @param T         Simulation horizon per trial
 * @param dt        Time step
 * @return Passivity test result
 */
PBC_PassivityTest pbc_test_passivity_nonlinear(PBC_System* sys,
    int n_trials, double T, double dt);

#endif /* PBC_PASSIVITY_H */

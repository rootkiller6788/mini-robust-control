#ifndef PU_APP_H
#define PU_APP_H

#include "pu_core.h"
#include "pu_kharitonov.h"
#include "pu_stability.h"
#include "pu_polytope.h"
#include "pu_lmi.h"

/* ============================================================================
 * Applications of Parametric Uncertainty Analysis
 *
 * Real-world engineering systems where parametric uncertainty is critical:
 *
 * 1. Aerospace: F-16 pitch control with uncertain aerodynamic coefficients
 *    (Sobel, Shapiro, Andry, 1994; Stevens, Lewis, Johnson, 2016)
 *
 * 2. Automotive: engine idle speed control with uncertain parameters
 *    (Butts, Sivashankar, Sun, 1999; Cook, Powell, 1988)
 *
 * 3. Power Systems: robust stability under uncertain load conditions
 *    (Kundur, 1994; Anderson, Fouad, 2003)
 *
 * 4. Robotics: manipulator control with uncertain link parameters
 *    (Spong, Hutchinson, Vidyasagar, 2006)
 *
 * 5. Chemical Process: robust temperature control with uncertain kinetics
 *    (Luyben, 1990; Seborg, Edgar, Mellichamp, 2010)
 * ============================================================================ */

/* --- DC Motor with Parametric Uncertainty --- */

/**
 * DC motor model with uncertain parameters:
 *
 *   J * domega/dt = K_t * i - B * omega - T_L
 *   L * di/dt = V - K_e * omega - R * i
 *
 * State: x = [omega; i]
 * Input: u = [V; T_L]
 *
 * Uncertain parameters: R (armature resistance), J (inertia), B (damping)
 *
 * @param R_nom       Nominal resistance (Ohm)
 * @param L           Inductance (H)
 * @param Kt_nom      Nominal torque constant (Nm/A)
 * @param Ke_nom      Nominal back-EMF constant (Vs/rad)
 * @param J_nom       Nominal inertia (kg*m^2)
 * @param B_nom       Nominal damping (Nm*s/rad)
 * @param delta_R     Relative uncertainty in R
 * @param delta_J     Relative uncertainty in J
 * @param delta_B     Relative uncertainty in B
 * @return            Uncertain state-space model
 */
PU_UncertainStateSpace pu_dc_motor_model(double R_nom, double L,
                                           double Kt_nom, double Ke_nom,
                                           double J_nom, double B_nom,
                                           double delta_R, double delta_J,
                                           double delta_B);

/**
 * Test robust stability of a DC motor speed control system with PI controller.
 *
 * @param motor       DC motor uncertain model
 * @param Kp          Proportional gain
 * @param Ki          Integral gain
 * @return            Robust stability status
 */
PU_StabilityStatus pu_dc_motor_robust_stability(const PU_UncertainStateSpace *motor,
                                                  double Kp, double Ki);

/* --- Mass-Spring-Damper with Uncertain Parameters --- */

/**
 * Mass-spring-damper system:
 *   m * d2x/dt2 + c * dx/dt + k * x = F
 *
 * State: x = [position; velocity]
 * Uncertain: m in [m_min, m_max], k in [k_min, k_max], c in [c_min, c_max]
 *
 * @param m_nom       Nominal mass (kg)
 * @param c_nom       Nominal damping (N*s/m)
 * @param k_nom       Nominal stiffness (N/m)
 * @param delta_m     Relative uncertainty in m
 * @param delta_c     Relative uncertainty in c
 * @param delta_k     Relative uncertainty in k
 * @return            Uncertain state-space model
 */
PU_UncertainStateSpace pu_msd_model(double m_nom, double c_nom, double k_nom,
                                      double delta_m, double delta_c, double delta_k);

/**
 * Design a robust state-feedback controller for MSD system.
 * Uses pole placement at the nominal system and then verifies robust stability.
 *
 * @param msd         MSD uncertain model
 * @param wn_des      Desired natural frequency
 * @param xi_des      Desired damping ratio
 * @param K_out       Output: feedback gain [k1, k2]
 * @return            Robust stability status
 */
PU_StabilityStatus pu_msd_robust_controller(const PU_UncertainStateSpace *msd,
                                              double wn_des, double xi_des,
                                              double *K_out);

/* --- Flight Control: F-16 Short-Period with Uncertainty --- */

/**
 * F-16 short-period longitudinal dynamics:
 *
 *   [alpha_dot; q_dot] = [Z_alpha/V, 1; M_alpha, M_q] * [alpha; q]
 *                          + [Z_delta/V; M_delta] * delta_e
 *
 * where alpha = angle of attack, q = pitch rate, delta_e = elevator deflection
 *
 * Uncertain parameters: Z_alpha, M_alpha, M_delta (aerodynamic coefficients)
 * vary with Mach number and altitude.
 *
 * @param V           Trim velocity (ft/s)
 * @param Z_alpha     Nominal Z force derivative w.r.t. alpha (1/s)
 * @param M_alpha     Nominal pitching moment derivative (1/s^2)
 * @param M_q         Nominal pitch damping (1/s)
 * @param Z_delta     Nominal Z force derivative w.r.t. elevator
 * @param M_delta     Nominal control effectiveness (1/s^2)
 * @param unc_level   Relative uncertainty level (e.g., 0.2 = 20%)
 * @return            Uncertain state-space model
 */
PU_UncertainStateSpace pu_f16_short_period_model(double V,
                                                   double Z_alpha, double M_alpha,
                                                   double M_q,
                                                   double Z_delta, double M_delta,
                                                   double unc_level);

/**
 * Design a pitch-rate damper for F-16 with robust stability guarantee.
 *
 * @param f16         F-16 uncertain model
 * @param Kq          Pitch-rate feedback gain
 * @param margin      Output: robust stability margin
 * @return            Robust stability status
 */
PU_StabilityStatus pu_f16_pitch_damper(const PU_UncertainStateSpace *f16,
                                         double Kq, double *margin);

/* --- Power System: Single-Machine Infinite-Bus --- */

/**
 * SMIB model with parametric uncertainty:
 *
 *   d(delta)/dt = omega - omega_0
 *   d(omega)/dt = (P_m - P_e - D*(omega-omega_0)) / M
 *   P_e = (E'*V/X) * sin(delta)
 *
 * Uncertain: inertia M, damping D, reactance X
 *
 * The linearized model at operating point delta_0:
 *   dx/dt = [        0           ,         1       ] * [delta; omega]
 *           [-(E'*V/(M*X))*cos(d0),    -D/M       ]
 *
 * @param E_prime     Internal voltage (pu)
 * @param V_inf       Infinite bus voltage (pu)
 * @param X_nom       Nominal reactance (pu)
 * @param M_nom       Nominal inertia (s)
 * @param D_nom       Nominal damping (pu)
 * @param P_m         Mechanical power (pu)
 * @param delta_R     Uncertainty in each parameter
 * @return            Uncertain state-space model
 */
PU_UncertainStateSpace pu_smib_model(double E_prime, double V_inf,
                                       double X_nom, double M_nom,
                                       double D_nom, double P_m,
                                       double delta_R);

/**
 * Analyze robust small-signal stability of SMIB system.
 * Returns the robust stability margin.
 */
double pu_smib_robust_margin(const PU_UncertainStateSpace *smib);

/* --- Generic Application Helpers --- */

/**
 * Compute worst-case gain margin under parametric uncertainty.
 * Gain margin is the factor by which the DC gain can be multiplied
 * before instability. Under uncertainty, we compute the worst case
 * over the uncertainty set.
 *
 * @param uss         Uncertain state-space
 * @param K_nominal   Nominal controller gain
 * @param wc_gain     Output: worst-case gain margin
 * @param wc_phase    Output: worst-case phase margin (degrees)
 */
void pu_worst_case_margins(const PU_UncertainStateSpace *uss,
                             double K_nominal,
                             double *wc_gain, double *wc_phase);

/**
 * Monte Carlo robust stability simulation.
 * Samples N parameter combinations from the uncertainty set,
 * checks stability of each, and reports the fraction that are stable.
 *
 * @param uss         Uncertain state-space
 * @param n_samples   Number of Monte Carlo samples
 * @param stable_frac Output: fraction of samples that are stable
 */
void pu_monte_carlo_stability(const PU_UncertainStateSpace *uss,
                                int n_samples, double *stable_frac);

/**
 * Print system eigenvalues with uncertainty bounds.
 * For each eigenvalue, computes the range over the uncertainty set.
 *
 * @param uss         Uncertain state-space
 * @param n_samples   Number of samples for eigenvalue range estimation
 */
void pu_eigenvalue_uncertainty(const PU_UncertainStateSpace *uss, int n_samples);

#endif /* PU_APP_H */

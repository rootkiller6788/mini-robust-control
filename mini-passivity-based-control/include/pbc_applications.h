#ifndef PBC_APPLICATIONS_H
#define PBC_APPLICATIONS_H

#include "pbc_core.h"
#include "pbc_el_systems.h"

/* ==============================================================
 * pbc_applications.h — Real-World Applications of PBC
 *
 * Passivity-based control has been successfully applied to:
 *
 * 1. ROBOT MANIPULATORS (since Takegaki & Arimoto 1981):
 *    - Industrial robots (KUKA, ABB, Fanuc)
 *    - Space robotics (NASA Robonaut, Canadarm)
 *    - Surgical robots (da Vinci system)
 *
 * 2. POWER ELECTRONICS:
 *    - DC-DC converters (boost, buck, buck-boost)
 *    - Motor drives (AC/DC)
 *    - Active power filters
 *
 * 3. ELECTROMECHANICAL SYSTEMS:
 *    - Magnetic levitation (maglev trains)
 *    - Electric vehicles (Tesla motor control)
 *    - Wind turbine generators
 *
 * 4. AEROSPACE:
 *    - Satellite attitude control
 *    - Quadrotor UAVs
 *    - Spacecraft docking
 *
 * 5. PROCESS CONTROL:
 *    - Chemical reactor temperature control
 *    - pH neutralization
 *    - Bioreactor control
 *
 * References:
 *   Kelly, Santibañez, Loria (2005) "Control of Robot Manipulators"
 *   Ortega et al. (1998) IEEE Control Systems (PBC tutorial)
 *   Escobar et al. (1999) IEEE TAC (PBC for power converters)
 *   Astolfi, Karagiannis, Ortega (2008) "Nonlinear and Adaptive Control"
 * ============================================================== */

/* --- DC-DC Converter types --- */
typedef enum {
    PBC_CONVERTER_BUCK       = 0,  /* Step-down */
    PBC_CONVERTER_BOOST      = 1,  /* Step-up */
    PBC_CONVERTER_BUCK_BOOST = 2,  /* Inverting step-up/down */
    PBC_CONVERTER_CUK        = 3,  /* Cuk converter */
    PBC_CONVERTER_SEPIC      = 4   /* Single-ended primary-inductor */
} PBC_ConverterType;

/* --- DC-DC converter model --- */
typedef struct {
    PBC_ConverterType type;
    double L;          /* Inductance (H) */
    double C;          /* Capacitance (F) */
    double R;          /* Load resistance (Ω) */
    double Vin;        /* Input voltage (V) */
    double Vref;       /* Desired output voltage (V) */
    double iL;         /* Inductor current (A) */
    double vC;         /* Capacitor voltage (V) */
    double duty;       /* Duty cycle (0 to 1) */
    double switching_freq; /* Switching frequency (Hz) */
} PBC_DCDC_Converter;

/* --- Motor control types --- */
typedef enum {
    PBC_MOTOR_DC_BRUSHED  = 0,
    PBC_MOTOR_PMSM        = 1,  /* Permanent Magnet Synchronous */
    PBC_MOTOR_INDUCTION   = 2,
    PBC_MOTOR_BLDC        = 3   /* Brushless DC */
} PBC_MotorType;

/* --- Motor model --- */
typedef struct {
    PBC_MotorType type;
    double R;          /* Armature resistance (Ω) */
    double L;          /* Armature inductance (H) */
    double J;          /* Rotor inertia (kg·m²) */
    double B;          /* Viscous friction (N·m·s/rad) */
    double Kt;         /* Torque constant (N·m/A) */
    double Ke;         /* Back-EMF constant (V·s/rad) */
    double V_supply;   /* Supply voltage (V) */
    double i;          /* Current (A) */
    double omega;      /* Angular velocity (rad/s) */
    double theta;      /* Angular position (rad) */
    double tau_load;   /* Load torque (N·m) */
} PBC_DC_Motor;

/* --- Quadrotor UAV model --- */
typedef struct {
    double mass;       /* Total mass (kg) */
    double Ixx;        /* Roll inertia (kg·m²) */
    double Iyy;        /* Pitch inertia (kg·m²) */
    double Izz;        /* Yaw inertia (kg·m²) */
    double arm_length; /* Distance from CoM to motor (m) */
    double kt;         /* Thrust coefficient (N·s²) */
    double km;         /* Torque coefficient (N·m·s²) */
    double g;          /* Gravity (m/s²) */
    /* State */
    double x, y, z;    /* Position (m) */
    double dx, dy, dz; /* Velocity (m/s) */
    double phi;        /* Roll angle (rad) */
    double theta;      /* Pitch angle (rad) */
    double psi;        /* Yaw angle (rad) */
    double p, q, r;    /* Angular rates (rad/s) */
} PBC_Quadrotor;

/* --- Maglev system model --- */
typedef struct {
    double mass;        /* Levitated mass (kg) */
    double R;           /* Coil resistance (Ω) */
    double L;           /* Coil inductance (H) */
    double K_force;     /* Force constant (N·A²/m²) */
    double gap_nominal; /* Nominal air gap (m) */
    double g;           /* Gravity (m/s²) */
    double gap;         /* Current air gap (m) */
    double gap_dot;     /* Gap velocity (m/s) */
    double i;           /* Coil current (A) */
    double V;           /* Coil voltage (V) */
} PBC_Maglev;

/* --- pH neutralization process --- */
typedef struct {
    double V;           /* Tank volume (L) */
    double F;           /* Flow rate (L/min) */
    double Ca;          /* Acid concentration (mol/L) */
    double Cb;          /* Base concentration (mol/L) */
    double Ka;          /* Acid dissociation constant */
    double Kw;          /* Water dissociation constant */
    double pH_setpoint; /* Desired pH */
    double pH;          /* Current pH */
    double x1;          /* State: [H⁺] - [OH⁻] */
    double x2;          /* State: [A⁻] */
    double u;           /* Control: base flow rate */
} PBC_pH_Process;

/* ==============================================================
 * DC-DC Converter PBC
 * ============================================================== */

/**
 * Initialize a DC-DC converter model.
 */
PBC_DCDC_Converter* pbc_dcdc_create(PBC_ConverterType type,
    double L, double C, double R, double Vin, double Vref);

/** Free converter model. */
void pbc_dcdc_free(PBC_DCDC_Converter* conv);

/**
 * Compute PBC control law for DC-DC converters.
 *
 * For boost converter (standard PBC):
 *   The system in pH form has:
 *     H(x) = ½ L i² + ½ C v²
 *
 * The PBC law shapes the energy to regulate output voltage:
 *   d = 1 - (Vin + Kp*(vC - Vref) + Kd*iL_disturbance) / vC
 *
 * @param conv     Converter model
 * @param Kp       Proportional gain
 * @param Kd       Damping injection gain
 * @return Duty cycle (0 to 1, clamped)
 */
double pbc_dcdc_control(PBC_DCDC_Converter* conv, double Kp, double Kd);

/**
 * Simulate DC-DC converter dynamics (one switching period).
 * Uses averaged model.
 */
void pbc_dcdc_step(PBC_DCDC_Converter* conv, double duty, double dt);

/* ==============================================================
 * Motor PBC Control
 * ============================================================== */

/**
 * Initialize a DC motor model.
 */
PBC_DC_Motor* pbc_motor_create(PBC_MotorType type, double R, double L,
    double J, double B, double Kt, double Ke, double V_supply);

/** Free motor model. */
void pbc_motor_free(PBC_DC_Motor* motor);

/**
 * PBC-based speed control for DC motor.
 *
 * Motor model (2-state):
 *   L di/dt = -R i - Ke ω + V
 *   J dω/dt = Kt i - B ω - τ_L
 *
 * Energy function: H = ½ L i² + ½ J ω²
 *
 * PBC speed control law:
 *   V = R i_ref + Ke ω_ref + Kp (ω - ω_ref) + Kd (i - i_ref)
 *
 * @param motor      Motor model
 * @param omega_ref  Desired speed (rad/s)
 * @param Kp         Speed proportional gain
 * @param Kd         Current damping gain
 * @return Control voltage (V)
 */
double pbc_motor_speed_control(PBC_DC_Motor* motor, double omega_ref,
                                double Kp, double Kd);

/**
 * PBC-based position control for DC motor (with energy shaping).
 *
 * Extended energy function: H_d = ½ L(i-i*)² + ½ J(ω-ω*)² + ½ K_p(θ-θ*)²
 *
 * @param motor      Motor model
 * @param theta_ref  Desired angle (rad)
 * @param Kp         Position gain
 * @param Kd         Velocity damping gain
 * @param Ki         Current gain
 * @return Control voltage (V)
 */
double pbc_motor_position_control(PBC_DC_Motor* motor, double theta_ref,
                                   double Kp, double Kd, double Ki);

/** Simulate motor dynamics one step forward. */
void pbc_motor_step(PBC_DC_Motor* motor, double V, double dt);

/* ==============================================================
 * Quadrotor PBC Control
 * ============================================================== */

/**
 * Initialize a quadrotor model.
 */
PBC_Quadrotor* pbc_quadrotor_create(double mass, double Ixx, double Iyy,
    double Izz, double arm_length, double kt, double km);

/** Free quadrotor model. */
void pbc_quadrotor_free(PBC_Quadrotor* quad);

/**
 * PBC attitude stabilization for quadrotor.
 *
 * The rotational dynamics are passive with kinetic energy as storage.
 * Control torques are computed via energy shaping:
 *   τ = -Kp (attitude_error) - Kd (angular_rate)
 *
 * @param quad      Quadrotor model
 * @param phi_ref   Desired roll (rad)
 * @param theta_ref Desired pitch (rad)
 * @param psi_ref   Desired yaw (rad)
 * @param Kp        Attitude gains [3]
 * @param Kd        Rate damping gains [3]
 * @param tau       Output torques [τ_phi, τ_theta, τ_psi]
 */
void pbc_quadrotor_attitude_control(PBC_Quadrotor* quad,
    double phi_ref, double theta_ref, double psi_ref,
    const double* Kp, const double* Kd, double* tau);

/**
 * Convert control torques to motor thrusts (quadrotor mixing matrix).
 *
 * For X-configuration:
 *   [F1]   [ 1/(4kt)   1/(4kt)   1/(4kt)  -1/(4km) ]
 *   [F2] = [ 1/(4kt)  -1/(4kt)  -1/(4kt)  -1/(4km) ] × [Fz]
 *   [F3]   [ 1/(4kt)  -1/(4kt)   1/(4kt)   1/(4km) ]   [τx]
 *   [F4]   [ 1/(4kt)   1/(4kt)  -1/(4kt)   1/(4km) ]   [τy]
 *                                                        [τz]
 */
void pbc_quadrotor_mixing(const PBC_Quadrotor* quad, double Fz,
    const double* tau, double* motor_thrusts);

/** Simulate quadrotor dynamics one step (Euler integration). */
void pbc_quadrotor_step(PBC_Quadrotor* quad, const double* motor_thrusts,
                         double dt);

/* ==============================================================
 * Maglev System PBC
 * ============================================================== */

/**
 * Initialize a maglev system model (e.g., Transrapid, Hyperloop concepts).
 *
 * Dynamics:
 *   m gap̈ = m g - K_force (i/gap)²
 *   L di/dt = -R i + V
 *
 * The system is inherently unstable (needs active control).
 * PBC provides a physically motivated stabilization approach.
 */
PBC_Maglev* pbc_maglev_create(double mass, double R, double L,
    double K_force, double gap_nominal);

/** Free maglev model. */
void pbc_maglev_free(PBC_Maglev* maglev);

/**
 * PBC control for maglev gap stabilization.
 *
 * Energy shaping with desired energy:
 *   H_d = ½ m gaṗ² + ½ Kp (gap - gap*)² + ½ L (i - i*)²
 *
 * The PBC control law computes voltage V to achieve:
 *   dH_d/dt = -Kd gaṗ² ≤ 0
 */
double pbc_maglev_control(PBC_Maglev* maglev, double Kp, double Kd);

/** Simulate maglev dynamics one step. */
void pbc_maglev_step(PBC_Maglev* maglev, double V, double dt);

/* ==============================================================
 * pH Process Control
 * ============================================================== */

/**
 * Initialize a pH neutralization process model.
 */
PBC_pH_Process* pbc_ph_process_create(double V, double F, double Ca,
    double Cb, double Ka, double pH_setpoint);

/** Free pH process model. */
void pbc_ph_process_free(PBC_pH_Process* proc);

/**
 * PBC-based pH control.
 * Uses the concept of passivity of reaction invariants.
 */
double pbc_ph_process_control(PBC_pH_Process* proc, double Kp, double Ki);

/** Simulate pH process one step. */
void pbc_ph_process_step(PBC_pH_Process* proc, double u, double dt);

/* ==============================================================
 * Utility: Control Performance Metrics
 * ============================================================== */

/**
 * Compute Integral of Absolute Error (IAE): ∫|e(t)| dt.
 */
double pbc_metric_IAE(const double* t, const double* error, int n_pts);

/**
 * Compute Integral of Squared Error (ISE): ∫ e²(t) dt.
 */
double pbc_metric_ISE(const double* t, const double* error, int n_pts);

/**
 * Compute Integral of Time-weighted Absolute Error (ITAE): ∫ t|e(t)| dt.
 */
double pbc_metric_ITAE(const double* t, const double* error, int n_pts);

/**
 * Compute settling time (2% criterion).
 */
double pbc_metric_settling_time(const double* t, const double* y,
    double y_ss, int n_pts, double band_pct);

/**
 * Compute overshoot percentage.
 */
double pbc_metric_overshoot(const double* y, double y_ss, int n_pts);

#endif /* PBC_APPLICATIONS_H */

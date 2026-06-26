#ifndef KH_APPLICATION_H
#define KH_APPLICATION_H
#include "kh_core.h"

/* ==============================================================
 * kh_application.h - Application-Level Robust Control Design
 *
 * Applications of the Kharitonov theorem to real-world
 * control problems with parametric uncertainty.
 *
 * Each application domain provides:
 *   1. Interval polynomial model of the uncertain system
 *   2. Nominal controller design
 *   3. Kharitonov-based robust stability verification
 *   4. Stability margin quantification
 *
 * Application domains:
 *   - DC motor with uncertain winding resistance / inductance
 *   - Quadrotor UAV with uncertain aerodynamic coefficients
 *   - Aircraft pitch control (F-35 inspired) with parameter variation
 *   - PID robust tuning under gain uncertainty
 *
 * References:
 *   Bhattacharyya et al. (1995). Robust Control.
 *   Dorf & Bishop (2011). Modern Control Systems, 12th ed.
 *   Beard & McLain (2012). Small Unmanned Aircraft.
 * ============================================================== */

/* ---- DC Motor Application ---- */
typedef struct {
    double R_nom;            /* nominal armature resistance (Ω) */
    double R_pct;            /* resistance uncertainty (±%) */
    double L_nom;            /* nominal inductance (H) */
    double L_pct;            /* inductance uncertainty (±%) */
    double J_nom;            /* nominal inertia (kg·m²) */
    double J_pct;            /* inertia uncertainty (±%) */
    double Kt;               /* torque constant (N·m/A) */
    double Ke;               /* back EMF constant (V·s/rad) */
    double B;                /* viscous friction (N·m·s/rad) */
} KH_DCMotorParams;

KH_IntervalPoly* kh_dc_motor_characteristic(const KH_DCMotorParams* params,
                                              double Kp, double Ki, double Kd);
double           kh_dc_motor_stability_margin(const KH_DCMotorParams* params,
                                               double Kp, double Ki, double Kd);

/* ---- Quadrotor UAV Application ---- */
typedef struct {
    double mass_nom;         /* nominal mass (kg) */
    double mass_pct;         /* mass uncertainty (±%) */
    double Ixx_nom;          /* nominal inertia about x (kg·m²) */
    double Ixx_pct;
    double Iyy_nom;
    double Iyy_pct;
    double Izz_nom;
    double Izz_pct;
    double l;                /* arm length (m) */
    double kf;               /* thrust coefficient */
    double km;               /* moment coefficient */
} KH_QuadrotorParams;

KH_IntervalPoly* kh_quadrotor_altitude_poly(const KH_QuadrotorParams* params,
                                              double Kp, double Kd);
KH_IntervalPoly* kh_quadrotor_attitude_poly(const KH_QuadrotorParams* params,
                                              double Kp, double Kd);

/* ---- Aircraft Pitch Control (F-35 inspired) ---- */
typedef struct {
    double M_alpha_nom;      /* pitching moment derivative */
    double M_alpha_pct;
    double M_q_nom;          /* pitch damping derivative */
    double M_q_pct;
    double Z_alpha_nom;      /* normal force derivative */
    double Z_alpha_pct;
    double M_delta_nom;      /* elevator effectiveness */
    double M_delta_pct;
} KH_AircraftPitchParams;

KH_IntervalPoly* kh_aircraft_pitch_poly(const KH_AircraftPitchParams* params);

/* ---- Generic PID robust tuning ---- */
typedef struct {
    double plant_a[5];       /* denominator coefficients a0..a4 */
    double plant_b[3];       /* numerator coefficients b0..b2 */
    int    den_degree;
    int    num_degree;
    double coeff_pct[8];     /* parameter uncertainty percentages */
} KH_PIDPlantParams;

KH_IntervalPoly* kh_pid_closed_loop_poly(const KH_PIDPlantParams* params,
                                           double Kp, double Ki, double Kd);

bool kh_pid_robust_tuning(const KH_PIDPlantParams* params,
                           double Kp_min, double Kp_max, int Kp_steps,
                           double Ki_min, double Ki_max, int Ki_steps,
                           double Kd_min, double Kd_max, int Kd_steps,
                           double* best_Kp, double* best_Ki, double* best_Kd,
                           double* best_margin);

/* ---- Robust stability analysis utilities ---- */
double kh_pole_sensitivity(const KH_IntervalPoly* ip, int pole_idx);
void   kh_worst_case_parameter_set(const KH_IntervalPoly* ip,
                                    double* worst_case_params);
bool   kh_is_robustly_stabilizable(const KH_IntervalPoly* plant);

/* ---- Application report generation ---- */
typedef struct {
    const char* application_name;
    KH_IntervalPoly* characteristic_poly;
    KH_StabilityReport* robust_report;
    double gain_margin_min;
    double phase_margin_min;
    double bandwidth;
    bool is_feasible;
} KH_AppReport;

KH_AppReport* kh_app_report_create(const char* name);
void          kh_app_report_free(KH_AppReport* ar);
void          kh_app_report_print(const KH_AppReport* ar);

#endif /* KH_APPLICATION_H */

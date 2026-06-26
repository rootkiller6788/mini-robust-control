#define _POSIX_C_SOURCE 200809L
#include "kh_application.h"
#include "kh_construction.h"
#include "kh_hurwitz.h"
#include "kh_interval.h"
#include "kh_verification.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

/* ==============================================================
 * kh_application.c - Application-Level Robust Control Design
 *
 * Real-world applications of the Kharitonov theorem.
 * Each application constructs an interval characteristic
 * polynomial from physical parameters, then uses the
 * Kharitonov theorem for robust stability verification.
 *
 * References:
 *   Bhattacharyya, Chapellat, Keel (1995). Robust Control.
 *   Dorf & Bishop (2011). Modern Control Systems.
 *   Beard & McLain (2012). Small Unmanned Aircraft.
 * ============================================================== */

/* ==============================================================
 * DC Motor Application (L7 Applications)
 *
 * Knowledge Point: DC motor transfer function with PID control
 *
 * Physical model (armature-controlled):
 *   Electrical:  L*di/dt + R*i + Ke*omega = V
 *   Mechanical:  J*domega/dt + B*omega = Kt*i
 *
 * Transfer function: G(s) = Kt / ((L*s+R)*(J*s+B) + Kt*Ke)
 *
 * Closed-loop polynomial with PID C(s) = Kp + Ki/s + Kd*s:
 *   phi(s)=L*J*s^3+(L*B+R*J+Kd*Kt)*s^2+(R*B+Kt*Ke+Kp*Kt)*s+Ki*Kt
 *
 * Reference: Dorf & Bishop (2011), Chapter 2.7, 10.4.
 * ============================================================== */

KH_IntervalPoly* kh_dc_motor_characteristic(const KH_DCMotorParams* params,
                                              double Kp, double Ki, double Kd) {
    if (!params) return NULL;
    KH_IntervalPoly* ip = kh_interval_poly_create(3);
    if (!ip) return NULL;

    double L_lo = params->L_nom * (1.0 - params->L_pct / 100.0);
    double L_hi = params->L_nom * (1.0 + params->L_pct / 100.0);
    double J_lo = params->J_nom * (1.0 - params->J_pct / 100.0);
    double J_hi = params->J_nom * (1.0 + params->J_pct / 100.0);
    double R_lo = params->R_nom * (1.0 - params->R_pct / 100.0);
    double R_hi = params->R_nom * (1.0 + params->R_pct / 100.0);

    KH_Interval L_iv = kh_interval_make(L_lo, L_hi);
    KH_Interval J_iv = kh_interval_make(J_lo, J_hi);
    KH_Interval a3 = kh_interval_mul(&L_iv, &J_iv);
    kh_interval_poly_set_coeff(ip, 3, a3.lo, a3.hi);

    double B = params->B, Kt = params->Kt;
    double LB_lo = L_lo * B, LB_hi = L_hi * B;
    double RJ_lo = R_lo * J_lo, RJ_hi = R_hi * J_hi;
    double KdKt = Kd * Kt;
    kh_interval_poly_set_coeff(ip, 2, LB_lo + RJ_lo + KdKt, LB_hi + RJ_hi + KdKt);

    double ReB_lo = R_lo * B, ReB_hi = R_hi * B;
    double KtKe = Kt * params->Ke, KpKt = Kp * Kt;
    kh_interval_poly_set_coeff(ip, 1, ReB_lo + KtKe + KpKt, ReB_hi + KtKe + KpKt);

    kh_interval_poly_set_coeff(ip, 0, Ki * Kt, Ki * Kt);
    return ip;
}

double kh_dc_motor_stability_margin(const KH_DCMotorParams* params,
                                     double Kp, double Ki, double Kd) {
    KH_IntervalPoly* ip = kh_dc_motor_characteristic(params, Kp, Ki, Kd);
    if (!ip) return -INFINITY;
    double margin = kh_stability_margin_family(ip);
    kh_interval_poly_free(ip);
    return margin;
}

/* ==============================================================
 * Quadrotor UAV Application (L7 Applications)
 *
 * Knowledge Point: Quadrotor altitude and attitude dynamics
 *
 * Altitude (linearized about hover):
 *   m * ddot(delta_z) + Kd*dot(delta_z) + Kp*delta_z = 0
 *   Characteristic: m*s^2 + Kd*s + Kp = 0
 *
 * Attitude (roll axis):
 *   Ixx * ddot_phi + Kd*l*kf*dot_phi + Kp*l*kf*phi = 0
 *   Characteristic: Ixx*s^2 + Kd*l*kf*s + Kp*l*kf = 0
 *
 * With mass and inertia uncertainty, the Kharitonov theorem
 * provides necessary and sufficient condition for robust stability.
 *
 * Reference: Beard & McLain (2012), Ch. 6-7.
 * ============================================================== */

KH_IntervalPoly* kh_quadrotor_altitude_poly(const KH_QuadrotorParams* params,
                                              double Kp, double Kd) {
    if (!params) return NULL;
    KH_IntervalPoly* ip = kh_interval_poly_create(2);
    if (!ip) return NULL;

    double m_lo = params->mass_nom * (1.0 - params->mass_pct / 100.0);
    double m_hi = params->mass_nom * (1.0 + params->mass_pct / 100.0);

    kh_interval_poly_set_coeff(ip, 2, m_lo, m_hi);
    kh_interval_poly_set_coeff(ip, 1, Kd, Kd);
    kh_interval_poly_set_coeff(ip, 0, Kp, Kp);
    return ip;
}

KH_IntervalPoly* kh_quadrotor_attitude_poly(const KH_QuadrotorParams* params,
                                              double Kp, double Kd) {
    if (!params) return NULL;
    KH_IntervalPoly* ip = kh_interval_poly_create(2);
    if (!ip) return NULL;

    double I_lo = params->Ixx_nom * (1.0 - params->Ixx_pct / 100.0);
    double I_hi = params->Ixx_nom * (1.0 + params->Ixx_pct / 100.0);
    double lkf = params->l * params->kf;

    kh_interval_poly_set_coeff(ip, 2, I_lo, I_hi);
    kh_interval_poly_set_coeff(ip, 1, Kd * lkf, Kd * lkf);
    kh_interval_poly_set_coeff(ip, 0, Kp * lkf, Kp * lkf);
    return ip;
}

/* ==============================================================
 * Aircraft Pitch Control (F-35 Inspired) (L7 Applications)
 *
 * Knowledge Point: Short-period longitudinal dynamics
 *
 * State-space model for short-period approximation:
 *   [dot_alpha] = [Z_alpha/V    1  ] [alpha] + [0        ] delta_e
 *   [dot_q    ]   [M_alpha      M_q] [q    ]   [M_delta_e]
 *
 * Closed-loop with state feedback u = -K*[alpha; q]:
 *   phi(s) = det(s*I - (A - B*K))
 *
 * For F-35 class aircraft:
 *   M_alpha ~ -4.0, M_q ~ -1.5, Z_alpha ~ -150, M_delta ~ -30
 *
 * Reference: Stevens & Lewis (2003), Aircraft Control and
 *   Simulation, 2nd ed., Ch. 2-4.
 * ============================================================== */

KH_IntervalPoly* kh_aircraft_pitch_poly(const KH_AircraftPitchParams* params) {
    if (!params) return NULL;

    double K_alpha = 2.0, K_q = 1.0, V = 200.0;
    KH_IntervalPoly* ip = kh_interval_poly_create(2);
    if (!ip) return NULL;

    double M_alpha_lo = params->M_alpha_nom * (1.0 - params->M_alpha_pct / 100.0);
    double M_alpha_hi = params->M_alpha_nom * (1.0 + params->M_alpha_pct / 100.0);
    double M_q_lo = params->M_q_nom * (1.0 - params->M_q_pct / 100.0);
    double M_q_hi = params->M_q_nom * (1.0 + params->M_q_pct / 100.0);
    double Z_alpha_lo = params->Z_alpha_nom * (1.0 - params->Z_alpha_pct / 100.0);
    double Z_alpha_hi = params->Z_alpha_nom * (1.0 + params->Z_alpha_pct / 100.0);
    double M_delta_lo = params->M_delta_nom * (1.0 - params->M_delta_pct / 100.0);
    double M_delta_hi = params->M_delta_nom * (1.0 + params->M_delta_pct / 100.0);

    /* Leading coefficient */
    kh_interval_poly_set_coeff(ip, 2, 1.0, 1.0);

    /* a_1 = -Z_alpha/V - M_q + M_delta*K_q
     * Compute interval bounds analytically */
    double ZV_lo = -Z_alpha_hi / V, ZV_hi = -Z_alpha_lo / V;
    double Mq_term_lo = -M_q_hi, Mq_term_hi = -M_q_lo;
    double MdKq_lo = M_delta_lo * K_q, MdKq_hi = M_delta_hi * K_q;
    kh_interval_poly_set_coeff(ip, 1,
        ZV_lo + Mq_term_lo + MdKq_lo,
        ZV_hi + Mq_term_hi + MdKq_hi);

    /* a_0 = Z_alpha*M_q/V - M_alpha - Z_alpha*M_delta*K_q/V + M_delta*K_alpha
     * Enumerate 2^4 = 16 corner combinations for tight bounds */
    double a0_lo = INFINITY, a0_hi = -INFINITY;
    double Zv[2] = {Z_alpha_lo, Z_alpha_hi};
    double Mqv[2] = {M_q_lo, M_q_hi};
    double Malv[2] = {M_alpha_lo, M_alpha_hi};
    double Mdv[2] = {M_delta_lo, M_delta_hi};

    for (int z = 0; z < 2; z++)
        for (int mq = 0; mq < 2; mq++)
            for (int ma = 0; ma < 2; ma++)
                for (int md = 0; md < 2; md++) {
                    double term1 = Zv[z] * Mqv[mq] / V;
                    double term2 = -Malv[ma];
                    double term3 = -Zv[z] * Mdv[md] * K_q / V;
                    double term4 = Mdv[md] * K_alpha;
                    double val = term1 + term2 + term3 + term4;
                    if (val < a0_lo) a0_lo = val;
                    if (val > a0_hi) a0_hi = val;
                }
    kh_interval_poly_set_coeff(ip, 0, a0_lo, a0_hi);
    return ip;
}

/* ==============================================================
 * Generic PID Robust Tuning (L7 Applications)
 *
 * Knowledge Point: Closed-loop polynomial under PID control
 * phi(s) = s*D(s) + (Kd*s^2 + Kp*s + Ki)*N(s)
 * With uncertain plant coefficients, phi(s) is an interval
 * polynomial. Grid search over PID gains maximizes worst-case
 * stability margin verified by the Kharitonov theorem.
 *
 * Reference: Bhattacharyya et al. (1995), Ch. 8-9.
 * ============================================================== */

KH_IntervalPoly* kh_pid_closed_loop_poly(const KH_PIDPlantParams* params,
                                           double Kp, double Ki, double Kd) {
    if (!params) return NULL;

    int cl_degree = params->den_degree + 1;
    int nmax = params->num_degree + 2;
    if (nmax > cl_degree) cl_degree = nmax;

    KH_IntervalPoly* ip = kh_interval_poly_create(cl_degree);
    if (!ip) return NULL;
    for (int i = 0; i <= cl_degree; i++)
        kh_interval_poly_set_coeff(ip, i, 0.0, 0.0);

    /* s*D(s) contribution */
    for (int i = 0; i <= params->den_degree; i++) {
        double pct = params->coeff_pct[i];
        double nom = params->plant_a[i];
        double lo = nom * (1.0 - pct / 100.0);
        double hi = nom * (1.0 + pct / 100.0);
        int idx = i + 1;
        if (idx <= cl_degree)
            kh_interval_poly_set_coeff(ip, idx,
                ip->coeffs[idx].lo + lo, ip->coeffs[idx].hi + hi);
    }

    /* (Kd*s^2 + Kp*s + Ki)*N(s) contribution */
    for (int j = 0; j <= params->num_degree; j++) {
        int pidx = params->den_degree + 1 + j;
        if (pidx >= 8) pidx = 7;
        double pct = params->coeff_pct[pidx];
        double nom = params->plant_b[j];
        double lo = nom * (1.0 - pct / 100.0);
        double hi = nom * (1.0 + pct / 100.0);

        if (j <= cl_degree)
            kh_interval_poly_set_coeff(ip, j,
                ip->coeffs[j].lo + Ki * lo, ip->coeffs[j].hi + Ki * hi);
        if (j + 1 <= cl_degree)
            kh_interval_poly_set_coeff(ip, j + 1,
                ip->coeffs[j + 1].lo + Kp * lo, ip->coeffs[j + 1].hi + Kp * hi);
        if (j + 2 <= cl_degree)
            kh_interval_poly_set_coeff(ip, j + 2,
                ip->coeffs[j + 2].lo + Kd * lo, ip->coeffs[j + 2].hi + Kd * hi);
    }
    return ip;
}

bool kh_pid_robust_tuning(const KH_PIDPlantParams* params,
                           double Kp_min, double Kp_max, int Kp_steps,
                           double Ki_min, double Ki_max, int Ki_steps,
                           double Kd_min, double Kd_max, int Kd_steps,
                           double* best_Kp, double* best_Ki, double* best_Kd,
                           double* best_margin) {
    if (!params || !best_Kp || !best_Ki || !best_Kd || !best_margin)
        return false;
    if (Kp_steps < 1 || Ki_steps < 1 || Kd_steps < 1) return false;

    double dKp = (Kp_max - Kp_min) / (double)(Kp_steps - 1);
    double dKi = (Ki_max - Ki_min) / (double)(Ki_steps - 1);
    double dKd = (Kd_max - Kd_min) / (double)(Kd_steps - 1);

    *best_margin = -INFINITY;
    *best_Kp = Kp_min; *best_Ki = Ki_min; *best_Kd = Kd_min;
    bool found = false;

    for (int ip_k = 0; ip_k < Kp_steps; ip_k++) {
        double Kp = Kp_min + dKp * ip_k;
        for (int ii_k = 0; ii_k < Ki_steps; ii_k++) {
            double Ki = Ki_min + dKi * ii_k;
            for (int id_k = 0; id_k < Kd_steps; id_k++) {
                double Kd = Kd_min + dKd * id_k;
                KH_IntervalPoly* ip = kh_pid_closed_loop_poly(params, Kp, Ki, Kd);
                if (!ip) continue;
                KH_StabilityReport* report = kh_verify_interval_stability(ip);
                bool stable = (report && report->overall == KH_STABLE);
                double margin = kh_stability_margin_family(ip);
                if (stable && margin > *best_margin) {
                    *best_margin = margin;
                    *best_Kp = Kp; *best_Ki = Ki; *best_Kd = Kd;
                    found = true;
                }
                kh_report_free(report);
                kh_interval_poly_free(ip);
            }
        }
    }
    return found;
}

/* ==============================================================
 * Robust Stability Analysis Utilities (L8 Advanced Topics)
 *
 * Knowledge Point: Pole sensitivity via Kharitonov margin analysis
 *
 * The sensitivity of the stability margin to per-coefficient
 * perturbation is estimated by individually perturbing each
 * coefficient interval and measuring the change in the worst-case
 * Kharitonov stability margin.
 *
 * Knowledge Point: Worst-case parameter identification
 *
 * The Kharitonov theorem guarantees that the extremal behavior
 * (worst stability margin) occurs at one of the four Kharitonov
 * vertices. This function identifies which vertex is critical.
 *
 * Knowledge Point: Robust stabilizability criterion
 *
 * A plant is robustly stabilizable iff its nominal polynomial
 * is Hurwitz AND the Kharitonov stability margin is strictly
 * positive. A positive margin means all four Kharitonov
 * polynomials remain Hurwitz.
 *
 * Reference: Barmish (1994), New Tools for Robustness, Ch. 14.
 *   Vidyasagar (1985), Control System Synthesis.
 * ============================================================== */

double kh_pole_sensitivity(const KH_IntervalPoly* ip, int pole_idx) {
    if (!ip || pole_idx < 0 || pole_idx >= ip->degree) return INFINITY;
    double base_margin = kh_stability_margin_family(ip);
    double max_sens = 0.0;

    for (int i = 0; i <= ip->degree; i++) {
        KH_IntervalPoly* ip_pert = kh_interval_poly_create(ip->degree);
        if (!ip_pert) continue;
        for (int j = 0; j <= ip->degree; j++) {
            if (j == i) {
                double w = ip->coeffs[j].hi - ip->coeffs[j].lo;
                double delta = w * 0.01 + 1e-9;
                kh_interval_poly_set_coeff(ip_pert, j,
                    ip->coeffs[j].lo - delta, ip->coeffs[j].hi + delta);
            } else {
                kh_interval_poly_set_coeff(ip_pert, j,
                    ip->coeffs[j].lo, ip->coeffs[j].hi);
            }
        }
        double pm = kh_stability_margin_family(ip_pert);
        double sens = fabs(base_margin - pm);
        if (sens > max_sens) max_sens = sens;
        kh_interval_poly_free(ip_pert);
    }
    return max_sens;
}

void kh_worst_case_parameter_set(const KH_IntervalPoly* ip,
                                  double* worst_case_params) {
    if (!ip || !worst_case_params) return;
    KH_KharitonovSet* ks = kh_kharitonov_construct(ip);
    if (!ks) return;

    const KH_Polynomial* kp[4] = {&ks->K1, &ks->K2, &ks->K3, &ks->K4};
    double worst_margin = INFINITY;
    int worst_idx = 0;

    for (int i = 0; i < 4; i++) {
        KH_StabilityResult r = kh_polynomial_stability(kp[i], NULL);
        if (r != KH_STABLE) { worst_idx = i; break; }
        KH_RouthTable* rt = kh_routh_create(kp[i]);
        if (!rt) continue;
        double min_col = INFINITY;
        for (int rr = 0; rr < rt->n_rows; rr++)
            if (rt->rows[rr][0] > 0.0 && rt->rows[rr][0] < min_col)
                min_col = rt->rows[rr][0];
        if (min_col < worst_margin) { worst_margin = min_col; worst_idx = i; }
        kh_routh_free(rt);
    }

    for (int i = 0; i <= ip->degree; i++)
        worst_case_params[i] = kp[worst_idx]->coeffs[i];
    kh_kharitonov_free(ks);
}

bool kh_is_robustly_stabilizable(const KH_IntervalPoly* plant) {
    if (!plant) return false;
    KH_Polynomial* nom = kh_poly_create(plant->degree);
    if (!nom) return false;
    for (int i = 0; i <= plant->degree; i++)
        nom->coeffs[i] = plant->nominal[i];
    KH_StabilityResult r = kh_polynomial_stability(nom, NULL);
    kh_poly_free(nom);
    if (r != KH_STABLE) return false;
    double margin = kh_stability_margin_family(plant);
    return margin > 0.0;
}

/* ==============================================================
 * Application Report (L7 Applications)
 * ============================================================== */

KH_AppReport* kh_app_report_create(const char* name) {
    KH_AppReport* ar = (KH_AppReport*)malloc(sizeof(KH_AppReport));
    if (!ar) return NULL;
    ar->application_name = name ? strdup(name) : strdup("Unnamed");
    ar->characteristic_poly = NULL;
    ar->robust_report = NULL;
    ar->gain_margin_min = INFINITY;
    ar->phase_margin_min = INFINITY;
    ar->bandwidth = 0.0;
    ar->is_feasible = false;
    return ar;
}

void kh_app_report_free(KH_AppReport* ar) {
    if (!ar) return;
    free((void*)ar->application_name);
    kh_interval_poly_free(ar->characteristic_poly);
    kh_report_free(ar->robust_report);
    free(ar);
}

void kh_app_report_print(const KH_AppReport* ar) {
    if (!ar) return;
    printf("============================================\n");
    printf("Application Report: %s\n", ar->application_name);
    printf("============================================\n");
    printf("Feasibility:       %s\n", ar->is_feasible ? "FEASIBLE" : "NOT FEASIBLE");
    printf("Gain margin (min): %.4f dB\n", ar->gain_margin_min);
    printf("Phase margin (min): %.4f deg\n", ar->phase_margin_min);
    printf("Bandwidth:         %.4f rad/s\n", ar->bandwidth);
    if (ar->robust_report) {
        printf("\nRobust Stability Report:\n");
        kh_report_print(ar->robust_report);
    }
    if (ar->characteristic_poly) {
        printf("\nCharacteristic Polynomial:\n");
        kh_interval_poly_print(ar->characteristic_poly);
    }
    printf("============================================\n");
}
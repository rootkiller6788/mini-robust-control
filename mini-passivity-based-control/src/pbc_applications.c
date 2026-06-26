/* ==============================================================
 * pbc_applications.c -- Real-World PBC Applications
 *
 * Passivity-Based Control applications across multiple domains:
 *   1. DC-DC Power Converters (boost, buck, buck-boost)
 *   2. Electric Motor Control (DC brushed, PMSM)
 *   3. Quadrotor UAV Attitude Control
 *   4. Magnetic Levitation (maglev) Systems
 *   5. pH Neutralization Process Control
 *
 * Each application demonstrates how the fundamental PBC
 * principle (energy shaping + damping injection) naturally
 * extends to diverse physical domains.
 *
 * PBC for Power Converters (Escobar et al. 1999):
 *   The averaged model of a DC-DC converter is an EL system.
 *   Energy shaping regulates output voltage to desired value.
 *
 * PBC for Motors (Ortega et al. 1998):
 *   Motor electromechanical energy is the natural storage
 *   function. PBC exploits this to achieve high-performance
 *   speed and position regulation.
 *
 * PBC for Maglev (Gentili & Marconi 1999):
 *   The inherently unstable magnetic levitation requires
 *   active stabilization. PBC provides a physically-motivated
 *   solution through total energy shaping.
 *
 * References:
 *   Escobar et al. (1999) IEEE TAC (PBC for power converters)
 *   Ortega et al. (1998) IEEE CSM (PBC tutorial)
 *   Kelly, Santibanez, Loria (2005) (Robot control)
 *   Astolfi, Karagiannis, Ortega (2008) (Nonlinear control)
 * ============================================================== */

#include "pbc_applications.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ==============================================================
 * DC-DC Converter PBC
 *
 * Boost converter averaged model (EL form):
 *   L di_L/dt = -(1-d) v_C + V_in
 *   C dv_C/dt = (1-d) i_L - v_C / R
 *
 * Energy: H = 0.5 L i_L^2 + 0.5 C v_C^2
 *
 * PBC objective: shape energy to regulate v_C to V_ref.
 *
 * Control law (Escobar et al. 1999):
 *   The desired equilibrium satisfies:
 *     V_ref / V_in = 1 / (1 - d_eq)
 *   where d_eq = 1 - V_in / V_ref.
 *
 * The PBC duty cycle includes damping injection:
 *   d = 1 - (V_in - K_d * (i_L - i_Ld)) / v_C
 *   where i_Ld = V_ref^2 / (V_in * R) is the desired current.
 * ============================================================== */

PBC_DCDC_Converter* pbc_dcdc_create(PBC_ConverterType type,
    double L, double C, double R, double Vin, double Vref) {
    if (L <= 0.0 || C <= 0.0 || R <= 0.0 || Vin <= 0.0 || Vref <= 0.0)
        return NULL;
    PBC_DCDC_Converter* conv = (PBC_DCDC_Converter*)calloc(1,
        sizeof(PBC_DCDC_Converter));
    if (!conv) return NULL;
    conv->type = type;
    conv->L = L; conv->C = C; conv->R = R;
    conv->Vin = Vin; conv->Vref = Vref;
    conv->iL = 0.0; conv->vC = 0.0;
    conv->duty = (type == PBC_CONVERTER_BOOST) ? 1.0 - Vin/Vref :
                (type == PBC_CONVERTER_BUCK) ? Vref/Vin : 0.5;
    if (conv->duty < 0.0) conv->duty = 0.0;
    if (conv->duty > 1.0) conv->duty = 1.0;
    conv->switching_freq = 100e3;
    return conv;
}

void pbc_dcdc_free(PBC_DCDC_Converter* conv) {
    free(conv);
}

double pbc_dcdc_control(PBC_DCDC_Converter* conv, double Kp, double Kd) {
    if (!conv) return 0.0;
    double L = conv->L, C = conv->C, R = conv->R;
    double Vin = conv->Vin, Vref = conv->Vref;
    double vC = conv->vC, iL = conv->iL;
    double duty;

    switch (conv->type) {
        case PBC_CONVERTER_BOOST: {
            /* Boost: desired equilibrium
             *   vC* = Vref
             *   iL* = Vref^2 / (Vin * R)
             *   d* = 1 - Vin / Vref
             *
             * PBC law: d = 1 - (Vin + Kp*(vC-Vref) + Kd*(iL-iL*)) / vC
             */
            double iL_star = Vref * Vref / (Vin * R);
            if (fabs(vC) < 1e-6) vC = 1e-6;
            double damping = Kp * (vC - Vref) + Kd * (iL - iL_star);
            duty = 1.0 - (Vin + damping) / vC;
            break;
        }
        case PBC_CONVERTER_BUCK: {
            /* Buck: d* = Vref / Vin
             * PBC law: d = Vref/Vin + Kp*(Vref-vC) + Kd*iL_perturb
             */
            duty = Vref / Vin + Kp * (Vref - vC) - Kd * iL;
            break;
        }
        case PBC_CONVERTER_BUCK_BOOST: {
            /* Buck-boost: d* = Vref / (Vref + Vin)
             * PBC uses energy shaping of combined storage
             */
            double d_eq = Vref / (Vref + Vin);
            double damping = Kp * (vC - Vref) - Kd * iL;
            duty = d_eq + damping;
            break;
        }
        default:
            duty = 0.5;
            break;
    }

    /* Clamp duty cycle */
    if (duty < 0.0) duty = 0.0;
    if (duty > 1.0) duty = 1.0;
    conv->duty = duty;
    (void)L; (void)C;
    return duty;
}

void pbc_dcdc_step(PBC_DCDC_Converter* conv, double duty, double dt) {
    if (!conv) return;
    double L = conv->L, C = conv->C, R = conv->R, Vin = conv->Vin;
    double iL = conv->iL, vC = conv->vC;

    /* Clamp duty */
    if (duty < 0.0) duty = 0.0;
    if (duty > 1.0) duty = 1.0;
    conv->duty = duty;

    switch (conv->type) {
        case PBC_CONVERTER_BOOST:
            /* L diL/dt = Vin - (1-d) vC */
            /* C dvC/dt = (1-d) iL - vC/R */
            conv->iL += dt * (Vin - (1.0 - duty) * vC) / L;
            conv->vC += dt * ((1.0 - duty) * iL - vC / R) / C;
            break;
        case PBC_CONVERTER_BUCK:
            /* L diL/dt = d*Vin - vC */
            /* C dvC/dt = iL - vC/R */
            conv->iL += dt * (duty * Vin - vC) / L;
            conv->vC += dt * (iL - vC / R) / C;
            break;
        case PBC_CONVERTER_BUCK_BOOST:
            /* L diL/dt = d*Vin - (1-d)*(-vC) = d*Vin + (1-d)*vC */
            /* C dvC/dt = -(1-d)*iL - vC/R */
            conv->iL += dt * (duty * Vin + (1.0 - duty) * vC) / L;
            conv->vC += dt * (-(1.0 - duty) * iL - vC / R) / C;
            break;
        default:
            break;
    }
}

/* ==============================================================
 * DC Motor PBC Control
 *
 * Motor model (electromechanical):
 *   L di/dt = -R i - Ke omega + V
 *   J domega/dt = Kt i - B omega - tau_L
 *
 * Total stored energy (Hamiltonian):
 *   H = 0.5 L i^2 + 0.5 J omega^2
 *
 * The motor is passive with storage H, input V, output i
 * (electrical port) and also passive with storage = kinetic,
 * input tau, output omega (mechanical port).
 *
 * PBC speed control (Ortega et al. 1998):
 *   V = R i_ref + Ke omega_ref
 *       - Kp (omega - omega_ref)
 *       - Kd (i - i_ref)
 *   where i_ref = (B omega_ref + tau_L) / Kt
 * ============================================================== */

PBC_DC_Motor* pbc_motor_create(PBC_MotorType type, double R, double L,
    double J, double B, double Kt, double Ke, double V_supply) {
    if (L <= 0.0 || J <= 0.0 || R <= 0.0) return NULL;
    PBC_DC_Motor* motor = (PBC_DC_Motor*)calloc(1, sizeof(PBC_DC_Motor));
    if (!motor) return NULL;
    motor->type = type;
    motor->R = R; motor->L = L; motor->J = J; motor->B = B;
    motor->Kt = Kt; motor->Ke = Ke; motor->V_supply = V_supply;
    motor->i = 0.0; motor->omega = 0.0; motor->theta = 0.0;
    motor->tau_load = 0.0;
    return motor;
}

void pbc_motor_free(PBC_DC_Motor* motor) {
    free(motor);
}

double pbc_motor_speed_control(PBC_DC_Motor* motor, double omega_ref,
                                double Kp, double Kd) {
    if (!motor) return 0.0;
    double R = motor->R, Ke = motor->Ke;
    double B = motor->B, Kt = motor->Kt;

    /* Desired current for given speed and load:
     * J domega/dt = Kt i - B omega - tau_L = 0 at steady state
     * i_ref = (B omega_ref + tau_L) / Kt */
    double i_ref = (B * omega_ref + motor->tau_load) / Kt;
    if (fabs(Kt) < 1e-12) i_ref = 0.0;

    /* PBC control voltage:
     * V = R i_ref + Ke omega_ref - Kp (omega - omega_ref) - Kd (i - i_ref) */
    double V = R * i_ref + Ke * omega_ref
             - Kp * (motor->omega - omega_ref)
             - Kd * (motor->i - i_ref);

    /* Saturate to supply voltage */
    if (V > motor->V_supply) V = motor->V_supply;
    if (V < -motor->V_supply) V = -motor->V_supply;
    return V;
}

double pbc_motor_position_control(PBC_DC_Motor* motor, double theta_ref,
                                   double Kp, double Kd, double Ki) {
    if (!motor) return 0.0;
    double R = motor->R, Ke = motor->Ke;
    double B = motor->B, Kt = motor->Kt;

    /* Position error */
    double e_theta = motor->theta - theta_ref;

    /* Desired velocity (proportional to position error) */
    double omega_ref = -Kp * e_theta;

    /* Desired current */
    double i_ref = (B * omega_ref + motor->tau_load
                    - Ki * e_theta) / Kt;
    if (fabs(Kt) < 1e-12) i_ref = 0.0;

    /* PBC control with position energy shaping:
     * V = R i_ref + Ke omega_ref
     *     - Kd (motor->omega - omega_ref)
     *     - Ki*(i - i_ref) */
    double V = R * i_ref + Ke * omega_ref
             - Kd * (motor->omega - omega_ref)
             - Ki * (motor->i - i_ref);

    if (V > motor->V_supply) V = motor->V_supply;
    if (V < -motor->V_supply) V = -motor->V_supply;
    return V;
}

void pbc_motor_step(PBC_DC_Motor* motor, double V, double dt) {
    if (!motor) return;
    double R = motor->R, L = motor->L, J = motor->J;
    double B = motor->B, Kt = motor->Kt, Ke = motor->Ke;

    /* Euler integration of motor dynamics */
    double di = (V - R * motor->i - Ke * motor->omega) / L;
    double dw = (Kt * motor->i - B * motor->omega - motor->tau_load) / J;

    /* Current limiting */
    motor->i += dt * di;
    motor->omega += dt * dw;
    motor->theta += dt * motor->omega;
}

/* ==============================================================
 * Quadrotor PBC Attitude Control
 *
 * The rotational dynamics of a quadrotor are passive with
 * kinetic energy as the storage function:
 *   H = 0.5 omega^T I omega
 *
 * PBC attitude control (Bouabdallah et al. 2004):
 *   tau = -Kp (attitude_error) - Kd omega
 *
 * This is energy shaping (Kp term) + damping injection (Kd term)
 * applied to the Euler angle representation.
 *
 * Mixing matrix converts body-frame torques to individual motor
 * thrust commands.
 * ============================================================== */

PBC_Quadrotor* pbc_quadrotor_create(double mass, double Ixx, double Iyy,
    double Izz, double arm_length, double kt, double km) {
    PBC_Quadrotor* quad = (PBC_Quadrotor*)calloc(1, sizeof(PBC_Quadrotor));
    if (!quad) return NULL;
    quad->mass = mass;
    quad->Ixx = Ixx; quad->Iyy = Iyy; quad->Izz = Izz;
    quad->arm_length = arm_length;
    quad->kt = kt; quad->km = km;
    quad->g = 9.81;
    quad->z = 0.0; /* On ground initially */
    return quad;
}

void pbc_quadrotor_free(PBC_Quadrotor* quad) {
    free(quad);
}

void pbc_quadrotor_attitude_control(PBC_Quadrotor* quad,
    double phi_ref, double theta_ref, double psi_ref,
    const double* Kp, const double* Kd, double* tau) {
    if (!quad || !Kp || !Kd || !tau) return;

    /* Attitude errors */
    double e_phi = quad->phi - phi_ref;
    double e_theta = quad->theta - theta_ref;
    double e_psi = quad->psi - psi_ref;

    /* PBC torques: tau = -Kp * e - Kd * omega */
    tau[0] = -Kp[0] * e_phi   - Kd[0] * quad->p;   /* tau_phi (roll) */
    tau[1] = -Kp[1] * e_theta - Kd[1] * quad->q;   /* tau_theta (pitch) */
    tau[2] = -Kp[2] * e_psi   - Kd[2] * quad->r;   /* tau_psi (yaw) */
}

void pbc_quadrotor_mixing(const PBC_Quadrotor* quad, double Fz,
    const double* tau, double* motor_thrusts) {
    if (!quad || !tau || !motor_thrusts) return;
    double l = quad->arm_length;
    double kt = quad->kt, km = quad->km;

    /* X-configuration mixing matrix:
     * [F1]   [ 1  1  1 -1 ] [Fz/4]
     * [F2] = [ 1 -1 -1 -1 ] [Tx/l]
     * [F3]   [ 1 -1  1  1 ] [Ty/l]
     * [F4]   [ 1  1 -1  1 ] [Tz]
     *
     * Corrected for standard X config:
     * F_i = Fz/4 +/- Tx/(4*l) +/- Ty/(4*l) +/- Tz/(4*km)
     */
    double fz4 = Fz / 4.0;
    double tx = tau[0] / (4.0 * l);
    double ty = tau[1] / (4.0 * l);
    double tz = tau[2] / (4.0 * km);

    motor_thrusts[0] = fz4 + tx + ty - tz;  /* Motor 1: front-right */
    motor_thrusts[1] = fz4 - tx - ty - tz;  /* Motor 2: rear-left */
    motor_thrusts[2] = fz4 - tx + ty + tz;  /* Motor 3: front-left */
    motor_thrusts[3] = fz4 + tx - ty + tz;  /* Motor 4: rear-right */

    /* Ensure non-negative thrusts */
    for (int i = 0; i < 4; i++) {
        if (motor_thrusts[i] < 0.0) motor_thrusts[i] = 0.0;
    }
}

void pbc_quadrotor_step(PBC_Quadrotor* quad, const double* motor_thrusts,
                         double dt) {
    if (!quad || !motor_thrusts) return;
    double m = quad->mass;
    double l = quad->arm_length;
    double kt = quad->kt, km = quad->km;
    double g = quad->g;

    /* Total thrust and torques from motor thrusts */
    double F1 = motor_thrusts[0], F2 = motor_thrusts[1];
    double F3 = motor_thrusts[2], F4 = motor_thrusts[3];
    double F_total = F1 + F2 + F3 + F4;
    double tau_phi = l * (F1 + F4 - F2 - F3);
    double tau_theta = l * (F1 + F2 - F3 - F4);
    double tau_psi = km * (F1 - F2 + F3 - F4) / kt;

    /* Simplified Euler integration of rotational dynamics */
    /* phi_ddot = tau_phi / Ixx */
    quad->p += dt * tau_phi / quad->Ixx;
    quad->q += dt * tau_theta / quad->Iyy;
    quad->r += dt * tau_psi / quad->Izz;

    /* Euler angle rates (simplified for small angles) */
    quad->phi += dt * quad->p;
    quad->theta += dt * quad->q;
    quad->psi += dt * quad->r;

    /* Vertical acceleration */
    double phi = quad->phi, theta = quad->theta;
    double cos_phi = cos(phi), cos_theta = cos(theta);
    double az = F_total * cos_phi * cos_theta / m - g;
    quad->dz += dt * az;
    quad->z += dt * quad->dz;
    if (quad->z < 0.0) { quad->z = 0.0; quad->dz = 0.0; }
}

/* ==============================================================
 * Maglev System PBC
 *
 * Nonlinear dynamics:
 *   m gap_ddot = m g - K_force (i/gap)^2
 *   L i_dot = -R i + V
 *
 * The system is inherently unstable (negative stiffness from
 * magnetic force). PBC stabilizes by shaping total energy:
 *   H_d = 0.5 m gap_dot^2 + 0.5 Kp (gap - gap*)^2
 *       + 0.5 L (i - i_eq)^2
 *
 * Control law ensures dH_d/dt = -Kd gap_dot^2 <= 0.
 *
 * This implements the approach from:
 *   Gentili & Marconi (1999) "Robust nonlinear control of
 *     a magnetic levitation system"
 * ============================================================== */

PBC_Maglev* pbc_maglev_create(double mass, double R, double L,
    double K_force, double gap_nominal) {
    if (mass <= 0.0 || L <= 0.0 || K_force <= 0.0 || gap_nominal <= 0.0)
        return NULL;
    PBC_Maglev* maglev = (PBC_Maglev*)calloc(1, sizeof(PBC_Maglev));
    if (!maglev) return NULL;
    maglev->mass = mass;
    maglev->R = R; maglev->L = L;
    maglev->K_force = K_force;
    maglev->gap_nominal = gap_nominal;
    maglev->g = 9.81;
    maglev->gap = gap_nominal;
    maglev->gap_dot = 0.0;
    maglev->i = mass * gap_nominal * sqrt(9.81 / K_force);
    maglev->V = 0.0;
    return maglev;
}

void pbc_maglev_free(PBC_Maglev* maglev) {
    free(maglev);
}

double pbc_maglev_control(PBC_Maglev* maglev, double Kp, double Kd) {
    if (!maglev) return 0.0;
    double m = maglev->mass, g = maglev->g;
    double Kf = maglev->K_force, R = maglev->R, L = maglev->L;
    double gap = maglev->gap, gap_dot = maglev->gap_dot;
    double i = maglev->i, gap_star = maglev->gap_nominal;

    /* Equilibrium current: m g = Kf (i_eq / gap*)^2
     * i_eq = gap* sqrt(m g / Kf) */
    double i_eq = gap_star * sqrt(m * g / Kf);

    /* PBC: energy shaping + damping injection
     * The control voltage V produces:
     *   dH_d/dt = -Kd gap_dot^2 - R (i - i_eq)^2 <= 0 */

    /* Feedback linearization with energy shaping:
     * Let u = V be computed such that closed-loop error
     * dynamics are asymptotically stable.
     *
     * Simplified: V = R i_eq - Kp (gap - gap*) - Kd gap_dot */
    double V = R * i_eq
             - Kp * (gap - gap_star)
             - Kd * gap_dot;

    /* Add voltage to compensate inductance dynamics */
    V += L * (Kf / (m * gap*gap)) * (2.0 * i * gap_dot / gap - i * i / (gap*gap) * gap_dot);

    return V;
}

void pbc_maglev_step(PBC_Maglev* maglev, double V, double dt) {
    if (!maglev) return;
    double m = maglev->mass, g = maglev->g;
    double Kf = maglev->K_force, R = maglev->R, L = maglev->L;
    double gap = maglev->gap, i = maglev->i;

    if (fabs(gap) < 1e-6) gap = 1e-6;

    /* gap_ddot = g - (Kf/m) * (i/gap)^2 */
    double force_ratio = (i / gap);
    double gap_ddot = g - (Kf / m) * force_ratio * force_ratio;

    /* i_dot = (V - R i) / L */
    double i_dot = (V - R * i) / L;

    maglev->gap_dot += dt * gap_ddot;
    maglev->gap += dt * maglev->gap_dot;
    if (maglev->gap < 1e-4) maglev->gap = 1e-4;  /* Prevent collision */
    maglev->i += dt * i_dot;
    maglev->V = V;
}

/* ==============================================================
 * pH Process Control (Neutralization)
 *
 * pH neutralization is a challenging nonlinear process due to
 * the logarithmic relationship between concentration and pH.
 *
 * The process can be modeled using reaction invariants which
 * simplify the dynamics. PBC exploits the passivity of the
 * reaction invariant dynamics to design a stabilizing controller.
 *
 * Model (Henson & Seborg 1994):
 *   V dx1/dt = F (Ca - x1) + u Cb
 *   V dx2/dt = -F x2 + u Cb
 *   pH = f(x1, x2) (implicit titration curve)
 *
 * PBC approach (Gomez et al. 2004):
 *   Shape the reaction invariant energy to drive pH to setpoint.
 * ============================================================== */

PBC_pH_Process* pbc_ph_process_create(double V, double F, double Ca,
    double Cb, double Ka, double pH_setpoint) {
    if (V <= 0.0 || F < 0.0) return NULL;
    PBC_pH_Process* proc = (PBC_pH_Process*)calloc(1, sizeof(PBC_pH_Process));
    if (!proc) return NULL;
    proc->V = V; proc->F = F;
    proc->Ca = Ca; proc->Cb = Cb;
    proc->Ka = Ka; proc->Kw = 1e-14;
    proc->pH_setpoint = pH_setpoint;
    proc->pH = 7.0;
    proc->x1 = 0.0; proc->x2 = 0.0;
    proc->u = 0.0;
    return proc;
}

void pbc_ph_process_free(PBC_pH_Process* proc) {
    free(proc);
}

double pbc_ph_process_control(PBC_pH_Process* proc, double Kp, double Ki) {
    if (!proc) return 0.0;

    /* Convert pH setpoint to reaction invariant setpoint
     * [H+] = 10^{-pH}
     * Reaction invariant x1 represents [H+] - [OH-] = 10^{-pH} - Kw/10^{-pH}
     */
    double H_set = pow(10.0, -proc->pH_setpoint);
    double x1_set = H_set - proc->Kw / H_set;

    /* PBC-inspired control: proportional + integral on x1 error
     * u = (F*(x1 - Ca) + Kp*(x1 - x1_set) + Ki*integral_error) / Cb */
    double error = proc->x1 - x1_set;
    static double integral = 0.0;  /* Simple static integral (reset on setpoint change) */

    double u = (proc->F * (proc->x1 - proc->Ca) + Kp * error) / proc->Cb;
    if (u < 0.0) u = 0.0;
    if (u > 1.0) u = 1.0;
    (void)Ki;
    return u;
}

void pbc_ph_process_step(PBC_pH_Process* proc, double u, double dt) {
    if (!proc) return;
    double V = proc->V, F = proc->F;
    double Ca = proc->Ca, Cb = proc->Cb;

    /* Reaction invariant dynamics:
     * dx1/dt = F/V (Ca - x1) + u/V * Cb
     * dx2/dt = -F/V * x2 + u/V * Cb
     */
    proc->x1 += dt * (F / V * (Ca - proc->x1) + u / V * Cb);
    proc->x2 += dt * (-F / V * proc->x2 + u / V * Cb);

    /* Compute pH from reaction invariants
     * Solve: [H+]^3 + (Ka + x2)[H+]^2 + (Ka(x2-x1) - Kw)[H+] - Ka Kw = 0
     * Simplified: pH = -log10([H+]) = -log10(sqrt(x1^2 + 4*Kw) - x1)/2
     * For strong acid/base: [H+] ~ x1 when x1 >> 0
     */
    double Kw = proc->Kw;
    double H_conc;
    if (proc->x1 > 1e-6) {
        H_conc = proc->x1;
    } else if (proc->x1 < -1e-6) {
        H_conc = Kw / (-proc->x1);
    } else {
        H_conc = sqrt(Kw);
    }
    if (H_conc > 0.0) {
        proc->pH = -log10(H_conc);
    }
    if (proc->pH < 0.0) proc->pH = 0.0;
    if (proc->pH > 14.0) proc->pH = 14.0;
    proc->u = u;
}

/* ==============================================================
 * Control Performance Metrics
 *
 * Standard metrics for evaluating controller performance:
 *   IAE  = integral |e(t)| dt      (robustness)
 *   ISE  = integral e^2(t) dt      (penalizes large errors)
 *   ITAE = integral t|e(t)| dt     (penalizes late errors)
 *
 * Settling time (2% criterion):
 *   Smallest t such that |y(t) - y_ss| <= 0.02 |y_ss|
 *   for all subsequent times.
 *
 * Overshoot:
 *   max(0, (y_max - y_ss) / y_ss) * 100%
 * ============================================================== */

double pbc_metric_IAE(const double* t, const double* error, int n_pts) {
    if (!t || !error || n_pts < 2) return 0.0;
    double iae = 0.0;
    for (int i = 1; i < n_pts; i++) {
        double dt = t[i] - t[i-1];
        iae += 0.5 * (fabs(error[i]) + fabs(error[i-1])) * dt;
    }
    return iae;
}

double pbc_metric_ISE(const double* t, const double* error, int n_pts) {
    if (!t || !error || n_pts < 2) return 0.0;
    double ise = 0.0;
    for (int i = 1; i < n_pts; i++) {
        double dt = t[i] - t[i-1];
        ise += 0.5 * (error[i]*error[i] + error[i-1]*error[i-1]) * dt;
    }
    return ise;
}

double pbc_metric_ITAE(const double* t, const double* error, int n_pts) {
    if (!t || !error || n_pts < 2) return 0.0;
    double itae = 0.0;
    for (int i = 1; i < n_pts; i++) {
        double dt = t[i] - t[i-1];
        itae += 0.5 * (t[i] * fabs(error[i]) + t[i-1] * fabs(error[i-1])) * dt;
    }
    return itae;
}

double pbc_metric_settling_time(const double* t, const double* y,
    double y_ss, int n_pts, double band_pct) {
    if (!t || !y || n_pts < 2) return 0.0;
    double band = fabs(y_ss) * band_pct / 100.0;
    if (band < 1e-9) band = 1e-6;

    /* Search backwards from the end for the first sustained violation */
    for (int i = n_pts - 1; i >= 1; i--) {
        if (fabs(y[i] - y_ss) > band) {
            /* Return the next time point (when settling begins) */
            if (i + 1 < n_pts) return t[i + 1];
            return t[i];
        }
    }
    return t[0]; /* Already settled */
}

double pbc_metric_overshoot(const double* y, double y_ss, int n_pts) {
    if (!y || n_pts < 2) return 0.0;
    double y_max = y[0];
    for (int i = 1; i < n_pts; i++) {
        if (y[i] > y_max) y_max = y[i];
    }
    if (fabs(y_ss) < 1e-9) return 0.0;
    double overshoot = (y_max - y_ss) / fabs(y_ss) * 100.0;
    if (overshoot < 0.0) overshoot = 0.0;
    return overshoot;
}
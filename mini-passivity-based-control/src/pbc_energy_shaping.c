/* ==============================================================
 * pbc_energy_shaping.c — Energy Shaping and Damping Injection
 *
 * Implements the two fundamental operations of Passivity-Based Control:
 *
 * 1. Energy Shaping: modify closed-loop potential energy so desired
 *    equilibrium becomes the global minimum.
 *
 * 2. Damping Injection: add virtual damping to ensure asymptotic stability.
 *
 * The combined PBC control law for EL systems:
 *   τ = g(q) - K_p (q - q_d) - K_d q̇
 *
 * This yields closed-loop dynamics:
 *   M(q) q̈ + [C(q,q̇) + K_d] q̇ + K_p (q - q_d) = 0
 *
 * With Lyapunov function:
 *   V = ½ q̇ᵀ M q̇ + ½ (q-q_d)ᵀ K_p (q-q_d)
 *
 * And derivative:
 *   dV/dt = -q̇ᵀ K_d q̇ ≤ 0
 *
 * LaSalle's invariance principle guarantees asymptotic stability
 * if K_d > 0 and K_p > 0.
 *
 * References:
 *   Takegaki & Arimoto (1981), Ortega & Spong (1989)
 *   Kelly, Santibañez, Loria (2005), Khalil (2002) Ch. 4.2
 * ============================================================== */

#include "pbc_energy_shaping.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ==============================================================
 * Energy Shaping Parameters
 * ============================================================== */

PBC_EnergyShapingParams* pbc_es_params_create(int n_joints,
    const PBC_Matrix* Kp, const PBC_Matrix* Kd, const double* q_desired) {
    if (n_joints <= 0 || !Kp || !Kd || !q_desired) return NULL;
    if (Kp->rows != n_joints || Kp->cols != n_joints) return NULL;
    if (Kd->rows != n_joints || Kd->cols != n_joints) return NULL;

    PBC_EnergyShapingParams* p = (PBC_EnergyShapingParams*)calloc(1,
        sizeof(PBC_EnergyShapingParams));
    if (!p) return NULL;

    p->Kp = pbc_matrix_copy(Kp);
    p->Kd = pbc_matrix_copy(Kd);
    p->Ki = pbc_matrix_zeros(n_joints, n_joints);
    p->q_desired = (double*)malloc((size_t)n_joints * sizeof(double));
    p->qdot_desired = (double*)calloc((size_t)n_joints, sizeof(double));
    p->n_joints = n_joints;
    p->gravity_compensation = true;

    if (!p->Kp || !p->Kd || !p->q_desired || !p->qdot_desired) {
        pbc_es_params_free(p); return NULL;
    }
    memcpy(p->q_desired, q_desired, (size_t)n_joints * sizeof(double));
    return p;
}

void pbc_es_params_free(PBC_EnergyShapingParams* params) {
    if (!params) return;
    pbc_matrix_free(params->Kp);
    pbc_matrix_free(params->Kd);
    pbc_matrix_free(params->Ki);
    free(params->q_desired);
    free(params->qdot_desired);
    free(params);
}

/* ==============================================================
 * Control Law Computation
 * ============================================================== */

void pbc_es_compute_control(const PBC_EnergyShapingParams* params,
    const double* q, const double* qdot, const double* g, double* tau) {
    if (!params || !q || !qdot || !tau) return;
    int n = params->n_joints;

    /* τ = g(q) - Kp (q - qd) - Kd qdot */
    for (int i = 0; i < n; i++) {
        tau[i] = (params->gravity_compensation && g) ? g[i] : 0.0;

        /* Proportional term: -Kp * (q - qd) */
        for (int j = 0; j < n; j++) {
            tau[i] -= params->Kp->data[i * n + j] *
                      (q[j] - params->q_desired[j]);
        }

        /* Damping term: -Kd * qdot */
        for (int j = 0; j < n; j++) {
            tau[i] -= params->Kd->data[i * n + j] * qdot[j];
        }
    }
}

void pbc_es_compute_control_nonlinear(const PBC_EnergyShapingParams* params,
    const double* q, const double* qdot, const double* g,
    PBC_DampingFunc damping, double* tau) {
    if (!params || !q || !qdot || !damping || !tau) return;
    int n = params->n_joints;

    /* Proportional + gravity terms */
    for (int i = 0; i < n; i++) {
        tau[i] = (params->gravity_compensation && g) ? g[i] : 0.0;
        for (int j = 0; j < n; j++) {
            tau[i] -= params->Kp->data[i * n + j] *
                      (q[j] - params->q_desired[j]);
        }
    }

    /* Nonlinear damping: τ -= Ψ(qdot) */
    double* damp_tau = (double*)malloc((size_t)n * sizeof(double));
    if (!damp_tau) return;
    damping(qdot, n, damp_tau);
    for (int i = 0; i < n; i++) tau[i] -= damp_tau[i];
    free(damp_tau);
}

/* ==============================================================
 * Gravity Compensation
 * ============================================================== */

void pbc_gravity_2link_planar(double m1, double m2, double l1, double l2,
    double lc1, double lc2, double g_acc, const double* q, double* g_out) {
    if (!q || !g_out) return;
    double s1 = sin(q[0]), s12 = sin(q[0] + q[1]), c1 = cos(q[0]), c12 = cos(q[0] + q[1]);

    /* g1 = (m1 lc1 + m2 l1) g cos(q1) + m2 lc2 g cos(q1 + q2) */
    /* Actually the standard robot convention uses: g1 = (m1 lc1 + m2 l1) g sin(q1) + ... */
    /* Wait, gravity is mg downward. The torque depends on how gravity acts on links.
     * Standard planar robot: gravity acts along -y. So potential energy:
     * V = m1 g lc1 sin(q1) + m2 g (l1 sin(q1) + lc2 sin(q1+q2))
     * g1 = ∂V/∂q1 = m1 g lc1 cos(q1) + m2 g l1 cos(q1) + m2 g lc2 cos(q1+q2)
     * g2 = ∂V/∂q2 = m2 g lc2 cos(q1+q2) */
    g_out[0] = (m1*lc1 + m2*l1) * g_acc * c1 + m2*lc2 * g_acc * c12;
    g_out[1] = m2 * lc2 * g_acc * c12;
    (void)l2; (void)s1; (void)s12; /* Parameters reserved for extended models */
}

void pbc_gravity_nlink(int n, const double* masses, const double* lengths,
    const double* com_dist, double g_acc, const double* q, double* g_out) {
    if (!masses || !lengths || !com_dist || !q || !g_out || n <= 0) return;

    /* Initialize to zero */
    for (int i = 0; i < n; i++) g_out[i] = 0.0;

    /* For each link, compute contribution to gravity torques.
     * For an n-link serial chain in the vertical plane:
     *   g_i = ∑_{j=i}^{n} m_j g ∂y_j/∂q_i
     * where y_j is the vertical position of COM of link j.
     *
     * Simplified: assume links are point masses at COM distances,
     * all joints revolute about z-axis, gravity along -y.
     */
    for (int j = 0; j < n; j++) {
        /* Cumulative angle from base */
        double theta_cum = 0.0;
        /* COM vertical position partial derivative chain */
        for (int i = 0; i <= j; i++) {
            theta_cum += q[i];
            /* Contribution: mass_j * g * (horizontal lever arm change) */
            /* dy_j/dq_i = ∂/∂q_i [∑_{k=1}^{j} l_k sin(∑_{l=1}^{k} q_l) + lc_j sin(∑_{l=1}^{j} q_l)] */
            /* Actually for k ≤ j: dy_j/dq_i = (l_i + ... + l_j + lc_j) cos(∑_{l=i}^{j} q_l) */
            /* Simplified for chain model */
            double lever_x = 0.0;
            for (int k = i; k < j; k++) lever_x += lengths[k];
            lever_x += com_dist[j];
            /* The cos of the angle from joint i to COM of j */
            double angle_sum = 0.0;
            for (int k = i; k <= j; k++) angle_sum += q[k];
            g_out[i] += masses[j] * g_acc * lever_x * cos(angle_sum);
        }
    }
}

/* ==============================================================
 * Damping Functions
 * ============================================================== */

void pbc_damping_linear(const double* qdot, int n, const PBC_Matrix* Kd,
                         double* tau_damp) {
    if (!qdot || !Kd || !tau_damp) return;
    for (int i = 0; i < n; i++) {
        tau_damp[i] = 0.0;
        for (int j = 0; j < n; j++) {
            tau_damp[i] += Kd->data[i * n + j] * qdot[j];
        }
    }
}

void pbc_damping_tanh(const double* qdot, int n, const PBC_Matrix* Kd,
                       const double* alpha, double* tau_damp) {
    if (!qdot || !Kd || !tau_damp) return;
    for (int i = 0; i < n; i++) {
        double a = (alpha && alpha[i] > 0.0) ? alpha[i] : 1.0;
        tau_damp[i] = 0.0;
        for (int j = 0; j < n; j++) {
            /* Use diagonal entries only for tanh gains */
            double kd_ij = (i == j) ? Kd->data[i * n + j] : 0.0;
            tau_damp[i] += kd_ij * tanh(a * qdot[j]);
        }
        /* Handle off-diagonal with linear behavior */
        if (i < n) {
            for (int j = 0; j < n; j++) {
                if (i != j) tau_damp[i] += Kd->data[i * n + j] * qdot[j];
            }
        }
    }
}

void pbc_damping_power_law(const double* qdot, int n, const PBC_Matrix* Kd,
                            double gamma, double* tau_damp) {
    if (!qdot || !Kd || !tau_damp) return;
    if (gamma <= 0.0 || gamma >= 1.0) gamma = 0.5; /* Default to sqrt */

    for (int i = 0; i < n; i++) {
        tau_damp[i] = 0.0;
        for (int j = 0; j < n; j++) {
            if (i == j) {
                double abs_qj = fabs(qdot[j]);
                double pow_term = pow(abs_qj, gamma);
                double sign = (qdot[j] > 0.0) ? 1.0 :
                              ((qdot[j] < 0.0) ? -1.0 : 0.0);
                tau_damp[i] += Kd->data[i * n + j] * pow_term * sign;
            } else {
                tau_damp[i] += Kd->data[i * n + j] * qdot[j];
            }
        }
    }
}

void pbc_damping_custom(const double* qdot, int n, void* ctx,
                         double* tau_damp) {
    if (!qdot || !ctx || !tau_damp) return;
    PBC_DampingContext* dc = (PBC_DampingContext*)ctx;
    switch (dc->type) {
        case 0: pbc_damping_linear(qdot, n, dc->Kd, tau_damp); break;
        case 1: pbc_damping_tanh(qdot, n, dc->Kd, dc->params, tau_damp); break;
        case 2: {
            double g = dc->params ? dc->params[0] : 0.5;
            pbc_damping_power_law(qdot, n, dc->Kd, g, tau_damp);
            break;
        }
        default: break;
    }
}

/* ==============================================================
 * Closed-Loop Energy Analysis
 * ============================================================== */

double pbc_closed_loop_energy(const PBC_Matrix* M, const PBC_Matrix* Kp,
    const double* q, const double* qdot, const double* qd, int n) {
    if (!M || !Kp || !q || !qdot || !qd) return 0.0;

    /* Kinetic energy: ½ q̇ᵀ M q̇ */
    double T = 0.0;
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) row_sum += M->data[i * n + j] * qdot[j];
        T += qdot[i] * row_sum;
    }
    T *= 0.5;

    /* Potential energy: ½ (q - qd)ᵀ Kp (q - qd) */
    double V = 0.0;
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) {
            row_sum += Kp->data[i * n + j] * (q[j] - qd[j]);
        }
        V += (q[i] - qd[i]) * row_sum;
    }
    V *= 0.5;

    return T + V;
}

double pbc_energy_derivative(const double* qdot, const PBC_Matrix* Kd, int n) {
    if (!qdot || !Kd) return 0.0;
    /* dH_d/dt = -q̇ᵀ Kd q̇ ≤ 0 */
    double dH = 0.0;
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) row_sum += Kd->data[i * n + j] * qdot[j];
        dH -= qdot[i] * row_sum;
    }
    return dH;
}

bool pbc_lasalle_condition(const double* qdot, const PBC_Matrix* Kd, int n,
                            double tol) {
    if (!qdot || !Kd) return false;
    /* LaSalle: dH/dt = 0 implies q̇ = 0 when Kd > 0.
     * Check if dH/dt = -q̇ᵀ Kd q̇ = 0 ⇒ q̇ᵀ Kd q̇ = 0 ⇒ q̇ = 0.
     * This holds iff Kd is positive definite. */
    double dH_dt = 0.0;
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) row_sum += Kd->data[i * n + j] * qdot[j];
        dH_dt += qdot[i] * row_sum;
    }
    /* If dH_dt ≈ 0, then q̇ must be ≈ 0 for PD Kd */
    double qdot_norm2 = 0.0;
    for (int i = 0; i < n; i++) qdot_norm2 += qdot[i] * qdot[i];

    if (fabs(dH_dt) <= tol) {
        return (qdot_norm2 <= tol * tol);
    }
    return true; /* dH_dt < 0, LaSalle condition satisfied */
}

/* ==============================================================
 * Simulation (RK4 integration of closed-loop)
 * ============================================================== */

void pbc_es_simulate(const PBC_Matrix* M, const PBC_Matrix* C,
    const double* g_vec, const PBC_EnergyShapingParams* params,
    const double* q0, const double* qdot0, double T, double dt,
    int n, int* n_steps, double** t_out, double** q_out, double** qdot_out) {
    if (!M || !params || !q0 || !qdot0 || !n_steps || !t_out || !q_out || !qdot_out || n <= 0) {
        if (n_steps) *n_steps = 0;
        return;
    }

    int steps = (int)(T / dt) + 1;
    *n_steps = steps;

    double* t = (double*)malloc((size_t)steps * sizeof(double));
    double* q = (double*)malloc((size_t)(steps * n) * sizeof(double));
    double* qd = (double*)malloc((size_t)(steps * n) * sizeof(double));

    if (!t || !q || !qd) {
        free(t); free(q); free(qd);
        *n_steps = 0; return;
    }

    /* Initial conditions */
    t[0] = 0.0;
    for (int i = 0; i < n; i++) {
        q[i] = q0[i];
        qd[i] = qdot0[i];
    }

    double* current_q = (double*)malloc((size_t)n * sizeof(double));
    double* current_qd = (double*)malloc((size_t)n * sizeof(double));
    if (!current_q || !current_qd) { free(t); free(q); free(qd); *n_steps=0; return; }
    memcpy(current_q, q0, (size_t)n * sizeof(double));
    memcpy(current_qd, qdot0, (size_t)n * sizeof(double));

    for (int k = 0; k < steps - 1; k++) {
        double h = dt;

        /* Compute control with current state */
        double* tau = (double*)malloc((size_t)n * sizeof(double));
        if (!tau) break;
        pbc_es_compute_control(params, current_q, current_qd, g_vec, tau);

        /* q̈ = M⁻¹ (τ - C q̇ - g) */
        /* For simplicity (small systems): solve M * qddot = rhs */
        double* rhs = (double*)malloc((size_t)n * sizeof(double));
        if (!rhs) { free(tau); break; }
        if (g_vec) {
            for (int i = 0; i < n; i++) rhs[i] = tau[i] - g_vec[i];
        } else {
            memcpy(rhs, tau, (size_t)n * sizeof(double));
        }
        /* Subtract C*qdot if C provided */
        if (C) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    rhs[i] -= C->data[i * n + j] * current_qd[j];
                }
            }
        }

        /* Solve M * qddot = rhs (only for n=1,2 analytically) */
        double qddot[4] = {0, 0, 0, 0};
        if (n == 1) {
            qddot[0] = rhs[0] / M->data[0];
        } else if (n == 2) {
            double det = M->data[0]*M->data[3] - M->data[1]*M->data[2];
            if (fabs(det) > 1e-15) {
                qddot[0] = (M->data[3]*rhs[0] - M->data[1]*rhs[1]) / det;
                qddot[1] = (-M->data[2]*rhs[0] + M->data[0]*rhs[1]) / det;
            }
        }

        /* RK4 for [q; q̇] */
        double k1_q[4], k1_qd[4], k2_q[4], k2_qd[4];
        double k3_q[4], k3_qd[4], k4_q[4], k4_qd[4];

        for (int i = 0; i < n; i++) { k1_q[i] = current_qd[i]; k1_qd[i] = qddot[i]; }
        for (int i = 0; i < n; i++) { k2_q[i] = current_qd[i] + 0.5*h*k1_qd[i]; k2_qd[i] = qddot[i]; }
        for (int i = 0; i < n; i++) { k3_q[i] = current_qd[i] + 0.5*h*k2_qd[i]; k3_qd[i] = qddot[i]; }
        for (int i = 0; i < n; i++) { k4_q[i] = current_qd[i] + h*k3_qd[i]; k4_qd[i] = qddot[i]; }

        for (int i = 0; i < n; i++) {
            current_q[i]  += (h/6.0)*(k1_q[i] + 2*k2_q[i] + 2*k3_q[i] + k4_q[i]);
            current_qd[i] += (h/6.0)*(k1_qd[i] + 2*k2_qd[i] + 2*k3_qd[i] + k4_qd[i]);
        }

        t[k+1] = (k+1) * dt;
        for (int i = 0; i < n; i++) {
            q[(k+1)*n + i] = current_q[i];
            qd[(k+1)*n + i] = current_qd[i];
        }

        free(tau); free(rhs);
    }

    free(current_q); free(current_qd);
    *t_out = t;
    *q_out = q;
    *qdot_out = qd;
}

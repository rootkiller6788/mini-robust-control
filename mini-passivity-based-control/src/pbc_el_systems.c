/* ==============================================================
 * pbc_el_systems.c -- Euler-Lagrange Mechanical Systems
 *
 * EL systems are the natural domain of PBC. The fundamental
 * passivity property dE/dt = q_dot^T tau makes them ideally
 * suited for energy-based control design.
 *
 * Standard form: M(q) q_ddot + C(q,qdot) q_dot + g(q) = tau
 *
 * Key properties exploited by PBC:
 *   1. M(q) = M^T(q) > 0    (positive definite)
 *   2. M_dot - 2C skew-symmetric (passivity of tau->q_dot map)
 *   3. g(q) = dV/dq          (gradient of potential)
 *
 * References: Goldstein (2002), Spong et al. (2006),
 *   Ortega et al. (1998), Kelly et al. (2005)
 * ============================================================== */

#include "pbc_el_systems.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

PBC_EL_System* pbc_el_create(int n_dof, const char* name) {
    if (n_dof <= 0) return NULL;
    PBC_EL_System* el = (PBC_EL_System*)calloc(1, sizeof(PBC_EL_System));
    if (!el) return NULL;
    el->n_dof = n_dof;
    el->M = pbc_matrix_zeros(n_dof, n_dof);
    el->C = pbc_matrix_zeros(n_dof, n_dof);
    el->g = (double*)calloc((size_t)n_dof, sizeof(double));
    el->q = (double*)calloc((size_t)n_dof, sizeof(double));
    el->qdot = (double*)calloc((size_t)n_dof, sizeof(double));
    el->qddot = (double*)calloc((size_t)n_dof, sizeof(double));
    el->T = 0.0; el->V = 0.0; el->E = 0.0;
    el->name = name ? _strdup(name) : _strdup("unnamed_el");
    if (!el->M || !el->C || !el->g || !el->q || !el->qdot || !el->qddot) {
        pbc_el_free(el); return NULL;
    }
    return el;
}

void pbc_el_free(PBC_EL_System* el) {
    if (!el) return;
    pbc_matrix_free(el->M); pbc_matrix_free(el->C);
    free(el->g); free(el->q); free(el->qdot); free(el->qddot);
    free(el->name); free(el);
}

void pbc_el_update(PBC_EL_System* el, const double* tau) {
    if (!el) return;
    int n = el->n_dof;
    el->T = 0.0;
    for (int i = 0; i < n; i++) {
        double row = 0.0;
        for (int j = 0; j < n; j++) row += el->M->data[i * n + j] * el->qdot[j];
        el->T += 0.5 * el->qdot[i] * row;
    }
    el->E = el->T + el->V;
    double* rhs = (double*)malloc((size_t)n * sizeof(double));
    if (!rhs) return;
    for (int i = 0; i < n; i++) {
        rhs[i] = (tau ? tau[i] : 0.0) - el->g[i];
        for (int j = 0; j < n; j++) rhs[i] -= el->C->data[i * n + j] * el->qdot[j];
    }
    if (n == 1) {
        el->qddot[0] = (fabs(el->M->data[0]) > 1e-15) ? rhs[0] / el->M->data[0] : 0.0;
    } else if (n == 2) {
        double a=el->M->data[0], b=el->M->data[1], c=el->M->data[2], d=el->M->data[3];
        double det = a*d - b*c;
        if (fabs(det) > 1e-15) {
            el->qddot[0] = (d*rhs[0] - b*rhs[1]) / det;
            el->qddot[1] = (-c*rhs[0] + a*rhs[1]) / det;
        }
    } else {
        for (int i = 0; i < n; i++) el->qddot[i] = 0.0;
    }
    free(rhs);
}
bool pbc_el_verify_properties(const PBC_EL_System* el,
                               const double* qdot_dot, double tol) {
    if (!el) return false;
    int n = el->n_dof;
    (void)qdot_dot;
    if (!pbc_matrix_is_symmetric(el->M, tol)) return false;
    if (!pbc_matrix_is_positive_definite(el->M)) return false;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double sum = el->C->data[i * n + j] + el->C->data[j * n + i];
            if (fabs(sum) > tol) return false;
        }
    }
    return true;
}

double pbc_el_passivity_check(const PBC_EL_System* el, const double* tau) {
    if (!el || !tau) return 0.0;
    int n = el->n_dof;
    double dE_dt = 0.0;
    for (int i = 0; i < n; i++) dE_dt += el->qdot[i] * tau[i];
    return dE_dt;
}

PBC_RobotModel* pbc_robot_create(int n_links, const PBC_LinkParams* links) {
    if (n_links <= 0 || !links) return NULL;
    PBC_RobotModel* robot = (PBC_RobotModel*)calloc(1, sizeof(PBC_RobotModel));
    if (!robot) return NULL;
    robot->n_links = n_links;
    robot->links = (PBC_LinkParams*)malloc((size_t)n_links * sizeof(PBC_LinkParams));
    if (!robot->links) { free(robot); return NULL; }
    memcpy(robot->links, links, (size_t)n_links * sizeof(PBC_LinkParams));
    robot->has_gravity = true;
    robot->gravity_acc = 9.81;
    robot->el.n_dof = n_links;
    robot->el.M = pbc_matrix_zeros(n_links, n_links);
    robot->el.C = pbc_matrix_zeros(n_links, n_links);
    robot->el.g = (double*)calloc((size_t)n_links, sizeof(double));
    robot->el.q = (double*)calloc((size_t)n_links, sizeof(double));
    robot->el.qdot = (double*)calloc((size_t)n_links, sizeof(double));
    robot->el.qddot = (double*)calloc((size_t)n_links, sizeof(double));
    robot->el.T = 0.0; robot->el.V = 0.0; robot->el.E = 0.0;
    robot->el.name = _strdup("robot");
    return robot;
}

void pbc_robot_free(PBC_RobotModel* robot) {
    if (!robot) return;
    free(robot->links);
    pbc_matrix_free(robot->el.M); pbc_matrix_free(robot->el.C);
    free(robot->el.g); free(robot->el.q);
    free(robot->el.qdot); free(robot->el.qddot);
    free(robot->el.name); free(robot);
}

void pbc_robot_inertia_matrix(const PBC_RobotModel* robot,
                               const double* q, PBC_Matrix* M_out) {
    if (!robot || !q || !M_out) return;
    int n = robot->n_links;
    if (M_out->rows != n || M_out->cols != n) return;
    for (int i = 0; i < n * n; i++) M_out->data[i] = 0.0;
    double* com_x = (double*)malloc((size_t)n * sizeof(double));
    double* com_y = (double*)malloc((size_t)n * sizeof(double));
    if (!com_x || !com_y) { free(com_x); free(com_y); return; }
    double cum_angle = 0.0, cum_x = 0.0, cum_y = 0.0;
    for (int k = 0; k < n; k++) {
        cum_angle += q[k];
        double joint_x = cum_x, joint_y = cum_y;
        com_x[k] = joint_x + robot->links[k].com_dist * cos(cum_angle);
        com_y[k] = joint_y + robot->links[k].com_dist * sin(cum_angle);
        cum_x += robot->links[k].length * cos(cum_angle);
        cum_y += robot->links[k].length * sin(cum_angle);
    }
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double Mij = 0.0;
            for (int k = j; k < n; k++) {
                double dx = 0.0, dy = 0.0;
                for (int m = i; m <= k; m++) {
                    double a = 0.0;
                    for (int p = i; p <= m; p++) a += q[p];
                    double vec_x, vec_y;
                    if (m < k) {
                        vec_x = robot->links[m].length * cos(a);
                        vec_y = robot->links[m].length * sin(a);
                    } else {
                        vec_x = robot->links[k].com_dist * cos(a);
                        vec_y = robot->links[k].com_dist * sin(a);
                    }
                    dx += -vec_y;
                    dy += vec_x;
                }
                Mij += robot->links[k].mass * (dx * dx + dy * dy);
            }
            M_out->data[i * n + j] = Mij;
            M_out->data[j * n + i] = Mij;
        }
    }
    free(com_x); free(com_y);
}

void pbc_robot_coriolis_matrix(const PBC_RobotModel* robot,
    const double* q, const double* qdot, PBC_Matrix* C_out) {
    if (!robot || !q || !qdot || !C_out) return;
    int n = robot->n_links;
    if (C_out->rows != n || C_out->cols != n) return;
    for (int i = 0; i < n * n; i++) C_out->data[i] = 0.0;
    double eps = 1e-6;
    for (int k = 0; k < n; k++) {
        PBC_Matrix* Mp = pbc_matrix_zeros(n, n);
        PBC_Matrix* Mm = pbc_matrix_zeros(n, n);
        if (!Mp || !Mm) { pbc_matrix_free(Mp); pbc_matrix_free(Mm); return; }
        double* qp = (double*)malloc((size_t)n * sizeof(double));
        if (!qp) { pbc_matrix_free(Mp); pbc_matrix_free(Mm); return; }
        memcpy(qp, q, (size_t)n * sizeof(double));
        qp[k] = q[k] + eps;
        pbc_robot_inertia_matrix(robot, qp, Mp);
        qp[k] = q[k] - eps;
        pbc_robot_inertia_matrix(robot, qp, Mm);
        free(qp);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double dMij = (Mp->data[i*n+j] - Mm->data[i*n+j]) / (2.0*eps);
                C_out->data[i*n+j] += 0.5 * dMij * qdot[k];
            }
        }
        pbc_matrix_free(Mp); pbc_matrix_free(Mm);
    }
}

void pbc_robot_gravity_vector(const PBC_RobotModel* robot,
                               const double* q, double* g_out) {
    if (!robot || !q || !g_out) return;
    int n = robot->n_links;
    double g0 = robot->gravity_acc;
    for (int i = 0; i < n; i++) g_out[i] = 0.0;
    for (int k = 0; k < n; k++) {
        for (int i = 0; i <= k; i++) {
            double lever = 0.0;
            for (int m = i; m < k; m++) lever += robot->links[m].length;
            lever += robot->links[k].com_dist;
            double angle = 0.0;
            for (int m = i; m <= k; m++) angle += q[m];
            g_out[i] += robot->links[k].mass * g0 * lever * cos(angle);
        }
    }
}

void pbc_robot_forward_dynamics(const PBC_RobotModel* robot,
    const double* q, const double* qdot, const double* tau, double* qddot) {
    if (!robot || !q || !qdot || !tau || !qddot) return;
    int n = robot->n_links;
    PBC_Matrix* M = pbc_matrix_zeros(n, n);
    PBC_Matrix* C = pbc_matrix_zeros(n, n);
    double* g = (double*)malloc((size_t)n * sizeof(double));
    if (!M || !C || !g) { pbc_matrix_free(M); pbc_matrix_free(C); free(g); return; }
    pbc_robot_inertia_matrix(robot, q, M);
    pbc_robot_coriolis_matrix(robot, q, qdot, C);
    pbc_robot_gravity_vector(robot, q, g);
    double* rhs = (double*)malloc((size_t)n * sizeof(double));
    if (!rhs) { pbc_matrix_free(M); pbc_matrix_free(C); free(g); return; }
    for (int i = 0; i < n; i++) {
        rhs[i] = tau[i] - g[i];
        for (int j = 0; j < n; j++) rhs[i] -= C->data[i*n+j] * qdot[j];
    }
    if (n == 1) {
        qddot[0] = (fabs(M->data[0]) > 1e-15) ? rhs[0]/M->data[0] : 0.0;
    } else if (n == 2) {
        double det = M->data[0]*M->data[3] - M->data[1]*M->data[2];
        if (fabs(det) > 1e-15) {
            qddot[0] = (M->data[3]*rhs[0] - M->data[1]*rhs[1]) / det;
            qddot[1] = (-M->data[2]*rhs[0] + M->data[0]*rhs[1]) / det;
        }
    } else {
        for (int i = 0; i < n; i++) qddot[i] = 0.0;
    }
    pbc_matrix_free(M); pbc_matrix_free(C); free(g); free(rhs);
}

void pbc_robot_inverse_dynamics(const PBC_RobotModel* robot,
    const double* q, const double* qdot, const double* qddot, double* tau) {
    if (!robot || !q || !qdot || !qddot || !tau) return;
    int n = robot->n_links;
    PBC_Matrix* M = pbc_matrix_zeros(n, n);
    PBC_Matrix* C = pbc_matrix_zeros(n, n);
    double* g = (double*)malloc((size_t)n * sizeof(double));
    if (!M || !C || !g) { pbc_matrix_free(M); pbc_matrix_free(C); free(g); return; }
    pbc_robot_inertia_matrix(robot, q, M);
    pbc_robot_coriolis_matrix(robot, q, qdot, C);
    pbc_robot_gravity_vector(robot, q, g);
    for (int i = 0; i < n; i++) {
        tau[i] = g[i];
        for (int j = 0; j < n; j++) {
            tau[i] += M->data[i*n+j] * qddot[j];
            tau[i] += C->data[i*n+j] * qdot[j];
        }
    }
    pbc_matrix_free(M); pbc_matrix_free(C); free(g);
}

PBC_RobotModel* pbc_robot_1dof_pendulum(double mass, double length,
                                         double com_dist, double inertia) {
    PBC_LinkParams link = {mass, length, com_dist, inertia, 0.0, 1.0};
    return pbc_robot_create(1, &link);
}

PBC_RobotModel* pbc_robot_2dof_planar(const PBC_LinkParams* link1,
                                       const PBC_LinkParams* link2) {
    PBC_LinkParams links[2] = {*link1, *link2};
    return pbc_robot_create(2, links);
}

PBC_RobotModel* pbc_robot_3dof_planar(const PBC_LinkParams* link1,
                                       const PBC_LinkParams* link2,
                                       const PBC_LinkParams* link3) {
    PBC_LinkParams links[3] = {*link1, *link2, *link3};
    return pbc_robot_create(3, links);
}

PBC_Trajectory* pbc_robot_simulate(PBC_RobotModel* robot,
    const double* q0, const double* qdot0, double T, double dt,
    PBC_TorqueFunc torque_func, void* user_data) {
    if (!robot || !q0 || !qdot0 || !torque_func) return NULL;
    if (T <= 0.0 || dt <= 0.0) return NULL;
    int n = robot->n_links;
    int steps = (int)(T / dt) + 1;
    PBC_Trajectory* traj = (PBC_Trajectory*)calloc(1, sizeof(PBC_Trajectory));
    if (!traj) return NULL;
    traj->n_steps = steps; traj->n_dof = n;
    traj->t = (double*)malloc((size_t)steps * sizeof(double));
    traj->q = (double*)calloc((size_t)(steps * n), sizeof(double));
    traj->qdot = (double*)calloc((size_t)(steps * n), sizeof(double));
    traj->tau = (double*)calloc((size_t)(steps * n), sizeof(double));
    if (!traj->t || !traj->q || !traj->qdot || !traj->tau) {
        pbc_trajectory_free(traj); return NULL;
    }
    double* q_c = (double*)malloc((size_t)n * sizeof(double));
    double* qd = (double*)malloc((size_t)n * sizeof(double));
    double* tau_v = (double*)malloc((size_t)n * sizeof(double));
    if (!q_c || !qd || !tau_v) {
        free(q_c); free(qd); free(tau_v);
        pbc_trajectory_free(traj); return NULL;
    }
    memcpy(q_c, q0, (size_t)n * sizeof(double));
    memcpy(qd, qdot0, (size_t)n * sizeof(double));
    for (int k = 0; k < steps; k++) {
        traj->t[k] = k * dt;
        for (int i = 0; i < n; i++) {
            traj->q[k*n+i] = q_c[i];
            traj->qdot[k*n+i] = qd[i];
        }
        torque_func(q_c, qd, traj->t[k], user_data, tau_v);
        for (int i = 0; i < n; i++) traj->tau[k*n+i] = tau_v[i];
        if (k < steps - 1) {
            double h = dt;
            double *k1 = (double*)malloc((size_t)(8*n) * sizeof(double));
            if (!k1) break;
            double *k1q=k1, *k1qd=k1+n, *k2q=k1+2*n, *k2qd=k1+3*n;
            double *k3q=k1+4*n, *k3qd=k1+5*n, *k4q=k1+6*n, *k4qd=k1+7*n;
            double *qt = (double*)malloc((size_t)(2*n) * sizeof(double));
            double *qdt = qt + n;
            if (!qt) { free(k1); break; }
            memcpy(k1q, qd, (size_t)n*sizeof(double));
            pbc_robot_forward_dynamics(robot, q_c, qd, tau_v, k1qd);
            for (int i=0;i<n;i++) {qt[i]=q_c[i]+0.5*h*k1q[i]; qdt[i]=qd[i]+0.5*h*k1qd[i];}
            memcpy(k2q, qdt, (size_t)n*sizeof(double));
            torque_func(qt, qdt, traj->t[k]+0.5*h, user_data, tau_v);
            pbc_robot_forward_dynamics(robot, qt, qdt, tau_v, k2qd);
            for (int i=0;i<n;i++) {qt[i]=q_c[i]+0.5*h*k2q[i]; qdt[i]=qd[i]+0.5*h*k2qd[i];}
            memcpy(k3q, qdt, (size_t)n*sizeof(double));
            torque_func(qt, qdt, traj->t[k]+0.5*h, user_data, tau_v);
            pbc_robot_forward_dynamics(robot, qt, qdt, tau_v, k3qd);
            for (int i=0;i<n;i++) {qt[i]=q_c[i]+h*k3q[i]; qdt[i]=qd[i]+h*k3qd[i];}
            memcpy(k4q, qdt, (size_t)n*sizeof(double));
            torque_func(qt, qdt, traj->t[k]+h, user_data, tau_v);
            pbc_robot_forward_dynamics(robot, qt, qdt, tau_v, k4qd);
            for (int i=0;i<n;i++) {
                q_c[i]+=(h/6.0)*(k1q[i]+2*k2q[i]+2*k3q[i]+k4q[i]);
                qd[i]+=(h/6.0)*(k1qd[i]+2*k2qd[i]+2*k3qd[i]+k4qd[i]);
            }
            free(k1); free(qt);
        }
    }
    free(q_c); free(qd); free(tau_v);
    return traj;
}

void pbc_trajectory_free(PBC_Trajectory* traj) {
    if (!traj) return;
    free(traj->t); free(traj->q); free(traj->qdot); free(traj->tau);
    free(traj);
}

void pbc_trajectory_stats(const PBC_Trajectory* traj) {
    if (!traj) { printf("Trajectory: (null)\n"); return; }
    printf("=== Trajectory: %d steps, %d DOF ===\n", traj->n_steps, traj->n_dof);
    printf("Time: %.3f to %.3f\n", traj->t[0], traj->t[traj->n_steps-1]);
    int n = traj->n_dof;
    for (int j = 0; j < n; j++) {
        double q_min = traj->q[j], q_max = traj->q[j];
        for (int k = 1; k < traj->n_steps; k++) {
            double val = traj->q[k * n + j];
            if (val < q_min) q_min = val;
            if (val > q_max) q_max = val;
        }
        printf("  Joint %d: q [%.3f, %.3f]\n", j, q_min, q_max);
    }
}
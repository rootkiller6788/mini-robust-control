/* ==============================================================
 * pbc_ida.c -- Interconnection and Damping Assignment PBC
 *
 * IDA-PBC systematically stabilizes nonlinear systems by
 * assigning a desired port-Hamiltonian closed-loop structure.
 *
 * Matching Equation: [J-R]grad_H + g u = [J_d-R_d]grad_H_d
 *
 * For fully-actuated mechanical systems (Ortega et al. 2002):
 *   tau = g(q) + [M M_d^{-1} - I] K_d (q - q*)
 *         + [M M_d^{-1} C - C] qdot - M M_d^{-1} K_v qdot
 *
 * References: Ortega et al. (2002), Ortega & Garcia-Canseco (2004)
 * ============================================================== */

#include "pbc_ida.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

PBC_IDA_Structure* pbc_ida_structure_create(int n,
    const PBC_Matrix* J_d, const PBC_Matrix* R_d, const double* x_star) {
    if (n <= 0 || !J_d || !R_d || !x_star) return NULL;
    if (J_d->rows != n || J_d->cols != n) return NULL;
    if (R_d->rows != n || R_d->cols != n) return NULL;
    PBC_IDA_Structure* s = (PBC_IDA_Structure*)calloc(1, sizeof(PBC_IDA_Structure));
    if (!s) return NULL;
    s->n = n;
    s->J_d = pbc_matrix_copy(J_d);
    s->R_d = pbc_matrix_copy(R_d);
    s->x_star = (double*)malloc((size_t)n * sizeof(double));
    if (!s->x_star) { pbc_ida_structure_free(s); return NULL; }
    memcpy(s->x_star, x_star, (size_t)n * sizeof(double));
    return s;
}

void pbc_ida_structure_free(PBC_IDA_Structure* s) {
    if (!s) return;
    pbc_matrix_free(s->J_d);
    pbc_matrix_free(s->R_d);
    free(s->x_star);
    free(s);
}

PBC_IDA_Hamiltonian* pbc_ida_hamiltonian_create(int n_dof,
    const PBC_Matrix* M_d, const PBC_Matrix* K_d, const double* q_star) {
    if (n_dof <= 0 || !M_d || !K_d || !q_star) return NULL;
    if (M_d->rows != n_dof || M_d->cols != n_dof) return NULL;
    if (K_d->rows != n_dof || K_d->cols != n_dof) return NULL;
    if (!pbc_matrix_is_positive_definite(M_d)) return NULL;
    PBC_IDA_Hamiltonian* Hd = (PBC_IDA_Hamiltonian*)calloc(1, sizeof(PBC_IDA_Hamiltonian));
    if (!Hd) return NULL;
    Hd->n_dof = n_dof;
    Hd->M_d = pbc_matrix_copy(M_d);
    Hd->K_d = pbc_matrix_copy(K_d);
    Hd->q_star = (double*)malloc((size_t)n_dof * sizeof(double));
    Hd->p_star = (double*)calloc((size_t)n_dof, sizeof(double));
    Hd->H_min = 0.0;
    if (!Hd->q_star || !Hd->p_star) { pbc_ida_hamiltonian_free(Hd); return NULL; }
    memcpy(Hd->q_star, q_star, (size_t)n_dof * sizeof(double));
    return Hd;
}

void pbc_ida_hamiltonian_free(PBC_IDA_Hamiltonian* Hd) {
    if (!Hd) return;
    pbc_matrix_free(Hd->M_d);
    pbc_matrix_free(Hd->K_d);
    free(Hd->q_star);
    free(Hd->p_star);
    free(Hd);
}

PBC_IDA_Controller* pbc_ida_controller_create(int n,
    const PBC_Matrix* M, const PBC_Matrix* C,
    const PBC_Matrix* J_d, const PBC_Matrix* R_d,
    const PBC_IDA_Hamiltonian* Hd) {
    if (n <= 0 || !M || !J_d || !R_d || !Hd) return NULL;
    if (M->rows != n || M->cols != n) return NULL;
    PBC_IDA_Controller* ctrl = (PBC_IDA_Controller*)calloc(1, sizeof(PBC_IDA_Controller));
    if (!ctrl) return NULL;
    ctrl->n = n;
    ctrl->M = pbc_matrix_copy(M);
    ctrl->C = C ? pbc_matrix_copy(C) : pbc_matrix_zeros(n, n);
    ctrl->g = (double*)calloc((size_t)n, sizeof(double));
    ctrl->structure.n = Hd->n_dof;
    ctrl->structure.J_d = pbc_matrix_copy(J_d);
    ctrl->structure.R_d = pbc_matrix_copy(R_d);
    ctrl->structure.x_star = (double*)malloc((size_t)n * sizeof(double));
    if (ctrl->structure.x_star && Hd->q_star) {
        memcpy(ctrl->structure.x_star, Hd->q_star, (size_t)n * sizeof(double));
    }
    ctrl->hamiltonian.M_d = pbc_matrix_copy(Hd->M_d);
    ctrl->hamiltonian.K_d = pbc_matrix_copy(Hd->K_d);
    ctrl->hamiltonian.n_dof = Hd->n_dof;
    ctrl->hamiltonian.q_star = (double*)malloc((size_t)n * sizeof(double));
    ctrl->hamiltonian.p_star = (double*)calloc((size_t)n, sizeof(double));
    if (ctrl->hamiltonian.q_star && Hd->q_star) {
        memcpy(ctrl->hamiltonian.q_star, Hd->q_star, (size_t)n * sizeof(double));
    }
    if (!ctrl->M || !ctrl->C || !ctrl->g ||
        !ctrl->structure.J_d || !ctrl->structure.R_d ||
        !ctrl->hamiltonian.M_d || !ctrl->hamiltonian.K_d) {
        pbc_ida_controller_free(ctrl); return NULL;
    }
    ctrl->is_configured = true;
    return ctrl;
}

void pbc_ida_controller_free(PBC_IDA_Controller* ctrl) {
    if (!ctrl) return;
    pbc_matrix_free(ctrl->M);
    pbc_matrix_free(ctrl->C);
    free(ctrl->g);
    pbc_matrix_free(ctrl->structure.J_d);
    pbc_matrix_free(ctrl->structure.R_d);
    free(ctrl->structure.x_star);
    pbc_matrix_free(ctrl->hamiltonian.M_d);
    pbc_matrix_free(ctrl->hamiltonian.K_d);
    free(ctrl->hamiltonian.q_star);
    free(ctrl->hamiltonian.p_star);
    free(ctrl);
}


PBC_IDA_Result pbc_ida_compute_control(PBC_IDA_Controller* ctrl,
    const double* q, const double* qdot, const double* g_vec, double* tau) {
    PBC_IDA_Result result;
    memset(&result, 0, sizeof(result));
    if (!ctrl || !q || !qdot || !tau) return result;
    int n = ctrl->n;
    double* Kd = ctrl->hamiltonian.K_d->data;
    double* Rd = ctrl->structure.R_d->data;
    double* qstar = ctrl->hamiltonian.q_star;
    for (int i = 0; i < n; i++) {
        tau[i] = (g_vec) ? g_vec[i] : 0.0;
        for (int j = 0; j < n; j++) {
            tau[i] -= Kd[i * n + j] * (q[j] - qstar[j]);
            tau[i] -= Rd[i * n + j] * qdot[j];
        }
    }
    result.H_d_value = 0.0;
    for (int i = 0; i < n; i++) {
        double row_k = 0.0, row_d = 0.0;
        for (int j = 0; j < n; j++) {
            row_k += Kd[i * n + j] * (q[j] - qstar[j]);
            row_d += Rd[i * n + j] * qdot[j];
        }
        result.H_d_value += 0.5 * (q[i] - qstar[i]) * row_k;
        result.H_d_dot -= qdot[i] * row_d;
    }
    for (int i = 0; i < n; i++) {
        double row = 0.0;
        for (int j = 0; j < n; j++) {
            row += ctrl->M->data[i * n + j] * qdot[j];
        }
        result.H_d_value += 0.5 * qdot[i] * row;
    }
    result.u_pbc = tau;
    result.solved = true;
    return result;
}

double pbc_ida_1dof_solve(double m, double c, double md, double kd,
                           double kv, double q, double qdot,
                           double qstar, double g_val) {
    (void)c; /* natural damping not used in simplified solution */
    if (m <= 0.0 || md <= 0.0) return g_val;
    double eq = q - qstar;
    double restoring = kv * qdot + kd * eq;
    double ratio = m / md;
    return g_val + (ratio - 1.0) * restoring - kv * qdot - kd * eq;
}

void pbc_ida_mechanical_ndof(const PBC_Matrix* M, const PBC_Matrix* C,
    const PBC_Matrix* Md, const PBC_Matrix* Kd, const PBC_Matrix* Kv,
    const double* q, const double* qdot, const double* qstar,
    const double* g_vec, double* tau, int n) {
    if (!M || !Kd || !Kv || !q || !qdot || !qstar || !tau || n <= 0) return;
    double Md_inv[9] = {0};
    int has_inv = 0;
    if (n == 1) {
        if (fabs(Md->data[0]) > 1e-15) { Md_inv[0] = 1.0 / Md->data[0]; has_inv = 1; }
    } else if (n == 2) {
        double det = Md->data[0]*Md->data[3] - Md->data[1]*Md->data[2];
        if (fabs(det) > 1e-15) {
            Md_inv[0] = Md->data[3] / det; Md_inv[1] = -Md->data[1] / det;
            Md_inv[2] = -Md->data[2] / det; Md_inv[3] = Md->data[0] / det;
            has_inv = 1;
        }
    } else if (n == 3) {
        double d0 = Md->data[4]*Md->data[8] - Md->data[5]*Md->data[7];
        double d1 = Md->data[3]*Md->data[8] - Md->data[5]*Md->data[6];
        double d2 = Md->data[3]*Md->data[7] - Md->data[4]*Md->data[6];
        double det = Md->data[0]*d0 - Md->data[1]*d1 + Md->data[2]*d2;
        if (fabs(det) > 1e-15) {
            Md_inv[0] = d0 / det;
            Md_inv[1] = -(Md->data[1]*Md->data[8] - Md->data[2]*Md->data[7]) / det;
            Md_inv[2] = (Md->data[1]*Md->data[5] - Md->data[2]*Md->data[4]) / det;
            Md_inv[3] = -d1 / det;
            Md_inv[4] = (Md->data[0]*Md->data[8] - Md->data[2]*Md->data[6]) / det;
            Md_inv[5] = -(Md->data[0]*Md->data[5] - Md->data[2]*Md->data[3]) / det;
            Md_inv[6] = d2 / det;
            Md_inv[7] = -(Md->data[0]*Md->data[7] - Md->data[1]*Md->data[6]) / det;
            Md_inv[8] = (Md->data[0]*Md->data[4] - Md->data[1]*Md->data[3]) / det;
            has_inv = 1;
        }
    }
    for (int i = 0; i < n; i++) {
        tau[i] = (g_vec) ? g_vec[i] : 0.0;
        for (int j = 0; j < n; j++) {
            tau[i] -= Kd->data[i * n + j] * (q[j] - qstar[j]);
            tau[i] -= Kv->data[i * n + j] * qdot[j];
            if (C) tau[i] -= C->data[i * n + j] * qdot[j];
        }
    }
    int md_is_identity = 1;
    for (int i = 0; i < n && md_is_identity; i++) {
        for (int j = 0; j < n && md_is_identity; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (fabs(Md->data[i * n + j] - expected) > 1e-12) md_is_identity = 0;
        }
    }
    if (!md_is_identity && has_inv) {
        double MMd_inv[9] = {0};
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < n; k++) sum += M->data[i * n + k] * Md_inv[k * n + j];
                MMd_inv[i * n + j] = sum;
            }
        }
        for (int i = 0; i < n; i++) {
            double restoring = 0.0;
            for (int j = 0; j < n; j++) {
                restoring += Kd->data[i * n + j] * (q[j] - qstar[j]);
                restoring += Kv->data[i * n + j] * qdot[j];
                if (C) restoring += C->data[i * n + j] * qdot[j];
            }
            for (int j = 0; j < n; j++) {
                double correction = (MMd_inv[i * n + j] - (i == j ? 1.0 : 0.0)) * restoring;
                tau[i] += correction;
            }
        }
    }
}

double pbc_ida_verify_matching(const PBC_IDA_Controller* ctrl,
    const double* q, const double* qdot, const double* u, double tol) {
    if (!ctrl || !q || !qdot || !u) return 1e10;
    int n = ctrl->n;
    (void)tol;
    if (n == 1) {
        double m = ctrl->M->data[0];
        double c_val = ctrl->C ? ctrl->C->data[0] : 0.0;
        double g_val = ctrl->g[0];
        double qddot = (fabs(m) > 1e-15) ? (u[0] - c_val * qdot[0] - g_val) / m : 0.0;
        double desired = ctrl->hamiltonian.M_d->data[0] * qddot
                       + ctrl->hamiltonian.K_d->data[0] * (q[0] - ctrl->hamiltonian.q_star[0])
                       + ctrl->structure.R_d->data[0] * qdot[0];
        return fabs(desired);
    }
    double residual = 0.0;
    for (int i = 0; i < n; i++) {
        double di = 0.0;
        for (int j = 0; j < n; j++) {
            di += ctrl->hamiltonian.K_d->data[i * n + j]
                 * (q[j] - ctrl->hamiltonian.q_star[j]);
            di += ctrl->structure.R_d->data[i * n + j] * qdot[j];
        }
        if (ctrl->g) di += ctrl->g[i];
        residual += di * di;
    }
    return sqrt(residual);
}

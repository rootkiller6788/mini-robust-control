/* ==============================================================
 * pbc_port_hamiltonian.c -- Port-Hamiltonian Systems Implementation
 *
 * The port-Hamiltonian (pH) framework unifies physical system
 * modeling through energy-based structures. Every pH system with
 * R(x) >= 0 is inherently passive with storage function H(x).
 *
 * Standard form (van der Schaft 2000):
 *   dx/dt = [J(x) - R(x)] grad_H(x) + g(x) u
 *        y = g^T(x) grad_H(x)
 *
 * Power balance (inherent passivity proof):
 *   dH/dt = grad_H^T [J - R] grad_H + grad_H^T g u
 *         = -grad_H^T R grad_H + y^T u  (J skew-symmetric)
 *         <= y^T u                       (R >= 0)
 *
 * This file implements:
 *   - pH system construction and dynamics
 *   - Canonical pH models (MSD, RLC, DC motor)
 *   - Interconnection structures (feedback, parallel, series)
 *   - Dirac structures and bond graphs
 *
 * References:
 *   Maschke & van der Schaft (1992), van der Schaft (2000),
 *   Duindam et al. (2009), Ortega et al. (2002)
 * ============================================================== */

#include "pbc_port_hamiltonian.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

PBC_PH_System* pbc_ph_create(int n, int m, const PBC_Matrix* J,
                              const PBC_Matrix* R, const PBC_Matrix* g) {
    if (n <= 0 || m <= 0 || !J || !R || !g) return NULL;
    if (J->rows != n || J->cols != n) return NULL;
    if (R->rows != n || R->cols != n) return NULL;
    if (g->rows != n || g->cols != m) return NULL;
    PBC_PH_System* ph = (PBC_PH_System*)calloc(1, sizeof(PBC_PH_System));
    if (!ph) return NULL;
    ph->n = n; ph->m = m;
    ph->J = pbc_matrix_copy(J);
    ph->R = pbc_matrix_copy(R);
    ph->g = pbc_matrix_copy(g);
    ph->x = (double*)calloc((size_t)n, sizeof(double));
    ph->dH = (double*)calloc((size_t)n, sizeof(double));
    ph->H = 0.0; ph->H_dot = 0.0;
    ph->is_conservative = true;
    int nr = R->rows * R->cols;
    for (int i = 0; i < nr; i++) {
        if (fabs(R->data[i]) > 1e-15) { ph->is_conservative = false; break; }
    }
    if (!ph->J || !ph->R || !ph->g || !ph->x || !ph->dH) {
        pbc_ph_free(ph); return NULL;
    }
    return ph;
}

void pbc_ph_free(PBC_PH_System* ph) {
    if (!ph) return;
    pbc_matrix_free(ph->J); pbc_matrix_free(ph->R); pbc_matrix_free(ph->g);
    free(ph->x); free(ph->dH); free(ph);
}

void pbc_ph_dynamics(const PBC_PH_System* ph, const double* dH,
                      const double* u, double* xdot) {
    if (!ph || !dH || !xdot) return;
    int n = ph->n, m = ph->m;
    for (int i = 0; i < n; i++) {
        double J_R_dH = 0.0;
        for (int j = 0; j < n; j++) {
            J_R_dH += (ph->J->data[i * n + j] - ph->R->data[i * n + j]) * dH[j];
        }
        xdot[i] = J_R_dH;
        if (u) {
            for (int j = 0; j < m; j++) {
                xdot[i] += ph->g->data[i * m + j] * u[j];
            }
        }
    }
}

void pbc_ph_output(const PBC_PH_System* ph, const double* dH, double* y) {
    if (!ph || !dH || !y) return;
    int n = ph->n, m = ph->m;
    for (int i = 0; i < m; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++) {
            y[i] += ph->g->data[j * m + i] * dH[j];
        }
    }
}

double pbc_ph_power_balance(const PBC_PH_System* ph, const double* dH,
                             const double* u) {
    if (!ph || !dH) return 0.0;
    int n = ph->n, m = ph->m;
    double dissipated = 0.0;
    for (int i = 0; i < n; i++) {
        double row = 0.0;
        for (int j = 0; j < n; j++) {
            row += ph->R->data[i * n + j] * dH[j];
        }
        dissipated += dH[i] * row;
    }
    double supplied = 0.0;
    if (u) {
        for (int i = 0; i < m; i++) {
            double yi = 0.0;
            for (int j = 0; j < n; j++) {
                yi += ph->g->data[j * m + i] * dH[j];
            }
            supplied += yi * u[i];
        }
    }
    return supplied - dissipated;
}

bool pbc_ph_is_passive(const PBC_PH_System* ph) {
    if (!ph || !ph->R) return false;
    return pbc_matrix_is_positive_semidefinite(ph->R);
}

PBC_PH_System* pbc_ph_mass_spring_damper(double m, double k, double b) {
    if (m <= 0.0 || k < 0.0 || b < 0.0) return NULL;
    PBC_Matrix* J = pbc_matrix_zeros(2, 2);
    PBC_Matrix* R_mat = pbc_matrix_zeros(2, 2);
    PBC_Matrix* g_mat = pbc_matrix_zeros(2, 1);
    if (!J || !R_mat || !g_mat) {
        pbc_matrix_free(J); pbc_matrix_free(R_mat);
        pbc_matrix_free(g_mat); return NULL;
    }
    pbc_matrix_set(J, 0, 1, 1.0);
    pbc_matrix_set(J, 1, 0, -1.0);
    pbc_matrix_set(R_mat, 1, 1, b);
    pbc_matrix_set(g_mat, 1, 0, 1.0);
    PBC_PH_System* ph = pbc_ph_create(2, 1, J, R_mat, g_mat);
    pbc_matrix_free(J); pbc_matrix_free(R_mat); pbc_matrix_free(g_mat);
    return ph;
}

PBC_PH_System* pbc_ph_rlc_circuit(double R_val, double L, double C, bool is_series) {
    if (L <= 0.0 || C <= 0.0 || R_val < 0.0) return NULL;
    PBC_Matrix* J = pbc_matrix_zeros(2, 2);
    PBC_Matrix* R_mat = pbc_matrix_zeros(2, 2);
    PBC_Matrix* g_mat = pbc_matrix_zeros(2, 1);
    if (!J || !R_mat || !g_mat) {
        pbc_matrix_free(J); pbc_matrix_free(R_mat);
        pbc_matrix_free(g_mat); return NULL;
    }
    pbc_matrix_set(J, 0, 1, 1.0);
    pbc_matrix_set(J, 1, 0, -1.0);
    if (is_series) {
        pbc_matrix_set(R_mat, 1, 1, R_val);
        pbc_matrix_set(g_mat, 1, 0, 1.0);
    } else {
        pbc_matrix_set(R_mat, 0, 0, 1.0 / R_val);
        pbc_matrix_set(g_mat, 0, 0, 1.0);
    }
    PBC_PH_System* ph = pbc_ph_create(2, 1, J, R_mat, g_mat);
    pbc_matrix_free(J); pbc_matrix_free(R_mat); pbc_matrix_free(g_mat);
    return ph;
}

PBC_PH_System* pbc_ph_dc_motor(double R_a, double L_a, double J_m,
                                double B_m, double K_t, double K_e) {
    if (L_a <= 0.0 || J_m <= 0.0) return NULL;
    PBC_Matrix* J = pbc_matrix_zeros(3, 3);
    PBC_Matrix* R_mat = pbc_matrix_zeros(3, 3);
    PBC_Matrix* g_mat = pbc_matrix_zeros(3, 2);
    if (!J || !R_mat || !g_mat) {
        pbc_matrix_free(J); pbc_matrix_free(R_mat);
        pbc_matrix_free(g_mat); return NULL;
    }
    pbc_matrix_set(J, 0, 1, -K_e);
    pbc_matrix_set(J, 1, 0, K_t);
    pbc_matrix_set(R_mat, 0, 0, R_a);
    pbc_matrix_set(R_mat, 1, 1, B_m);
    pbc_matrix_set(g_mat, 0, 0, 1.0);
    pbc_matrix_set(g_mat, 1, 1, -1.0);
    PBC_PH_System* ph = pbc_ph_create(3, 2, J, R_mat, g_mat);
    pbc_matrix_free(J); pbc_matrix_free(R_mat); pbc_matrix_free(g_mat);
    return ph;
}

PBC_PH_System* pbc_ph_feedback_interconnect(PBC_PH_System* ph1,
                                             PBC_PH_System* ph2) {
    if (!ph1 || !ph2) return NULL;
    int n1 = ph1->n, m1 = ph1->m, n2 = ph2->n, m2 = ph2->m;
    int N = n1 + n2, M = m1 + m2;
    PBC_Matrix* Jc = pbc_matrix_zeros(N, N);
    PBC_Matrix* Rc = pbc_matrix_zeros(N, N);
    PBC_Matrix* gc = pbc_matrix_zeros(N, M);
    if (!Jc || !Rc || !gc) {
        pbc_matrix_free(Jc); pbc_matrix_free(Rc);
        pbc_matrix_free(gc); return NULL;
    }
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n1; j++) {
            Jc->data[i * N + j] = ph1->J->data[i * n1 + j];
            Rc->data[i * N + j] = ph1->R->data[i * n1 + j];
        }
        for (int j = 0; j < m1; j++) {
            gc->data[i * M + j] = ph1->g->data[i * m1 + j];
        }
    }
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < n2; j++) {
            Jc->data[(n1 + i) * N + (n1 + j)] = ph2->J->data[i * n2 + j];
            Rc->data[(n1 + i) * N + (n1 + j)] = ph2->R->data[i * n2 + j];
        }
        for (int j = 0; j < m2; j++) {
            gc->data[(n1 + i) * M + (m1 + j)] = ph2->g->data[i * m2 + j];
        }
    }
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            double val = 0.0;
            for (int k = 0; k < m1; k++) {
                val -= ph1->g->data[i * m1 + k] * ph2->g->data[j * m2 + k];
            }
            Jc->data[i * N + (n1 + j)] = val;
        }
    }
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < n1; j++) {
            double val = 0.0;
            for (int k = 0; k < m1; k++) {
                val += ph2->g->data[i * m2 + k] * ph1->g->data[j * m1 + k];
            }
            Jc->data[(n1 + i) * N + j] = val;
        }
    }
    PBC_PH_System* combined = pbc_ph_create(N, M, Jc, Rc, gc);
    pbc_matrix_free(Jc); pbc_matrix_free(Rc); pbc_matrix_free(gc);
    return combined;
}

PBC_PH_System* pbc_ph_parallel_interconnect(PBC_PH_System* ph1,
                                              PBC_PH_System* ph2) {
    if (!ph1 || !ph2) return NULL;
    if (ph1->m != ph2->m) return NULL;
    int n1 = ph1->n, n2 = ph2->n, m = ph1->m, N = n1 + n2;
    PBC_Matrix* J = pbc_matrix_zeros(N, N);
    PBC_Matrix* R_mat = pbc_matrix_zeros(N, N);
    PBC_Matrix* g_mat = pbc_matrix_zeros(N, m);
    if (!J || !R_mat || !g_mat) {
        pbc_matrix_free(J); pbc_matrix_free(R_mat);
        pbc_matrix_free(g_mat); return NULL;
    }
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n1; j++) {
            J->data[i * N + j] = ph1->J->data[i * n1 + j];
            R_mat->data[i * N + j] = ph1->R->data[i * n1 + j];
        }
        for (int j = 0; j < m; j++) {
            g_mat->data[i * m + j] = ph1->g->data[i * m + j];
        }
    }
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < n2; j++) {
            J->data[(n1 + i) * N + (n1 + j)] = ph2->J->data[i * n2 + j];
            R_mat->data[(n1 + i) * N + (n1 + j)] = ph2->R->data[i * n2 + j];
        }
        for (int j = 0; j < m; j++) {
            g_mat->data[(n1 + i) * m + j] = ph2->g->data[i * m + j];
        }
    }
    PBC_PH_System* combined = pbc_ph_create(N, m, J, R_mat, g_mat);
    pbc_matrix_free(J); pbc_matrix_free(R_mat); pbc_matrix_free(g_mat);
    return combined;
}

PBC_PH_System* pbc_ph_series_interconnect(PBC_PH_System* ph1,
                                           PBC_PH_System* ph2) {
    if (!ph1 || !ph2) return NULL;
    if (ph1->m != ph2->m) return NULL;
    int n1 = ph1->n, n2 = ph2->n, m = ph1->m, N = n1 + n2;
    PBC_Matrix* J = pbc_matrix_zeros(N, N);
    PBC_Matrix* R_mat = pbc_matrix_zeros(N, N);
    PBC_Matrix* g_mat = pbc_matrix_zeros(N, m);
    if (!J || !R_mat || !g_mat) {
        pbc_matrix_free(J); pbc_matrix_free(R_mat);
        pbc_matrix_free(g_mat); return NULL;
    }
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n1; j++) {
            J->data[i * N + j] = ph1->J->data[i * n1 + j];
            R_mat->data[i * N + j] = ph1->R->data[i * n1 + j];
        }
        for (int j = 0; j < m; j++) {
            g_mat->data[i * m + j] = ph1->g->data[i * m + j];
        }
    }
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < n2; j++) {
            J->data[(n1 + i) * N + (n1 + j)] = ph2->J->data[i * n2 + j];
            R_mat->data[(n1 + i) * N + (n1 + j)] = ph2->R->data[i * n2 + j];
        }
    }
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < n1; j++) {
            double val = 0.0;
            for (int k = 0; k < m; k++) {
                val += ph2->g->data[i * m + k] * ph1->g->data[j * m + k];
            }
            J->data[(n1 + i) * N + j] = val;
        }
    }
    PBC_PH_System* combined = pbc_ph_create(N, m, J, R_mat, g_mat);
    pbc_matrix_free(J); pbc_matrix_free(R_mat); pbc_matrix_free(g_mat);
    return combined;
}

PBC_DiracStructure* pbc_dirac_create(const PBC_Matrix* D,
                                      int dim_flow, int dim_effort) {
    if (!D || dim_flow <= 0 || dim_effort <= 0) return NULL;
    if (D->rows != dim_effort || D->cols != dim_flow) return NULL;
    PBC_DiracStructure* ds = (PBC_DiracStructure*)calloc(1,
        sizeof(PBC_DiracStructure));
    if (!ds) return NULL;
    ds->D = pbc_matrix_copy(D);
    ds->num_ports = dim_flow;
    ds->dim_flow = dim_flow;
    ds->dim_effort = dim_effort;
    return ds;
}

void pbc_dirac_free(PBC_DiracStructure* ds) {
    if (!ds) return;
    pbc_matrix_free(ds->D);
    free(ds);
}

bool pbc_dirac_is_valid(const PBC_DiracStructure* ds, double tol) {
    if (!ds || !ds->D) return false;
    if (ds->dim_flow != ds->dim_effort) return false;
    return pbc_matrix_is_skew_symmetric(ds->D, tol);
}

PBC_BG_Element* pbc_bg_element_create(int id, PBC_BG_ElementType type,
                                       const char* name, double param) {
    PBC_BG_Element* e = (PBC_BG_Element*)calloc(1, sizeof(PBC_BG_Element));
    if (!e) return NULL;
    e->id = id;
    e->type = type;
    e->name = name ? _strdup(name) : _strdup("unnamed");
    e->parameter = param;
    e->num_ports = 1;
    e->effort = (double*)calloc(1, sizeof(double));
    e->flow = (double*)calloc(1, sizeof(double));
    e->power = 0.0;
    return e;
}

void pbc_bg_element_free(PBC_BG_Element* e) {
    if (!e) return;
    free(e->name);
    free(e->effort);
    free(e->flow);
    free(e);
}

void pbc_bg_element_constitutive(PBC_BG_Element* e) {
    if (!e) return;
    switch (e->type) {
        case PBC_BG_RESISTOR:
            e->effort[0] = e->parameter * e->flow[0];
            break;
        case PBC_BG_CAPACITOR:
            break;
        case PBC_BG_INERTANCE:
            break;
        case PBC_BG_TRANSFORMER:
            if (e->num_ports >= 2) {
                e->effort[0] = e->parameter * e->effort[1];
                e->flow[1] = e->parameter * e->flow[0];
            }
            break;
        case PBC_BG_GYRATOR:
            if (e->num_ports >= 2) {
                e->effort[0] = e->parameter * e->flow[1];
                e->effort[1] = e->parameter * e->flow[0];
            }
            break;
        case PBC_BG_JUNCTION_0:
            break;
        case PBC_BG_JUNCTION_1:
            break;
        default: break;
    }
    e->power = e->effort[0] * e->flow[0];
}

PBC_PH_System* pbc_bg_to_ph(const PBC_BondGraph* bg) {
    if (!bg || bg->num_elements <= 0) return NULL;
    int n_states = 0, n_sources = 0;
    for (int i = 0; i < bg->num_elements; i++) {
        PBC_BG_ElementType t = bg->elements[i].type;
        if (t == PBC_BG_CAPACITOR || t == PBC_BG_INERTANCE) n_states++;
        if (t == PBC_BG_SOURCE_EFFORT || t == PBC_BG_SOURCE_FLOW) n_sources++;
    }
    if (n_states == 0) return NULL;
    int m_ports = (n_sources > 0) ? n_sources : 1;
    PBC_Matrix* J = pbc_matrix_zeros(n_states, n_states);
    PBC_Matrix* R_mat = pbc_matrix_zeros(n_states, n_states);
    PBC_Matrix* g_mat = pbc_matrix_zeros(n_states, m_ports);
    if (!J || !R_mat || !g_mat) {
        pbc_matrix_free(J); pbc_matrix_free(R_mat);
        pbc_matrix_free(g_mat); return NULL;
    }
    int si = 0, srci = 0, c_idx = -1, i_idx = -1;
    for (int k = 0; k < bg->num_elements; k++) {
        PBC_BG_Element* el = &bg->elements[k];
        if (el->type == PBC_BG_CAPACITOR) { c_idx = si; si++; }
        else if (el->type == PBC_BG_INERTANCE) { i_idx = si; si++; }
        else if (el->type == PBC_BG_RESISTOR) {
            if (i_idx >= 0)
                R_mat->data[i_idx * n_states + i_idx] = el->parameter;
            else if (c_idx >= 0)
                R_mat->data[c_idx * n_states + c_idx] = el->parameter;
        } else if (el->type == PBC_BG_SOURCE_EFFORT) {
            if (i_idx >= 0 && srci < m_ports)
                g_mat->data[i_idx * m_ports + srci++] = 1.0;
        } else if (el->type == PBC_BG_SOURCE_FLOW) {
            if (c_idx >= 0 && srci < m_ports)
                g_mat->data[c_idx * m_ports + srci++] = 1.0;
        }
    }
    if (c_idx >= 0 && i_idx >= 0) {
        J->data[c_idx * n_states + i_idx] = -1.0;
        J->data[i_idx * n_states + c_idx] = 1.0;
    }
    PBC_PH_System* ph = pbc_ph_create(n_states, m_ports, J, R_mat, g_mat);
    pbc_matrix_free(J); pbc_matrix_free(R_mat); pbc_matrix_free(g_mat);
    return ph;
}

void pbc_bg_print(const PBC_BondGraph* bg) {
    if (!bg) { printf("BondGraph: (null)\n"); return; }
    printf("=== Bond Graph: %d elements ===\n", bg->num_elements);
    const char* type_names[] = {
        "R", "C", "I", "TF", "GY", "Se", "Sf", "0-junc", "1-junc"
    };
    for (int i = 0; i < bg->num_elements; i++) {
        PBC_BG_Element* e = &bg->elements[i];
        const char* tn = (e->type <= PBC_BG_JUNCTION_1) ?
                         type_names[e->type] : "?";
        printf("  [%d] %s \"%s\" param=%.4f power=%.4f\n",
               e->id, tn, e->name ? e->name : "", e->parameter, e->power);
    }
}

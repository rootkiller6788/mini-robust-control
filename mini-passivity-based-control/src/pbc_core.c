/* ==============================================================
 * pbc_core.c — Passivity-Based Control Core Implementation
 *
 * Implements foundational data structures and operations for
 * passivity-based control theory:
 *   - Matrix/Vector algebra (the mathematical substrate)
 *   - Storage function construction and evaluation
 *   - Generic passive system descriptor
 *   - Dissipation inequality verification
 *   - Numerical integration (Euler, RK4)
 *
 * All operations are numerically grounded and handle edge cases
 * (NULL pointers, dimension mismatches, singular matrices).
 *
 * References:
 *   Willems (1972), Hill & Moylan (1976), Khalil (2002)
 *   van der Schaft (2000), Ortega et al. (1998)
 * ============================================================== */

#include "pbc_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ==============================================================
 * Internal helpers
 * ============================================================== */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double pbc_max(double a, double b) { return (a > b) ? a : b; }
static double pbc_min(double a, double b) { return (a < b) ? a : b; }

static bool pbc_index_valid(int row, int col, int rows, int cols) {
    return (row >= 0 && row < rows && col >= 0 && col < cols);
}

/* ==============================================================
 * Matrix Operations
 * ============================================================== */

PBC_Matrix* pbc_matrix_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0) return NULL;
    PBC_Matrix* m = (PBC_Matrix*)malloc(sizeof(PBC_Matrix));
    if (!m) return NULL;
    m->data = (double*)calloc((size_t)(rows * cols), sizeof(double));
    if (!m->data) {
        free(m);
        return NULL;
    }
    m->rows = rows;
    m->cols = cols;
    return m;
}

void pbc_matrix_free(PBC_Matrix* m) {
    if (!m) return;
    free(m->data);
    free(m);
}

PBC_Matrix* pbc_matrix_copy(const PBC_Matrix* m) {
    if (!m) return NULL;
    PBC_Matrix* copy = pbc_matrix_create(m->rows, m->cols);
    if (!copy) return NULL;
    int n = m->rows * m->cols;
    for (int i = 0; i < n; i++) {
        copy->data[i] = m->data[i];
    }
    return copy;
}

PBC_Matrix* pbc_matrix_identity(int n) {
    PBC_Matrix* m = pbc_matrix_zeros(n, n);
    if (!m) return NULL;
    for (int i = 0; i < n; i++) {
        m->data[i * n + i] = 1.0;
    }
    return m;
}

PBC_Matrix* pbc_matrix_zeros(int rows, int cols) {
    return pbc_matrix_create(rows, cols);
}

void pbc_matrix_set(PBC_Matrix* m, int row, int col, double val) {
    if (!m || !pbc_index_valid(row, col, m->rows, m->cols)) return;
    m->data[row * m->cols + col] = val;
}

double pbc_matrix_get(const PBC_Matrix* m, int row, int col) {
    if (!m || !pbc_index_valid(row, col, m->rows, m->cols)) return 0.0;
    return m->data[row * m->cols + col];
}

void pbc_matrix_fill(PBC_Matrix* m, double val) {
    if (!m) return;
    int n = m->rows * m->cols;
    for (int i = 0; i < n; i++) m->data[i] = val;
}

PBC_Matrix* pbc_matrix_add(const PBC_Matrix* a, const PBC_Matrix* b) {
    if (!a || !b || a->rows != b->rows || a->cols != b->cols) return NULL;
    PBC_Matrix* c = pbc_matrix_create(a->rows, a->cols);
    if (!c) return NULL;
    int n = a->rows * a->cols;
    for (int i = 0; i < n; i++) c->data[i] = a->data[i] + b->data[i];
    return c;
}

PBC_Matrix* pbc_matrix_sub(const PBC_Matrix* a, const PBC_Matrix* b) {
    if (!a || !b || a->rows != b->rows || a->cols != b->cols) return NULL;
    PBC_Matrix* c = pbc_matrix_create(a->rows, a->cols);
    if (!c) return NULL;
    int n = a->rows * a->cols;
    for (int i = 0; i < n; i++) c->data[i] = a->data[i] - b->data[i];
    return c;
}

PBC_Matrix* pbc_matrix_mul(const PBC_Matrix* a, const PBC_Matrix* b) {
    if (!a || !b || a->cols != b->rows) return NULL;
    PBC_Matrix* c = pbc_matrix_create(a->rows, b->cols);
    if (!c) return NULL;
    for (int i = 0; i < a->rows; i++) {
        for (int k = 0; k < a->cols; k++) {
            double aik = a->data[i * a->cols + k];
            if (fabs(aik) < 1e-30) continue;
            for (int j = 0; j < b->cols; j++) {
                c->data[i * c->cols + j] += aik * b->data[k * b->cols + j];
            }
        }
    }
    return c;
}

PBC_Matrix* pbc_matrix_scale(const PBC_Matrix* a, double alpha) {
    if (!a) return NULL;
    PBC_Matrix* c = pbc_matrix_create(a->rows, a->cols);
    if (!c) return NULL;
    int n = a->rows * a->cols;
    for (int i = 0; i < n; i++) c->data[i] = a->data[i] * alpha;
    return c;
}

PBC_Matrix* pbc_matrix_transpose(const PBC_Matrix* a) {
    if (!a) return NULL;
    PBC_Matrix* t = pbc_matrix_create(a->cols, a->rows);
    if (!t) return NULL;
    for (int i = 0; i < a->rows; i++) {
        for (int j = 0; j < a->cols; j++) {
            t->data[j * t->cols + i] = a->data[i * a->cols + j];
        }
    }
    return t;
}

double pbc_matrix_norm_frobenius(const PBC_Matrix* a) {
    if (!a) return NAN;
    double sum = 0.0;
    int n = a->rows * a->cols;
    for (int i = 0; i < n; i++) sum += a->data[i] * a->data[i];
    return sqrt(sum);
}

double pbc_matrix_trace(const PBC_Matrix* a) {
    if (!a || a->rows != a->cols) return NAN;
    double tr = 0.0;
    for (int i = 0; i < a->rows; i++) tr += a->data[i * a->cols + i];
    return tr;
}

double pbc_matrix_det(const PBC_Matrix* a) {
    if (!a || a->rows != a->cols) return NAN;
    if (a->rows == 1) return a->data[0];
    if (a->rows == 2) {
        return a->data[0] * a->data[3] - a->data[1] * a->data[2];
    }
    if (a->rows == 3) {
        double d = a->data[0] * (a->data[4] * a->data[8] - a->data[5] * a->data[7])
                 - a->data[1] * (a->data[3] * a->data[8] - a->data[5] * a->data[6])
                 + a->data[2] * (a->data[3] * a->data[7] - a->data[4] * a->data[6]);
        return d;
    }
    return NAN; /* Only 1x1, 2x2, 3x3 supported directly */
}

bool pbc_matrix_is_symmetric(const PBC_Matrix* a, double tol) {
    if (!a || a->rows != a->cols) return false;
    for (int i = 0; i < a->rows; i++) {
        for (int j = i + 1; j < a->cols; j++) {
            double diff = a->data[i * a->cols + j] - a->data[j * a->cols + i];
            if (fabs(diff) > tol) return false;
        }
    }
    return true;
}

bool pbc_matrix_is_positive_definite(const PBC_Matrix* a) {
    if (!a || a->rows != a->cols || a->rows > 3) return false;

    /* Check symmetry first */
    if (!pbc_matrix_is_symmetric(a, 1e-10)) return false;

    /* Sylvester's criterion: all leading principal minors > 0 */
    if (a->rows >= 1 && a->data[0] <= 0.0) return false;
    if (a->rows >= 2) {
        double d2 = a->data[0]*a->data[3] - a->data[1]*a->data[2];
        if (d2 <= 0.0) return false;
    }
    if (a->rows >= 3) {
        double d3 = pbc_matrix_det(a);
        if (d3 <= 0.0) return false;
    }
    return true;
}

bool pbc_matrix_is_positive_semidefinite(const PBC_Matrix* a) {
    if (!a || a->rows != a->cols || a->rows > 3) return false;
    if (!pbc_matrix_is_symmetric(a, 1e-10)) return false;
    if (a->rows >= 1 && a->data[0] < -1e-12) return false;
    if (a->rows >= 2) {
        double d2 = a->data[0]*a->data[3] - a->data[1]*a->data[2];
        if (d2 < -1e-12) return false;
    }
    if (a->rows >= 3) {
        double d3 = pbc_matrix_det(a);
        if (d3 < -1e-12) return false;
    }
    return true;
}

bool pbc_matrix_is_skew_symmetric(const PBC_Matrix* a, double tol) {
    if (!a || a->rows != a->cols) return false;
    for (int i = 0; i < a->rows; i++) {
        /* Diagonal must be zero */
        if (fabs(a->data[i * a->cols + i]) > tol) return false;
        for (int j = i + 1; j < a->cols; j++) {
            double sum = a->data[i * a->cols + j] + a->data[j * a->cols + i];
            if (fabs(sum) > tol) return false;
        }
    }
    return true;
}

void pbc_matrix_print(const PBC_Matrix* m, const char* label) {
    if (!m) { printf("%s: (null)\n", label ? label : "matrix"); return; }
    printf("%s (%d×%d):\n", label ? label : "matrix", m->rows, m->cols);
    for (int i = 0; i < m->rows; i++) {
        printf("  ");
        for (int j = 0; j < m->cols; j++) {
            printf("% 10.4f ", m->data[i * m->cols + j]);
        }
        printf("\n");
    }
}

/* ==============================================================
 * Vector Operations
 * ============================================================== */

PBC_Vector* pbc_vector_create(int size) {
    if (size <= 0) return NULL;
    PBC_Vector* v = (PBC_Vector*)malloc(sizeof(PBC_Vector));
    if (!v) return NULL;
    v->data = (double*)calloc((size_t)size, sizeof(double));
    if (!v->data) {
        free(v);
        return NULL;
    }
    v->size = size;
    return v;
}

void pbc_vector_free(PBC_Vector* v) {
    if (!v) return;
    free(v->data);
    free(v);
}

void pbc_vector_set(PBC_Vector* v, int idx, double val) {
    if (!v || idx < 0 || idx >= v->size) return;
    v->data[idx] = val;
}

double pbc_vector_get(const PBC_Vector* v, int idx) {
    if (!v || idx < 0 || idx >= v->size) return 0.0;
    return v->data[idx];
}

double pbc_vector_norm(const PBC_Vector* v) {
    if (!v) return NAN;
    double sum = 0.0;
    for (int i = 0; i < v->size; i++) sum += v->data[i] * v->data[i];
    return sqrt(sum);
}

double pbc_vector_dot(const PBC_Vector* a, const PBC_Vector* b) {
    if (!a || !b || a->size != b->size) return 0.0;
    double d = 0.0;
    for (int i = 0; i < a->size; i++) d += a->data[i] * b->data[i];
    return d;
}

void pbc_vector_print(const PBC_Vector* v, const char* label) {
    if (!v) { printf("%s: (null)\n", label ? label : "vector"); return; }
    printf("%s (%d): [ ", label ? label : "vector", v->size);
    for (int i = 0; i < v->size; i++) printf("% 8.4f ", v->data[i]);
    printf("]\n");
}

/* ==============================================================
 * Storage Function
 * ============================================================== */

PBC_StorageFunction* pbc_storage_quadratic_init(int n, const PBC_Matrix* P) {
    if (n <= 0 || !P || P->rows != n || P->cols != n) return NULL;
    PBC_StorageFunction* sf = (PBC_StorageFunction*)malloc(
        sizeof(PBC_StorageFunction));
    if (!sf) return NULL;
    sf->type = PBC_STORAGE_QUADRATIC;
    sf->P = pbc_matrix_copy(P);
    sf->gamma = 0.0;
    sf->min_value = 0.0;
    sf->is_positive_definite = pbc_matrix_is_positive_definite(P);
    sf->is_radially_unbounded = sf->is_positive_definite;
    sf->n_states = n;
    return sf;
}

PBC_StorageFunction* pbc_storage_mechanical_init(int n,
    const PBC_Matrix* M, const PBC_Matrix* K) {
    if (n <= 0 || !M || !K || M->rows != n || K->rows != n) return NULL;

    /* For mechanical systems, storage = kinetic + potential energy.
     * Construct block-diagonal matrix P = blkdiag(M, K) for state [q̇; q].
     * S(q,q̇) = ½ q̇ᵀ M q̇ + ½ qᵀ K q = ½ [q̇ᵀ qᵀ] [M 0; 0 K] [q̇; q] */
    PBC_Matrix* P = pbc_matrix_zeros(2 * n, 2 * n);
    if (!P) { return NULL; }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            P->data[i * (2*n) + j] = M->data[i * n + j];
            P->data[(n+i) * (2*n) + (n+j)] = K->data[i * n + j];
        }
    }

    PBC_StorageFunction* sf = (PBC_StorageFunction*)malloc(
        sizeof(PBC_StorageFunction));
    if (!sf) { pbc_matrix_free(P); return NULL; }
    sf->type = PBC_STORAGE_TOTAL_ENERGY;
    sf->P = P;
    sf->gamma = 0.0;
    sf->min_value = 0.0;
    sf->is_positive_definite = pbc_matrix_is_positive_definite(M) &&
                               pbc_matrix_is_positive_semidefinite(K);
    sf->is_radially_unbounded = sf->is_positive_definite;
    sf->n_states = 2 * n;
    return sf;
}

void pbc_storage_free(PBC_StorageFunction* sf) {
    if (!sf) return;
    pbc_matrix_free(sf->P);
    free(sf);
}

double pbc_storage_eval(const PBC_StorageFunction* sf, const double* x) {
    if (!sf || !sf->P || !x) return 0.0;
    int n = sf->P->rows;
    double val = 0.0;
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) {
            row_sum += sf->P->data[i * n + j] * x[j];
        }
        val += x[i] * row_sum;
    }
    return 0.5 * val;
}

void pbc_storage_gradient(const PBC_StorageFunction* sf,
                           const double* x, double* grad) {
    if (!sf || !sf->P || !x || !grad) return;
    int n = sf->P->rows;
    for (int i = 0; i < n; i++) {
        grad[i] = 0.0;
        for (int j = 0; j < n; j++) {
            grad[i] += sf->P->data[i * n + j] * x[j];
        }
    }
}

void pbc_storage_hessian(const PBC_StorageFunction* sf,
                          const double* x, PBC_Matrix* H) {
    if (!sf || !sf->P || !H) return;
    (void)x; /* Hessian is constant P for quadratic */
    int n = sf->P->rows;
    for (int i = 0; i < n * n; i++) H->data[i] = sf->P->data[i];
}

/* ==============================================================
 * Passive System
 * ============================================================== */

PBC_System* pbc_system_create(const char* name, int n, int m, int p,
                               PBC_StateFunc f, PBC_OutputFunc h,
                               PBC_StorageEval S, void* data) {
    if (n <= 0 || m <= 0 || p <= 0) return NULL;
    PBC_System* sys = (PBC_System*)calloc(1, sizeof(PBC_System));
    if (!sys) return NULL;
    sys->name = name ? _strdup(name) : _strdup("unnamed");
    sys->n = n;
    sys->m = m;
    sys->p = p;
    sys->f = f;
    sys->h = h;
    sys->S = S;
    sys->data = data;
    sys->x = (double*)calloc((size_t)n, sizeof(double));
    if (!sys->x) { free(sys->name); free(sys); return NULL; }
    sys->t = 0.0;
    sys->is_initialized = true;
    return sys;
}

void pbc_system_free(PBC_System* sys) {
    if (!sys) return;
    free(sys->name);
    free(sys->x);
    free(sys);
}

void pbc_system_set_state(PBC_System* sys, const double* x0) {
    if (!sys || !x0) return;
    for (int i = 0; i < sys->n; i++) sys->x[i] = x0[i];
}

void pbc_system_f(const PBC_System* sys, const double* x,
                   const double* u, double* xdot) {
    if (!sys || !sys->f || !x || !xdot) return;
    sys->f(sys, x, u, xdot);
}

void pbc_system_h(const PBC_System* sys, const double* x, double* y) {
    if (!sys || !sys->h || !x || !y) return;
    sys->h(sys, x, y);
}

double pbc_supply_rate(const double* u, const double* y, int dim) {
    if (!u || !y || dim <= 0) return 0.0;
    double w = 0.0;
    for (int i = 0; i < dim; i++) w += u[i] * y[i];
    return w;
}

/* ==============================================================
 * Dissipation Inequality Verification
 *
 * Checks: S(x(T)) - S(x(0)) ≤ ∫₀ᵀ uᵀ(t) y(t) dt
 *
 * For a passive system, the "dissipation" value (S(T) - S(0) - ∫uᵀy)
 * must be ≤ 0.
 *
 * Uses trapezoidal integration for the supply rate integral.
 * ============================================================== */

PBC_Dissipation pbc_verify_dissipation(const PBC_System* sys,
    const double* t_span, const double* u_traj, int n_pts) {
    PBC_Dissipation result;
    memset(&result, 0, sizeof(result));

    if (!sys || !t_span || !u_traj || n_pts < 2) {
        result.satisfied = false;
        return result;
    }

    /* Save original state */
    double* x_save = (double*)malloc((size_t)sys->n * sizeof(double));
    if (!x_save) { result.satisfied = false; return result; }
    memcpy(x_save, sys->x, (size_t)sys->n * sizeof(double));

    /* Storage at initial time */
    result.S_initial = sys->S ? sys->S(sys, sys->x) : 0.0;

    /* Integrate supply rate ∫ uᵀ y dt using trapezoidal rule */
    double* y = (double*)malloc((size_t)sys->p * sizeof(double));
    if (!y) { free(x_save); result.satisfied = false; return result; }

    double integral = 0.0;
    for (int k = 0; k < n_pts; k++) {
        /* Set state and compute output */
        sys->h(sys, sys->x, y);
        const double* uk = &u_traj[k * sys->m];
        double w_k = pbc_supply_rate(uk, y, sys->m);

        if (k == 0 || k == n_pts - 1) {
            integral += 0.5 * w_k * (t_span[1] - t_span[0]);
        } else {
            integral += w_k * (t_span[1] - t_span[0]);
        }

        /* Step forward (if not last point) */
        if (k < n_pts - 1) {
            double dt = t_span[k+1] - t_span[k];
            pbc_step_euler((PBC_System*)sys, uk, dt);
        }
    }

    result.supply_integral = integral;
    result.S_current = sys->S ? sys->S(sys, sys->x) : 0.0;
    result.dissipation = result.S_current - result.S_initial - integral;
    result.satisfied = (result.dissipation <= 1e-8);

    /* Restore state */
    memcpy(sys->x, x_save, (size_t)sys->n * sizeof(double));

    free(x_save);
    free(y);
    return result;
}

/* ==============================================================
 * Numerical Integration
 * ============================================================== */

void pbc_step_euler(PBC_System* sys, const double* u, double dt) {
    if (!sys || !u || dt <= 0.0) return;
    double* xdot = (double*)malloc((size_t)sys->n * sizeof(double));
    if (!xdot) return;
    sys->f(sys, sys->x, u, xdot);
    for (int i = 0; i < sys->n; i++) {
        sys->x[i] += dt * xdot[i];
    }
    sys->t += dt;
    free(xdot);
}

void pbc_step_rk4(PBC_System* sys, const double* u, double dt) {
    if (!sys || !u || dt <= 0.0) return;
    int n = sys->n;
    double *k1, *k2, *k3, *k4, *xtmp;
    k1 = (double*)malloc((size_t)(5 * n) * sizeof(double));
    if (!k1) return;
    k2 = k1 + n; k3 = k1 + 2*n; k4 = k1 + 3*n; xtmp = k1 + 4*n;

    double* x_orig = (double*)malloc((size_t)n * sizeof(double));
    if (!x_orig) { free(k1); return; }
    memcpy(x_orig, sys->x, (size_t)n * sizeof(double));

    /* k1 */
    sys->f(sys, sys->x, u, k1);
    /* k2 */
    for (int i = 0; i < n; i++) xtmp[i] = sys->x[i] + 0.5*dt*k1[i];
    sys->f(sys, xtmp, u, k2);
    /* k3 */
    for (int i = 0; i < n; i++) xtmp[i] = sys->x[i] + 0.5*dt*k2[i];
    sys->f(sys, xtmp, u, k3);
    /* k4 */
    for (int i = 0; i < n; i++) xtmp[i] = sys->x[i] + dt*k3[i];
    sys->f(sys, xtmp, u, k4);
    /* Update */
    for (int i = 0; i < n; i++) {
        sys->x[i] += (dt/6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
    sys->t += dt;
    free(x_orig);
    free(k1);
}

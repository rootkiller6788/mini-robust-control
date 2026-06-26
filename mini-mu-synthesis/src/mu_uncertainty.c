/* ====================================================================
 * mu_uncertainty.c — Uncertainty Modeling for u-Synthesis
 *
 * Models various types of uncertainty and builds the structured
 * uncertainty representation Delta for real-world systems.
 *
 * Uncertainty Types:
 *   1. Parametric: physical parameters within bounds (real scalars)
 *   2. Unmodeled dynamics: neglected high-frequency dynamics (complex full)
 *   3. LTI uncertainty: neglected linear dynamics (complex full)
 *
 * Knowledge Coverage:
 *   L1: Uncertainty weight template and weighting function creation
 *   L3: Parametric uncertainty modeling via LFT extraction
 *   L3: Uncertainty structure assembly from block descriptions
 *   L3: Input/output multiplicative uncertainty modeling
 *   L3: Additive uncertainty modeling
 *   L3: Coprime factor uncertainty (McFarlane & Glover 1990)
 *   L3: Diagonal input uncertainty for MIMO systems
 *   L7: Boeing 747 lateral-directional uncertainty model
 *   L7: Tesla electric vehicle lateral uncertainty model
 *   L7: SpaceX Falcon 9 rocket attitude uncertainty model
 *
 * Reference:
 *   Skogestad & Postlethwaite (2005), Ch. 7-8
 *   Zhou, Doyle & Glover (1996), Ch. 9
 *   McFarlane & Glover (1990), IEEE TAC
 * ==================================================================== */

#include "mu_core.h"
#include "mu_uncertainty.h"
#include "mu_matrix.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================
 * L1: Uncertainty Weight Function
 *
 * W_u(s) = gain * (s + omega_b * r0) / (s * rinf + omega_b)
 *
 * This weight represents low uncertainty at low frequencies
 * (good low-frequency model) and high uncertainty at high
 * frequencies (unmodeled dynamics).
 *
 * State-space realization:
 *   A = -omega_b / rinf
 *   B = 1
 *   C = gain * (r0 - rinf) * omega_b / (rinf * rinf)
 *   D = gain * r0 / rinf
 *
 * Complexity: O(1)
 * Reference: Skogestad & Postlethwaite (2005), §7.4
 * ============================ */

MUSystem* mu_unc_weight_create(const MUUncertaintyWeight *w) {
    if (!w || w->rinf < MU_EPSILON) return NULL;

    MUSystem *W = mu_system_create(1, 1, 1);
    if (!W) return NULL;

    double a_val = -w->omega_b / w->rinf;
    double b_val = 1.0;
    double c_val = w->gain * (w->r0 - w->rinf) * w->omega_b / (w->rinf * w->rinf);
    double d_val = w->gain * w->r0 / w->rinf;

    mu_real_mat_set(W->A, 0, 0, a_val);
    mu_real_mat_set(W->B, 0, 0, b_val);
    mu_real_mat_set(W->C, 0, 0, c_val);
    mu_real_mat_set(W->D, 0, 0, d_val);

    return W;
}

/* ============================
 * L3: Build Uncertainty Structure from Block Descriptions
 *
 * Assemble Delta = diag(block_1, block_2, ..., block_n)
 *
 * Each block is characterized by type, size, and bound.
 * The total dimension is sum of block dimensions.
 * ============================ */

MUUncertaintyStructure* mu_unc_build_structure(
    const MUUncertaintyBlock *blocks, int num_blocks) {
    if (!blocks || num_blocks <= 0) return NULL;

    MUUncertaintyStructure *s = mu_unc_struct_create(num_blocks);
    if (!s) return NULL;

    int total = 0;
    for (int i = 0; i < num_blocks; i++) {
        s->blocks[i] = blocks[i];
        total += blocks[i].size;
        if (blocks[i].type == MU_UNC_REPEATED_SCALAR ||
            blocks[i].type == MU_UNC_REAL_SCALAR) {
            /* Repeated scalar counts as 1 uncertainty parameter */
        }
    }
    s->total_dim = total;
    return s;
}

/* ============================
 * L3: Input Multiplicative Uncertainty
 *
 * G_true = G_nom * (I + W_I * Delta_I), ||Delta_I||_inf <= 1
 *
 * LFT form:
 *   [ z ]   [ 0       W_I * G_nom ] [ w ]
 *   [ y ] = [ G_nom   W_I * G_nom ] [ u ]
 *
 * Complexity: O(n^3)
 * ============================ */

MUSystem* mu_unc_input_multiplicative(const MUSystem *G_nom,
                                       const MUSystem *W_I) {
    if (!G_nom) return NULL;

    int nG = G_nom->n, mG = G_nom->m, pG = G_nom->p;
    int nW = (W_I ? W_I->n : 0);

    int n = nG + nW;
    int m = mG * 2;   /* uncertainty input + control input */
    int p = pG * 2;   /* performance output + measurement output */

    MUSystem *P = mu_system_create(n, m, p);
    if (!P) return NULL;

    /* Place G_nom in the interconnection */
    for (int i = 0; i < nG; i++)
        for (int j = 0; j < nG; j++)
            mu_real_mat_set(P->A, i, j,
                mu_real_mat_get(G_nom->A, i, j));

    /* B1 = G_nom->B (disturbance from uncertainty) */
    /* B2 = G_nom->B (control input) */
    for (int i = 0; i < nG; i++) {
        for (int j = 0; j < mG; j++) {
            mu_real_mat_set(P->B, i, j,
                mu_real_mat_get(G_nom->B, i, j));
            mu_real_mat_set(P->B, i, mG + j,
                mu_real_mat_get(G_nom->B, i, j));
        }
    }

    /* C1 = G_nom->C (uncertainty output to performance) */
    /* C2 = G_nom->C (measurement output) */
    for (int i = 0; i < pG; i++)
        for (int j = 0; j < nG; j++) {
            mu_real_mat_set(P->C, i, j,
                mu_real_mat_get(G_nom->C, i, j));
            mu_real_mat_set(P->C, pG + i, j,
                mu_real_mat_get(G_nom->C, i, j));
        }

    /* D matrix */
    for (int i = 0; i < pG; i++) {
        for (int j = 0; j < mG; j++) {
            /* D11 (uncertainty): 0 initially, modified by W_I */
            /* D12 (control): G_nom->D */
            /* D21 (measurement from uncertainty): G_nom->D */
            /* D22 (measurement from control): G_nom->D */
            mu_real_mat_set(P->D, i, mG + j,
                mu_real_mat_get(G_nom->D, i, j));
            mu_real_mat_set(P->D, pG + i, j,
                mu_real_mat_get(G_nom->D, i, j));
            mu_real_mat_set(P->D, pG + i, mG + j,
                mu_real_mat_get(G_nom->D, i, j));
        }
    }

    /* Place W_I dynamics */
    if (W_I && nW > 0) {
        for (int i = 0; i < nW; i++)
            for (int j = 0; j < nW; j++)
                mu_real_mat_set(P->A, nG + i, nG + j,
                    mu_real_mat_get(W_I->A, i, j));
        /* W_I input from control signal */
        for (int i = 0; i < nW; i++)
            for (int j = 0; j < mG; j++)
                mu_real_mat_set(P->B, nG + i, j,
                    mu_real_mat_get(W_I->B, i, j));
        /* W_I output adds to performance */
        for (int i = 0; i < pG; i++)
            for (int j = 0; j < nW; j++)
                mu_real_mat_set(P->A, i, nG + j, 0.0);
    }

    return P;
}

/* ============================
 * L3: Output Multiplicative Uncertainty
 *
 * G_true = (I + W_O * Delta_O) * G_nom
 *
 * LFT form:
 *   [ z ]   [ 0       G_nom * W_O ] [ w ]
 *   [ y ] = [ G_nom   G_nom * W_O ] [ u ]
 * ============================ */

MUSystem* mu_unc_output_multiplicative(const MUSystem *G_nom,
                                        const MUSystem *W_O) {
    if (!G_nom) return NULL;

    int nG = G_nom->n, mG = G_nom->m, pG = G_nom->p;
    int nW = (W_O ? W_O->n : 0);

    int n = nG + nW;
    int m = mG * 2;
    int p = pG * 2;

    MUSystem *P = mu_system_create(n, m, p);
    if (!P) return NULL;

    /* Similar structure to input multiplicative but with output weighting */
    for (int i = 0; i < nG; i++)
        for (int j = 0; j < nG; j++)
            mu_real_mat_set(P->A, i, j,
                mu_real_mat_get(G_nom->A, i, j));

    for (int i = 0; i < nG; i++) {
        for (int j = 0; j < mG; j++) {
            mu_real_mat_set(P->B, i, j, 0.0);
            mu_real_mat_set(P->B, i, mG + j,
                mu_real_mat_get(G_nom->B, i, j));
        }
    }

    for (int i = 0; i < pG; i++)
        for (int j = 0; j < nG; j++) {
            mu_real_mat_set(P->C, i, j, 0.0);
            mu_real_mat_set(P->C, pG + i, j,
                mu_real_mat_get(G_nom->C, i, j));
        }

    for (int i = 0; i < pG; i++) {
        for (int j = 0; j < mG; j++) {
            mu_real_mat_set(P->D, i, mG + j, 0.0);
            mu_real_mat_set(P->D, pG + i, j, 0.0);
            mu_real_mat_set(P->D, pG + i, mG + j,
                mu_real_mat_get(G_nom->D, i, j));
        }
    }

    /* Couple W_O */
    if (W_O && nW > 0) {
        for (int i = 0; i < nW; i++)
            for (int j = 0; j < nW; j++)
                mu_real_mat_set(P->A, nG + i, nG + j,
                    mu_real_mat_get(W_O->A, i, j));
        for (int i = 0; i < nW; i++)
            for (int j = 0; j < mG; j++)
                mu_real_mat_set(P->B, nG + i, j, 0.0);
        for (int i = 0; i < pG; i++)
            for (int j = 0; j < nW; j++)
                mu_real_mat_set(P->C, i, nG + j,
                    mu_real_mat_get(W_O->C, i, j));
        for (int i = 0; i < nW; i++)
            for (int j = 0; j < mG; j++)
                mu_real_mat_set(P->B, nG + i, mG + j,
                    mu_real_mat_get(W_O->D, i, j));
    }

    return P;
}

/* ============================
 * L3: Additive Uncertainty
 *
 * G_true = G_nom + W_A * Delta_A, ||Delta_A||_inf <= 1
 *
 * LFT form:
 *   [ z ]   [ 0      W_A ] [ w ]
 *   [ y ] = [ I      G_nom ] [ u ]
 * ============================ */

MUSystem* mu_unc_additive(const MUSystem *G_nom, const MUSystem *W_A) {
    if (!G_nom || !W_A) return NULL;

    int nG = G_nom->n, mG = G_nom->m, pG = G_nom->p;
    int nW = W_A->n, mW = W_A->m;

    int n = nG + nW;
    int m = mW + mG;
    int p = pG + pG;

    MUSystem *P = mu_system_create(n, m, p);
    if (!P) return NULL;

    /* Place G_nom and W_A in block-diagonal form */
    for (int i = 0; i < nG; i++)
        for (int j = 0; j < nG; j++)
            mu_real_mat_set(P->A, i, j,
                mu_real_mat_get(G_nom->A, i, j));
    for (int i = 0; i < nW; i++)
        for (int j = 0; j < nW; j++)
            mu_real_mat_set(P->A, nG + i, nG + j,
                mu_real_mat_get(W_A->A, i, j));

    for (int i = 0; i < nG; i++)
        for (int j = 0; j < mG; j++)
            mu_real_mat_set(P->B, i, mW + j,
                mu_real_mat_get(G_nom->B, i, j));
    for (int i = 0; i < nW; i++)
        for (int j = 0; j < mW; j++)
            mu_real_mat_set(P->B, nG + i, j,
                mu_real_mat_get(W_A->B, i, j));

    for (int i = 0; i < pG; i++) {
        for (int j = 0; j < nG; j++)
            mu_real_mat_set(P->C, pG + i, j,
                mu_real_mat_get(G_nom->C, i, j));
        for (int j = 0; j < nW; j++)
            mu_real_mat_set(P->C, i, nG + j,
                mu_real_mat_get(W_A->C, i, j));
    }

    for (int i = 0; i < pG; i++) {
        for (int j = 0; j < mG; j++)
            mu_real_mat_set(P->D, pG + i, mW + j,
                mu_real_mat_get(G_nom->D, i, j));
        for (int j = 0; j < mW; j++)
            mu_real_mat_set(P->D, i, j,
                mu_real_mat_get(W_A->D, i, j));
    }

    return P;
}

/* ============================
 * L3: Coprime Factor Uncertainty
 *
 * G = (M + Delta_M)^{-1} (N + Delta_N), ||[Delta_N, Delta_M]||_inf <= epsilon
 *
 * Normalized coprime factorization: G = M^{-1} N
 *
 * The robust stabilization problem has an explicit solution:
 *   gamma_min = (1 + lambda_max(XZ))^{1/2}
 * where X and Z are solutions to the normalized coprime factor
 * robust stabilization Riccati equations.
 *
 * Complexity: O(n^3)
 * Reference: McFarlane & Glover (1990), IEEE TAC
 * ============================ */

MUSystem* mu_unc_coprime_factor(const MUSystem *G_nom, double epsilon) {
    if (!G_nom || epsilon <= 0.0) return NULL;

    /* Build the extended plant for coprime factor uncertainty */
    int n = G_nom->n, m = G_nom->m, p = G_nom->p;

    MUSystem *P = mu_system_create(n, m + p, p + m);
    if (!P) return NULL;

    /* P = [ A,  [B, -H*C'], [B; 0]; [C; 0], [D, I; I, 0], [D, I; 0, 0] ] */
    /* Simplified: return a copy of G with epsilon weighting */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            mu_real_mat_set(P->A, i, j,
                mu_real_mat_get(G_nom->A, i, j));
        for (int j = 0; j < m; j++)
            mu_real_mat_set(P->B, i, j,
                mu_real_mat_get(G_nom->B, i, j) * epsilon);
    }
    for (int i = 0; i < p; i++)
        for (int j = 0; j < n; j++)
            mu_real_mat_set(P->C, i, j,
                mu_real_mat_get(G_nom->C, i, j));
    for (int i = 0; i < p; i++)
        for (int j = 0; j < m; j++)
            mu_real_mat_set(P->D, i, j,
                mu_real_mat_get(G_nom->D, i, j) * epsilon);

    return P;
}

/* ============================
 * L3: Diagonal Input Uncertainty
 *
 * For MIMO systems with m independent actuators, each with
 * independent uncertainty Delta_I = diag(delta_1, ..., delta_m).
 *
 * Complexity: O(n^3)
 * ============================ */

MUSystem* mu_unc_diagonal_input(const MUSystem *G_nom,
                                 const double *weights, int m) {
    if (!G_nom || !weights || m <= 0) return NULL;

    int nG = G_nom->n, pG = G_nom->p;

    int n = nG;
    int p = pG + m;

    MUSystem *P = mu_system_create(n, m * 2, p);
    if (!P) return NULL;

    /* Plant dynamics */
    for (int i = 0; i < nG; i++)
        for (int j = 0; j < nG; j++)
            mu_real_mat_set(P->A, i, j,
                mu_real_mat_get(G_nom->A, i, j));

    /* B: uncertainty inputs (weighted) + control inputs */
    for (int i = 0; i < nG; i++)
        for (int j = 0; j < m; j++) {
            mu_real_mat_set(P->B, i, j,
                mu_real_mat_get(G_nom->B, i, j) * weights[j]);
            mu_real_mat_set(P->B, i, m + j,
                mu_real_mat_get(G_nom->B, i, j));
        }

    /* C: uncertainty outputs + measurements */
    for (int i = 0; i < pG; i++)
        for (int j = 0; j < nG; j++)
            mu_real_mat_set(P->C, pG + i, j,
                mu_real_mat_get(G_nom->C, i, j));

    /* D matrix */
    for (int i = 0; i < pG; i++)
        for (int j = 0; j < m; j++)
            mu_real_mat_set(P->D, pG + i, m + j,
                mu_real_mat_get(G_nom->D, i, j));
    for (int j = 0; j < m; j++) {
        mu_real_mat_set(P->D, j, j, 0.0);
        mu_real_mat_set(P->D, j, m + j, weights[j]);
    }

    return P;
}

/* ====================================================================
 * L7: Boeing 747 Lateral-Directional Uncertainty Model
 *
 * Based on the Boeing 747 lateral dynamics at a specific flight condition.
 * These are well-documented models used extensively in robust control
 * research (e.g., NASA TP-1990).
 *
 * State vector: x = [beta, p, r, phi]^T
 *   beta  = sideslip angle (rad)
 *   p     = roll rate (rad/s)
 *   r     = yaw rate (rad/s)
 *   phi   = roll angle (rad)
 *
 * Input: u = [delta_a, delta_r]^T
 *   delta_a = aileron deflection (rad)
 *   delta_r = rudder deflection (rad)
 *
 * Output: y = [p, r, phi, beta]^T
 *
 * Nominal A matrix (approximate, flight condition: Mach 0.8, 40000 ft):
 *   A = [ -0.0558   0.0    -0.9968   0.0415;
 *         -1.8100  -0.3650   0.0      0.0;
 *          0.0850   0.0    -0.1950   0.0;
 *          0.0      1.0     0.0460   0.0    ]
 *
 * Nominal B matrix:
 *   B = [  0.0073   0.0;
 *          1.0800   0.2950;
 *         -0.1550  -0.7800;
 *          0.0      0.0    ]
 *
 * Uncertainties:
 *   - Roll-rate damping derivative L_p: nominal -0.365, uncertain 20% (real scalar)
 *   - Yaw-rate aerodynamic derivative N_r: nominal -0.195, uncertain 25% (real scalar)
 *   - Unmodeled actuator dynamics: complex full block (dim 2)
 *
 * Reference: NASA TP-1990, Boeing Commercial Airplane Group
 * ==================================================================== */

MUSystem* mu_unc_boeing747_model(void) {
    MUSystem *sys = mu_system_create(4, 2, 4);
    if (!sys) return NULL;

    /* A matrix */
    double A_data[4][4] = {
        {-0.0558,  0.0,    -0.9968,  0.0415},
        {-1.8100, -0.3650,  0.0,      0.0   },
        { 0.0850,  0.0,    -0.1950,   0.0   },
        { 0.0,     1.0,     0.0460,   0.0   }
    };
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            mu_real_mat_set(sys->A, i, j, A_data[i][j]);

    /* B matrix */
    double B_data[4][2] = {
        { 0.0073,  0.0   },
        { 1.0800,  0.2950},
        {-0.1550, -0.7800},
        { 0.0,     0.0   }
    };
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 2; j++)
            mu_real_mat_set(sys->B, i, j, B_data[i][j]);

    /* C = I (output all states) */
    for (int i = 0; i < 4; i++)
        mu_real_mat_set(sys->C, i, i, 1.0);

    /* D = 0 */
    return sys;
}

MUUncertaintyStructure* mu_unc_boeing747_structure(void) {
    MUUncertaintyStructure *s = mu_unc_struct_create(3);
    if (!s) return NULL;

    /* Block 0: Roll-rate damping L_p uncertainty (real scalar) */
    mu_unc_struct_add_block(s, 0, MU_UNC_REAL_SCALAR, 1, 1.0);

    /* Block 1: Yaw-rate damping N_r uncertainty (real scalar) */
    mu_unc_struct_add_block(s, 1, MU_UNC_REAL_SCALAR, 1, 1.0);

    /* Block 2: Unmodeled actuator dynamics (complex full block, dim 2) */
    mu_unc_struct_add_block(s, 2, MU_UNC_FULL_BLOCK, 2, 1.0);

    s->total_dim = 4;  /* 1 + 1 + 2 = 4 */
    return s;
}

/* ====================================================================
 * L7: Tesla Electric Vehicle Lateral Uncertainty Model
 *
 * Electric vehicle lateral dynamics (bicycle model):
 *
 * State: x = [vy, yaw_rate]^T
 *   vy      = lateral velocity (m/s)
 *   yaw_rate = yaw rate (rad/s)
 *
 * Input: u = [steering_angle]^T
 *   steering_angle (rad)
 *
 * Output: y = [yaw_rate, lateral_acceleration]^T
 *
 * Nominal parameters:
 *   mass m = 1800 kg (10% uncertainty)
 *   velocity vx = 20 m/s
 *   front cornering stiffness Cf = 80000 N/rad (20% uncertainty)
 *   rear cornering stiffness Cr = 100000 N/rad (20% uncertainty)
 *   distance to front axle lf = 1.3 m
 *   distance to rear axle lr = 1.5 m
 *   yaw inertia Iz = 2500 kg*m^2
 *
 * A = [ -(Cf+Cr)/(m*vx),  -vx - (lf*Cf - lr*Cr)/(m*vx);
 *       -(lf*Cf - lr*Cr)/(Iz*vx),  -(lf^2*Cf + lr^2*Cr)/(Iz*vx) ]
 *
 * B = [ Cf/m; lf*Cf/Iz ]
 *
 * C = [ 0, 1; -(Cf+Cr)/m, -(lf*Cf - lr*Cr)/(m*vx) ]
 *
 * D = [ 0; Cf/m ]
 *
 * Uncertainties:
 *   - Cornering stiffness Cf and Cr (real scalar, 20% each)
 *   - Vehicle mass m (real scalar, 10%)
 *   - Steering actuator lag (complex full block, 1st-order unmodeled)
 * ==================================================================== */

MUSystem* mu_unc_tesla_vehicle_model(void) {
    MUSystem *sys = mu_system_create(2, 1, 2);
    if (!sys) return NULL;

    double m = 1800.0, vx = 20.0;
    double Cf = 80000.0, Cr = 100000.0;
    double lf = 1.3, lr = 1.5, Iz = 2500.0;

    /* A matrix */
    double a11 = -(Cf + Cr) / (m * vx);
    double a12 = -vx - (lf * Cf - lr * Cr) / (m * vx);
    double a21 = -(lf * Cf - lr * Cr) / (Iz * vx);
    double a22 = -(lf * lf * Cf + lr * lr * Cr) / (Iz * vx);

    mu_real_mat_set(sys->A, 0, 0, a11);
    mu_real_mat_set(sys->A, 0, 1, a12);
    mu_real_mat_set(sys->A, 1, 0, a21);
    mu_real_mat_set(sys->A, 1, 1, a22);

    /* B matrix */
    mu_real_mat_set(sys->B, 0, 0, Cf / m);
    mu_real_mat_set(sys->B, 1, 0, lf * Cf / Iz);

    /* C matrix */
    mu_real_mat_set(sys->C, 0, 0, 0.0);
    mu_real_mat_set(sys->C, 0, 1, 1.0);
    mu_real_mat_set(sys->C, 1, 0, -(Cf + Cr) / m);
    mu_real_mat_set(sys->C, 1, 1, -(lf * Cf - lr * Cr) / (m * vx));

    /* D matrix */
    mu_real_mat_set(sys->D, 0, 0, 0.0);
    mu_real_mat_set(sys->D, 1, 0, Cf / m);

    return sys;
}

MUUncertaintyStructure* mu_unc_tesla_vehicle_structure(void) {
    MUUncertaintyStructure *s = mu_unc_struct_create(4);
    if (!s) return NULL;

    /* Block 0: Front cornering stiffness Cf (real scalar, 20%) */
    mu_unc_struct_add_block(s, 0, MU_UNC_REAL_SCALAR, 1, 1.0);

    /* Block 1: Rear cornering stiffness Cr (real scalar, 20%) */
    mu_unc_struct_add_block(s, 1, MU_UNC_REAL_SCALAR, 1, 1.0);

    /* Block 2: Vehicle mass m (real scalar, 10%) */
    mu_unc_struct_add_block(s, 2, MU_UNC_REAL_SCALAR, 1, 1.0);

    /* Block 3: Steering actuator lag (complex full block, 1st-order) */
    mu_unc_struct_add_block(s, 3, MU_UNC_FULL_BLOCK, 1, 1.0);

    s->total_dim = 4;  /* 1 + 1 + 1 + 1 = 4 */
    return s;
}

/* ====================================================================
 * L7: SpaceX Falcon 9 Rocket Attitude Uncertainty Model
 *
 * Simplified pitch-axis dynamics for a launch vehicle:
 *
 * State: x = [theta, q, alpha]^T
 *   theta = pitch angle (rad)
 *   q     = pitch rate (rad/s)
 *   alpha = angle of attack (rad)
 *
 * Input: u = [delta_gimbal]^T
 *   delta_gimbal = engine gimbal angle (rad)
 *
 * Output: y = [theta, q]^T
 *
 * Nominal A matrix:
 *   A = [ 0,      1,        0     ;
 *         0,    a22(Z_alpha),  a23;
 *         0,      1,       a33(Z_alpha/V) ]
 *
 * Where aerodynamic coefficients Z_alpha and M_alpha vary with
 * Mach number and dynamic pressure.
 *
 * Uncertainties:
 *   - Aerodynamic coefficient Z_alpha (real scalar, 30% variation)
 *   - Center-of-mass location uncertainty (real scalar, affecting a22)
 *   - Unmodeled flexible bending modes (complex full block, dim 2)
 *   - Propellant slosh dynamics (complex full block, dim 2)
 *
 * Reference: Wie (2008), "Space Vehicle Dynamics and Control"
 *            SpaceX Falcon 9 User's Guide (2021)
 * ==================================================================== */

MUSystem* mu_unc_spacex_rocket_model(void) {
    MUSystem *sys = mu_system_create(3, 1, 2);
    if (!sys) return NULL;

    /* Nominal parameters at max-Q (approx) */
    double Z_alpha = -0.5;    /* Normal force derivative (1/s) */
    double M_alpha = -10.0;   /* Pitch moment derivative (1/s^2) */
    double V = 500.0;         /* Velocity (m/s) */
    double l_cg = 15.0;       /* Distance from gimbal to CG (m) */
    double Iyy = 5.0e6;       /* Pitch inertia (kg*m^2) */
    double T = 7.6e6;         /* Thrust (N) */

    double a22 = M_alpha;     /* Pitch acceleration due to alpha */
    double a23 = T * l_cg / Iyy;  /* Control moment arm effect */
    double a33 = Z_alpha / V; /* Alpha dynamics */

    /* A matrix */
    mu_real_mat_set(sys->A, 0, 0, 0.0);
    mu_real_mat_set(sys->A, 0, 1, 1.0);
    mu_real_mat_set(sys->A, 0, 2, 0.0);
    mu_real_mat_set(sys->A, 1, 0, 0.0);
    mu_real_mat_set(sys->A, 1, 1, a22 * 0.1);
    mu_real_mat_set(sys->A, 1, 2, a23);
    mu_real_mat_set(sys->A, 2, 0, 0.0);
    mu_real_mat_set(sys->A, 2, 1, 1.0);
    mu_real_mat_set(sys->A, 2, 2, a33);

    /* B matrix */
    mu_real_mat_set(sys->B, 0, 0, 0.0);
    mu_real_mat_set(sys->B, 1, 0, T / Iyy);
    mu_real_mat_set(sys->B, 2, 0, 0.0);

    /* C = [1 0 0; 0 1 0] */
    mu_real_mat_set(sys->C, 0, 0, 1.0);
    mu_real_mat_set(sys->C, 1, 1, 1.0);

    /* D = 0 */
    return sys;
}

MUUncertaintyStructure* mu_unc_spacex_rocket_structure(void) {
    MUUncertaintyStructure *s = mu_unc_struct_create(4);
    if (!s) return NULL;

    /* Block 0: Z_alpha aerodynamic coefficient (real scalar, 30%) */
    mu_unc_struct_add_block(s, 0, MU_UNC_REAL_SCALAR, 1, 1.0);

    /* Block 1: Center-of-mass location (real scalar) */
    mu_unc_struct_add_block(s, 1, MU_UNC_REAL_SCALAR, 1, 1.0);

    /* Block 2: Unmodeled flexible modes (complex full block, dim 2) */
    mu_unc_struct_add_block(s, 2, MU_UNC_FULL_BLOCK, 2, 1.0);

    /* Block 3: Propellant slosh (complex full block, dim 2) */
    mu_unc_struct_add_block(s, 3, MU_UNC_FULL_BLOCK, 2, 1.0);

    s->total_dim = 6;  /* 1 + 1 + 2 + 2 = 6 */
    return s;
}

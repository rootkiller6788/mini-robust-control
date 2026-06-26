#ifndef MU_UNCERTAINTY_H
#define MU_UNCERTAINTY_H

#include "mu_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ====================================================================
 * mu_uncertainty.h — Uncertainty Modeling for u-Synthesis
 *
 * Physical systems have various types of uncertainty:
 *
 * 1. Parametric uncertainty: physical parameters (mass, stiffness,
 *    damping) are known only within bounds. Modeled as REAL repeated
 *    scalar blocks in Delta.
 *
 * 2. Unmodeled dynamics: high-frequency dynamics neglected in the
 *    nominal model. Modeled as COMPLEX full blocks with frequency-
 *    dependent weighting functions.
 *
 * 3. LTI uncertainty: neglected linear dynamics (e.g., actuator lags).
 *    Modeled as COMPLEX full blocks.
 *
 * The uncertainty structure Delta is built by combining these blocks:
 *   Delta = diag(delta_1 I_{r1}, ..., delta_S I_{rS}, Delta_1, ..., Delta_F)
 *
 * Reference: Skogestad & Postlethwaite (2005), Ch. 7-8
 *            Zhou, Doyle & Glover (1996), Ch. 9
 * ==================================================================== */

/*
 * L1: Uncertainty Weight Template
 *
 * Uncertainty weighting functions shape the frequency-dependent
 * magnitude of uncertainty. Common forms:
 *
 *   W_u(s) = (s + omega_b * r_0) / (s * r_inf + omega_b)
 *
 * where r_0 = relative uncertainty at DC
 *       r_inf = relative uncertainty at high frequency
 *       omega_b = corner frequency
 *
 * This represents: low uncertainty at low freq (good model),
 * high uncertainty at high freq (unmodeled dynamics).
 */
typedef struct {
    double r0;         /* DC relative uncertainty */
    double rinf;       /* high-frequency relative uncertainty */
    double omega_b;    /* corner frequency (rad/s) */
    double gain;       /* overall scaling */
} MUUncertaintyWeight;

MUSystem* mu_unc_weight_create(const MUUncertaintyWeight *w);

/*
 * L3: Parametric Uncertainty Modeling
 *
 * For a parameter p = p_nom * (1 + r_p * delta) where |delta| <= 1:
 *
 * Extract the uncertain parameter from the plant using LFT:
 *   G(s, p) = F_u(P(s), Delta_p)
 *
 * This function constructs the LFT representation of a system
 * with uncertain parameters.
 *
 * Example: uncertain mass m = m0 * (1 + w_m * delta)
 * The parameter depends linearly on delta through a feedback structure.
 *
 * Complexity: O(n^3)
 */
typedef struct {
    double nominal;
    double relative_uncertainty;  /* r_p */
    int    state_index;           /* which state equation contains this param */
} MUParametricUncertainty;

/*
 * L3: Building the Uncertainty Structure
 *
 * Assemble the complete Delta structure from individual block descriptions.
 *
 * For a system with:
 *   - 2 uncertain masses (real repeated scalars, size 1 each)
 *   - 1 unmodeled actuator dynamics (complex full block, size 2)
 *
 * Delta = diag(delta_m1, delta_m2, Delta_act)
 *        where delta_m1, delta_m2 in R, |delta| <= 1
 *              Delta_act in C^{2x2}, ||Delta_act||_inf <= 1
 */
MUUncertaintyStructure* mu_unc_build_structure(
    const MUUncertaintyBlock *blocks, int num_blocks);

/*
 * L3: Input/Output Multiplicative Uncertainty
 *
 * Models uncertainty at the plant input or output:
 *
 * Input multiplicative:  G_true = G_nom * (I + W_I * Delta_I)
 * Output multiplicative: G_true = (I + W_O * Delta_O) * G_nom
 *
 * where W_I, W_O are weighting functions and |Delta|_inf <= 1.
 *
 * The LFT representation pulls out Delta:
 *   [ z ]   [ 0       W_I * G_nom ] [ w ]
 *   [ y ] = [ G_nom   W_I * G_nom ] [ u ]   (input mult)
 *
 * Complexity: O(n^3) for LFT construction
 */
MUSystem* mu_unc_input_multiplicative(const MUSystem *G_nom,
                                       const MUSystem *W_I);
MUSystem* mu_unc_output_multiplicative(const MUSystem *G_nom,
                                        const MUSystem *W_O);

/*
 * L3: Additive Uncertainty
 *
 * G_true = G_nom + W_A * Delta_A
 *
 * where ||Delta_A||_inf <= 1. The LFT form:
 *
 *   [ z ]   [ 0      W_A ] [ w ]
 *   [ y ] = [ I      G_nom ] [ u ]   (additive)
 */
MUSystem* mu_unc_additive(const MUSystem *G_nom, const MUSystem *W_A);

/*
 * L3: Feedback Uncertainty (Coprime Factor)
 *
 * G_true = (M + Delta_M)^{-1} (N + Delta_N)
 *
 * where ||[Delta_N, Delta_M]||_inf <= epsilon.
 * This is the gap metric / coprime factor uncertainty.
 *
 * Reference: McFarlane & Glover (1990), "Robust controller design
 *            using normalized coprime factor plant descriptions"
 */
MUSystem* mu_unc_coprime_factor(const MUSystem *G_nom, double epsilon);

/*
 * L3: Diagonal Input Uncertainty
 *
 * For systems with multiple independent actuators, each with
 * independent uncertainty:
 *
 *   Delta_I = diag(delta_1, delta_2, ..., delta_m)
 *
 * Each delta_i is a SISO uncertainty (scalar complex or real).
 *
 * With weighting W_I = diag(w_1, w_2, ..., w_m), the LFT form:
 *
 *   [ z ]   [ 0      W_I ] [ w ]
 *   [ y ] = [ G_nom  G_nom * W_I ] [ u ]
 */
MUSystem* mu_unc_diagonal_input(const MUSystem *G_nom,
                                 const double *weights, int m);

/*
 * L7: Boeing 747 Lateral-Directional Uncertainty Model
 *
 * Real-world model from Boeing 747 lateral dynamics with:
 *   - Uncertain roll-rate damping coefficient (real scalar)
 *   - Uncertain yaw-rate aerodynamic derivative (real scalar)
 *   - Unmodeled actuator dynamics (complex full block, dim 2)
 *
 * State: [beta, p, r, phi]  (sideslip, roll rate, yaw rate, roll angle)
 * Input: [delta_a, delta_r]  (aileron, rudder)
 * Output: [p, r, phi, beta]
 *
 * Reference: "Robust Control of the Boeing 747", NASA TP-1990
 *            Boeing Commercial Airplane Group, lateral-directional model
 */
MUSystem* mu_unc_boeing747_model(void);
MUUncertaintyStructure* mu_unc_boeing747_structure(void);

/*
 * L7: Tesla Vehicle Lateral Uncertainty Model
 *
 * Electric vehicle lateral dynamics with:
 *   - Uncertain cornering stiffness (real scalar, 20% variation)
 *   - Uncertain vehicle mass (real scalar, 10% variation)
 *   - Unmodeled steering actuator lag (complex full block)
 *
 * State: [vy, yaw_rate]  (lateral velocity, yaw rate)
 * Input: [steering_angle]
 * Output: [yaw_rate, lateral_acceleration]
 */
MUSystem* mu_unc_tesla_vehicle_model(void);
MUUncertaintyStructure* mu_unc_tesla_vehicle_structure(void);

/*
 * L7: SpaceX Rocket Attitude Uncertainty Model
 *
 * Launch vehicle attitude dynamics with:
 *   - Uncertain aerodynamic coefficients (real scalars)
 *   - Uncertain center-of-mass location (real scalar)
 *   - Unmodeled flexible modes (complex full block)
 *   - Slosh dynamics uncertainty (complex full block)
 *
 * Reference: "Launch Vehicle Attitude Control", SpaceX Falcon 9
 *            Wie (2008), "Space Vehicle Dynamics and Control"
 */
MUSystem* mu_unc_spacex_rocket_model(void);
MUUncertaintyStructure* mu_unc_spacex_rocket_structure(void);

#ifdef __cplusplus
}
#endif

#endif /* MU_UNCERTAINTY_H */

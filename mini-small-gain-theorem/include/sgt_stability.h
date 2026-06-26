#ifndef SGT_STABILITY_H
#define SGT_STABILITY_H

#include "sgt_core.h"

/* ============================================================================
 * Small Gain Theorem — Stability Theory
 *
 * This header covers the core stability concepts related to the small gain
 * theorem: input-output stability, finite-gain stability, interconnection
 * analysis, and the statement of the theorem itself.
 *
 * The small gain theorem connects the open-loop gains of subsystems to
 * the closed-loop stability of their feedback interconnection. It is
 * fundamentally an input-output (operator) theory result that does not
 * require state-space models — only that each subsystem defines a
 * finite-gain stable operator on L2.
 *
 * Theorem (Small Gain, Zames 1966):
 *   Let H1, H2 : L2e -> L2e be causal operators with finite L2 gains
 *   gamma1 = ||H1|| and gamma2 = ||H2||. If gamma1 * gamma2 < 1, then
 *   the feedback interconnection of H1 and H2 is finite-gain L2 stable.
 *
 * For LTI systems: gamma_i = ||Hi||_inf (H-infinity norm).
 * ============================================================================ */

/* --- Interconnection Structure --- */

/**
 * Feedback interconnection of two systems.
 *
 * Standard configuration (positive feedback on u1):
 *   u1 = w1 + y2    (w1 is external input to channel 1)
 *   u2 = w2 + y1    (w2 is external input to channel 2)
 *   y1 = H1(u1)     (forward system output)
 *   y2 = H2(u2)     (feedback system output)
 *
 * For negative feedback, negate the output of H2 before interconnection.
 */
typedef struct {
    SGTLTISystem *H1;
    SGTLTISystem *H2;
    SGTInterconnectionType topology;
    char *name;

    double gamma1;
    double gamma2;
    double gamma_product;
    SGTVerifyResult verification;

    SGTLTISystem *closed_loop;
    SGTStabilityStatus cl_stability;
    bool closed_loop_computed;
} SGTInterconnection;

/**
 * Input-output signal pair for time-domain simulation.
 */
typedef struct {
    double *time;
    double *input;
    double *output;
    int length;
    double dt;
} SGTSignal;

/**
 * Bounded-input bounded-output (BIBO) test data.
 * Stores the result of injecting a bounded test signal
 * and measuring the output-to-input energy ratio.
 */
typedef struct {
    double input_energy;
    double output_energy;
    double gain_estimate;
    bool is_bounded;
} SGTBIBOResult;

/* --- Input-Output Stability --- */

/**
 * Finite-gain Lp stability definition.
 * A system H is finite-gain Lp stable if there exists gamma >= 0
 * and beta >= 0 such that for all u in Lp and all T >= 0:
 *   ||(H(u))_T||_p <= gamma * ||u_T||_p + beta
 * The smallest such gamma is the induced Lp gain of H.
 */
typedef struct {
    double gain;
    double bias;
    SGTNormType norm_type;
    bool is_finite_gain_stable;
} SGTFiniteGainStability;

/* --- Stability Verification API --- */

SGTInterconnection* sgt_interconnection_create(const char *name,
                                                SGTLTISystem *H1,
                                                SGTLTISystem *H2,
                                                SGTInterconnectionType topology);
void sgt_interconnection_free(SGTInterconnection *ic);
void sgt_interconnection_verify(SGTInterconnection *ic, int n_freq_points);
bool sgt_scaled_small_gain(SGTInterconnection *ic, int max_d_iter);
void sgt_build_closed_loop(SGTInterconnection *ic);
SGTStabilityStatus sgt_check_eigenvalue_stability(const SGTMatrix *A);
SGTBIBOResult sgt_time_domain_verify(const SGTLTISystem *sys,
                                      const SGTVector *x0,
                                      double duration, double dt);

/* --- Passivity and Small Gain --- */

double sgt_passivity_index(const SGTLTISystem *sys, int n_freq_points);
SGTLTISystem* sgt_scattering_transform(const SGTLTISystem *H);

/* --- Stability Margins --- */

double sgt_stability_margin(const SGTInterconnection *ic,
                             const SGTStructuredUncertainty *delta_struct);
double sgt_gap_metric(const SGTLTISystem *P1, const SGTLTISystem *P2);
double sgt_nu_gap_metric(const SGTLTISystem *P1, const SGTLTISystem *P2);

#endif /* SGT_STABILITY_H */

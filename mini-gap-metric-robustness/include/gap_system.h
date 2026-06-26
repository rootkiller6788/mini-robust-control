/* gap_system.h ? State-Space System Representation for Robust Control
 *
 * Defines linear time-invariant (LTI) systems in state-space form,
 * transfer function evaluation, system interconnections, and the
 * graph symbol representation central to gap metric theory.
 *
 * Knowledge coverage:
 *   L1: State-space system, transfer function, graph symbol (definitions)
 *   L2: Controllability, observability, stability, minimal realization
 *   L3: State-space (A,B,C,D), McMillan degree, frequency response
 *   L5: System interconnection (series, parallel, feedback)
 *
 * Ref: K. Zhou, J.C. Doyle, K. Glover, "Robust and Optimal Control" (1996)
 *      T. Kailath, "Linear Systems" (1980)
 *      D.C. McFarlane & K. Glover, "Robust Controller Design Using
 *      Normalized Coprime Factor Plant Descriptions" (1990)
 */

#ifndef GAP_SYSTEM_H
#define GAP_SYSTEM_H

#include "gap_core.h"

/* ================================================================
 * L1: Core Definitions - LTI System and Graph Symbol
 * ================================================================ */

/** GapSystem - Linear Time-Invariant system in state-space form.
 *
 * Continuous-time: dx/dt = A x + B u,  y = C x + D u
 * Discrete-time:   x[k+1] = A x[k] + B u[k],  y[k] = C x[k] + D u[k]
 *
 * Transfer function: G(s) = C (sI - A)^{-1} B + D
 */
typedef struct {
    GapMatrix *A;
    GapMatrix *B;
    GapMatrix *C;
    GapMatrix *D;
    int        n;
    int        m;
    int        p;
    bool       is_discrete;
    double     ts;
} GapSystem;

/** GapFreqResp - Frequency response G(jw) or G(e^{jw}) at a point. */
typedef struct {
    int          p;
    int          m;
    double       omega;
    GapComplex  *resp;
} GapFreqResp;

/** GapGraphSymbol - Normalized coprime factor graph symbol.
 *
 * For plant G with NRCF (N,M): G = N M^{-1}, N*N + M*M = I.
 * Graph symbol G_sym = [M; N] encodes closed-loop geometry.
 * The gap metric is defined via this graph symbol.
 *
 * Ref: Georgiou & Smith, "Optimal robustness in the gap metric" (1990)
 */
typedef struct {
    GapSystem *M_sys;
    GapSystem *N_sys;
    bool       is_normalized;
    int        n_g;
} GapGraphSymbol;

/** GapFeedbackLoop - Standard feedback interconnection.
 *  u1 = r1 - y2, u2 = r2 + y1; y1 = G1 u1, y2 = G2 u2
 */
typedef struct {
    GapSystem *plant;
    GapSystem *controller;
    bool       is_stable;
    double     stab_margin;
} GapFeedbackLoop;

/* ================================================================
 * L2: Core Concepts - Controllability, Observability, Stability
 * ================================================================ */

/** Controllability matrix: Ctr = [B, AB, A^2 B, ..., A^{n-1} B].
 *  Ref: R.E. Kalman (1960) */
GapMatrix* gap_sys_controllability_matrix(const GapSystem *sys);

/** Observability matrix: Obs = [C; CA; CA^2; ...; CA^{n-1}]. */
GapMatrix* gap_sys_observability_matrix(const GapSystem *sys);

/** Check if (A,B) is controllable via Gramian/PBH test. */
bool gap_sys_is_controllable(const GapSystem *sys);

/** Check if (A,C) is observable. */
bool gap_sys_is_observable(const GapSystem *sys);

/** Check minimality (controllable + observable). */
bool gap_sys_is_minimal(const GapSystem *sys);

/** Check continuous-time: Re(eig(A)) < 0; discrete: |eig(A)| < 1. */
bool gap_sys_is_stable(const GapSystem *sys);

/** Stability margin = max{ alpha : A + alpha*I has eig with Re=0 }. */
double gap_sys_stability_margin(const GapSystem *sys);

/* ================================================================
 * L3: Mathematical Structures - System Operations
 * ================================================================ */

/* ---- System Lifecycle ---- */
GapSystem* gap_sys_create(int n, int m, int p, bool is_discrete, double ts);
GapSystem* gap_sys_create_from(GapMatrix *A, GapMatrix *B,
                                GapMatrix *C, GapMatrix *D,
                                bool is_discrete, double ts);
GapSystem* gap_sys_clone(const GapSystem *sys);
void gap_sys_free(GapSystem *sys);
void gap_sys_print(const GapSystem *sys, const char *name);

/* ---- Frequency Response ---- */
GapFreqResp* gap_sys_freqresp(const GapSystem *sys, double omega);
GapFreqResp* gap_sys_freqresp_grid(const GapSystem *sys,
                                    const double *omega, int n_omega);
void gap_freqresp_free(GapFreqResp *fr);

/** H-infinity norm via bisection on Hamiltonian eigenvalues.
 *  ||G||inf = sup_omega sigma_max(G(j omega)). */
double gap_sys_hinf_norm(const GapSystem *sys);

/** H2 norm via controllability Gramian.
 *  ||G||_2^2 = trace(C P C^T), A P + P A^T + B B^T = 0. */
double gap_sys_h2_norm(const GapSystem *sys);

/* ---- Gramians ---- */
GapMatrix* gap_sys_controllability_gramian(const GapSystem *sys);
GapMatrix* gap_sys_observability_gramian(const GapSystem *sys);

/** Hankel singular values: sqrt(eigenvalues(P Q)). */
double* gap_sys_hankel_singular_values(const GapSystem *sys, int *n_sv);

/* ---- System Interconnection ---- */
GapSystem* gap_sys_series(const GapSystem *G1, const GapSystem *G2);
GapSystem* gap_sys_parallel(const GapSystem *G1, const GapSystem *G2);
GapSystem* gap_sys_lft(const GapSystem *G, const GapSystem *K);
GapSystem* gap_sys_generalized_plant(const GapSystem *G,
                                      const GapMatrix *W1,
                                      const GapMatrix *W2);

/* ---- Pole Placement ---- */
/** State feedback u = -Kx for pole placement (Ackermann). */
GapMatrix* gap_sys_acker(const GapSystem *sys, const double *poles, int np);

/** Observer gain L for pole placement (dual Ackermann). */
GapMatrix* gap_sys_observer_gain(const GapSystem *sys,
                                  const double *poles, int np);

/* ---- Lyapunov Equation ---- */
/** Solve A X + X A^T + Q = 0 via Bartels-Stewart algorithm. */
GapMatrix* gap_sys_lyapunov(const GapMatrix *A, const GapMatrix *Q);

#endif /* GAP_SYSTEM_H */

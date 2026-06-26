#ifndef SSV_SYNTHESIS_H
#define SSV_SYNTHESIS_H

#include "ssv_core.h"
#include "ssv_uncertainty.h"
#include "ssv_bounds.h"
#include <complex.h>
#include <stdbool.h>

/* ============================================================================
 * mu-Synthesis (D-K Iteration)
 *
 * mu-synthesis designs a controller K that minimizes the structured singular
 * value mu of the closed-loop system, thereby achieving robust performance.
 *
 * The standard approach is D-K iteration (Doyle, 1985):
 *   1. Fix D-scaling -> Solve H-infinity problem for K   (K-step)
 *   2. Fix K           -> Optimize D-scaling at each freq  (D-step)
 *   3. Iterate until mu < 1 or convergence
 *
 * While D-K iteration does not guarantee convergence to the global optimum,
 * it works well in practice for many engineering problems.
 *
 * References:
 *   Doyle, "Structured uncertainty in control system design" (CDC, 1985)
 *   Balas, Doyle, Glover, Packard, Smith — mu-Analysis and Synthesis Toolbox (1998)
 *   Zhou, Doyle & Glover — Robust and Optimal Control, Ch. 17 (1996)
 * ============================================================================ */

/* --- State-Space System Representation --- */

/** Continuous-time state-space system:
 *    dx/dt = A * x + B * u
 *    y     = C * x + D * u
 */
typedef struct {
    int n_states;          /* number of states (order of A) */
    int n_inputs;          /* number of inputs */
    int n_outputs;         /* number of outputs */
    SSVRealMatrix *A;      /* system matrix, n x n */
    SSVRealMatrix *B;      /* input matrix, n x m */
    SSVRealMatrix *C;      /* output matrix, p x n */
    SSVRealMatrix *D;      /* feedthrough matrix, p x m */
} SSVStateSpace;

/** Create a state-space system.
 *  @param n_states state dimension
 *  @param n_inputs input dimension
 *  @param n_outputs output dimension
 */
SSVStateSpace* ssv_ss_create(int n_states, int n_inputs, int n_outputs);

/** Free state-space system. */
void ssv_ss_free(SSVStateSpace *sys);

/** Compute frequency response: G(jw) = C * (jw*I - A)^{-1} * B + D
 *  @param sys state-space system
 *  @param w frequency in rad/s
 *  @return complex transfer matrix at given frequency
 */
SSVComplexMatrix* ssv_ss_freqresp(const SSVStateSpace *sys, double w);

/** Compute frequency response at multiple frequencies. */
SSVComplexMatrix** ssv_ss_freqresp_grid(const SSVStateSpace *sys,
                                          const double *w, int n_freq);

/** Print state-space system matrices. */
void ssv_ss_print(const SSVStateSpace *sys);

/* --- Closed-Loop Interconnection --- */

/** Generalized plant P (for H∞ / μ synthesis):
 *
 *   [ z ]   [ P11  P12 ] [ w ]
 *   [ y ] = [ P21  P22 ] [ u ]
 *
 *   where w = exogenous inputs (disturbances, references, noise)
 *         z = regulated outputs (errors, weighted outputs)
 *         u = control inputs
 *         y = measured outputs
 *
 *   Closed-loop with K: z = Fl(P, K) * w
 *   where Fl(P,K) = P11 + P12 * K * (I - P22*K)^{-1} * P21
 */
typedef struct {
    int n_w;               /* number of exogenous inputs */
    int n_u;               /* number of control inputs */
    int n_z;               /* number of regulated outputs */
    int n_y;               /* number of measured outputs */
    int n_states;           /* number of states of P */
    SSVRealMatrix *A;      /* n x n */
    SSVRealMatrix *B1;     /* n x n_w (exogenous input matrix) */
    SSVRealMatrix *B2;     /* n x n_u (control input matrix) */
    SSVRealMatrix *C1;     /* n_z x n (regulated output matrix) */
    SSVRealMatrix *C2;     /* n_y x n (measured output matrix) */
    SSVRealMatrix *D11;    /* n_z x n_w */
    SSVRealMatrix *D12;    /* n_z x n_u */
    SSVRealMatrix *D21;    /* n_y x n_w */
    SSVRealMatrix *D22;    /* n_y x n_u */
} SSVGeneralizedPlant;

/** Create a generalized plant. */
SSVGeneralizedPlant* ssv_gplant_create(int n_states, int n_w, int n_u, int n_z, int n_y);

/** Free generalized plant. */
void ssv_gplant_free(SSVGeneralizedPlant *P);

/** Set a sub-matrix of the generalized plant by name.
 *  submatrix_name: "A", "B1", "B2", "C1", "C2", "D11", "D12", "D21", "D22"
 */
void ssv_gplant_set_matrix(SSVGeneralizedPlant *P, const char *submatrix_name,
                            const SSVRealMatrix *M);

/** Form closed-loop system: Fl(P, K) = lower LFT of P and K.
 *  @return state-space of closed-loop system
 */
SSVStateSpace* ssv_closed_loop(const SSVGeneralizedPlant *P, const SSVStateSpace *K);

/** Form the mu-analysis interconnection: augment M with uncertainty channels.
 *  Given P and K, extract M = Fl(P, K) for mu-analysis.
 *
 *  The uncertainty channels are separated from performance channels:
 *    M = [ M11  M12 ]
 *        [ M21  M22 ]
 *  where M11 connects to Delta (uncertainty), M22 connects to performance.
 */
SSVComplexMatrix* ssv_form_m_interconnection(const SSVGeneralizedPlant *P,
                                               const SSVStateSpace *K,
                                               double w);

/* --- D-K Iteration --- */

/** D-K iteration configuration. */
typedef struct {
    int max_iterations;          /* maximum D-K iterations (default 10) */
    double mu_target;            /* target mu value (default 1.0) */
    double tol_dk;               /* D-K convergence tolerance (default 0.01) */
    SSVDScaleOptions dscale_opts; /* D-optimization options */
    double *freq_grid;           /* frequency grid points (rad/s) */
    int n_freq;                  /* number of frequency points */
    bool verbose;                /* print iteration progress */
} SSVDKIterOptions;

/** Default D-K iteration options. */
SSVDKIterOptions ssv_dk_options_default(void);

/** Result of one D-K iteration. */
typedef struct {
    int iteration;                      /* iteration number (0-based) */
    double mu_peak;                     /* peak mu value across frequency */
    double mu_peak_freq;                /* frequency of peak mu */
    SSVStateSpace *controller;          /* current controller K */
    double *mu_curve;                   /* mu values across frequency grid */
    int n_freq;                         /* number of frequency points */
    double hinf_norm;                   /* H-infinity norm achieved in K-step */
} SSVDKIterResult;

/** Free D-K iteration result. */
void ssv_dk_result_free(SSVDKIterResult *result);

/** Perform one step of D-K iteration:
 *   K-step: Design H-infinity controller for scaled plant
 *   D-step: Optimize D-scaling at each frequency
 *
 *  Simplified H-infinity design uses a Riccati-based approach
 *  (Doyle-Glover-Khargonekar-Francis 2-Riccati solution, simplified).
 *
 *  @param P generalized plant
 *  @param ustruct uncertainty structure
 *  @param K_prev previous controller (NULL for initial), shape preserved
 *  @param options D-K iteration options
 *  @param iteration iteration number
 *  @return iteration result
 */
SSVDKIterResult* ssv_dk_step(const SSVGeneralizedPlant *P,
                               const SSVUncertaintyStructure *ustruct,
                               const SSVStateSpace *K_prev,
                               const SSVDKIterOptions *options,
                               int iteration);

/** Run full D-K synthesis.
 *  Iterates the D-K steps until convergence or max_iterations.
 *
 *  @param P generalized plant
 *  @param ustruct uncertainty structure
 *  @param options D-K iteration options
 *  @param n_results output: number of iterations performed
 *  @return array of iteration results (length n_results)
 */
SSVDKIterResult** ssv_dk_synthesize(const SSVGeneralizedPlant *P,
                                      const SSVUncertaintyStructure *ustruct,
                                      const SSVDKIterOptions *options,
                                      int *n_results);

/** Extract the controller from a D-K result. */
SSVStateSpace* ssv_dk_get_controller(const SSVDKIterResult *result);

/* --- H-infinity Norm Computation --- */

/** Compute the H-infinity norm of a stable continuous-time system
 *  via the Hamiltonian matrix and its eigenvalues.
 *
 *  H∞ norm: ||G||_∞ = sup_w sigma_max(G(jw))
 *
 *  Algorithm via bisection on gamma:
 *    Check if Hamiltonian H_gamma has eigenvalues on imaginary axis.
 *    H_gamma = [ A - B*R^{-1}*D^T*C,  -gamma*B*R^{-1}*B^T;
 *                gamma*C^T*S^{-1}*C,   -(A - B*R^{-1}*D^T*C)^T ]
 *    where R = D^T*D - gamma^2*I, S = D*D^T - gamma^2*I
 *
 *  @param sys stable state-space system
 *  @param tol bisection tolerance
 *  @return H-infinity norm
 */
double ssv_hinf_norm(const SSVStateSpace *sys, double tol);

/** Compute the H-infinity norm via frequency grid (brute-force check).
 *  Used as a fallback or for verification.
 */
double ssv_hinf_norm_grid(const SSVStateSpace *sys,
                           const double *freq_grid, int n_freq);

#endif /* SSV_SYNTHESIS_H */

#ifndef SGT_CORE_H
#define SGT_CORE_H

#include <stdbool.h>
#include <stddef.h>

/* ============================================================================
 * Small Gain Theorem — Core Definitions
 *
 * The Small Gain Theorem (Zames, 1966) is a cornerstone of robust control
 * theory. It provides a sufficient condition for the stability of a feedback
 * interconnection of two systems: if the product of their induced L2 gains
 * is strictly less than unity, the closed-loop system is stable.
 *
 * Reference works:
 *   - G. Zames, "On the input-output stability of time-varying nonlinear
 *     feedback systems, Parts I and II," IEEE Trans. AC, 1966.
 *   - C.A. Desoer and M. Vidyasagar, "Feedback Systems: Input-Output
 *     Properties," Academic Press, 1975 (reprinted SIAM, 2009).
 *   - J.C. Doyle, B.A. Francis, A.R. Tannenbaum, "Feedback Control Theory,"
 *     Macmillan, 1992.
 *   - K. Zhou, J.C. Doyle, K. Glover, "Robust and Optimal Control,"
 *     Prentice Hall, 1996.
 *   - A. van der Schaft, "L2-Gain and Passivity Techniques in Nonlinear
 *     Control," Springer, 2000 (3rd ed. 2017).
 * ============================================================================ */

/* --- Fundamental Enumerations --- */

/** System stability status classification. */
typedef enum {
    SGT_STABLE            = 0,
    SGT_MARGINALLY_STABLE = 1,
    SGT_UNSTABLE          = 2,
    SGT_UNDETERMINED      = 3
} SGTStabilityStatus;

/** Types of system norms used in gain computation. */
typedef enum {
    SGT_L2_GAIN     = 0,
    SGT_LINF_GAIN   = 1,
    SGT_L1_GAIN     = 2,
    SGT_HINF_NORM   = 3,
    SGT_H2_NORM     = 4,
    SGT_HANKEL_NORM = 5
} SGTNormType;

/** Interconnection topology types. */
typedef enum {
    SGT_SIMPLE_FEEDBACK = 0,
    SGT_CASCADE         = 1,
    SGT_PARALLEL        = 2,
    SGT_MULTI_LOOP      = 3,
    SGT_MIMO_FEEDBACK   = 4,
    SGT_NETWORKED       = 5
} SGTInterconnectionType;

/** Uncertainty model types. */
typedef enum {
    SGT_UNSTRUCTURED            = 0,
    SGT_ADDITIVE                = 1,
    SGT_MULTIPLICATIVE_INPUT    = 2,
    SGT_MULTIPLICATIVE_OUTPUT   = 3,
    SGT_COPRIME_FACTOR          = 4,
    SGT_PARAMETRIC              = 5,
    SGT_MIXED                   = 6
} SGTUncertaintyType;

/** Small gain verification result detail. */
typedef enum {
    SGT_PASS_FIRM         = 0,
    SGT_PASS_ADEQUATE     = 1,
    SGT_PASS_TIGHT        = 2,
    SGT_FAIL_BARELY       = 3,
    SGT_FAIL_SIGNIFICANT  = 4,
    SGT_VERIFY_INFEASIBLE = 5
} SGTVerifyResult;

/* --- Numeric Types --- */

/** Complex number for transfer function evaluation. */
typedef struct {
    double re;
    double im;
} SGTComplex;

/** Real matrix stored in row-major order. */
typedef struct {
    double *data;
    int rows;
    int cols;
    int stride;
} SGTMatrix;

/** Real vector. */
typedef struct {
    double *elements;
    int dim;
} SGTVector;

/** Frequency grid for H∞ norm evaluation. */
typedef struct {
    double *omega;
    double *magnitude;
    int n_points;
    double omega_min;
    double omega_max;
} SGTFrequencyGrid;

/* --- State-Space System Model --- */

/**
 * Continuous-time LTI state-space model:
 *     ẋ = A x + B u
 *     y = C x + D u
 */
typedef struct {
    char *name;
    int nx, nu, ny;
    double *A, *B, *C, *D;
    SGTStabilityStatus stability;
    double hinf_norm;
    double h2_norm;
    bool norms_computed;
} SGTLTISystem;

/* --- Transfer Function Model --- */

/** Rational transfer function: G(s) = num(s)/den(s).
 *  Polynomials in descending power order. */
typedef struct {
    double *num_coeffs;
    double *den_coeffs;
    int num_order;
    int den_order;
    double dc_gain;
} SGTTransferFunction;

/* --- Uncertainty Model --- */

/** Uncertainty block: norm-bounded operator Δ with ||Δ||∞ ≤ bound. */
typedef struct {
    SGTUncertaintyType type;
    int input_dim;
    int output_dim;
    double norm_bound;
    double *data;
    SGTLTISystem *dynamics;
    bool is_active;
} SGTUncertaintyBlock;

/** Structured uncertainty for μ-analysis. */
typedef struct {
    SGTUncertaintyBlock *blocks;
    int n_blocks;
    int total_dim;
    double mu_lower;
    double mu_upper;
} SGTStructuredUncertainty;

/* --- Core API: System Construction --- */

SGTLTISystem* sgt_lti_create(const char *name, int nx, int nu, int ny);
void sgt_lti_free(SGTLTISystem *sys);
void sgt_lti_set_A(SGTLTISystem *sys, const double *A_data);
void sgt_lti_set_B(SGTLTISystem *sys, const double *B_data);
void sgt_lti_set_C(SGTLTISystem *sys, const double *C_data);
void sgt_lti_set_D(SGTLTISystem *sys, const double *D_data);
SGTLTISystem* sgt_lti_first_order(const char *name, double gain, double tau);
SGTLTISystem* sgt_lti_second_order(const char *name, double gain,
                                    double wn, double zeta);
void sgt_lti_print(const SGTLTISystem *sys);

/* --- Core API: Matrix Operations --- */

SGTMatrix* sgt_matrix_create(int rows, int cols);
void sgt_matrix_free(SGTMatrix *mat);
double sgt_matrix_get(const SGTMatrix *mat, int i, int j);
void sgt_matrix_set(SGTMatrix *mat, int i, int j, double value);
void sgt_matrix_multiply(const SGTMatrix *A, const SGTMatrix *B, SGTMatrix *C);
void sgt_matrix_transpose(const SGTMatrix *A, SGTMatrix *B);
bool sgt_matrix_eig2(const SGTMatrix *A,
                     double *re1, double *im1,
                     double *re2, double *im2);
double sgt_matrix_sv_max(const SGTMatrix *A, int max_iters, double tol);
double sgt_matrix_frobenius_norm(const SGTMatrix *A);
SGTMatrix* sgt_lyapunov_2x2(const SGTMatrix *A, const SGTMatrix *Q);

/* --- Core API: Vector Operations --- */

SGTVector* sgt_vector_create(int dim);
void sgt_vector_free(SGTVector *v);
double sgt_vector_norm(const SGTVector *v);
double sgt_vector_dot(const SGTVector *a, const SGTVector *b);

/* --- Core API: Frequency Grid --- */

SGTFrequencyGrid* sgt_freq_grid_create_log(double omega_min, double omega_max,
                                            int n_points);
void sgt_freq_grid_free(SGTFrequencyGrid *grid);
void sgt_freq_evaluate(const SGTLTISystem *sys, SGTFrequencyGrid *grid);

/* --- Core API: Transfer Function --- */

SGTTransferFunction* sgt_tf_create(const double *num, int num_ord,
                                    const double *den, int den_ord);
void sgt_tf_free(SGTTransferFunction *tf);
SGTComplex sgt_tf_evaluate(const SGTTransferFunction *tf,
                            double sigma, double omega);
SGTLTISystem* sgt_tf_to_ss(const SGTTransferFunction *tf);

/* --- Utility Functions --- */

SGTVerifyResult sgt_check_small_gain(double gamma1, double gamma2);
void sgt_feedback_interconnect(const SGTLTISystem *G1, const SGTLTISystem *G2,
                                SGTMatrix *cl_A);
double sgt_rayleigh_iteration(const SGTMatrix *A, double shift, int max_iter);

#endif /* SGT_CORE_H */

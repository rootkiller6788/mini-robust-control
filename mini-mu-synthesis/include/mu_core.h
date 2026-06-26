#ifndef MU_CORE_H
#define MU_CORE_H

#include <stdbool.h>
#include <stddef.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ====================================================================
 * u-Synthesis Core Definitions
 * Based on: Doyle (1982), Packard & Doyle (1993), Zhou, Doyle & Glover (1996)
 *
 * Structured singular value u: the fundamental tool for analyzing
 * robustness of MIMO systems with structured uncertainty.
 *
 *   u_Delta(M) = 1 / min{ sigma_bar(Delta) : Delta in Delta, det(I - M*Delta) = 0 }
 *
 * unless no Delta makes I - M*Delta singular, in which case u_Delta(M) = 0.
 *
 * Key property (Small-u Theorem):
 *   The closed-loop system F_u(M, Delta) is robustly stable for all
 *   Delta in Delta with ||Delta||_inf <= 1 iff u_Delta(M) < 1 for all omega.
 *
 * Companion source: mu_core.c
 * ==================================================================== */

#define MU_EPSILON      1e-12
#define MU_MAX_DIM      128
#define MU_MAX_ITER     200
#define MU_PI           3.14159265358979323846
#define MU_TWO_PI       6.28318530717958647692
#define MU_TOL           1e-8
#define MU_POWER_TOL     1e-6

/* --- L1: Core Data Types --- */

typedef struct {
    double re;
    double im;
} MUComplex;

typedef struct {
    MUComplex *data;
    int        rows;
    int        cols;
    int        capacity;
} MUMatrix;

typedef struct {
    double *data;
    int     rows;
    int     cols;
    int     capacity;
} MURealMatrix;

/* --- L1: Uncertainty Block Structure --- */
/*
 * Delta = { diag(delta_1 I_{r1}, ..., delta_S I_{rS}, Delta_1, ..., Delta_F) }
 *
 * where delta_i in C (repeated scalar blocks) and Delta_j in C^{mj x mj} (full blocks).
 *
 * Reference: Packard & Doyle (1993), "The complex structured singular value"
 */

typedef enum {
    MU_UNC_REPEATED_SCALAR = 0,
    MU_UNC_FULL_BLOCK      = 1,
    MU_UNC_REAL_SCALAR     = 2,
} MUUncBlockType;

typedef struct {
    MUUncBlockType  type;
    int             size;
    double          bound;
    bool            is_real;
} MUUncertaintyBlock;

typedef struct {
    MUUncertaintyBlock *blocks;
    int                 num_blocks;
    int                 total_dim;
} MUUncertaintyStructure;

/* --- L1: State-Space System --- */
/*
 * Continuous:  dx/dt = A x + B u,  y = C x + D u
 * Discrete:    x[k+1] = A x[k] + B u[k],  y[k] = C x[k] + D u[k]
 *
 * Reference: Zhou, Doyle & Glover (1996), Ch. 3
 */

typedef struct {
    MURealMatrix *A;
    MURealMatrix *B;
    MURealMatrix *C;
    MURealMatrix *D;
    int n;  /* state dimension */
    int m;  /* input dimension */
    int p;  /* output dimension */
    bool is_discrete;
    double sample_time;
} MUSystem;

/* --- L1: u Analysis Result --- */

typedef struct {
    double mu_lower;
    double mu_upper;
    double frequency;
    int    iterations;
    bool   converged;
} MUMuResult;

/* --- L1: D-Scaling Matrices --- */
/*
 * Upper bound: u(M) <= inf_D sigma_bar(D M D^{-1})
 *
 * D in D_Delta = { diag(d1 I_{r1}, ..., dS I_{rS}, D1, ..., DF) :
 *                   di > 0, Dj = Dj* > 0 }
 *
 * D commutes with Delta: D Delta = Delta D for all Delta in Delta.
 */

typedef struct {
    MUComplex  *data;
    int         dim;
    MUUncertaintyStructure *structure;
} MUDScale;

/* --- L1: DK-Iteration State --- */
/*
 * D-step: Fix K, fit D-scales to minimize u upper bound
 * K-step: Fix D(s), solve H-infinity problem for shaped plant
 *
 * Reference: Doyle (1985), Balas et al. (1991)
 */

typedef struct {
    int         iteration;
    double      mu_peak;
    double      mu_peak_prev;
    double      gamma;
    double      gamma_prev;
    bool        converged;
    const char *stop_reason;
    int         max_iterations;
    double      tolerance;
} MUDKState;

/* --- L2: Frequency Grid --- */

typedef struct {
    double *omega;
    int     num_points;
    double  w_min;
    double  w_max;
} MUFrequencyGrid;

/* --- L1: Interconnection Structure --- */
/*
 * Generalized plant P:
 *   [ z ]   [ P11  P12 ] [ w ]
 *   [ y ] = [ P21  P22 ] [ u ]
 */

typedef struct {
    MUSystem *plant;
} MUInterconnection;

/* --- Complex Number Utilities (inline) --- */

static inline MUComplex mu_complex(double re, double im) {
    MUComplex c = { re, im };
    return c;
}

static inline MUComplex mu_cadd(MUComplex a, MUComplex b) {
    return mu_complex(a.re + b.re, a.im + b.im);
}

static inline MUComplex mu_csub(MUComplex a, MUComplex b) {
    return mu_complex(a.re - b.re, a.im - b.im);
}

static inline MUComplex mu_cmul(MUComplex a, MUComplex b) {
    return mu_complex(a.re * b.re - a.im * b.im,
                      a.re * b.im + a.im * b.re);
}

static inline MUComplex mu_cdiv(MUComplex a, MUComplex b) {
    double den = b.re * b.re + b.im * b.im;
    if (den < MU_EPSILON) return mu_complex(0.0, 0.0);
    return mu_complex((a.re * b.re + a.im * b.im) / den,
                      (a.im * b.re - a.re * b.im) / den);
}

static inline double mu_cabs(MUComplex z) {
    return sqrt(z.re * z.re + z.im * z.im);
}

static inline double mu_cabs2(MUComplex z) {
    return z.re * z.re + z.im * z.im;
}

static inline MUComplex mu_conj(MUComplex z) {
    return mu_complex(z.re, -z.im);
}

static inline MUComplex mu_cexp(MUComplex z) {
    double e = exp(z.re);
    return mu_complex(e * cos(z.im), e * sin(z.im));
}

static inline MUComplex mu_cneg(MUComplex z) {
    return mu_complex(-z.re, -z.im);
}

static inline MUComplex mu_cinv(MUComplex z) {
    double d = z.re * z.re + z.im * z.im;
    if (d < MU_EPSILON) return mu_complex(0.0, 0.0);
    return mu_complex(z.re / d, -z.im / d);
}

/* --- Matrix Allocation/Free --- */

MUMatrix*     mu_mat_create(int rows, int cols);
MURealMatrix* mu_mat_create_real(int rows, int cols);
void          mu_mat_free(MUMatrix *m);
void          mu_mat_free_real(MURealMatrix *m);

void      mu_mat_set(MUMatrix *m, int i, int j, MUComplex val);
MUComplex mu_mat_get(const MUMatrix *m, int i, int j);

void   mu_real_mat_set(MURealMatrix *m, int i, int j, double val);
double mu_real_mat_get(const MURealMatrix *m, int i, int j);

/* --- Matrix Operations --- */

MUMatrix* mu_mat_copy(const MUMatrix *src);
MUMatrix* mu_mat_identity(int n);
MUMatrix* mu_mat_add(const MUMatrix *a, const MUMatrix *b);
MUMatrix* mu_mat_sub(const MUMatrix *a, const MUMatrix *b);
MUMatrix* mu_mat_mul(const MUMatrix *a, const MUMatrix *b);
MUMatrix* mu_mat_scale(const MUMatrix *a, MUComplex alpha);
MUMatrix* mu_mat_transpose(const MUMatrix *a);
MUMatrix* mu_mat_conjugate_transpose(const MUMatrix *a);
MUMatrix* mu_mat_inverse(const MUMatrix *a);

/* --- Matrix Decompositions --- */

void   mu_mat_svd(const MUMatrix *a, MUMatrix **u, double *s, MUMatrix **v);
int    mu_mat_eig(const MUMatrix *a, MUComplex *eigenvalues);
double mu_mat_spectral_norm(const MUMatrix *a);
double mu_mat_frobenius_norm(const MUMatrix *a);
double mu_mat_residual(const MUMatrix *a);

/* --- System Operations --- */

MUMatrix* mu_system_freqresp(const MUSystem *sys, double omega);
MUSystem* mu_system_create(int n, int m, int p);
void      mu_system_free(MUSystem *sys);
bool      mu_system_is_stable(const MUSystem *sys);

/* --- Utility --- */

MUFrequencyGrid* mu_frequency_grid_create(double w_min, double w_max, int num_points);
void             mu_frequency_grid_free(MUFrequencyGrid *grid);

MUUncertaintyStructure* mu_unc_struct_create(int num_blocks);
void                    mu_unc_struct_free(MUUncertaintyStructure *s);
void                    mu_unc_struct_add_block(MUUncertaintyStructure *s,
                            int idx, MUUncBlockType type, int size, double bound);

#ifdef __cplusplus
}
#endif

#endif /* MU_CORE_H */

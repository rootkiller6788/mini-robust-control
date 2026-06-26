#ifndef HINF_TYPES_H
#define HINF_TYPES_H
#include <stddef.h>
#include <stdint.h>
#define HINF_MAX_STATES 64
#define HINF_MAX_ITER 2000
#define HINF_EPS 1e-12
#define HINF_EPS_SQRT 1e-6
#define HINF_MIN_GAMMA 1e-8
#define HINF_MAX_GAMMA 1e8
#define HINF_METHOD_DGKF 0
#define HINF_METHOD_LMI 1
#define HINF_METHOD_LOOPSHAPE 2
#define HINF_METHOD_MIXED 3
#define HINF_REDUCE_NONE 0
#define HINF_REDUCE_BALANCED 1
#define HINF_REDUCE_HANKEL 2
#define HINF_REDUCE_MODAL 3
typedef struct { double norm, peak_frequency; int converged, iterations; double lower_bound, upper_bound; } HinfNorm;
typedef struct { double wb, A, M, wbt; int order_w1, order_w3; double rho; } HinfWeights;
typedef struct { double *data; int rows, cols, ld, owner; } HinfMatrix;
typedef struct { HinfMatrix A,B,C,D; int n,m,p,ct,valid; } HinfSS;
typedef struct { HinfMatrix A,B1,B2,C1,C2,D11,D12,D21,D22; int n,nw,nu,nz,ny,valid; } HinfGenPlant;
typedef struct { HinfMatrix Ak,Bk,Ck,Dk; int n,ny,nu; double gamma; int valid; } HinfController;
typedef struct { HinfWeights weights; double gamma,gamma_min,gamma_max,gamma_tol; int gamma_auto,method,reduction,verify; } HinfSpec;
typedef struct { HinfMatrix H; int n; } HinfHamiltonian;
typedef struct { HinfMatrix X; int n,stabilizing,converged,iterations; double residual; } HinfCARE;
typedef struct { HinfMatrix F0,*Fk; int m,k; double *x; int feasible; } HinfLMI;
typedef struct { int bounded_real; double gamma; double *eigvals; int n_eigvals; } HinfBRL;
typedef struct { double mu_ub,mu_lb; int iterations,converged; double *freq,*mu_values; int nfreq; } HinfMu;
typedef struct { HinfSS W1,W2,Gs; double emax,epsilon; } HinfLoopShape;
typedef struct { HinfSS plant,controller,closed_loop; double gamma_achieved; } HinfClosedLoop;
HinfMatrix hinf_matrix_alloc(int rows, int cols);
void hinf_matrix_free(HinfMatrix *M);
HinfMatrix hinf_matrix_view(double *data, int rows, int cols, int ld);
void hinf_matrix_copy(HinfMatrix *dst, const HinfMatrix *src);
void hinf_matrix_set(HinfMatrix *M, int i, int j, double val);
double hinf_matrix_get(const HinfMatrix *M, int i, int j);
void hinf_matrix_print(const HinfMatrix *M, const char *name);
HinfSS hinf_ss_alloc(int n, int m, int p, int ct);
void hinf_ss_free(HinfSS *sys);
HinfGenPlant hinf_genplant_create(int n, int nw, int nu, int nz, int ny);
void hinf_genplant_free(HinfGenPlant *P);
HinfController hinf_controller_alloc(int n, int ny, int nu);
void hinf_controller_free(HinfController *K);
HinfSpec hinf_spec_default(void);
int hinf_matrix_check(const HinfMatrix *M);
void hinf_matrix_set_identity(HinfMatrix *M);
void hinf_matrix_fill(HinfMatrix *M, double val);
void hinf_matrix_set_diag(HinfMatrix *M, double val);

/* ---- Matrix arithmetic ---- */
void hinf_mat_scale(HinfMatrix *C, double alpha, const HinfMatrix *A);
void hinf_mat_add(HinfMatrix *C, const HinfMatrix *A, const HinfMatrix *B);
void hinf_mat_sub(HinfMatrix *C, const HinfMatrix *A, const HinfMatrix *B);
void hinf_mat_mul_elem(HinfMatrix *C, const HinfMatrix *A, const HinfMatrix *B);
void hinf_mat_mul(HinfMatrix *C, const HinfMatrix *A, const HinfMatrix *B);
double hinf_mat_norm_fro(const HinfMatrix *A);
double hinf_mat_norm_inf(const HinfMatrix *A);
double hinf_mat_norm_maxabs(const HinfMatrix *A);
double hinf_mat_trace(const HinfMatrix *A);
double hinf_mat_det(HinfMatrix *A);
int hinf_lu_decomp(HinfMatrix *A, int *ipiv);
int hinf_lu_solve(HinfMatrix *A, int *ipiv, double *b);
int hinf_mat_inv(HinfMatrix *A);
void hinf_mat_transpose(HinfMatrix *B, const HinfMatrix *A);
int hinf_mat_is_equal(const HinfMatrix *A, const HinfMatrix *B, double tol);
void hinf_mat_symmetrize(HinfMatrix *A);
int hinf_cholesky(HinfMatrix *A);
#endif

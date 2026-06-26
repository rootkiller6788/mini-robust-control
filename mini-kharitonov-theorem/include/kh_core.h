#ifndef KH_CORE_H
#define KH_CORE_H
#include <stdbool.h>
#include <stddef.h>

/* ==============================================================
 * kh_core.h - Kharitonov Theorem Core Types
 *
 * Kharitonov (1978) proved that an interval polynomial family
 * is Hurwitz-stable iff four specific "Kharitonov polynomials"
 * are Hurwitz-stable. This reduces checking infinitely many
 * polynomials to checking exactly four.
 *
 * Interval polynomial:
 *   p(s; a) = a_n s^n + a_{n-1} s^{n-1} + ... + a_0
 *   where a_i ∈ [a_i^-, a_i^+], a_n > 0
 *
 * The four Kharitonov polynomials K1, K2, K3, K4 are formed by
 * taking alternating min/max coefficient values.
 *
 * Theorem (Kharitonov, 1978):
 *   The interval family is Hurwitz-stable ⇔ K1, K2, K3, K4 are
 *   all Hurwitz-stable.
 *
 * References:
 *   Kharitonov, V.L. (1978). "Asymptotic stability of an
 *     equilibrium position of a family of systems of linear
 *     differential equations." Differential Equations,
 *     14(11):2086-2088.
 *   Bhattacharyya, Chapellat, Keel (1995). Robust Control:
 *     The Parametric Approach. Prentice Hall.
 *   Barmish, B.R. (1994). New Tools for Robustness of Linear
 *     Systems. Macmillan.
 * ============================================================== */

/* ---- Stability result types ---- */
typedef enum {
    KH_STABLE        = 0,  /* All roots in open LHP */
    KH_UNSTABLE      = 1,  /* At least one root in open RHP */
    KH_MARGINAL      = 2,  /* Roots on jω-axis, none in RHP */
    KH_UNCERTAIN     = 3,  /* Cannot determine (numerical) */
    KH_INCONCLUSIVE  = 4   /* Test inconclusive */
} KH_StabilityResult;

/* ---- Polynomial type flags ---- */
typedef enum {
    KH_POLY_MONIC           = 0x01,  /* Leading coeff = 1 */
    KH_POLY_REAL_COEFFS     = 0x02,  /* All real coefficients */
    KH_POLY_HURWITZ         = 0x04,  /* All roots in LHP */
    KH_POLY_INTERVAL        = 0x08   /* Interval polynomial */
} KH_PolyFlags;

/* ---- A single real interval [lo, hi] ---- */
typedef struct {
    double lo;   /* lower bound */
    double hi;   /* upper bound */
} KH_Interval;

/* ---- Interval polynomial structure ---- */
typedef struct {
    int          degree;       /* degree n ≥ 1 */
    KH_Interval* coeffs;       /* length degree+1 */
    double*      nominal;      /* nominal (center) polynomial coefficients */
    bool         is_monic;     /* a_n ≡ 1 */
} KH_IntervalPoly;

/* ---- Single polynomial (fixed coefficients) ---- */
typedef struct {
    int     degree;
    double* coeffs;    /* coeffs[0]=a_0 ... coeffs[degree]=a_n */
} KH_Polynomial;

/* ---- The four Kharitonov polynomials ---- */
typedef struct {
    KH_Polynomial  K1;   /* min min max max min min ... */
    KH_Polynomial  K2;   /* max max min min max max ... */
    KH_Polynomial  K3;   /* max min min max max min ... */
    KH_Polynomial  K4;   /* min max max min min max ... */
    int            degree;
} KH_KharitonovSet;

/* ---- Comprehensive stability report ---- */
typedef struct {
    KH_StabilityResult  overall;           /* family stability status */
    KH_StabilityResult  k1_result;         /* K1 polynomial stability */
    KH_StabilityResult  k2_result;         /* K2 polynomial stability */
    KH_StabilityResult  k3_result;         /* K3 polynomial stability */
    KH_StabilityResult  k4_result;         /* K4 polynomial stability */
    int                 n_checked;         /* number of roots checked */
    int                 n_stable_roots;    /* roots with negative real part */
    int                 n_unstable_roots;  /* roots with positive real part */
    double              worst_case_real;   /* max real part among roots */
    double              stability_margin;  /* -max(Re(λ)) */
    const char*         message;           /* human-readable summary */
} KH_StabilityReport;

/* ---- Routh-Hurwitz table structure ---- */
typedef struct {
    double** rows;        /* Routh table rows */
    int      n_rows;      /* number of rows = degree+1 */
    int      n_cols;      /* number of columns for table */
    int      degree;      /* polynomial degree */
    int      sign_changes; /* number of sign changes in first column */
} KH_RouthTable;

/* ---- Value set for frequency analysis ---- */
typedef struct {
    double    omega;         /* frequency point (rad/s) */
    double*   real_part;     /* Re(p(jω)) values for sampled corners */
    double*   imag_part;     /* Im(p(jω)) values for sampled corners */
    double    real_min;      /* min Re(p(jω)) */
    double    real_max;      /* max Re(p(jω)) */
    double    imag_min;      /* min Im(p(jω)) */
    double    imag_max;      /* max Im(p(jω)) */
    int       n_samples;     /* number of sampled values */
    bool      zero_excluded; /* does rectangle exclude zero? */
} KH_ValueSet;

/* ---- Core API ---- */

/* Polynomial lifecycle */
KH_Polynomial*  kh_poly_create(int degree);
void            kh_poly_free(KH_Polynomial* p);
KH_Polynomial*  kh_poly_copy(const KH_Polynomial* p);
KH_Polynomial*  kh_poly_derivative(const KH_Polynomial* p);
void            kh_poly_set_coeff(KH_Polynomial* p, int idx, double val);
double          kh_poly_eval(const KH_Polynomial* p, double x);
double          kh_poly_eval_complex_real(const KH_Polynomial* p, double omega);
double          kh_poly_eval_complex_imag(const KH_Polynomial* p, double omega);
void            kh_poly_scale(KH_Polynomial* p, double s);
void            kh_poly_print(const KH_Polynomial* p);

/* Interval polynomial lifecycle */
KH_IntervalPoly* kh_interval_poly_create(int degree);
void             kh_interval_poly_free(KH_IntervalPoly* ip);
void             kh_interval_poly_set_coeff(KH_IntervalPoly* ip, int idx, double lo, double hi);
void             kh_interval_poly_set_symmetric(KH_IntervalPoly* ip, int idx, double nominal, double pct);
void             kh_interval_poly_print(const KH_IntervalPoly* ip);

/* Interval operations */
KH_Interval  kh_interval_make(double lo, double hi);
double       kh_interval_center(const KH_Interval* iv);
double       kh_interval_radius(const KH_Interval* iv);
double       kh_interval_width(const KH_Interval* iv);
bool         kh_interval_contains(const KH_Interval* iv, double x);
bool         kh_interval_is_valid(const KH_Interval* iv);
KH_Interval  kh_interval_add(const KH_Interval* a, const KH_Interval* b);
KH_Interval  kh_interval_sub(const KH_Interval* a, const KH_Interval* b);
KH_Interval  kh_interval_mul(const KH_Interval* a, const KH_Interval* b);
KH_Interval  kh_interval_div(const KH_Interval* a, const KH_Interval* b);
KH_Interval  kh_interval_neg(const KH_Interval* a);
KH_Interval  kh_interval_intersect(const KH_Interval* a, const KH_Interval* b);
KH_Interval  kh_interval_hull(const KH_Interval* a, const KH_Interval* b);

/* Stability report lifecycle */
KH_StabilityReport* kh_report_create(void);
void                kh_report_free(KH_StabilityReport* r);
void                kh_report_print(const KH_StabilityReport* r);

#endif /* KH_CORE_H */

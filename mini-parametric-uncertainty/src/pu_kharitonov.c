/* ============================================================================
 * pu_kharitonov.c — Kharitonov Theorem Implementation
 *
 * V.L. Kharitonov (1978): An interval polynomial family
 *   P(s,q) = sum_{i=0}^{n} [a_i^-, a_i^+] s^i
 * is Hurwitz stable iff the four Kharitonov corner polynomials are Hurwitz.
 *
 * Theoretical foundation (L4: Fundamental Law):
 *   The four Kharitonov polynomials K1, K2, K3, K4 are constructed by
 *   selecting coefficient bounds according to a modulo-4 pattern:
 *
 *   For coefficient index i with mod4 = i % 4:
 *     K1: {0:low, 1:low, 2:high, 3:high}  (--++ pattern)
 *     K2: {0:high, 1:low, 2:low, 3:high}  (+--+ pattern)
 *     K3: {0:high, 1:high, 2:low, 3:low}   (++-- pattern)
 *     K4: {0:low, 1:high, 2:high, 3:low}   (-++- pattern)
 *
 * Prerequisites for Kharitonov's theorem:
 *   - Independent interval coefficients
 *   - Fixed degree (no degree drop)
 *   - Continuous-time Hurwitz stability
 *   - Real coefficients
 *
 * References:
 *   Kharitonov, V.L. (1978). "Asymptotic stability of an equilibrium
 *     position of a family of systems of linear differential equations."
 *     Differentsial'nye Uravneniya, 14(11):2086-2088.
 *   Barmish, B.R. (1994). "New Tools for Robustness of Linear Systems."
 *     Macmillan, Ch. 4-5.
 *   Bhattacharyya, Chapellat, Keel (1995). "Robust Control: The Parametric
 *     Approach." Prentice Hall, Ch. 3.
 * ============================================================================ */

#include "pu_kharitonov.h"
#include "pu_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Kharitonov pattern lookup: which bound to use for each coefficient index.
 * pattern[k][mod4] = 0 -> use lower bound, 1 -> use upper bound. */
static const int KH_PATTERN[4][4] = {
    {0, 0, 1, 1},  /* K1: --++ */
    {1, 0, 0, 1},  /* K2: +--+ */
    {1, 1, 0, 0},  /* K3: ++-- */
    {0, 1, 1, 0}   /* K4: -++- */
};

/**
 * Build a single Kharitonov polynomial from the interval bounds.
 *
 * Kharitonov's theorem relies on the observation that the Hurwitz
 * stability of an interval polynomial depends only on the four
 * extremal (corner) polynomials constructed from alternating
 * coefficient bounds.
 *
 * This is a consequence of the convexity of the Hurwitz region
 * in the coefficient space for even and odd parts of the polynomial.
 *
 * @param ip      Source interval polynomial
 * @param k_idx   Kharitonov index: 0=K1, 1=K2, 2=K3, 3=K4
 * @param kp      Output Kharitonov polynomial (caller provides storage)
 */
static void build_single_kharitonov(const PU_IntervalPolynomial *ip,
                                     int k_idx, PU_KharitonovPolynomial *kp)
{
    int n = ip->degree + 1;
    kp->degree = ip->degree;
    kp->coeff = (double *)malloc((size_t)n * sizeof(double));
    kp->kharitonov_index = k_idx;
    kp->is_hurwitz = false;
    for (int i = 0; i < n; i++) {
        int mod4 = i % 4;
        kp->coeff[i] = KH_PATTERN[k_idx][mod4]
                       ? ip->coeff_upper[i]
                       : ip->coeff_lower[i];
    }
}

int pu_kharitonov_build(const PU_IntervalPolynomial *ip, PU_KharitonovPolynomial kp[4])
{
    if (!ip || ip->degree < 1 || !kp) return 0;
    for (int k = 0; k < 4; k++) {
        build_single_kharitonov(ip, k, &kp[k]);
    }
    return 4;
}

PU_KharitonovPolynomial pu_kharitonov_build_one(const PU_IntervalPolynomial *ip, int index)
{
    PU_KharitonovPolynomial kp;
    if (!ip || index < 0 || index > 3) {
        kp.degree = -1; kp.coeff = NULL;
        kp.kharitonov_index = -1; kp.is_hurwitz = false;
        return kp;
    }
    build_single_kharitonov(ip, index, &kp);
    return kp;
}

void pu_kharitonov_free(PU_KharitonovPolynomial *kp)
{
    if (!kp) return;
    free(kp->coeff);
    kp->coeff = NULL;
    kp->degree = 0;
}

void pu_kharitonov_free_all(PU_KharitonovPolynomial kp[4], int n)
{
    for (int i = 0; i < n && i < 4; i++) {
        pu_kharitonov_free(&kp[i]);
    }
}

/* ============================================================================
 * Kharitonov Stability Test (L4: Fundamental Law)
 *
 * The Kharitonov test is necessary and sufficient for interval polynomials
 * with independent coefficients. It reduces an infinite family of polynomials
 * to just four fixed-coefficient polynomials.
 *
 * Key insight: The Hurwitz region is a convex set in even-odd coefficient
 * space, and the image of an axis-aligned box under the even-odd decomposition
 * is also convex, so extremal values occur at corners.
 * ============================================================================ */

PU_StabilityStatus pu_kharitonov_test(const PU_IntervalPolynomial *ip)
{
    if (!ip || ip->degree < 1) return PU_INDETERMINATE;
    /* Necessary condition: leading coefficient interval must not contain 0 */
    if (ip->coeff_lower[ip->degree] <= 0.0)
        return PU_UNSTABLE;
    PU_KharitonovPolynomial kp[4];
    pu_kharitonov_build(ip, kp);
    for (int k = 0; k < 4; k++) {
        kp[k].is_hurwitz = pu_is_hurwitz_stable(kp[k].coeff, kp[k].degree);
        if (!kp[k].is_hurwitz) {
            pu_kharitonov_free_all(kp, 4);
            return PU_UNSTABLE;
        }
    }
    pu_kharitonov_free_all(kp, 4);
    return PU_STABLE;
}

PU_StabilityStatus pu_kharitonov_test_verbose(const PU_IntervalPolynomial *ip)
{
    printf("=== Kharitonov Robust Stability Test ===\n");
    pu_interval_poly_print(ip);
    printf("\n");
    if (!ip || ip->degree < 1) {
        printf("Result: INDETERMINATE (invalid input)\n");
        return PU_INDETERMINATE;
    }
    if (ip->coeff_lower[ip->degree] <= 0.0) {
        printf("Result: UNSTABLE (leading coefficient interval contains 0)\n");
        return PU_UNSTABLE;
    }
    PU_KharitonovPolynomial kp[4];
    pu_kharitonov_build(ip, kp);
    const char *knames[4] = {"K1(--++)", "K2(+--+)", "K3(++--)", "K4(-++-)"};
    bool all_stable = true;
    for (int k = 0; k < 4; k++) {
        kp[k].is_hurwitz = pu_is_hurwitz_stable(kp[k].coeff, kp[k].degree);
        printf("%s: ", knames[k]);
        for (int i = kp[k].degree; i >= 0; i--)
            printf("%+.4f*s^%d ", kp[k].coeff[i], i);
        printf(" -> %s\n", kp[k].is_hurwitz ? "STABLE" : "UNSTABLE");
        if (!kp[k].is_hurwitz) all_stable = false;
    }
    printf("\nResult: %s\n", all_stable ? "ROBUSTLY STABLE" : "NOT ROBUSTLY STABLE");
    pu_kharitonov_free_all(kp, 4);
    return all_stable ? PU_STABLE : PU_UNSTABLE;
}

/* ============================================================================
 * Stability Margin Computation (L5: Algorithm)
 *
 * The stability margin rho^* is the maximum factor by which the uncertainty
 * intervals can be uniformly expanded while preserving stability.
 *
 * For interval polynomial:
 *   a_i in [a_i^0 - rho*w_i, a_i^0 + rho*w_i]
 * where a_i^0 is the nominal coefficient and w_i = (upper-lower)/2.
 *
 * We use bisection to find the critical rho^*.
 * ============================================================================ */

double pu_kharitonov_stability_margin(const PU_IntervalPolynomial *ip)
{
    if (!ip || ip->degree < 1) return 0.0;
    /* Compute nominal and half-width for each coefficient */
    int n = ip->degree + 1;
    double *nom = (double *)malloc((size_t)n * sizeof(double));
    double *hw  = (double *)malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        nom[i] = 0.5 * (ip->coeff_lower[i] + ip->coeff_upper[i]);
        hw[i]  = 0.5 * (ip->coeff_upper[i] - ip->coeff_lower[i]);
    }
    /* Bisection: find rho in [0, rho_max] where stability transitions */
    double rho_lo = 0.0, rho_hi = 10.0;
    /* Find an upper bound by doubling until unstable */
    int found_hi = 0;
    for (int t = 0; t < 10; t++) {
        /* Check stability at rho_hi */
        PU_IntervalPolynomial test = pu_interval_poly_create(ip->degree);
        for (int i = 0; i < n; i++) {
            double lo_i = nom[i] - rho_hi * hw[i];
            double hi_i = nom[i] + rho_hi * hw[i];
            pu_interval_poly_set_coeff(&test, i, lo_i, hi_i);
        }
        if (ip->is_monic) pu_interval_poly_set_monic(&test, true);
        PU_StabilityStatus s = pu_kharitonov_test(&test);
        pu_interval_poly_free(&test);
        if (s != PU_STABLE) { found_hi = 1; break; }
        rho_hi *= 2.0;
    }
    if (!found_hi) { free(nom); free(hw); return PU_INF; }
    /* Bisection */
    for (int iter = 0; iter < 50; iter++) {
        double rho_mid = 0.5 * (rho_lo + rho_hi);
        if (rho_hi - rho_lo < 1e-6) break;
        PU_IntervalPolynomial test = pu_interval_poly_create(ip->degree);
        for (int i = 0; i < n; i++) {
            double lo_i = nom[i] - rho_mid * hw[i];
            double hi_i = nom[i] + rho_mid * hw[i];
            pu_interval_poly_set_coeff(&test, i, lo_i, hi_i);
        }
        if (ip->is_monic) pu_interval_poly_set_monic(&test, true);
        PU_StabilityStatus s = pu_kharitonov_test(&test);
        pu_interval_poly_free(&test);
        if (s == PU_STABLE) rho_lo = rho_mid;
        else rho_hi = rho_mid;
    }
    free(nom); free(hw);
    return rho_lo;
}

/* ============================================================================
 * Critical Perturbation (L5: Algorithm)
 *
 * Finds the coefficient values on the stability boundary by expanding
 * the uncertainty intervals until instability occurs.
 *
 * This is important for:
 *   - Identifying the "weakest link" in the system
 *   - Guiding robust controller redesign
 *   - Sensitivity analysis for parameter tolerances
 * ============================================================================ */

int pu_kharitonov_critical_perturbation(const PU_IntervalPolynomial *ip, double *crit_vec)
{
    if (!ip || !crit_vec || ip->degree < 1) return 0;
    int n = ip->degree + 1;
    double *nom = (double *)malloc((size_t)n * sizeof(double));
    double *hw  = (double *)malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        nom[i] = 0.5 * (ip->coeff_lower[i] + ip->coeff_upper[i]);
        hw[i]  = 0.5 * (ip->coeff_upper[i] - ip->coeff_lower[i]);
    }
    double rho_crit = pu_kharitonov_stability_margin(ip);
    for (int i = 0; i < n; i++) {
        crit_vec[i] = nom[i] + rho_crit * hw[i];
    }
    free(nom); free(hw);
    return 1;
}

/* ============================================================================
 * Monte Carlo Verification (L5: Algorithm — Statistical Validation)
 *
 * Empirically verifies Kharitonov's theorem by random sampling:
 *   - If all 4 Kharitonov polynomials are stable, all random samples should be stable.
 *   - If any Kharitonov polynomial is unstable, some samples might still be stable.
 *
 * This demonstrates the "necessary and sufficient" nature of the theorem.
 * A false positive would indicate a bug; a false negative would be impossible
 * for a correct implementation of the theorem.
 * ============================================================================ */

void pu_kharitonov_monte_carlo_verify(const PU_IntervalPolynomial *ip,
                                       int n_samples,
                                       int *false_pos, int *false_neg)
{
    *false_pos = 0;
    *false_neg = 0;
    if (!ip || n_samples <= 0) return;
    PU_StabilityStatus kharitonov_result = pu_kharitonov_test(ip);
    for (int s = 0; s < n_samples; s++) {
        double *sample_coeff = (double *)malloc((size_t)(ip->degree + 1) * sizeof(double));
        pu_interval_poly_sample(ip, sample_coeff);
        bool sample_stable = pu_is_hurwitz_stable(sample_coeff, ip->degree);
        if (kharitonov_result == PU_STABLE && !sample_stable) (*false_pos)++;
        if (kharitonov_result == PU_UNSTABLE && sample_stable) (*false_neg)++;
        free(sample_coeff);
    }
}

void pu_kharitonov_print(const PU_KharitonovPolynomial *kp)
{
    if (!kp) return;
    const char *knames[4] = {"K1", "K2", "K3", "K4"};
    printf("Kharitonov %s (deg=%d): ", knames[kp->kharitonov_index], kp->degree);
    for (int i = kp->degree; i >= 0; i--) {
        printf("%+.4g*s^%d ", kp->coeff[i], i);
    }
    printf("\n  Hurwitz: %s\n", kp->is_hurwitz ? "yes" : "no");
}

/* ============================================================================
 * Discrete-Time Kharitonov Analog (L4/L5: Jury-Cohn Criterion)
 *
 * For discrete-time Schur stability, the analog of Kharitonov's theorem
 * is more complex. A necessary and sufficient condition requires checking
 * the "Kharitonov box" of 4 corner polynomials in the Cohn transformed space.
 *
 * We construct the corner polynomials using the same pattern but then
 * apply the Jury stability test instead of Routh-Hurwitz.
 *
 * Reference: Jury (1964), "Theory and Application of the z-Transform Method"
 *            Ackermann (1993), "Robust Control: The Parameter Space Approach"
 * ============================================================================ */

int pu_kharitonov_build_discrete(const PU_IntervalPolynomial *ip,
                                  PU_KharitonovPolynomial kp_dt[4])
{
    /* Discrete-time corner polynomials use the same coefficient selection
     * pattern as continuous-time. The difference is in the stability test:
     * we use Jury/Schur instead of Routh-Hurwitz. */
    return pu_kharitonov_build(ip, kp_dt);
}

PU_StabilityStatus pu_kharitonov_test_discrete(const PU_IntervalPolynomial *ip)
{
    if (!ip || ip->degree < 1) return PU_INDETERMINATE;
    PU_KharitonovPolynomial kp[4];
    pu_kharitonov_build(ip, kp);
    for (int k = 0; k < 4; k++) {
        kp[k].is_hurwitz = pu_is_schur_stable(kp[k].coeff, kp[k].degree);
        if (!kp[k].is_hurwitz) {
            pu_kharitonov_free_all(kp, 4);
            return PU_UNSTABLE;
        }
    }
    pu_kharitonov_free_all(kp, 4);
    return PU_STABLE;
}

/* ============================================================================
 * Bialas Counterexample Handling (L8: Advanced Topic)
 *
 * Bialas (1983) showed that Kharitonov's theorem does NOT extend naively
 * to interval matrices. Checking just the vertex matrices is not sufficient
 * for stability of the entire interval matrix family.
 *
 * Bialas provided a necessary and sufficient condition that requires
 * checking all 2^(n^2) vertex matrices for stability. This is feasible
 * only for very small matrices (n ≤ 5 or 6).
 *
 * Reference: Bialas, S. (1983). "A necessary and sufficient condition for
 *   the stability of interval matrices." Int. J. Control, 37(4):717-722.
 * ============================================================================ */

/** Extract a vertex matrix from an interval matrix.
 *  The bits of vertex_idx determine whether to use lower (0) or upper (1)
 *  bound for each element. */
static void extract_vertex_matrix(const PU_IntervalMatrix *im,
                                   unsigned long long vertex_idx, double **out)
{
    int idx = 0;
    for (int i = 0; i < im->rows; i++) {
        for (int j = 0; j < im->cols; j++) {
            out[i][j] = (vertex_idx & (1ULL << idx)) ? im->upper[i][j] : im->lower[i][j];
            idx++;
        }
    }
}

bool pu_bialas_counterexample_check(const PU_IntervalMatrix *im)
{
    if (!im || im->rows != im->cols) return true;
    int n = im->rows;
    /* Check a necessary condition: are all diagonal intervals negative? */
    for (int i = 0; i < n; i++) {
        if (im->upper[i][i] >= 0.0) return true;
    }
    /* Check if the matrix is symmetric — Bialas counterexample relies on
     * non-symmetric structure. Symmetric interval matrices have stronger
     * stability properties. */
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (fabs(im->lower[i][j] - im->lower[j][i]) > PU_EPS ||
                fabs(im->upper[i][j] - im->upper[j][i]) > PU_EPS)
                return true;
        }
    }
    return false;
}

/** Compute eigenvalues of a matrix and check if all have negative real parts.
 *  Uses the QR algorithm from pu_core. */
static bool matrix_is_hurwitz(double **M, int n)
{
    double *re = (double *)malloc((size_t)n * sizeof(double));
    double *im = (double *)malloc((size_t)n * sizeof(double));
    pu_eigen_qr(M, n, re, im, PU_MAX_ITER_QR, PU_QR_TOL);
    bool stable = true;
    for (int i = 0; i < n; i++) {
        if (re[i] >= -PU_EPS) { stable = false; break; }
    }
    free(re); free(im);
    return stable;
}

PU_StabilityStatus pu_bialas_vertex_check(const PU_IntervalMatrix *im)
{
    if (!im || im->rows != im->cols || im->rows > 6) return PU_INDETERMINATE;
    int n = im->rows;
    int n_elements = n * n;
    /* 2^(n^2) vertex matrices — feasible only for n ≤ 5 (max 2^25 ≈ 33M)
     * or n = 6 (2^36 ≈ 68B — too many). We limit to n ≤ 5 in practice. */
    if (n > 5) return PU_INDETERMINATE;
    unsigned long long n_vertices = 1ULL << n_elements;
    double **V = pu_matrix_alloc(n, n);
    bool all_stable = true;
    /* Check each vertex matrix */
    for (unsigned long long v = 0; v < n_vertices && all_stable; v++) {
        extract_vertex_matrix(im, v, V);
        if (!matrix_is_hurwitz(V, n)) {
            all_stable = false;
        }
    }
    pu_matrix_free(V, n);
    return all_stable ? PU_STABLE : PU_UNSTABLE;
}
#include "kh_construction.h"
#include "kh_core.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* ==============================================================
 * kh_construction.c - Kharitonov Polynomial Construction
 *
 * Implements the construction of the four Kharitonov polynomials
 * from an interval polynomial. This is the central operation
 * of the Kharitonov theorem.
 *
 * Theorem (Kharitonov, 1978):
 *   The interval polynomial family is Hurwitz-stable
 *   iff the four Kharitonov polynomials K1, K2, K3, K4
 *   are all Hurwitz-stable.
 *
 * The pattern rule for each Kharitonov polynomial:
 *
 *   K1: coefficients follow (min, min, max, max) repeating
 *       a_0^-, a_1^-, a_2^+, a_3^+, a_4^-, a_5^-, ...
 *
 *   K2: coefficients follow (max, max, min, min) repeating
 *       a_0^+, a_1^+, a_2^-, a_3^-, a_4^+, a_5^+, ...
 *
 *   K3: coefficients follow (max, min, min, max) repeating
 *       a_0^+, a_1^-, a_2^-, a_3^+, a_4^+, a_5^-, ...
 *
 *   K4: coefficients follow (min, max, max, min) repeating
 *       a_0^-, a_1^+, a_2^+, a_3^-, a_4^-, a_5^+, ...
 *
 * Equivalently, using the parity of floor(i/2):
 *   K1: use a_i^- if floor(i/2) is even, else a_i^+
 *   K2: use a_i^+ if floor(i/2) is even, else a_i^-
 *   K3: use a_i^+ if floor(i/2) even, a_i^- if floor(i/2) odd
 *        when i even; reversed when i odd
 *   K4: inverse of K3
 *
 * ============================================================== */

/* --------------------------------------------------------------
 * Knowledge Point: Kharitonov K1 coefficient rule
 * (L1 Definitions / L5 Algorithms)
 *
 * K1: --++--++... (min,min,max,max repeating)
 * Translation: floor(i/2) even → a_i^-, else → a_i^+
 * ------------------------------------------------------------ */
double kh_kharitonov_coeff_K1(const KH_IntervalPoly* ip, int i) {
    if (!ip || i < 0 || i > ip->degree) return 0.0;
    int block = i / 2;
    return (block % 2 == 0) ? ip->coeffs[i].lo : ip->coeffs[i].hi;
}

/* K2: ++--++--... (max,max,min,min repeating) */
double kh_kharitonov_coeff_K2(const KH_IntervalPoly* ip, int i) {
    if (!ip || i < 0 || i > ip->degree) return 0.0;
    int block = i / 2;
    return (block % 2 == 0) ? ip->coeffs[i].hi : ip->coeffs[i].lo;
}

/* --------------------------------------------------------------
 * Knowledge Point: Kharitonov K3 coefficient rule
 * (L1 Definitions / L5 Algorithms)
 *
 * K3: +--++--+... pattern
 * Even indices i: use max if floor(i/2) even, min if odd
 * Odd indices i: use min if floor(i/2) even, max if odd
 *
 * More compactly:
 *   K3[i] = a_i^+ if (i + floor(i/2)) is even, else a_i^-
 * ------------------------------------------------------------ */
double kh_kharitonov_coeff_K3(const KH_IntervalPoly* ip, int i) {
    if (!ip || i < 0 || i > ip->degree) return 0.0;
    /* Even i: K3 uses max if block even, min if block odd
     * Odd  i: K3 uses min if block even, max if block odd */
    int block = i / 2;
    if (i % 2 == 0) {
        return (block % 2 == 0) ? ip->coeffs[i].hi : ip->coeffs[i].lo;
    } else {
        return (block % 2 == 0) ? ip->coeffs[i].lo : ip->coeffs[i].hi;
    }
}

/* K4: -++--++-... pattern (inverse of K3) */
double kh_kharitonov_coeff_K4(const KH_IntervalPoly* ip, int i) {
    if (!ip || i < 0 || i > ip->degree) return 0.0;
    int block = i / 2;
    if (i % 2 == 0) {
        return (block % 2 == 0) ? ip->coeffs[i].lo : ip->coeffs[i].hi;
    } else {
        return (block % 2 == 0) ? ip->coeffs[i].hi : ip->coeffs[i].lo;
    }
}

/* --------------------------------------------------------------
 * Knowledge Point: Construction of all four Kharitonov polynomials
 * (L5 Algorithms)
 *
 * This is the main construction function. Returns a set of
 * four Kharitonov polynomials K1, K2, K3, K4. Each polynomial
 * has the same degree as the original interval polynomial.
 * ------------------------------------------------------------ */
KH_KharitonovSet* kh_kharitonov_construct(const KH_IntervalPoly* ip) {
    if (!ip) return NULL;
    KH_KharitonovSet* ks = (KH_KharitonovSet*)malloc(sizeof(KH_KharitonovSet));
    if (!ks) return NULL;
    ks->degree = ip->degree;

    /* Initialize each Kharitonov polynomial */
    ks->K1.degree = ip->degree;
    ks->K2.degree = ip->degree;
    ks->K3.degree = ip->degree;
    ks->K4.degree = ip->degree;

    ks->K1.coeffs = (double*)calloc((size_t)(ip->degree + 1), sizeof(double));
    ks->K2.coeffs = (double*)calloc((size_t)(ip->degree + 1), sizeof(double));
    ks->K3.coeffs = (double*)calloc((size_t)(ip->degree + 1), sizeof(double));
    ks->K4.coeffs = (double*)calloc((size_t)(ip->degree + 1), sizeof(double));

    if (!ks->K1.coeffs || !ks->K2.coeffs || !ks->K3.coeffs || !ks->K4.coeffs) {
        free(ks->K1.coeffs); free(ks->K2.coeffs);
        free(ks->K3.coeffs); free(ks->K4.coeffs);
        free(ks); return NULL;
    }

    /* Fill coefficients for each Kharitonov polynomial */
    for (int i = 0; i <= ip->degree; i++) {
        ks->K1.coeffs[i] = kh_kharitonov_coeff_K1(ip, i);
        ks->K2.coeffs[i] = kh_kharitonov_coeff_K2(ip, i);
        ks->K3.coeffs[i] = kh_kharitonov_coeff_K3(ip, i);
        ks->K4.coeffs[i] = kh_kharitonov_coeff_K4(ip, i);
    }

    return ks;
}

void kh_kharitonov_free(KH_KharitonovSet* ks) {
    if (!ks) return;
    free(ks->K1.coeffs);
    free(ks->K2.coeffs);
    free(ks->K3.coeffs);
    free(ks->K4.coeffs);
    free(ks);
}

const KH_Polynomial* kh_kharitonov_get(const KH_KharitonovSet* ks, int idx) {
    if (!ks) return NULL;
    switch (idx) {
        case 0: return &ks->K1;
        case 1: return &ks->K2;
        case 2: return &ks->K3;
        case 3: return &ks->K4;
        default: return NULL;
    }
}

void kh_kharitonov_print(const KH_KharitonovSet* ks) {
    if (!ks) return;
    printf("Kharitonov Polynomials (degree %d):\n", ks->degree);
    printf("K1(s): "); kh_poly_print(&ks->K1);
    printf("K2(s): "); kh_poly_print(&ks->K2);
    printf("K3(s): "); kh_poly_print(&ks->K3);
    printf("K4(s): "); kh_poly_print(&ks->K4);
}

/* --------------------------------------------------------------
 * Knowledge Point: Generalized Kharitonov pattern generator
 * (L5 Algorithms)
 *
 * Given a pattern type and separate min/max arrays, generates
 * the corresponding Kharitonov polynomial coefficients.
 * This is useful for verification and testing.
 * ------------------------------------------------------------ */
void kh_kharitonov_pattern(KH_KharitonovPattern pattern, int degree,
                           const double* mins, const double* maxs,
                           double* coeffs_out) {
    if (!mins || !maxs || !coeffs_out || degree < 0) return;

    for (int i = 0; i <= degree; i++) {
        int block = i / 2;
        switch (pattern) {
            case KH_PAT_K1:
                coeffs_out[i] = (block % 2 == 0) ? mins[i] : maxs[i];
                break;
            case KH_PAT_K2:
                coeffs_out[i] = (block % 2 == 0) ? maxs[i] : mins[i];
                break;
            case KH_PAT_K3:
                if (i % 2 == 0)
                    coeffs_out[i] = (block % 2 == 0) ? maxs[i] : mins[i];
                else
                    coeffs_out[i] = (block % 2 == 0) ? mins[i] : maxs[i];
                break;
            case KH_PAT_K4:
                if (i % 2 == 0)
                    coeffs_out[i] = (block % 2 == 0) ? mins[i] : maxs[i];
                else
                    coeffs_out[i] = (block % 2 == 0) ? maxs[i] : mins[i];
                break;
        }
    }
}

/* --------------------------------------------------------------
 * Knowledge Point: Corner polynomial enumeration (L5 Algorithms)
 *
 * An interval polynomial with n+1 coefficients (degree n) has
 * 2^(n+1) "corner polynomials" corresponding to all min/max
 * choices for each coefficient.
 *
 * The Kharitonov theorem reduces checking all 2^(n+1) to
 * checking just 4 specific corner polynomials. This function
 * enumerates all corners for verification purposes.
 *
 * Complexity: O(2^(n+1) * n)
 * ------------------------------------------------------------ */
KH_CornerPolys* kh_corner_polys_create(const KH_IntervalPoly* ip) {
    if (!ip) return NULL;
    int n_coeffs = ip->degree + 1;
    int n_corners = 1 << n_coeffs;
    if (n_corners > 4096) n_corners = 4096;

    KH_CornerPolys* cp = (KH_CornerPolys*)malloc(sizeof(KH_CornerPolys));
    if (!cp) return NULL;
    cp->degree = ip->degree;
    cp->n_polys = n_corners;
    cp->polys = (KH_Polynomial*)malloc((size_t)n_corners * sizeof(KH_Polynomial));
    if (!cp->polys) { free(cp); return NULL; }

    for (int corner = 0; corner < n_corners; corner++) {
        cp->polys[corner].degree = ip->degree;
        cp->polys[corner].coeffs = (double*)malloc((size_t)n_coeffs * sizeof(double));
        if (!cp->polys[corner].coeffs) {
            for (int k = 0; k < corner; k++) free(cp->polys[k].coeffs);
            free(cp->polys);
            free(cp);
            return NULL;
        }
        for (int i = 0; i < n_coeffs; i++) {
            cp->polys[corner].coeffs[i] = (corner & (1 << i))
                                          ? ip->coeffs[i].hi
                                          : ip->coeffs[i].lo;
        }
    }
    return cp;
}

void kh_corner_polys_free(KH_CornerPolys* cp) {
    if (!cp) return;
    for (int i = 0; i < cp->n_polys; i++) {
        free(cp->polys[i].coeffs);
    }
    free(cp->polys);
    free(cp);
}

/* --------------------------------------------------------------
 * Knowledge Point: Exhaustive corner verification
 * (L5 Algorithms)
 *
 * Checks all corner polynomials for Hurwitz stability.
 * If any corner is unstable, the family may be unstable.
 * But Kharitonov's theorem says only K1-K4 need to be checked.
 *
 * This function serves as a validation: checking all corners
 * should give the same result as checking K1-K4 for an
 * interval polynomial.
 * ------------------------------------------------------------ */
bool kh_check_all_corners_hurwitz(const KH_IntervalPoly* ip) {
    if (!ip) return false;
    KH_CornerPolys* cp = kh_corner_polys_create(ip);
    if (!cp) return false;

    extern KH_StabilityResult kh_polynomial_stability(const KH_Polynomial*, void*);
    for (int i = 0; i < cp->n_polys; i++) {
        KH_StabilityResult r = kh_polynomial_stability(&cp->polys[i], NULL);
        if (r != KH_STABLE) {
            kh_corner_polys_free(cp);
            return false;
        }
    }
    kh_corner_polys_free(cp);
    return true;
}

/* --------------------------------------------------------------
 * Knowledge Point: Edge polynomial (L8 Advanced Topics)
 *
 * The Edge Theorem (Bartlett, Hollot, Lin, 1988) states that
 * for a polytope of polynomials, the entire family is
 * Hurwitz-stable iff all exposed edges are Hurwitz-stable.
 *
 * An edge polynomial between vertices A and B is:
 *   p(s; lambda) = (1-lambda)*p_A(s) + lambda*p_B(s)
 *   for lambda in [0, 1]
 *
 * The Kharitonov theorem is a special case of the Edge Theorem
 * for interval (hyper-rectangle) polynomial families.
 * ------------------------------------------------------------ */
KH_Polynomial* kh_edge_polynomial(const KH_IntervalPoly* ip,
                                   int vertex_a, int vertex_b,
                                   double lambda) {
    if (!ip || lambda < 0.0 || lambda > 1.0) return NULL;
    KH_Polynomial* ep = kh_poly_create(ip->degree);
    if (!ep) return NULL;

    for (int i = 0; i <= ip->degree; i++) {
        double ca = (vertex_a & (1 << i)) ? ip->coeffs[i].hi : ip->coeffs[i].lo;
        double cb = (vertex_b & (1 << i)) ? ip->coeffs[i].hi : ip->coeffs[i].lo;
        ep->coeffs[i] = (1.0 - lambda) * ca + lambda * cb;
    }
    return ep;
}

/* --------------------------------------------------------------
 * Knowledge Point: Kharitonov polynomials in even/odd subfamilies
 * (L5 Algorithms)
 *
 * The even part of the interval polynomial gives a subfamily
 * for which the Kharitonov construction also applies.
 * This is used in the Hermite-Biehler proof of the
 * Kharitonov theorem.
 * ------------------------------------------------------------ */
KH_KharitonovSet* kh_kharitonov_even_part(const KH_IntervalPoly* ip) {
    if (!ip) return NULL;
    int even_deg = ip->degree / 2;
    KH_IntervalPoly* even_ip = kh_interval_poly_create(even_deg);
    if (!even_ip) return NULL;

    for (int i = 0; i <= even_deg; i++) {
        int orig_idx = 2 * i;
        if (orig_idx <= ip->degree) {
            even_ip->coeffs[i] = ip->coeffs[orig_idx];
            even_ip->nominal[i] = ip->nominal[orig_idx];
        }
    }

    KH_KharitonovSet* ks = kh_kharitonov_construct(even_ip);
    kh_interval_poly_free(even_ip);
    return ks;
}

KH_KharitonovSet* kh_kharitonov_odd_part(const KH_IntervalPoly* ip) {
    if (!ip) return NULL;
    int odd_deg = (ip->degree - 1) / 2;
    KH_IntervalPoly* odd_ip = kh_interval_poly_create(odd_deg);
    if (!odd_ip) return NULL;

    for (int i = 0; i <= odd_deg; i++) {
        int orig_idx = 2 * i + 1;
        if (orig_idx <= ip->degree) {
            odd_ip->coeffs[i] = ip->coeffs[orig_idx];
            odd_ip->nominal[i] = ip->nominal[orig_idx];
        }
    }

    KH_KharitonovSet* ks = kh_kharitonov_construct(odd_ip);
    kh_interval_poly_free(odd_ip);
    return ks;
}

/* --------------------------------------------------------------
 * Knowledge Point: Convex hull verification at frequency omega
 * (L8 Advanced Topics)
 *
 * Verifies that the convex hull of the four Kharitonov
 * polynomial values at frequency omega contains all 2^(n+1)
 * corner polynomial values.
 *
 * This is a computational verification of the Kharitonov
 * rectangle property: the value set is an axis-aligned
 * rectangle whose vertices are determined by K1..K4.
 * ------------------------------------------------------------ */
bool kh_verify_kharitonov_convex_hull(const KH_IntervalPoly* ip, double omega) {
    if (!ip) return false;
    int n_coeffs = ip->degree + 1;
    int n_corners = 1 << n_coeffs;
    if (n_corners > 512) n_corners = 512;

    /* Compute Kharitonov bounds */
    KH_KharitonovSet* ks = kh_kharitonov_construct(ip);
    if (!ks) return false;

    double k_values[4][2]; /* [i][0]=real, [i][1]=imag */
    const KH_Polynomial* kp[4] = {&ks->K1, &ks->K2, &ks->K3, &ks->K4};
    double real_min = INFINITY, real_max = -INFINITY;
    double imag_min = INFINITY, imag_max = -INFINITY;

    for (int i = 0; i < 4; i++) {
        k_values[i][0] = kh_poly_eval_complex_real(kp[i], omega);
        k_values[i][1] = kh_poly_eval_complex_imag(kp[i], omega);
        if (k_values[i][0] < real_min) real_min = k_values[i][0];
        if (k_values[i][0] > real_max) real_max = k_values[i][0];
        if (k_values[i][1] < imag_min) imag_min = k_values[i][1];
        if (k_values[i][1] > imag_max) imag_max = k_values[i][1];
    }

    /* Check all corner polynomials are within bounds */
    for (int corner = 0; corner < n_corners; corner++) {
        double re = 0.0, im = 0.0;
        double ow_re = 1.0, ow_im = 0.0;

        for (int i = 0; i <= ip->degree; i++) {
            double coef = (corner & (1 << i)) ? ip->coeffs[i].hi
                                               : ip->coeffs[i].lo;
            if (i % 4 == 0) re += coef * ow_re;
            else if (i % 4 == 1) im += coef * ow_im;
            else if (i % 4 == 2) re -= coef * ow_re;
            else im -= coef * ow_im;

            double next_re = ow_re * omega;
            double next_im = ow_im * omega;
            ow_re = next_re;
            ow_im = next_im;
        }

        /* Check containment with tolerance */
        if (re < real_min - 1e-8 || re > real_max + 1e-8 ||
            im < imag_min - 1e-8 || im > imag_max + 1e-8) {
            kh_kharitonov_free(ks);
            return false;
        }
    }

    kh_kharitonov_free(ks);
    return true;
}

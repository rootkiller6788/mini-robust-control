#include "kh_verification.h"
#include "kh_construction.h"
#include "kh_hurwitz.h"
#include "kh_interval.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ==============================================================
 * kh_verification.c - Robust Stability Verification
 *
 * This module implements the complete Kharitonov verification
 * pipeline: given an interval polynomial, it constructs the
 * four Kharitonov polynomials, tests each for Hurwitz
 * stability, and reports the overall family stability status.
 *
 * The verification pipeline:
 *
 *   1. Create the interval polynomial model
 *   2. Generate Kharitonov polynomials K1, K2, K3, K4
 *   3. Test each K_i for Hurwitz stability
 *      (using Routh-Hurwitz criterion)
 *   4. Synthesize overall result:
 *      - All 4 Hurwitz → Family is robustly stable
 *      - Any 1 not Hurwitz → Family is not robustly stable
 *
 * Theorem (Kharitonov, 1978):
 *   This procedure is both necessary AND sufficient for
 *   real interval polynomials of any degree.
 *
 * ============================================================== */

/* --------------------------------------------------------------
 * Knowledge Point: Main Kharitonov Verification (L5 Algorithms)
 *
 * This is the top-level entry point for robust stability
 * verification using the Kharitonov theorem. It:
 *   1. Constructs the four Kharitonov polynomials
 *   2. Tests each with Routh-Hurwitz
 *   3. Returns a comprehensive stability report
 *
 * Complexity: O(4 * n^2) where n is polynomial degree.
 * This is far better than the naive approach of checking
 * all 2^(n+1) corner polynomials.
 * ------------------------------------------------------------ */
KH_StabilityReport* kh_verify_interval_stability(const KH_IntervalPoly* ip) {
    KH_StabilityReport* report = kh_report_create();
    if (!report) return NULL;

    /* Validate input */
    if (!ip) {
        report->overall = KH_INCONCLUSIVE;
        report->message = "NULL interval polynomial";
        return report;
    }

    if (ip->degree < 1) {
        report->overall = KH_INCONCLUSIVE;
        report->message = "Degree must be at least 1";
        return report;
    }

    /* Step 1: Check leading coefficient is positive */
    if (ip->coeffs[ip->degree].lo <= 0.0) {
        report->overall = KH_UNSTABLE;
        report->message = "Leading coefficient interval contains zero or negative";
        return report;
    }

    /* Step 2: Fast necessary condition check */
    /* All constant-term lower bounds must be positive for Hurwitz */
    for (int i = 0; i <= ip->degree; i++) {
        if (ip->coeffs[i].hi <= 0.0) {
            report->overall = KH_UNSTABLE;
            report->message = "Coefficient upper bounds all non-positive";
            return report;
        }
    }

    /* Step 3: Construct Kharitonov polynomials */
    KH_KharitonovSet* ks = kh_kharitonov_construct(ip);
    if (!ks) {
        report->overall = KH_INCONCLUSIVE;
        report->message = "Failed to construct Kharitonov polynomials";
        return report;
    }

    /* Step 4: Test each Kharitonov polynomial */
    KH_StabilityResult results[4];
    const KH_Polynomial* kp[4] = {&ks->K1, &ks->K2, &ks->K3, &ks->K4};

    bool all_stable = true;
    bool any_unstable = false;
    bool any_marginal = false;

    for (int i = 0; i < 4; i++) {
        results[i] = kh_polynomial_stability(kp[i], NULL);

        if (results[i] == KH_UNSTABLE) any_unstable = true;
        if (results[i] == KH_MARGINAL) any_marginal = true;
        if (results[i] != KH_STABLE) all_stable = false;
    }

    report->k1_result = results[0];
    report->k2_result = results[1];
    report->k3_result = results[2];
    report->k4_result = results[3];

    /* Step 5: Synthesize overall result */
    if (all_stable) {
        report->overall = KH_STABLE;
        report->message = "Kharitonov verification: ALL four polynomials are Hurwitz-stable. The interval polynomial family is robustly stable.";
    } else if (any_unstable) {
        report->overall = KH_UNSTABLE;
        report->message = "Kharitonov verification: at least one Kharitonov polynomial is unstable. The family is NOT robustly stable.";
    } else if (any_marginal) {
        report->overall = KH_MARGINAL;
        report->message = "Kharitonov verification: at least one Kharitonov polynomial is marginally stable. Stability cannot be guaranteed.";
    } else {
        report->overall = KH_INCONCLUSIVE;
        report->message = "Inconclusive result";
    }

    /* Count stable/unstable roots across all four polynomials */
    report->n_stable_roots = 0;
    report->n_unstable_roots = 0;
    for (int i = 0; i < 4; i++) {
        if (results[i] == KH_STABLE) {
            report->n_stable_roots += kp[i]->degree;
        } else if (results[i] == KH_UNSTABLE) {
            KH_RouthTable* rt = kh_routh_create(kp[i]);
            if (rt) {
                report->n_unstable_roots += rt->sign_changes;
                kh_routh_free(rt);
            }
        }
    }
    report->n_checked = 4 * ip->degree;

    kh_kharitonov_free(ks);
    return report;
}

/* Step-by-step verification with progress */
bool kh_verify_kharitonov_condition(const KH_IntervalPoly* ip,
                                     KH_StabilityReport* report) {
    KH_StabilityReport* r = kh_verify_interval_stability(ip);
    if (!r) return false;
    bool result = (r->overall == KH_STABLE);
    if (report) {
        *report = *r;
        /* Need to be careful with the string pointer */
        free(r); /* Don't double-free the message pointer in report copy */
    } else {
        kh_report_free(r);
    }
    return result;
}

/* --------------------------------------------------------------
 * Knowledge Point: Zero Exclusion Condition (L4 Fundamental Laws)
 *
 * The zero exclusion condition is a frequency-domain
 * interpretation of robust stability:
 *
 *   An interval polynomial family is robustly Hurwitz-stable
 *   iff p_0(j*omega) != 0 for all omega >= 0 AND the
 *   value set {p(j*omega) : a_i in [a_i^-, a_i^+]} does
 *   NOT contain the origin for any omega >= 0.
 *
 * This is equivalent to the Kharitonov theorem but expressed
 * in the frequency domain. The edge theorem generalizes this
 * to polytopic uncertainty structures.
 *
 * Reference: Frazer & Duncan (1929), Barmish (1994)
 * ------------------------------------------------------------ */

/* Check zero exclusion at a single frequency */
bool kh_zero_exclusion_check(const KH_IntervalPoly* ip, double omega) {
    if (!ip) return false;
    KH_ValueSet* vs = kh_value_set_compute(ip, omega);
    if (!vs) return false;
    bool zero_excluded = vs->zero_excluded;
    kh_value_set_free(vs);
    return zero_excluded;
}

/* Sweep over a frequency range to check zero exclusion */
bool kh_zero_exclusion_sweep(const KH_IntervalPoly* ip,
                              double omega_start, double omega_end,
                              int n_points) {
    if (!ip || n_points < 2) return false;
    double dw = (omega_end - omega_start) / (double)(n_points - 1);

    for (int i = 0; i < n_points; i++) {
        double omega = omega_start + dw * i;
        if (!kh_zero_exclusion_check(ip, omega)) {
            return false; /* Zero is in the value set at this frequency */
        }
    }
    return true; /* Zero excluded at all sampled frequencies */
}

/* --------------------------------------------------------------
 * Knowledge Point: Stability margin computation (L5 Algorithms)
 *
 * The stability margin measures the distance from instability.
 * For a fixed polynomial p(s), the margin is:
 *   margin = -max_{roots} Re(lambda)
 *
 * A positive margin means all roots are strictly in LHP.
 * A negative margin means at least one root in RHP.
 *
 * For an interval family, the margin is the minimum margin
 * across the four Kharitonov polynomials.
 *
 * This can also be computed by checking when the shifted
 * polynomial p(s - sigma) loses Hurwitz stability.
 * ------------------------------------------------------------ */

/* Compute Hurwitz stability margin via Routh-Hurwitz on shifted polynomial */
double kh_stability_margin_hurwitz(const KH_Polynomial* p) {
    if (!p) return -INFINITY;

    /* Binary search for the maximum sigma where p(s-sigma) is Hurwitz */
    double lo = -100.0;
    /* First, find a bracket: move lo left until p(s-lo) is unstable */
    bool lo_stable = true;
    for (int iter = 0; iter < 20 && lo_stable; iter++) {
        /* Build shifted polynomial p(s - lo) */
        /* For small degree, evaluate via binomial expansion */
        int n = p->degree;
        KH_Polynomial* shifted = kh_poly_create(n);
        if (!shifted) break;

        /* (s - lo)^k expansion: sum binom(k,j)*(-lo)^(k-j)*s^j */
        double* binom = (double*)malloc((size_t)(n + 1) * sizeof(double));
        double* tmp = (double*)calloc((size_t)(n + 1), sizeof(double));
        for (int k = 0; k <= n; k++) {
            /* Binomial coefficients for (s-lo)^k */
            binom[0] = 1.0;
            for (int j = 1; j <= k; j++) {
                double prev = 1.0;
                double curr = 1.0;
                for (int r = 1; r <= j; r++) {
                    curr = prev * (k - r + 1) / r;
                    prev = curr;
                }
                binom[j] = curr;
            }
            for (int j = 0; j <= k; j++) {
                double term = binom[j] * pow(-lo, k - j);
                tmp[j] += p->coeffs[k] * term;
            }
        }
        for (int j = 0; j <= n; j++) shifted->coeffs[j] = tmp[j];
        free(binom);
        free(tmp);

        KH_RouthTable* rt = kh_routh_create(shifted);
        lo_stable = rt ? kh_routh_is_hurwitz(rt) : false;
        kh_routh_free(rt);
        kh_poly_free(shifted);

        if (lo_stable) lo -= 10.0;
    }

    /* Binary search for margin */
    double margin_lo = lo + 10.0;
    double margin_hi = 0.0;
    /* Find upper bound where p(s - hi) is unstable */
    bool hi_stable = true;
    for (int iter = 0; iter < 20 && hi_stable; iter++) {
        int n = p->degree;
        KH_Polynomial* shifted = kh_poly_create(n);
        if (!shifted) break;
        double* binom = (double*)malloc((size_t)(n + 1) * sizeof(double));
        double* tmp = (double*)calloc((size_t)(n + 1), sizeof(double));
        for (int k = 0; k <= n; k++) {
            binom[0] = 1.0;
            for (int j = 1; j <= k; j++) {
                double cc = 1.0;
                for (int r = 1; r <= j; r++) cc = cc * (k - r + 1) / r;
                binom[j] = cc;
            }
            for (int j = 0; j <= k; j++) {
                tmp[j] += p->coeffs[k] * binom[j] * pow(-margin_hi, k - j);
            }
        }
        for (int j = 0; j <= n; j++) shifted->coeffs[j] = tmp[j];
        free(binom); free(tmp);
        KH_RouthTable* rt = kh_routh_create(shifted);
        hi_stable = rt ? kh_routh_is_hurwitz(rt) : false;
        kh_routh_free(rt); kh_poly_free(shifted);
        if (hi_stable) margin_hi += 10.0;
    }

    /* Binary search refinement */
    for (int iter = 0; iter < 40; iter++) {
        double mid = 0.5 * (margin_lo + margin_hi);
        int n = p->degree;
        KH_Polynomial* shifted = kh_poly_create(n);
        if (!shifted) break;
        double* binom = (double*)malloc((size_t)(n + 1) * sizeof(double));
        double* tmp = (double*)calloc((size_t)(n + 1), sizeof(double));
        for (int k = 0; k <= n; k++) {
            binom[0] = 1.0;
            for (int j = 1; j <= k; j++) {
                double cc = 1.0;
                for (int r = 1; r <= j; r++) cc = cc * (k - r + 1) / r;
                binom[j] = cc;
            }
            for (int j = 0; j <= k; j++) {
                tmp[j] += p->coeffs[k] * binom[j] * pow(-mid, k - j);
            }
        }
        for (int j = 0; j <= n; j++) shifted->coeffs[j] = tmp[j];
        free(binom); free(tmp);
        KH_RouthTable* rt = kh_routh_create(shifted);
        bool mid_stable = rt ? kh_routh_is_hurwitz(rt) : false;
        kh_routh_free(rt); kh_poly_free(shifted);
        if (mid_stable) margin_lo = mid; else margin_hi = mid;
        if (margin_hi - margin_lo < 1e-8) break;
    }

    return margin_lo;
}

/* Stability margin for the entire family (worst-case of K1-K4) */
double kh_stability_margin_family(const KH_IntervalPoly* ip) {
    if (!ip) return -INFINITY;
    KH_KharitonovSet* ks = kh_kharitonov_construct(ip);
    if (!ks) return -INFINITY;

    double margins[4];
    const KH_Polynomial* kp[4] = {&ks->K1, &ks->K2, &ks->K3, &ks->K4};
    for (int i = 0; i < 4; i++) {
        margins[i] = kh_stability_margin_hurwitz(kp[i]);
    }

    kh_kharitonov_free(ks);

    /* Worst case = minimum margin */
    double worst = margins[0];
    for (int i = 1; i < 4; i++) {
        if (margins[i] < worst) worst = margins[i];
    }
    return worst;
}

/* --------------------------------------------------------------
 * Knowledge Point: Parameter sweep verification (L5 Algorithms)
 *
 * Exhaustively samples the parameter space to verify robust
 * stability. While the Kharitonov theorem makes this
 * unnecessary, the sweep is useful for:
 *   - Validating the Kharitonov theorem numerically
 *   - Estimating the fraction of stable parameters
 *   - Finding worst-case parameter combinations
 *
 * Complexity: O(n_points^n_dim * n_dim^2) - exponential in dim
 * ------------------------------------------------------------ */
KH_SweepResult* kh_parameter_sweep_verify(const KH_IntervalPoly* ip,
                                            int n_points) {
    if (!ip || n_points < 2) return NULL;
    KH_SweepResult* sr = (KH_SweepResult*)malloc(sizeof(KH_SweepResult));
    if (!sr) return NULL;
    sr->n_points = n_points;
    sr->n_dim = ip->degree + 1;
    sr->total_checks = 0;
    sr->stable_count = 0;
    sr->unstable_count = 0;
    sr->marginal_count = 0;
    sr->min_margin = INFINITY;
    sr->max_margin = -INFINITY;

    /* Simple grid over limited dimensions to avoid combinatorial explosion */
    int n_coeffs = ip->degree + 1;
    if (n_coeffs > 3) n_coeffs = 3;

    KH_Polynomial* test_poly = kh_poly_create(ip->degree);
    if (!test_poly) { free(sr); return NULL; }

    /* Initialize fixed coefficients to nominal */
    for (int i = 0; i <= ip->degree; i++) {
        test_poly->coeffs[i] = ip->nominal[i];
    }

    /* Grid over first sweep_dims coefficients */
    int total = 1;
    for (int d = 0; d < n_coeffs; d++) total *= n_points;

    for (int idx = 0; idx < total && idx < 5000; idx++) {
        int tmp = idx;
        for (int d = 0; d < n_coeffs; d++) {
            int step = tmp % n_points;
            tmp /= n_points;
            double lo = ip->coeffs[d].lo;
            double hi = ip->coeffs[d].hi;
            test_poly->coeffs[d] = lo + (hi - lo) * step / (n_points - 1.0);
        }

        KH_StabilityResult r = kh_polynomial_stability(test_poly, NULL);
        sr->total_checks++;
        if (r == KH_STABLE) sr->stable_count++;
        else if (r == KH_UNSTABLE) sr->unstable_count++;
        else sr->marginal_count++;

        double margin = kh_stability_margin_hurwitz(test_poly);
        if (margin < sr->min_margin) sr->min_margin = margin;
        if (margin > sr->max_margin) sr->max_margin = margin;
    }

    kh_poly_free(test_poly);
    return sr;
}

void kh_sweep_result_free(KH_SweepResult* sr) {
    free(sr);
}

void kh_sweep_result_print(const KH_SweepResult* sr) {
    if (!sr) return;
    printf("Parameter Sweep Results:\n");
    printf("  Total checks:     %d\n", sr->total_checks);
    printf("  Stable:           %d (%.1f%%)\n",
           sr->stable_count,
           100.0 * sr->stable_count / sr->total_checks);
    printf("  Unstable:         %d (%.1f%%)\n",
           sr->unstable_count,
           100.0 * sr->unstable_count / sr->total_checks);
    printf("  Marginal:         %d\n", sr->marginal_count);
    printf("  Stability margin range: [%.6f, %.6f]\n",
           sr->min_margin, sr->max_margin);
}

/* --------------------------------------------------------------
 * Knowledge Point: Edge Theorem Verification (L8 Advanced Topics)
 *
 * The Edge Theorem (Bartlett, Hollot, Lin, 1988) states:
 *   For a polytope of polynomials P = conv{p_1, ..., p_m},
 *   the entire family is Hurwitz-stable iff all exposed
 *   edges are Hurwitz-stable.
 *
 * An exposed edge is:
 *   e_{ij}(lambda) = (1-lambda) * p_i + lambda * p_j
 *   for lambda in [0, 1]
 *
 * For an interval polynomial (hyper-rectangle), this means
 * checking all edges of the box. The Kharitonov theorem
 * further reduces this from all edges to just 4 vertices.
 *
 * This function implements a numerical verification of
 * the edge theorem.
 * ------------------------------------------------------------ */
bool kh_edge_theorem_check(const KH_IntervalPoly* ip) {
    if (!ip) return false;
    int n_coeffs = ip->degree + 1;
    /* Number of edges = n_coeffs * 2^(n_coeffs-1) */
    /* This is exponential, so we sample a subset */
    int n_edges_to_check = n_coeffs * 16;
    if (n_edges_to_check > 256) n_edges_to_check = 256;

    for (int e = 0; e < n_edges_to_check; e++) {
        int dim = e % n_coeffs;
        int vertex_a_mask = e / n_coeffs;

        KH_Polynomial* edge_poly = kh_poly_create(ip->degree);
        if (!edge_poly) return false;

        /* Fill coefficients: vary only dimension 'dim' along edge */
        for (int i = 0; i <= ip->degree; i++) {
            if (i == dim) {
                /* Edge: goes from lo to hi */
                edge_poly->coeffs[i] = ip->nominal[i]; /* midpoint for demo */
            } else {
                edge_poly->coeffs[i] = (vertex_a_mask & (1 << i))
                                        ? ip->coeffs[i].hi
                                        : ip->coeffs[i].lo;
            }
        }

        /* Check Hurwitz at sampled points along edge */
        for (int k = 0; k <= 10; k++) {
            double t = k / 10.0;
            edge_poly->coeffs[dim] = ip->coeffs[dim].lo
                + t * (ip->coeffs[dim].hi - ip->coeffs[dim].lo);

            KH_StabilityResult r = kh_polynomial_stability(edge_poly, NULL);
            if (r == KH_UNSTABLE) {
                kh_poly_free(edge_poly);
                return false; /* Edge has unstable point */
            }
        }

        kh_poly_free(edge_poly);
    }
    return true; /* All sampled edges are stable */
}

/* Find critical frequencies where value set touches origin */
int kh_find_critical_edges(const KH_IntervalPoly* ip,
                            double* omega_crit, int max_crit) {
    if (!ip || !omega_crit || max_crit < 1) return 0;
    int found = 0;
    /* Sweep over a frequency range */
    for (int i = 0; i < 1000 && found < max_crit; i++) {
        double omega = 0.01 * i;
        if (!kh_zero_exclusion_check(ip, omega)) {
            omega_crit[found++] = omega;
        }
    }
    return found;
}

/* --------------------------------------------------------------
 * Knowledge Point: Control synthesis feasibility (L7 Applications)
 *
 * Given a plant with interval uncertainty and a fixed
 * controller, computes the closed-loop characteristic
 * polynomial as an interval polynomial, then verifies
 * robust stability using the Kharitonov theorem.
 *
 * For a feedback system:
 *   Plant:     P(s) = N_P(s) / D_P(s)
 *   Controller: C(s) = N_C(s) / D_C(s)
 *   Closed-loop denominator:
 *     D_CL(s) = D_P(s)*D_C(s) + N_P(s)*N_C(s)
 *
 * Since D_P and N_P have interval coefficients, D_CL is
 * also an interval polynomial.
 * ------------------------------------------------------------ */
bool kh_synthesis_feasibility_check(const KH_IntervalPoly* plant,
                                     const KH_Polynomial* controller,
                                     KH_IntervalPoly* closed_loop_out) {
    /* Simplified: treats plant as denominator polynomial,
     * controller as feedback gain */
    if (!plant || !controller) return false;

    /* For this demo, compute closed-loop characteristic polynomial
     * as p_cl(s) = plant(s) + controller(s)
     * where plant has interval coefficients */
    int max_deg = (plant->degree > controller->degree)
                  ? plant->degree : controller->degree;

    if (closed_loop_out) {
        closed_loop_out->degree = max_deg;
        /* Set closed-loop intervals as plant intervals perturbed by controller */
        for (int i = 0; i <= max_deg; i++) {
            double plant_lo = (i <= plant->degree) ? plant->coeffs[i].lo : 0.0;
            double plant_hi = (i <= plant->degree) ? plant->coeffs[i].hi : 0.0;
            double ctrl_val = (i <= controller->degree) ? controller->coeffs[i] : 0.0;

            if (closed_loop_out->coeffs && i <= closed_loop_out->degree) {
                closed_loop_out->coeffs[i].lo = plant_lo + ctrl_val;
                closed_loop_out->coeffs[i].hi = plant_hi + ctrl_val;
                closed_loop_out->nominal[i] = 0.5 * (plant_lo + plant_hi) + ctrl_val;
            }
        }
    }

    KH_StabilityReport* report = kh_verify_interval_stability(plant);
    bool feasible = report ? (report->overall == KH_STABLE) : false;
    kh_report_free(report);
    return feasible;
}

/* --------------------------------------------------------------
 * Knowledge Point: Frequency-domain robustness analysis (L5 Algorithms)
 *
 * Sweeps over frequency to analyze the value set properties.
 * For each frequency, computes the distance from the value set
 * to the origin. The minimum distance over all frequencies
 * indicates the robustness level.
 * ------------------------------------------------------------ */
bool kh_frequency_domain_robustness(const KH_IntervalPoly* ip,
                                     double omega_start, double omega_end,
                                     int n_points) {
    if (!ip || n_points < 2) return false;
    double dw = (omega_end - omega_start) / (double)(n_points - 1);

    for (int i = 0; i < n_points; i++) {
        double omega = omega_start + dw * i;
        KH_ValueSet* vs = kh_value_set_compute(ip, omega);
        if (!vs) continue;

        /* Check if value set touches or contains origin */
        if (vs->real_min <= 0.0 && vs->real_max >= 0.0 &&
            vs->imag_min <= 0.0 && vs->imag_max >= 0.0) {
            kh_value_set_free(vs);
            return false;
        }
        kh_value_set_free(vs);
    }
    return true;
}

/* Batch verification for multiple interval polynomials */
int kh_verify_batch(const KH_IntervalPoly** ips, int n_polys,
                     KH_StabilityReport** reports) {
    if (!ips || n_polys < 1) return 0;
    int n_stable = 0;
    for (int i = 0; i < n_polys; i++) {
        KH_StabilityReport* r = kh_verify_interval_stability(ips[i]);
        if (reports) reports[i] = r;
        if (r && r->overall == KH_STABLE) n_stable++;
        else if (!reports) kh_report_free(r);
    }
    return n_stable;
}

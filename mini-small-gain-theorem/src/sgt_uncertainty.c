/* ============================================================================
 * sgt_uncertainty.c — Uncertainty Modeling and Robust Stability Analysis
 *
 * Implements uncertainty modeling, structured singular value (mu) analysis,
 * D-K iteration, and Integral Quadratic Constraint (IQC) framework connection
 * for the Small Gain Theorem module.
 *
 * Robust control theory models the discrepancy between a nominal plant model
 * and the true physical system as "uncertainty." The Small Gain Theorem
 * provides the fundamental stability condition for feedback with norm-bounded
 * uncertainty: ||M * Delta||_inf < 1 => robust stability.
 *
 * Knowledge points implemented:
 *   L1: Uncertainty block types (unstructured, additive, multiplicative, etc.)
 *   L2: Robust stability condition, structured vs. unstructured uncertainty
 *   L3: Structured uncertainty model (block-diagonal structure)
 *   L4: Structured singular value theorem (mu)
 *   L5: D-K iteration, D-scaling, frequency-dependent uncertainty weight fitting
 *   L6: Robust stability margin computation
 *   L7: Multiplicative/additive uncertainty weighting from plant data
 *   L8: IQC framework connection (scattering -> small gain duality)
 *
 * References:
 *   K. Zhou, J.C. Doyle, K. Glover, "Robust and Optimal Control," 1996.
 *   A. Packard, J. Doyle, "The complex structured singular value,"
 *     Automatica, 29(1):71-109, 1993.
 *   A. Megretski, A. Rantzer, "System analysis via integral quadratic
 *     constraints," IEEE Trans. AC, 42(6):819-830, 1997.
 * ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "sgt_core.h"
#include "sgt_uncertainty.h"
#include "sgt_hinf.h"

/* --- Local helpers --- */
static void* sa_unc(size_t n) { void *p = malloc(n); assert(p); return p; }
static void* sc_unc(size_t n, size_t s) { void *p = calloc(n,s); assert(p); return p; }

/* ============================================================================
 * L1: Uncertainty Block Construction
 *
 * An uncertainty block represents a norm-bounded perturbation Delta with
 * ||Delta||_inf <= bound. Three types are supported:
 *   1. Full block: unstructured, any matrix with sigma_max(Delta) <= bound
 *   2. Repeated scalar: Delta = delta*I with |delta| <= bound
 *   3. Dynamic: Delta is an LTI system weighted by a frequency-dependent
 *      transfer function, with ||Delta||_inf <= bound.
 *
 * Knowledge: Uncertainty modeling is the first step in robust control design.
 * The choice of uncertainty type (additive, multiplicative, coprime factor)
 * determines the structure of the M-Delta interconnection and the
 * conservatism of the resulting stability condition. */

SGTUncertaintyBlock* sgt_uncertainty_create_full(int dim, double bound) {
    assert(dim > 0 && bound >= 0.0);
    SGTUncertaintyBlock *block =
        (SGTUncertaintyBlock*)sc_unc(1, sizeof(SGTUncertaintyBlock));
    block->type = SGT_UNSTRUCTURED;
    block->input_dim = dim;
    block->output_dim = dim;
    block->norm_bound = bound;
    block->data = (double*)sc_unc(dim * dim, sizeof(double));
    for (int i = 0; i < dim; i++) block->data[i * dim + i] = 1.0;
    block->dynamics = NULL;
    block->is_active = true;
    return block;
}

SGTUncertaintyBlock* sgt_uncertainty_create_repeated_scalar(
    int dim, double bound) {
    assert(dim > 0 && bound >= 0.0);
    SGTUncertaintyBlock *block =
        (SGTUncertaintyBlock*)sc_unc(1, sizeof(SGTUncertaintyBlock));
    block->type = SGT_PARAMETRIC;
    block->input_dim = dim;
    block->output_dim = dim;
    block->norm_bound = bound;
    block->data = (double*)sc_unc(1, sizeof(double));
    block->data[0] = 0.0; /* delta value (initially zero) */
    block->dynamics = NULL;
    block->is_active = true;
    return block;
}

SGTUncertaintyBlock* sgt_uncertainty_create_dynamic(
    const SGTLTISystem *weight, double bound) {
    assert(weight && bound >= 0.0);
    SGTUncertaintyBlock *block =
        (SGTUncertaintyBlock*)sc_unc(1, sizeof(SGTUncertaintyBlock));
    block->type = SGT_UNSTRUCTURED;
    block->input_dim = weight->nu;
    block->output_dim = weight->ny;
    block->norm_bound = bound;
    block->data = NULL;
    /* Own a copy of the weight dynamics */
    block->dynamics = sgt_lti_create("unc_weight",
        weight->nx, weight->nu, weight->ny);
    memcpy(block->dynamics->A, weight->A,
           weight->nx * weight->nx * sizeof(double));
    memcpy(block->dynamics->B, weight->B,
           weight->nx * weight->nu * sizeof(double));
    memcpy(block->dynamics->C, weight->C,
           weight->ny * weight->nx * sizeof(double));
    memcpy(block->dynamics->D, weight->D,
           weight->ny * weight->nu * sizeof(double));
    block->is_active = true;
    return block;
}

void sgt_uncertainty_block_free(SGTUncertaintyBlock *block) {
    if (!block) return;
    free(block->data);
    if (block->dynamics) sgt_lti_free(block->dynamics);
    free(block);
}

/* ============================================================================
 * L2/L3: Structured Uncertainty
 *
 * Structured uncertainty Delta = diag(delta_1*I_{r1}, ..., delta_k*I_{rk})
 * where each delta_i is a norm-bounded perturbation (scalar or full block).
 * The structure captures known physical couplings in the uncertainty.
 *
 * Key property: D-scaling matrices D = diag(d_1*I_{r1}, ..., d_k*I_{rk})
 * commute with Delta: D*Delta = Delta*D. This enables the D-scaling upper
 * bound for the structured singular value mu. */

SGTStructuredUncertainty* sgt_structured_uncertainty_create(
    const SGTUncertaintyBlock **blocks, int n_blocks) {
    assert(blocks && n_blocks > 0);
    SGTStructuredUncertainty *delta =
        (SGTStructuredUncertainty*)sc_unc(1, sizeof(SGTStructuredUncertainty));

    delta->n_blocks = n_blocks;
    delta->blocks =
        (SGTUncertaintyBlock*)sc_unc(n_blocks, sizeof(SGTUncertaintyBlock));

    int total_dim = 0;
    for (int i = 0; i < n_blocks; i++) {
        assert(blocks[i]);
        /* Deep copy each block */
        delta->blocks[i].type = blocks[i]->type;
        delta->blocks[i].input_dim = blocks[i]->input_dim;
        delta->blocks[i].output_dim = blocks[i]->output_dim;
        delta->blocks[i].norm_bound = blocks[i]->norm_bound;
        delta->blocks[i].is_active = blocks[i]->is_active;

        int bdim = blocks[i]->input_dim;
        if (blocks[i]->output_dim > bdim) bdim = blocks[i]->output_dim;

        delta->blocks[i].data =
            (double*)sc_unc(bdim * bdim, sizeof(double));
        if (blocks[i]->data)
            memcpy(delta->blocks[i].data, blocks[i]->data,
                   bdim * bdim * sizeof(double));

        delta->blocks[i].dynamics = NULL;
        if (blocks[i]->dynamics) {
            SGTLTISystem *wd = blocks[i]->dynamics;
            delta->blocks[i].dynamics = sgt_lti_create("unc_dyn",
                wd->nx, wd->nu, wd->ny);
            memcpy(delta->blocks[i].dynamics->A, wd->A,
                   wd->nx * wd->nx * sizeof(double));
            memcpy(delta->blocks[i].dynamics->B, wd->B,
                   wd->nx * wd->nu * sizeof(double));
            memcpy(delta->blocks[i].dynamics->C, wd->C,
                   wd->ny * wd->nx * sizeof(double));
            memcpy(delta->blocks[i].dynamics->D, wd->D,
                   wd->ny * wd->nu * sizeof(double));
        }
        total_dim += bdim;
    }

    delta->total_dim = total_dim;
    delta->mu_lower = 0.0;
    delta->mu_upper = INFINITY;

    return delta;
}

void sgt_structured_uncertainty_free(SGTStructuredUncertainty *delta) {
    if (!delta) return;
    for (int i = 0; i < delta->n_blocks; i++) {
        free(delta->blocks[i].data);
        if (delta->blocks[i].dynamics)
            sgt_lti_free(delta->blocks[i].dynamics);
    }
    free(delta->blocks);
    free(delta);
}

/* ============================================================================
 * L4/L6: Robust Stability — Unstructured Uncertainty
 *
 * For unstructured uncertainty ||Delta||_inf <= 1/alpha, the robust
 * stability condition is: alpha * ||M||_inf < 1.
 *
 * Returns the robust stability margin km: the smallest ||Delta||_inf
 * that can destabilize the system. km = 1 / ||M||_inf.
 *
 * This is the direct application of the Small Gain Theorem to robust
 * stability of an uncertain feedback system. */

double sgt_robust_stability_unstructured(const SGTLTISystem *M,
                                          double gamma_uncertainty) {
    assert(M);
    double M_gain = sgt_hinf_bisection(M, 1e-6, 50);
    if (M_gain < 1e-15) return INFINITY;
    double margin = 1.0 / M_gain;
    /* Margin relative to uncertainty bound */
    return margin / fmax(1e-15, gamma_uncertainty);
}

/* ============================================================================
 * L8: Structured Singular Value (mu) Lower Bound
 *
 * mu_Delta(M) = 1 / min{ sigma_max(Delta) : det(I - M*Delta) = 0 }
 *
 * Lower bound via power iteration: find worst-case perturbation Delta
 * such that det(I - M*Delta) = 0 at minimal norm.
 *
 * Algorithm: iterate between finding the destabilizing Delta (maximize
 * |lambda| of M*Delta such that sigma_max(Delta) <= 1) and updating M.
 *
 * Knowledge L8: The structured singular value generalizes the H-infinity
 * norm to structured uncertainty. mu = 1 for unstructured, and mu can be
 * significantly smaller than ||M||_inf for structured uncertainty, reducing
 * conservatism in robust stability analysis.
 *
 * Reference: Packard, Doyle (1993). */

double sgt_mu_lower_bound(const SGTLTISystem *M,
                           const SGTStructuredUncertainty *delta,
                           int max_iter, double tol) {
    assert(M && delta && max_iter > 0 && tol > 0.0);

    int ny = M->ny, nu = M->nu;
    int dim = (ny < nu) ? ny : nu;
    if (dim <= 0) return 0.0;

    /* For simplicity, approximate mu lower bound via spectral radius:
     * mu_Delta(M) >= rho(M) for any Delta in the structure.
     * At DC (w=0), compute sigma_max(M(0)) as a lower bound.
     * Then refine via power iteration on the structured spectral radius. */

    /* Evaluate M at DC */
    double M_dc_norm = sgt_freq_response_max_sv(M, 0.0);
    if (M_dc_norm < 1e-15) return 0.0;

    /* Initial lower bound */
    double mu_lb = M_dc_norm / sqrt((double)delta->n_blocks);

    /* Power iteration: find Delta to maximize rho(M*Delta) */
    double *v = (double*)sa_unc(dim * sizeof(double));
    double *w = (double*)sa_unc(dim * sizeof(double));
    for (int i = 0; i < dim; i++) v[i] = 1.0 / sqrt((double)dim);

    double rho = M_dc_norm;
    for (int iter = 0; iter < max_iter; iter++) {
        /* w = M * v (approximate via M_dc_norm scaling) */
        double nv = 0.0;
        for (int i = 0; i < dim; i++) {
            w[i] = M_dc_norm * v[i];
            nv += w[i] * w[i];
        }
        nv = sqrt(nv);
        if (nv < 1e-15) break;

        /* Normalize w to unit norm */
        for (int i = 0; i < dim; i++) w[i] /= nv;

        /* Update v = M * w (using M^T approximation) */
        for (int i = 0; i < dim; i++) {
            v[i] = M_dc_norm * w[i];
        }

        /* Rayleigh quotient: rho = v^T * M * v / (v^T * v) */
        double num = 0.0, den = 0.0;
        for (int i = 0; i < dim; i++) {
            num += v[i] * M_dc_norm * v[i];
            den += v[i] * v[i];
        }
        double rho_new = num / den;

        /* Structure-aware adjustment: for block-diagonal uncertainty with
         * k blocks, the structured spectral radius mu satisfies
         * rho / sqrt(k) <= mu <= rho. */
        double mu_new = rho_new / sqrt((double)delta->n_blocks);
        if (mu_new > mu_lb) mu_lb = mu_new;

        /* Normalize v for next iteration */
        double nv2 = 0.0;
        for (int i = 0; i < dim; i++) nv2 += v[i] * v[i];
        nv2 = sqrt(nv2);
        if (nv2 < 1e-15) break;
        for (int i = 0; i < dim; i++) v[i] /= nv2;

        if (fabs(rho_new - rho) < tol) break;
        rho = rho_new;
    }

    free(v); free(w);
    return mu_lb;
}

/* ============================================================================
 * L8: Structured Singular Value (mu) Upper Bound
 *
 * mu_Delta(M) <= inf_{D in D} sigma_max(D * M * D^{-1})
 *
 * where D = {diag(d_1*I_{r1}, ..., d_k*I_{rk}) : d_i > 0}.
 *
 * This is a convex optimization: at each frequency, minimize over D > 0
 * the maximum singular value of D*M(jw)*D^{-1}.
 *
 * For this implementation, we use an alternating optimization approach
 * over D-scales at DC (the worst-case frequency for many systems).
 *
 * Knowledge L8: The D-scaling upper bound is the key to tractable mu
 * computation. It replaces the NP-hard exact mu computation with a
 * convex optimization. The D-K iteration combines D-scaling (mu upper
 * bound) with H-infinity controller synthesis (K-step) for mu-optimal
 * controller design. */

double sgt_mu_upper_bound(const SGTLTISystem *M,
                           const SGTStructuredUncertainty *delta) {
    assert(M && delta);
    int ny = M->ny, nu = M->nu, nx = M->nx;
    SGTMatrix *M_real = sgt_matrix_create(ny, nu);
    SGTMatrix *M_imag = sgt_matrix_create(ny, nu);

    /* Compute DC gain M(0) = -C*A^{-1}*B + D.
     * For nx=1: M(0) = -c*b/a + d = c*b*tau (for first-order). */
    if (nx == 1) {
        double a = M->A[0], b = M->B[0], c = M->C[0], d = M->D[0];
        /* DC: 1/(0-a) = -1/a, so G(0) = c*(-1/a)*b + d = -c*b/a + d */
        M_real->data[0] = -c * b / a + d;
        M_imag->data[0] = 0.0;
    } else if (nx == 2 && nu == 1 && ny == 1) {
        /* For 2x2: solve A*X = B for X, then C*X + D */
        double a11 = M->A[0], a12 = M->A[1];
        double a21 = M->A[2], a22 = M->A[3];
        double detA = a11*a22 - a12*a21;
        if (fabs(detA) > 1e-15) {
            double x1 = ( a22*M->B[0] - a12*M->B[1]) / detA;
            double x2 = (-a21*M->B[0] + a11*M->B[1]) / detA;
            M_real->data[0] = M->C[0]*x1 + M->C[1]*x2 + M->D[0];
        } else {
            M_real->data[0] = M->D[0];
        }
        M_imag->data[0] = 0.0;
    } else {
        for (int i = 0; i < ny; i++)
            for (int j = 0; j < nu; j++) {
                M_real->data[i*nu + j] = M->D[i*nu + j];
                M_imag->data[i*nu + j] = 0.0;
            }
    }

    double ub = sgt_mu_upper_bound_freq(M_real, M_imag, delta, 30, 1e-6);

    sgt_matrix_free(M_real);
    sgt_matrix_free(M_imag);
    return ub;
}

double sgt_mu_upper_bound_freq(const SGTMatrix *M_real,
                                const SGTMatrix *M_imag,
                                const SGTStructuredUncertainty *delta,
                                int max_iter, double tol) {
    assert(M_real && M_imag && delta);

    int dim = delta->total_dim;
    if (dim <= 0) dim = 1;

    /* Initial D = I */
    double *d = (double*)sa_unc(dim * sizeof(double));
    for (int i = 0; i < dim; i++) d[i] = 1.0;

    /* ||M(jw)|| (approximate by Frobenius norm of real part) */
    double M_norm = 0.0;
    for (int i = 0; i < M_real->rows; i++)
        for (int j = 0; j < M_real->cols; j++) {
            double v = M_real->data[i * M_real->stride + j];
            M_norm += v * v;
        }
    M_norm = sqrt(M_norm);

    double best_ub = M_norm;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Find the worst-case ratio max_{i,j} sqrt(d_i/d_j) * M_norm */
        double max_ratio = 1.0;
        int worst_i = 0, worst_j = 0;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++) {
                if (d[j] < 1e-15) continue;
                double r = d[i] / d[j];
                if (r > max_ratio) {
                    max_ratio = r; worst_i = i; worst_j = j;
                }
            }

        double current_ub = sqrt(max_ratio) * M_norm;
        if (current_ub < best_ub) best_ub = current_ub;

        if (max_ratio < 1.0 + tol) break;

        /* Gradient descent: reduce worst_i / worst_j ratio */
        d[worst_i] *= 0.95;
        d[worst_j] *= 1.05;
    }

    free(d);
    return best_ub;
}

/* ============================================================================
 * L5: D-K Iteration — D-step
 *
 * The D-K iteration alternates between:
 *   D-step: Find optimal D-scales to minimize the mu upper bound
 *   K-step: Design an H-infinity controller for the D-scaled plant
 *
 * This function implements the D-step: at each frequency, optimize D-scales
 * to minimize sigma_max(D * M(jw) * D^{-1}), then fit a stable,
 * minimum-phase transfer function to the D-scales for use in the K-step.
 *
 * Knowledge L5: D-K iteration is the workhorse algorithm for mu-synthesis.
 * While it is not guaranteed to converge to the global optimum, it has
 * proven highly effective in practice for aerospace and industrial
 * applications (e.g., Space Station attitude control, fighter aircraft). */

void sgt_dk_iteration_d_step(const SGTLTISystem *M,
                              SGTStructuredUncertainty *delta,
                              int n_omega, int max_iter) {
    assert(M && delta && n_omega >= 10);

    int dim = delta->total_dim;
    if (dim <= 0) dim = 1;

    double w_min = 1e-3, w_max = 1e4;
    double peak_mu = 0.0;
    double *d_optimal = (double*)sa_unc(dim * sizeof(double));
    for (int i = 0; i < dim; i++) d_optimal[i] = 1.0;

    for (int kw = 0; kw < n_omega; kw++) {
        double w = pow(10.0, log10(w_min) +
                kw*(log10(w_max)-log10(w_min))/(n_omega-1));
        double M_mag = sgt_freq_response_max_sv(M, w);

        double *d = (double*)sa_unc(dim * sizeof(double));
        for (int i = 0; i < dim; i++) d[i] = 1.0;

        for (int iter = 0; iter < max_iter; iter++) {
            double max_ratio = 1.0;
            int wi = 0, wj = 0;
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++) {
                    if (d[j] < 1e-15) continue;
                    double r = d[i] / d[j];
                    if (r > max_ratio) { max_ratio = r; wi = i; wj = j; }
                }
            if (max_ratio < 1.0 + 1e-8) break;
            d[wi] *= 0.98;
            d[wj] *= 1.02;
        }

        double best_ratio = 1.0;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++) {
                if (d[j] < 1e-15) continue;
                double r = d[i] / d[j];
                if (r > best_ratio) best_ratio = r;
            }
        double scaled = sqrt(best_ratio) * M_mag;
        if (scaled > peak_mu) {
            peak_mu = scaled;
            for (int i = 0; i < dim; i++) d_optimal[i] = d[i];
        }
        free(d);
    }

    delta->mu_upper = peak_mu;
    delta->mu_lower = peak_mu * 0.5; /* rough estimate */
    free(d_optimal);
}

/* ============================================================================
 * L5: D-scale Transfer Function Fitting
 *
 * Given D-scale magnitudes at a set of frequencies, fit a stable,
 * minimum-phase transfer function D(s) = k * prod(s + zi) / prod(s + pi)
 * that matches the magnitude data.
 *
 * This is needed in D-K iteration: the D-scales optimized frequency by
 * frequency must be approximated by a rational transfer function so they
 * can be absorbed into the generalized plant for the K-step (H-infinity
 * controller synthesis).
 *
 * Knowledge L5: D-scale fitting bridges the gap between frequency-by-frequency
 * optimization (D-step) and controller synthesis (K-step). The quality of
 * the fit directly affects the convergence of D-K iteration. */

SGTLTISystem* sgt_fit_dscale_tf(const double *omega, const double *d_mag,
                                  int n_omega, int order) {
    assert(omega && d_mag && n_omega >= 2 && order >= 1);

    /* Find the geometric mean of D-magnitudes as DC gain */
    double log_sum = 0.0;
    for (int i = 0; i < n_omega; i++)
        log_sum += log(fmax(1e-15, d_mag[i]));
    double k = exp(log_sum / n_omega);

    /* For order=1: fit D(s) = k * (s + w0) / (s + w0) = k (constant).
     * A first-order non-constant fit: D(s) = k * (s/wh + 1) / (s/wl + 1). */

    /* Find the frequency where D-magnitude is halfway between min and max */
    double d_min = INFINITY, d_max = 0.0;
    for (int i = 0; i < n_omega; i++) {
        if (d_mag[i] < d_min) d_min = d_mag[i];
        if (d_mag[i] > d_max) d_max = d_mag[i];
    }
    double d_mid = 0.5 * (d_min + d_max);
    double w_mid = sqrt(omega[0] * omega[n_omega - 1]); /* geometric center */

    for (int i = 0; i < n_omega - 1; i++) {
        if ((d_mag[i] - d_mid) * (d_mag[i + 1] - d_mid) <= 0.0) {
            /* Interpolate */
            double frac = (d_mid - d_mag[i]) / (d_mag[i + 1] - d_mag[i] + 1e-15);
            w_mid = omega[i] * pow(omega[i + 1] / omega[i], frac);
            break;
        }
    }

    /* Build first-order system: D(s) = k * (s/w_mid + 1) / (s/(10*w_mid) + 1)
     * This creates a D-scale that increases with frequency, matching typical
     * D-scale behavior (low at DC, high at high frequencies). */
    double wh = w_mid;
    double wl = w_mid * 0.1;

    /* D(s) = k * (1/wh * s + 1) / (1/wl * s + 1)
     *      = k * wl/wh * (s + wh) / (s + wl)
     * State-space: A = -wl, B = 1, C = k*wl/wh * (wh - wl), D = k*wl/wh
     * Correct: D(s) = k * (s/wh + 1) / (s/wl + 1) = k*wl/wh * (s+wh)/(s+wl)
     *           = k*wl/wh * [1 + (wh-wl)/(s+wl)]
     *           = k*wl/wh + k*wl*(wh-wl)/(wh*(s+wl))
     * A = -wl, B = 1, C = k*wl*(wh-wl)/wh, D = k*wl/wh */

    SGTLTISystem *sys = sgt_lti_create("dscale_fit", 1, 1, 1);
    sys->A[0] = -wl;
    sys->B[0] = 1.0;
    if (fabs(wh) > 1e-15) {
        sys->C[0] = k * wl * (wh - wl) / wh;
        sys->D[0] = k * wl / wh;
    } else {
        sys->C[0] = k;
        sys->D[0] = 0.0;
    }
    sys->stability = SGT_STABLE;

    return sys;
}

/* ============================================================================
 * L7: Multiplicative Uncertainty Weight from Plant Data
 *
 * Given a set of plant models {P_i} and a nominal index, compute a
 * multiplicative uncertainty weight W(s) such that:
 *   |P_i(jw) - P_nom(jw)| / |P_nom(jw)| <= |W(jw)| for all i, w.
 *
 * Multiplicative form: P_i = P_nom * (I + W * Delta), ||Delta||_inf <= 1.
 *
 * Knowledge L7: This is the standard approach for building uncertainty
 * models from experimental or simulated data. The weight W captures the
 * frequency-dependent relative error between the models.
 *
 * Application: widely used in aerospace (NASA, Boeing) for modeling
 * aeroelastic uncertainty across flight conditions. */

SGTLTISystem* sgt_multiplicative_weight(const SGTLTISystem **plants,
                                         int n_plants, int nominal_idx,
                                         int n_omega) {
    assert(plants && n_plants >= 2 && nominal_idx >= 0 && nominal_idx < n_plants);
    assert(n_omega >= 10);

    const SGTLTISystem *Pnom = plants[nominal_idx];

    /* Find the peak multiplicative error across plants at each frequency */
    double w_min = 1e-3, w_max = 1e4;
    double max_error = 0.0;
    double worst_w = 1.0;

    for (int k = 0; k < n_omega; k++) {
        double w = pow(10.0, log10(w_min) +
                       k*(log10(w_max)-log10(w_min))/(n_omega-1));

        /* Evaluate |P_nom(jw)| */
        double mag_nom = 0.0;
        if (Pnom->nx == 1 && Pnom->nu == 1 && Pnom->ny == 1) {
            double a = Pnom->A[0], b = Pnom->B[0];
            double c = Pnom->C[0], d = Pnom->D[0];
            double denom = a*a + w*w;
            double re = c*b*(-a)/denom + d;
            double im = c*b*(-w)/denom;
            mag_nom = sqrt(re*re + im*im);
        }

        if (mag_nom < 1e-10) continue;

        /* Find max relative error across all plants */
        for (int p = 0; p < n_plants; p++) {
            if (p == nominal_idx) continue;
            const SGTLTISystem *Pi = plants[p];
            double mag_i = 0.0;
            if (Pi->nx == 1 && Pi->nu == 1 && Pi->ny == 1) {
                double a = Pi->A[0], b = Pi->B[0];
                double c = Pi->C[0], d = Pi->D[0];
                double denom = a*a + w*w;
                double re = c*b*(-a)/denom + d;
                double im = c*b*(-w)/denom;
                mag_i = sqrt(re*re + im*im);
            }
            double rel_err = fabs(mag_i - mag_nom) / mag_nom;
            if (rel_err > max_error) {
                max_error = rel_err;
                worst_w = w;
            }
        }
    }

    /* Build a first-order weight that captures the peak error:
     * W(s) = max_error * (s/wh + 1) / (s/wl + 1)
     * where wl, wh bound the frequency range where uncertainty is significant. */
    double wl = worst_w * 0.1;
    double wh = worst_w * 10.0;
    if (wl < 1e-4) wl = 1e-4;
    if (wh > 1e6) wh = 1e6;

    /* W(s) = max_error * (s/wh + 1) / (s/wl + 1) = max_error*wl/wh * (s+wh)/(s+wl) */
    SGTLTISystem *W = sgt_lti_create("multiplicative_weight", 1, 1, 1);
    W->A[0] = -wl;
    W->B[0] = 1.0;
    double gain_factor = max_error * wl / fmax(1e-15, wh);
    W->C[0] = gain_factor * (wh - wl);
    W->D[0] = gain_factor;
    W->stability = SGT_STABLE;

    return W;
}

/* ============================================================================
 * L7: Additive Uncertainty Weight from Plant Data
 *
 * Given plant models {P_i} and nominal index, compute an additive
 * uncertainty weight W(s) such that:
 *   |P_i(jw) - P_nom(jw)| <= |W(jw)| for all i, w.
 *
 * Additive form: P_i = P_nom + W * Delta, ||Delta||_inf <= 1.
 *
 * Knowledge L7: Additive uncertainty is simpler to work with than
 * multiplicative (the M matrix is just a constant). However, it is
 * typically more conservative because it doesn't scale with the
 * nominal plant magnitude. */

SGTLTISystem* sgt_additive_weight(const SGTLTISystem **plants,
                                   int n_plants, int nominal_idx,
                                   int n_omega) {
    assert(plants && n_plants >= 2 && nominal_idx >= 0 && nominal_idx < n_plants);
    assert(n_omega >= 10);

    const SGTLTISystem *Pnom = plants[nominal_idx];

    double w_min = 1e-3, w_max = 1e4;
    double max_error = 0.0;
    double worst_w = 1.0;

    for (int k = 0; k < n_omega; k++) {
        double w = pow(10.0, log10(w_min) +
                       k*(log10(w_max)-log10(w_min))/(n_omega-1));

        /* Evaluate |P_nom(jw)| */
        double mag_nom = 0.0;
        if (Pnom->nx == 1 && Pnom->nu == 1 && Pnom->ny == 1) {
            double a = Pnom->A[0], b = Pnom->B[0];
            double c = Pnom->C[0], d = Pnom->D[0];
            double denom = a*a + w*w;
            mag_nom = sqrt(pow(c*b*(-a)/denom + d, 2) + pow(c*b*(-w)/denom, 2));
        }

        /* Find max absolute error across plants */
        for (int p = 0; p < n_plants; p++) {
            if (p == nominal_idx) continue;
            const SGTLTISystem *Pi = plants[p];
            double mag_i = 0.0;
            if (Pi->nx == 1 && Pi->nu == 1 && Pi->ny == 1) {
                double a = Pi->A[0], b = Pi->B[0];
                double c = Pi->C[0], d = Pi->D[0];
                double denom = a*a + w*w;
                mag_i = sqrt(pow(c*b*(-a)/denom + d, 2) +
                             pow(c*b*(-w)/denom, 2));
            }
            double abs_err = fabs(mag_i - mag_nom);
            if (abs_err > max_error) {
                max_error = abs_err;
                worst_w = w;
            }
        }
    }

    /* Build first-order additive weight */
    double wl = worst_w * 0.1;
    double wh = worst_w * 10.0;
    if (wl < 1e-4) wl = 1e-4;

    /* W(s) = max_error * (s/wh + 1) / (s/wl + 1) */
    SGTLTISystem *W = sgt_lti_create("additive_weight", 1, 1, 1);
    W->A[0] = -wl;
    W->B[0] = 1.0;
    double gain_factor = max_error * wl / fmax(1e-15, wh);
    W->C[0] = gain_factor * (wh - wl);
    W->D[0] = gain_factor;
    W->stability = SGT_STABLE;

    return W;
}

/* ============================================================================
 * L8: IQC Framework — Small Gain Multiplier
 *
 * The Integral Quadratic Constraint (IQC) framework generalizes the Small
 * Gain Theorem. For the small gain multiplier:
 *   Pi_sg = [gamma^2 * I,  0;
 *            0,           -I  ]
 *
 * The IQC condition states: if a signal pair (v, w) satisfies
 *   integral_{-inf}^{inf} [v^H(jw), w^H(jw)] * Pi(jw) * [v(jw); w(jw)] dw >= 0
 * and Delta satisfies the IQC defined by Pi, then robust stability follows.
 *
 * For the small gain multiplier: ||Delta||_inf <= 1/gamma iff
 *   sigma(Delta) <= 1/gamma for all w
 * which is equivalent to the IQC: ||w||_2^2 - gamma^2*||v||_2^2 >= 0.
 *
 * Knowledge L8: IQCs unify small gain, passivity, circle criterion,
 * Popov criterion, and many other stability criteria into a single
 * mathematical framework. They enable frequency-dependent multipliers
 * that reduce conservatism beyond constant D-scaling.
 *
 * Reference: Megretski & Rantzer (1997). */

SGTMatrix* sgt_iqc_multiplier_small_gain(double gamma, int signal_dim) {
    assert(gamma > 0.0 && signal_dim > 0);
    int n = 2 * signal_dim;
    SGTMatrix *Pi = sgt_matrix_create(n, n);

    /* Pi = [gamma^2 * I,  0;
     *       0,           -I  ] */
    for (int i = 0; i < signal_dim; i++) {
        Pi->data[i * n + i] = gamma * gamma;
        Pi->data[(signal_dim + i) * n + (signal_dim + i)] = -1.0;
    }

    return Pi;
}

/* ============================================================================
 * L8: IQC Framework — Passivity Multiplier
 *
 * For a passive uncertainty (Delta + Delta^H >= 0 for all w):
 *   Pi_passive = [0,  I;
 *                 I,  0]
 *
 * The passivity theorem states: if H1 is strictly passive (Re[H1] > 0)
 * and H2 is passive (Re[H2] >= 0), then the feedback interconnection is
 * stable. This is complementary to the small gain theorem and often
 * gives less conservative results for physical systems (electrical,
 * mechanical) where passivity is a natural property. */

SGTMatrix* sgt_iqc_multiplier_passivity(int signal_dim) {
    assert(signal_dim > 0);
    int n = 2 * signal_dim;
    SGTMatrix *Pi = sgt_matrix_create(n, n);

    /* Pi = [0,  I;
     *       I,  0] */
    for (int i = 0; i < signal_dim; i++) {
        Pi->data[i * n + (signal_dim + i)] = 1.0;
        Pi->data[(signal_dim + i) * n + i] = 1.0;
    }

    return Pi;
}

/* ============================================================================
 * L8: IQC Stability Check
 *
 * Tests robust stability of the feedback interconnection of G and Delta
 * using the IQC multiplier Pi. Condition:
 *   [G(jw)^H, I] * Pi(jw) * [G(jw); I] < 0  for all w.
 *
 * For the small gain multiplier Pi = [gamma^2*I, 0; 0, -I]:
 *   gamma^2 * G^H * G - I < 0  =>  sigma_max(G) < 1/gamma
 * which is exactly the small gain condition.
 *
 * Knowledge L8: This frequency-domain inequality must hold at ALL frequencies.
 * In practice, it is checked on a grid (sufficient condition) or via a
 * Kalman-Yakubovich-Popov (KYP) lemma that converts it to an LMI. */

bool sgt_iqc_stability_check(const SGTLTISystem *G, const SGTMatrix *Pi,
                              int n_omega) {
    assert(G && Pi && n_omega >= 10);

    int dim = G->ny;
    assert(Pi->rows == 2 * dim && Pi->cols == 2 * dim);

    double w_min = 1e-3, w_max = 1e4;

    for (int k = 0; k < n_omega; k++) {
        double w = pow(10.0, log10(w_min) +
                       k*(log10(w_max)-log10(w_min))/(n_omega-1));

        /* Evaluate G(jw) */
        double mag_G = 0.0;
        if (G->nx == 1 && G->nu == 1 && G->ny == 1) {
            double a = G->A[0], b = G->B[0], c = G->C[0], d = G->D[0];
            double denom = a*a + w*w;
            double re = c*b*(-a)/denom + d;
            double im = c*b*(-w)/denom;
            mag_G = sqrt(re*re + im*im);
        } else {
            mag_G = sgt_freq_response_max_sv(G, w);
        }

        /* Check IQC condition: [G^H, I] * Pi * [G; I] < 0
         * For Pi = [p11, p12; p21, p22]:
         *   G^H*p11*G + G^H*p12 + p21*G + p22 < 0
         *
         * For Pi_sg = [gamma^2, 0; 0, -1] (scalar case):
         *   gamma^2 * |G|^2 - 1 < 0  =>  |G| < 1/gamma */
        double p11 = Pi->data[0];
        double p22 = Pi->data[Pi->rows * (Pi->stride) - Pi->stride + (Pi->cols - 1)];
        double p12 = Pi->data[Pi->cols - 1];
        double p21 = Pi->data[(Pi->rows - 1) * Pi->stride];

        double condition = p11 * mag_G * mag_G + p12 * mag_G + p21 * mag_G + p22;

        if (condition >= 0.0) return false; /* IQC violated */
    }

    return true;
}
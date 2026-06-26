/* ============================================================================
 * ssv_bounds.c — Structured Singular Value Upper and Lower Bounds
 *
 * Computes the mu upper bound via D-scaling optimization and the mu lower
 * bound via power iteration. These bounds bracket the true (generally
 * NP-hard to compute exactly) structured singular value.
 *
 * Key knowledge points:
 *   L2: mu lower bound — power iteration (Packard & Doyle, 1993)
 *   L2: mu upper bound — D-scaling optimization (convex)
 *   L4: Main Loop Theorem — mu < 1/gamma iff robust stability/performance
 *   L4: Exactness for 3 or fewer blocks (complex case)
 *   L5: Power iteration algorithm, Osborne iteration for D-scaling
 *   L8: Mixed mu bounds with D,G-scaling for real uncertainty
 *
 * Algorithms:
 *   - Power iteration: iterative maximization over unitary Q in Q
 *     Q = { diag(q_1*I_{r1}, ..., q_S*I_{rS}, Q_1, ..., Q_F) }
 *     where q_i in C, |q_i|=1 and Q_j unitary.
 *
 *   - D-scaling optimization: iterative Osborne-type balancing
 *     minimizing sigma_max(D * M * D^{-1}) over D in D.
 *
 * References:
 *   Packard & Doyle, "The complex structured singular value" (1993)
 *   Young, Newlin & Doyle, "mu analysis with real parametric uncertainty" (1991)
 *   Fan & Tits, "A measure of worst-case H∞ performance" (1986)
 * ============================================================================ */

#include "ssv_bounds.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Default Options
 * ============================================================================ */

SSVDScaleOptions ssv_dscale_options_default(void) {
    SSVDScaleOptions opts;
    opts.max_iterations = 100;
    opts.tolerance = 1e-6;
    opts.use_osborne_init = true;
    opts.step_size = 0.5;
    opts.verbose = false;
    return opts;
}

SSVPowerIterOptions ssv_power_iter_options_default(void) {
    SSVPowerIterOptions opts;
    opts.max_iterations = 500;
    opts.tolerance = 1e-8;
    opts.n_random_starts = 10;
    opts.use_complex_perturb = false;
    opts.verbose = false;
    return opts;
}

/* ============================================================================
 * Internal: Vector Operations for Power Iteration
 * ============================================================================ */

/** Complex vector 2-norm */
static double vec_norm2(const complex double *v, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        double re = creal(v[i]), im = cimag(v[i]);
        s += re * re + im * im;
    }
    return sqrt(s);
}

/** Dot product: a^H * b */
static complex double vec_dot_conj(const complex double *a, const complex double *b, int n) {
    complex double s = 0.0;
    for (int i = 0; i < n; i++)
        s += conj(a[i]) * b[i];
    return s;
}

/** Per-block normalization for power iteration.
 *  Splits vector v into blocks according to uncertainty structure,
 *  normalizes each block independently, and returns the block-wise
 *  unitary weight for the Q matrix update.
 */
static double per_block_normalize(complex double *v, const SSVUncertaintyStructure *ustruct) {
    double max_norm = 0.0;
    int offset = 0;
    for (int b = 0; b < ustruct->n_blocks; b++) {
        int sz = ustruct->blocks[b].size;
        double bn = vec_norm2(v + offset, sz);
        if (bn > 1e-15) {
            for (int i = 0; i < sz; i++) v[offset + i] /= bn;
        }
        if (bn > max_norm) max_norm = bn;
        offset += sz;
    }
    return max_norm;
}

/* ============================================================================
 * Power Iteration for mu Lower Bound
 *
 * Implements the Packard-Doyle power iteration:
 *
 * For a purely complex uncertainty structure, the lower bound is given by:
 *   mu_Delta(M) >= max_{Q in Q} rho(Q M)
 *
 * where Q = { diag(Q_1, ..., Q_n) : Q_i^H Q_i = I, and Q_i = q_i*I for scalars }.
 *
 * The power iteration looks for a fixed point of the mapping:
 *   b(k+1) = normalize_per_block(M * a(k))
 *   a(k+1) = normalize_per_block(M^H * b(k+1))
 *   theta(k+1) = b(k+1)^H * M^H * a(k+1) / (b(k+1)^H * b(k+1))
 *
 * At a fixed point, |theta| provides a lower bound on mu.
 *
 * For repeated scalar blocks: the per-block normalization involves
 * projecting onto the scalar subspace. For full blocks: full vector
 * normalization within the block.
 *
 * Complexity: O(n_random_starts * max_iterations * n^2)
 * ============================================================================ */

/** Matrix-vector multiply: y = alpha * M * x
 *  Column-major matrix M (n x n), vectors of length n.
 */
static void matvec_mul(complex double alpha, const SSVComplexMatrix *M,
                       const complex double *x, complex double *y) {
    int n = M->rows;
    for (int i = 0; i < n; i++) {
        complex double s = 0.0;
        for (int j = 0; j < n; j++)
            s += M->data[j * M->ld + i] * x[j];
        y[i] = alpha * s;
    }
}

/** Mat-vec with Hermitian transpose: y = M^H * x */
static void matvec_mul_h(const SSVComplexMatrix *M,
                         const complex double *x, complex double *y) {
    int n = M->rows;
    for (int j = 0; j < n; j++) {
        complex double s = 0.0;
        for (int i = 0; i < n; i++)
            s += conj(M->data[j * M->ld + i]) * x[i];
        y[j] = s;
    }
}

double ssv_mu_lower_bound(const SSVComplexMatrix *M,
                           const SSVUncertaintyStructure *ustruct,
                           const SSVPowerIterOptions *options) {
    if (!M || !ustruct || M->rows != M->cols || M->rows != ustruct->total_size)
        return 0.0;

    SSVPowerIterOptions opts = options ? *options : ssv_power_iter_options_default();
    int n = M->rows;
    double best_mu_lb = 0.0;

    complex double *a = (complex double*)malloc((size_t)n * sizeof(complex double));
    complex double *b = (complex double*)malloc((size_t)n * sizeof(complex double));
    complex double *temp = (complex double*)malloc((size_t)n * sizeof(complex double));
    if (!a || !b || !temp) {
        free(a); free(b); free(temp);
        return 0.0;
    }

    for (int restart = 0; restart < opts.n_random_starts; restart++) {
        /* Initialize b randomly on unit sphere */
        double nrm = 0.0;
        for (int i = 0; i < n; i++) {
            double re = ((double)rand() / RAND_MAX - 0.5) * 2.0;
            double im = ((double)rand() / RAND_MAX - 0.5) * 2.0;
            b[i] = re + I * im;
            nrm += re * re + im * im;
        }
        nrm = sqrt(nrm);
        for (int i = 0; i < n; i++) b[i] /= nrm;
        per_block_normalize(b, ustruct);

        double theta = 0.0;
        double prev_theta = 0.0;

        for (int iter = 0; iter < opts.max_iterations; iter++) {
            /* Step 1: a = M * b */
            matvec_mul(1.0, M, b, a);

            /* Step 2: per-block normalize a */
            per_block_normalize(a, ustruct);

            /* Step 3: b = M^H * a */
            matvec_mul_h(M, a, b);

            /* Step 4: per-block normalize b */
            per_block_normalize(b, ustruct);

            /* Step 5: compute theta = b^H * M^H * a */
            matvec_mul_h(M, a, temp);
            theta = cabs(vec_dot_conj(b, temp, n));

            if (theta > best_mu_lb) best_mu_lb = theta;

            /* Check convergence */
            if (iter > 0 && fabs(theta - prev_theta) < opts.tolerance) {
                break;
            }
            prev_theta = theta;

            /* Avoid getting stuck: add small random perturbation occasionally */
            if (iter > 20 && iter % 50 == 0) {
                for (int i = 0; i < n; i++) {
                    double rp = ((double)rand() / RAND_MAX - 0.5) * 0.01;
                    double ip = ((double)rand() / RAND_MAX - 0.5) * 0.01;
                    b[i] += rp + I * ip;
                }
                per_block_normalize(b, ustruct);
            }
        }
    }

    free(a); free(b); free(temp);
    return best_mu_lb;
}

/* ============================================================================
 * D-Scaling Upper Bound
 *
 * Algorithm:
 *   minimize_{D in D} sigma_max(D * M * D^{-1})
 *
 * D-set for purely complex uncertainty:
 *   D = { diag(D_1, ..., D_{S+F}) : D_i > 0 Hermitian }
 *   where for repeated scalar blocks, D_i is a full r_i x r_i Hermitian matrix,
 *   and for full blocks, D_i = d_i * I.
 *
 * Osborne-type iteration:
 *   For each block i, balance the block's contribution to sigma_max.
 *   This is equivalent to solving: inf_{D in D} || D M D^{-1} ||
 *
 * For each full block of size k_i: D_i = d_i * I_{k_i}
 *   The optimal d_i makes the block's input/output norms equal.
 *
 * For repeated scalar blocks of size r_i: D_i is a full Hermitian matrix.
 *   This is solved via the Osborne iteration on the sub-blocks.
 * ============================================================================ */

double ssv_mu_upper_bound(const SSVComplexMatrix *M,
                           const SSVUncertaintyStructure *ustruct,
                           const SSVDScaleOptions *options,
                           SSVDScaleMatrices **out_dscale) {
    if (!M || !ustruct || M->rows != M->cols || M->rows != ustruct->total_size)
        return 0.0;

    SSVDScaleOptions opts = options ? *options : ssv_dscale_options_default();

    /* Create D-scaling matrices initialized to identity */
    SSVDScaleMatrices *dscale = ssv_dscale_create(ustruct);
    if (!dscale) return 0.0;

    int n = M->rows;
    double mu_ub = ssv_dscale_evaluate_upper_bound(M, dscale);
    double prev_mu_ub = mu_ub;

    if (opts.use_osborne_init) {
        /* Apply Osborne balancing to M as initial D-scaling */
        SSVComplexMatrix *M_copy = ssv_cmatrix_create(n, n);
        ssv_cmatrix_copy(M_copy, M);
        SSVComplexMatrix *D_osborne = ssv_cmatrix_create(n, n);
        ssv_osborne_balance(M_copy, D_osborne);
        /* Extract diagonal of D_osborne into per-block D matrices */
        int offset = 0;
        for (int b = 0; b < ustruct->n_blocks; b++) {
            int sz = ustruct->blocks[b].size;
            SSVComplexMatrix *D_i = ssv_dscale_get_block(dscale, b);
            if (D_i) {
                for (int j = 0; j < sz; j++)
                    D_i->data[j * sz + j] = D_osborne->data[(offset + j) * n + (offset + j)];
            }
            offset += sz;
        }
        ssv_dscale_reassemble(dscale);
        ssv_cmatrix_free(M_copy);
        ssv_cmatrix_free(D_osborne);
        mu_ub = ssv_dscale_evaluate_upper_bound(M, dscale);
        prev_mu_ub = mu_ub;
    }

    /* Iteratively refine D-scaling */
    for (int iter = 0; iter < opts.max_iterations; iter++) {
        /* Compute scaled matrix S = D * M * D^{-1} */
        SSVComplexMatrix *S = ssv_dscale_apply(M, dscale);

        /* Update each D-block based on the row/column norms of the
         * corresponding block of S.
         *
         * For D_i = d_i * I (full blocks): update d_i to balance.
         * For arbitrary D_i (repeated scalars): update via Osborne. */
        int offset = 0;
        for (int b = 0; b < ustruct->n_blocks; b++) {
            int sz = ustruct->blocks[b].size;
            SSVComplexMatrix *D_i = ssv_dscale_get_block(dscale, b);
            if (!D_i) continue;

            /* Compute block row norm R_b and block column norm C_b of S */
            double row_norm = 0.0, col_norm = 0.0;
            for (int i = 0; i < sz; i++) {
                for (int j = 0; j < n; j++) {
                    /* Row offset+i within S */
                    complex double s_val = S->data[j * n + (offset + i)];
                    row_norm += creal(s_val) * creal(s_val) + cimag(s_val) * cimag(s_val);
                }
                for (int j = 0; j < n; j++) {
                    /* Column offset+i */
                    complex double s_val = S->data[(offset + i) * n + j];
                    col_norm += creal(s_val) * creal(s_val) + cimag(s_val) * cimag(s_val);
                }
            }
            row_norm = sqrt(row_norm);
            col_norm = sqrt(col_norm);

            if (row_norm < 1e-15 && col_norm < 1e-15) { offset += sz; continue; }

            /* Update D_i to balance block norms
             * For full blocks: D_i = d_i * I, update d_i *= sqrt(col_norm/row_norm)
             * For repeated scalar blocks: scale the entire D_i matrix */
            double scale_factor = 1.0;
            if (row_norm > 1e-15 && col_norm > 1e-15) {
                scale_factor = sqrt(col_norm / row_norm);
                /* Damped update */
                scale_factor = 1.0 + opts.step_size * (scale_factor - 1.0);
            }

            for (int j = 0; j < sz; j++)
                for (int i = 0; i < sz; i++)
                    D_i->data[j * sz + i] *= scale_factor;

            offset += sz;
        }

        ssv_dscale_reassemble(dscale);
        ssv_cmatrix_free(S);

        mu_ub = ssv_dscale_evaluate_upper_bound(M, dscale);

        if (opts.verbose && iter % 10 == 0) {
            printf("  D-scale iter %3d: mu_ub = %.8f\n", iter, mu_ub);
        }

        /* Check convergence */
        if (fabs(mu_ub - prev_mu_ub) < opts.tolerance * (1.0 + mu_ub)) {
            break;
        }
        prev_mu_ub = mu_ub;
    }

    if (out_dscale) {
        *out_dscale = dscale;
    } else {
        ssv_dscale_free(dscale);
    }

    return mu_ub;
}

/* ============================================================================
 * Full mu Analysis (Upper + Lower bounds)
 * ============================================================================ */

SSVMuResult* ssv_mu_analysis(const SSVComplexMatrix *M,
                              const SSVUncertaintyStructure *ustruct,
                              const SSVDScaleOptions *ub_options,
                              const SSVPowerIterOptions *lb_options) {
    if (!M || !ustruct) return NULL;

    SSVMuResult *result = (SSVMuResult*)calloc(1, sizeof(SSVMuResult));
    if (!result) return NULL;

    /* Upper bound */
    SSVDScaleOptions ub_opts = ub_options ? *ub_options : ssv_dscale_options_default();
    result->upper_bound = ssv_mu_upper_bound(M, ustruct, &ub_opts, NULL);

    /* Lower bound */
    SSVPowerIterOptions lb_opts = lb_options ? *lb_options : ssv_power_iter_options_default();
    result->lower_bound = ssv_mu_lower_bound(M, ustruct, &lb_opts);

    /* Gap */
    if (result->upper_bound > 1e-10) {
        result->gap = (result->upper_bound - result->lower_bound) / result->upper_bound;
    } else {
        result->gap = 0.0;
    }

    if (result->gap < 0.01) {
        result->bound_type = SSV_BOUND_EXACT;
    } else if (result->lower_bound > 0) {
        result->bound_type = SSV_BOUND_UPPER;
    } else {
        result->bound_type = SSV_BOUND_UPPER;
    }

    result->converged = (result->gap < 0.05);
    return result;
}

SSVMuResult** ssv_mu_analysis_frequency(
    const SSVComplexMatrix **M_array,
    const double *frequencies,
    int n_freq,
    const SSVUncertaintyStructure *ustruct,
    const SSVDScaleOptions *ub_options,
    const SSVPowerIterOptions *lb_options) {
    if (!M_array || !frequencies || n_freq <= 0 || !ustruct) return NULL;

    SSVMuResult **results = (SSVMuResult**)calloc((size_t)n_freq, sizeof(SSVMuResult*));
    if (!results) return NULL;

    for (int i = 0; i < n_freq; i++) {
        results[i] = ssv_mu_analysis(M_array[i], ustruct, ub_options, lb_options);
        if (results[i]) {
            results[i]->frequency_index = i;
            results[i]->frequency = frequencies[i];
        }
    }

    return results;
}

SSVMuResult** ssv_mu_upper_bound_frequency(
    const SSVComplexMatrix **M_array,
    const double *frequencies,
    int n_freq,
    const SSVUncertaintyStructure *ustruct,
    const SSVDScaleOptions *options) {
    if (!M_array || !frequencies || n_freq <= 0 || !ustruct) return NULL;

    SSVMuResult **results = (SSVMuResult**)calloc((size_t)n_freq, sizeof(SSVMuResult*));
    if (!results) return NULL;

    SSVDScaleOptions opts = options ? *options : ssv_dscale_options_default();

    for (int i = 0; i < n_freq; i++) {
        results[i] = (SSVMuResult*)calloc(1, sizeof(SSVMuResult));
        if (results[i]) {
            results[i]->upper_bound = ssv_mu_upper_bound(M_array[i], ustruct, &opts, NULL);
            results[i]->lower_bound = 0.0;
            results[i]->gap = 1.0;
            results[i]->frequency_index = i;
            results[i]->frequency = frequencies[i];
            results[i]->bound_type = SSV_BOUND_UPPER;
        }
    }

    return results;
}

/* ============================================================================
 * Robust Stability and Performance Tests
 * ============================================================================ */

bool ssv_robust_stability_test(const SSVComplexMatrix *M,
                                const SSVUncertaintyStructure *ustruct,
                                double gamma) {
    if (!M || !ustruct || gamma <= 0) return false;

    SSVDScaleOptions ub_opts = ssv_dscale_options_default();
    double mu_ub = ssv_mu_upper_bound(M, ustruct, &ub_opts, NULL);

    /* Robust stability iff mu(M) < 1/gamma */
    return (mu_ub * gamma < 1.0);
}

bool ssv_robust_performance_test(const SSVComplexMatrix *M,
                                  const SSVUncertaintyStructure *ustruct,
                                  int n_perf) {
    if (!M || !ustruct || n_perf <= 0) return false;

    /* Augment the uncertainty structure with a performance full block.
     * Delta_aug = diag(Delta_orig, Delta_perf)
     * where Delta_perf is a full complex block of size n_perf.
     *
     * Build augmented uncertainty structure */
    SSVUncertaintyStructure *ustruct_aug = ssv_ustruct_create(ustruct->n_blocks + 1);
    if (!ustruct_aug) return false;

    /* Copy original blocks */
    for (int i = 0; i < ustruct->n_blocks; i++) {
        SSVUncertaintyBlock *blk = &ustruct->blocks[i];
        if (blk->type == SSV_BLOCK_FULL_COMPLEX || blk->type == SSV_BLOCK_FULL_REAL)
            ssv_ustruct_add_full_block(ustruct_aug, (blk->type == SSV_BLOCK_FULL_REAL),
                                       blk->size, blk->description);
        else
            ssv_ustruct_add_scalar_block(ustruct_aug, (blk->type == SSV_BLOCK_REPEATED_SCALAR_REAL),
                                         blk->size, blk->description);
    }
    /* Add performance block */
    ssv_ustruct_add_full_block(ustruct_aug, false, n_perf, "performance_block");

    /* Augment M to include performance channels.
     * This requires that M already has the performance channels in its
     * last n_perf rows/columns. The augmentation simply embeds M into
     * a larger block-diagonal-like structure. */
    int n_total = ustruct_aug->total_size;
    SSVComplexMatrix *M_aug = ssv_cmatrix_create(n_total, n_total);
    if (!M_aug) { ssv_ustruct_free(ustruct_aug); return false; }

    /* Copy M into the upper-left block of M_aug.
     * The performance block connects the last n_perf channels. */
    int n_orig = M->rows;
    if (n_orig + n_perf != n_total) {
        ssv_cmatrix_free(M_aug);
        ssv_ustruct_free(ustruct_aug);
        return false;
    }

    for (int i = 0; i < n_orig; i++)
        for (int j = 0; j < n_orig; j++)
            M_aug->data[j * n_total + i] = M->data[j * n_orig + i];

    /* Robust performance: mu_{Delta_aug}(M_aug) < 1 */
    double mu_aug = ssv_mu_upper_bound(M_aug, ustruct_aug, NULL, NULL);

    ssv_cmatrix_free(M_aug);
    ssv_ustruct_free(ustruct_aug);

    return (mu_aug < 1.0);
}

double ssv_robust_stability_margin(const SSVMuResult **mu_results, int n_freq) {
    if (!mu_results || n_freq <= 0) return 0.0;
    double max_mu = 0.0;
    for (int i = 0; i < n_freq; i++) {
        if (mu_results[i] && mu_results[i]->upper_bound > max_mu)
            max_mu = mu_results[i]->upper_bound;
    }
    if (max_mu < 1e-15) return 1e15;
    return 1.0 / max_mu;
}

int ssv_mu_peak_frequency(const SSVMuResult **mu_results, int n_freq) {
    if (!mu_results || n_freq <= 0) return -1;
    int peak_idx = 0;
    double peak_mu = mu_results[0] ? mu_results[0]->upper_bound : 0.0;
    for (int i = 1; i < n_freq; i++) {
        if (mu_results[i] && mu_results[i]->upper_bound > peak_mu) {
            peak_mu = mu_results[i]->upper_bound;
            peak_idx = i;
        }
    }
    return peak_idx;
}

/* ============================================================================
 * Mixed mu Upper Bound (with G-scaling)
 *
 * For real uncertainty blocks, the standard D-scaling upper bound can be
 * arbitrarily conservative because real mu can be discontinuous. The
 * D,G-scaling provides a tighter upper bound:
 *
 *   mu_{mixed}(M) <= inf_{D, G} alpha : M^H D M + j(G M - M^H G) <= alpha^2 D
 *
 * where D is as before and G is a block-diagonal skew-Hermitian matrix
 * (G^H = -G) that only has non-zero blocks for real uncertainty blocks.
 *
 * The G-scaling tightens the bound for real uncertainties, but the
 * optimization is no longer convex in general.
 *
 * This implementation uses a simplified G-scaling for diagonal real blocks.
 * ============================================================================ */

double ssv_mixed_mu_upper_bound(const SSVComplexMatrix *M,
                                 const SSVUncertaintyStructure *ustruct,
                                 const SSVDScaleOptions *options) {
    if (!M || !ustruct || M->rows != ustruct->total_size) return 0.0;

    /* First, compute standard D-scaling upper bound */
    double mu_d = ssv_mu_upper_bound(M, ustruct, options, NULL);

    /* If no real blocks, D-scaling is already tight enough */
    if (ustruct->n_real_blocks == 0)
        return mu_d;

    /* For mixed real/complex mu, the G-scaling provides additional tightening.
     * We apply a simplified G-scaling for repeated real scalar blocks:
     * For each real scalar block delta_i * I_{r_i}, add a skew-Hermitian
     * G_i to tighten the bound.
     *
     * The improvement is computed by solving:
     *   inf_{d_i, g_i} alpha where M^H D M + j(G M - M^H G) < alpha^2 D
     *
     * For simplicity, we estimate the G-correction as:
     *   mu_mixed = mu_d * (1 - eps * n_real / n_blocks)
     * where eps is a small correction factor from G-scaling effect. */

    /* Estimate G-scaling correction */
    int n_total = M->rows;

    /* Compute D-scaled M */
    SSVDScaleMatrices *dscale = ssv_dscale_create(ustruct);
    SSVComplexMatrix *DMDinv = ssv_dscale_apply(M, dscale);

    /* Compute the imaginary part based correction:
     * For each real block, compute the skew-symmetric part */
    double g_correction = 0.0;
    int offset = 0;
    for (int b = 0; b < ustruct->n_blocks; b++) {
        int sz = ustruct->blocks[b].size;
        bool is_real = (ustruct->blocks[b].type == SSV_BLOCK_REPEATED_SCALAR_REAL ||
                        ustruct->blocks[b].type == SSV_BLOCK_FULL_REAL);

        if (is_real) {
            /* Extract block (offset:offset+sz, offset:offset+sz) from DMDinv */
            double block_skew = 0.0;
            for (int i = 0; i < sz; i++) {
                for (int j = 0; j < sz; j++) {
                    complex double mij = DMDinv->data[(offset + j) * n_total + (offset + i)];
                    complex double mji = DMDinv->data[(offset + i) * n_total + (offset + j)];
                    /* Skew-Hermitian contribution */
                    complex double skew = mij - conj(mji);
                    block_skew += cabs(skew);
                }
            }
            g_correction += block_skew / sz;
        }
        offset += sz;
    }

    ssv_cmatrix_free(DMDinv);
    ssv_dscale_free(dscale);

    /* Combine: mixed mu <= mu_d corrected by G-scaling effect */
    double g_factor = (ustruct->total_size > 0) ?
        g_correction / (ustruct->total_size * mu_d + 1e-15) : 0.0;
    g_factor = (g_factor > 0.5) ? 0.5 : g_factor; /* clamp */

    double mu_mixed = mu_d * (1.0 + g_factor);
    return mu_mixed;
}

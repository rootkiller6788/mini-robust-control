/* ============================================================================
 * ssv_uncertainty.c — Structured Uncertainty Modeling
 *
 * Implements the block-diagonal uncertainty structure that defines the
 * "structure" in the structured singular value (mu). Each block represents
 * a different source and type of uncertainty.
 *
 * Key knowledge points:
 *   L1: Uncertainty block types (repeated scalar real/complex, full block)
 *   L2: Structured vs unstructured uncertainty, block-diagonal form
 *   L3: D-scaling matrices that commute with Delta, commutation condition
 *   L5: D-scaling optimization (Osborne-based) for mu upper bound
 *
 * The critical insight: mu reduces conservatism over the small-gain theorem
 * by exploiting the structure of uncertainty. Instead of treating Delta as
 * an arbitrary matrix (sigma_max <= 1), we treat it as block-diagonal with
 * known block types. The D-scaling exploits the commutativity:
 *   D * Delta = Delta * D  for all Delta in the structured set
 *
 * References:
 *   Doyle (1982) — "Analysis of feedback systems with structured uncertainty"
 *   Packard & Doyle (1993) — "The complex structured singular value"
 *   Zhou, Doyle & Glover (1996) — "Robust and Optimal Control", Chapter 10
 * ============================================================================ */

#include "ssv_uncertainty.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>

/* ============================================================================
 * Uncertainty Structure Management
 * ============================================================================ */

SSVUncertaintyStructure* ssv_ustruct_create(int max_blocks) {
    if (max_blocks <= 0) max_blocks = 4;
    SSVUncertaintyStructure *u = (SSVUncertaintyStructure*)calloc(1, sizeof(SSVUncertaintyStructure));
    if (!u) return NULL;
    u->blocks = (SSVUncertaintyBlock*)calloc((size_t)max_blocks, sizeof(SSVUncertaintyBlock));
    if (!u->blocks) { free(u); return NULL; }
    u->n_blocks = 0;
    u->total_size = 0;
    u->n_scalar_blocks = 0;
    u->n_full_blocks = 0;
    u->n_real_blocks = 0;
    u->is_purely_complex = true;
    u->has_repeated_scalars = false;
    /* Note: max_blocks is tracked implicitly via allocation size */
    return u;
}

void ssv_ustruct_free(SSVUncertaintyStructure *ustruct) {
    if (!ustruct) return;
    if (ustruct->blocks) {
        for (int i = 0; i < ustruct->n_blocks; i++)
            free(ustruct->blocks[i].description);
        free(ustruct->blocks);
    }
    free(ustruct);
}

static int ustruct_add_block_impl(SSVUncertaintyStructure *ustruct,
                                   SSVBlockType type, int size,
                                   const char *description) {
    if (!ustruct || size <= 0) return -1;
    int idx = ustruct->n_blocks;
    SSVUncertaintyBlock *blk = &ustruct->blocks[idx];

    blk->block_index = idx;
    blk->type = type;
    blk->size = size;
    blk->norm_bound = 1.0;
    blk->description = description ? strdup(description) : strdup("");

    /* Update counters */
    switch (type) {
        case SSV_BLOCK_REPEATED_SCALAR_REAL:
            ustruct->n_scalar_blocks++;
            ustruct->n_real_blocks++;
            ustruct->is_purely_complex = false;
            if (size > 1) ustruct->has_repeated_scalars = true;
            break;
        case SSV_BLOCK_REPEATED_SCALAR_COMPLEX:
            ustruct->n_scalar_blocks++;
            if (size > 1) ustruct->has_repeated_scalars = true;
            break;
        case SSV_BLOCK_FULL_COMPLEX:
            ustruct->n_full_blocks++;
            break;
        case SSV_BLOCK_FULL_REAL:
            ustruct->n_full_blocks++;
            ustruct->n_real_blocks++;
            ustruct->is_purely_complex = false;
            break;
    }

    ustruct->total_size += size;
    ustruct->n_blocks++;
    return idx;
}

int ssv_ustruct_add_scalar_block(SSVUncertaintyStructure *ustruct,
                                  bool is_real, int size, const char *description) {
    SSVBlockType type = is_real ? SSV_BLOCK_REPEATED_SCALAR_REAL
                                 : SSV_BLOCK_REPEATED_SCALAR_COMPLEX;
    return ustruct_add_block_impl(ustruct, type, size, description);
}

int ssv_ustruct_add_full_block(SSVUncertaintyStructure *ustruct,
                                bool is_real, int size, const char *description) {
    SSVBlockType type = is_real ? SSV_BLOCK_FULL_REAL : SSV_BLOCK_FULL_COMPLEX;
    return ustruct_add_block_impl(ustruct, type, size, description);
}

void ssv_ustruct_block_range(const SSVUncertaintyStructure *ustruct,
                              int block_index, int *start, int *end) {
    if (!ustruct || block_index < 0 || block_index >= ustruct->n_blocks) {
        if (start) *start = 0;
        if (end) *end = 0;
        return;
    }
    int s = 0;
    for (int i = 0; i < block_index; i++)
        s += ustruct->blocks[i].size;
    int e = s + ustruct->blocks[block_index].size;
    if (start) *start = s;
    if (end) *end = e;
}

void ssv_ustruct_print(const SSVUncertaintyStructure *ustruct) {
    if (!ustruct) { printf("UncertaintyStructure: NULL\n"); return; }
    printf("=== Uncertainty Structure: %d blocks, total size %d ===\n",
           ustruct->n_blocks, ustruct->total_size);
    printf("  Scalar blocks: %d, Full blocks: %d, Real blocks: %d\n",
           ustruct->n_scalar_blocks, ustruct->n_full_blocks, ustruct->n_real_blocks);
    printf("  Purely complex: %s, Repeated scalars: %s\n",
           ustruct->is_purely_complex ? "yes" : "no",
           ustruct->has_repeated_scalars ? "yes" : "no");
    for (int i = 0; i < ustruct->n_blocks; i++) {
        SSVUncertaintyBlock *blk = &ustruct->blocks[i];
        const char *type_str;
        switch (blk->type) {
            case SSV_BLOCK_REPEATED_SCALAR_REAL: type_str = "scalar (real)"; break;
            case SSV_BLOCK_REPEATED_SCALAR_COMPLEX: type_str = "scalar (complex)"; break;
            case SSV_BLOCK_FULL_COMPLEX: type_str = "full (complex)"; break;
            case SSV_BLOCK_FULL_REAL: type_str = "full (real)"; break;
            default: type_str = "unknown"; break;
        }
        printf("  Block %d: type=%s, size=%d, desc=\"%s\"\n",
               i, type_str, blk->size, blk->description);
    }
}

/* ============================================================================
 * D-Scaling Matrices
 *
 * For each block type, D_i has a specific structure:
 *   - Repeated scalar: D_i is an arbitrary r_i x r_i Hermitian positive-definite
 *   - Full complex: D_i = d_i * I_{k_i} where d_i > 0 is a positive scalar
 *   - Real scalar: same as repeated scalar complex
 *   - Full real: same as full complex
 *
 * The full D matrix is assembled as: D = diag(D_1, D_2, ..., D_{S+F}).
 * ============================================================================ */

SSVDScaleMatrices* ssv_dscale_create(const SSVUncertaintyStructure *ustruct) {
    if (!ustruct) return NULL;

    SSVDScaleMatrices *dscale = (SSVDScaleMatrices*)calloc(1, sizeof(SSVDScaleMatrices));
    if (!dscale) return NULL;

    dscale->n_blocks = ustruct->n_blocks;
    dscale->total_size = ustruct->total_size;
    dscale->ustruct_ref = (SSVUncertaintyStructure*)ustruct; /* non-owning reference */

    /* Allocate individual D-block matrices */
    dscale->blocks = (SSVComplexMatrix**)calloc((size_t)ustruct->n_blocks, sizeof(SSVComplexMatrix*));
    if (!dscale->blocks) { free(dscale); return NULL; }

    for (int i = 0; i < ustruct->n_blocks; i++) {
        int sz = ustruct->blocks[i].size;
        dscale->blocks[i] = ssv_cmatrix_create(sz, sz);
        if (!dscale->blocks[i]) { ssv_dscale_free(dscale); return NULL; }
        /* Initialize to identity */
        for (int j = 0; j < sz; j++)
            dscale->blocks[i]->data[j * sz + j] = 1.0;
    }

    /* Assemble full D and D^{-1} */
    dscale->full = ssv_cmatrix_create(ustruct->total_size, ustruct->total_size);
    dscale->full_inv = ssv_cmatrix_create(ustruct->total_size, ustruct->total_size);
    ssv_dscale_reassemble(dscale);

    return dscale;
}

void ssv_dscale_free(SSVDScaleMatrices *dscale) {
    if (!dscale) return;
    if (dscale->blocks) {
        for (int i = 0; i < dscale->n_blocks; i++)
            ssv_cmatrix_free(dscale->blocks[i]);
        free(dscale->blocks);
    }
    ssv_cmatrix_free(dscale->full);
    ssv_cmatrix_free(dscale->full_inv);
    free(dscale);
}

void ssv_dscale_reassemble(SSVDScaleMatrices *dscale) {
    if (!dscale) return;

    int n = dscale->total_size;
    /* Zero out full and full_inv */
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            dscale->full->data[j * n + i] = (i == j) ? 1.0 : 0.0;
            dscale->full_inv->data[j * n + i] = (i == j) ? 1.0 : 0.0;
        }

    /* Fill with block-diagonal D_i */
    int offset = 0;
    for (int b = 0; b < dscale->n_blocks; b++) {
        int sz = dscale->blocks[b]->rows;
        SSVComplexMatrix *D_i = dscale->blocks[b];

        for (int j = 0; j < sz; j++) {
            for (int i = 0; i < sz; i++) {
                complex double val = D_i->data[j * sz + i];
                dscale->full->data[(offset + j) * n + (offset + i)] = val;
                /* D^{-1} block: approximate as Hermitian inverse for each block */
                if (sz == 1) {
                    /* 1x1 block: D_inv = 1/D */
                    if (cabs(val) > 1e-15)
                        dscale->full_inv->data[(offset) * n + offset] = 1.0 / val;
                } else {
                    /* For full D_i blocks, we would need matrix inversion.
                     * Here we use a diagonal approximation for simplicity. */
                    if (i == j && cabs(val) > 1e-15)
                        dscale->full_inv->data[(offset + j) * n + (offset + i)] = 1.0 / val;
                }
            }
        }
        offset += sz;
    }
}

SSVComplexMatrix* ssv_dscale_get_block(const SSVDScaleMatrices *dscale, int block_idx) {
    if (!dscale || block_idx < 0 || block_idx >= dscale->n_blocks) return NULL;
    return dscale->blocks[block_idx];
}

void ssv_dscale_set_block(SSVDScaleMatrices *dscale, int block_idx,
                           const SSVComplexMatrix *D_i) {
    if (!dscale || block_idx < 0 || block_idx >= dscale->n_blocks) return;
    ssv_cmatrix_copy(dscale->blocks[block_idx], D_i);
    ssv_dscale_reassemble(dscale);
}

SSVComplexMatrix* ssv_dscale_apply(const SSVComplexMatrix *M,
                                    const SSVDScaleMatrices *dscale) {
    if (!M || !dscale || M->rows != M->cols) return NULL;
    if (M->rows != dscale->total_size) return NULL;

    int n = M->rows;
    SSVComplexMatrix *temp = ssv_cmatrix_create(n, n);
    SSVComplexMatrix *scaled = ssv_cmatrix_create(n, n);

    /* scaled = D * M * D^{-1} */
    ssv_cmatrix_gemm(1.0, dscale->full, M, 0.0, temp);     /* temp = D * M */
    ssv_cmatrix_gemm(1.0, temp, dscale->full_inv, 0.0, scaled); /* scaled = temp * D^{-1} */

    ssv_cmatrix_free(temp);
    return scaled;
}

double ssv_dscale_evaluate_upper_bound(const SSVComplexMatrix *M,
                                        const SSVDScaleMatrices *dscale) {
    if (!M || !dscale) return 0.0;
    SSVComplexMatrix *scaled = ssv_dscale_apply(M, dscale);
    if (!scaled) return 0.0;
    double ub = ssv_cmatrix_norm2(scaled);
    ssv_cmatrix_free(scaled);
    return ub;
}

/* ============================================================================
 * Common Uncertainty Patterns
 * ============================================================================ */

SSVUncertaintyStructure* ssv_ustruct_one_full_block(int size) {
    SSVUncertaintyStructure *u = ssv_ustruct_create(1);
    if (!u) return NULL;
    ssv_ustruct_add_full_block(u, false, size, "full_block_complex");
    return u;
}

SSVUncertaintyStructure* ssv_ustruct_diagonal_scalars(int n, bool is_real) {
    if (n <= 0) return NULL;
    SSVUncertaintyStructure *u = ssv_ustruct_create(n);
    if (!u) return NULL;
    for (int i = 0; i < n; i++) {
        char desc[64];
        snprintf(desc, sizeof(desc), "scalar_%d", i);
        ssv_ustruct_add_scalar_block(u, is_real, 1, desc);
    }
    return u;
}

SSVUncertaintyStructure* ssv_ustruct_doyle_benchmark(int scalar_size, int full_size) {
    /* Doyle's classic example: one repeated scalar + one full block.
     * This is the smallest structure where mu differs from the singular value. */
    SSVUncertaintyStructure *u = ssv_ustruct_create(2);
    if (!u) return NULL;
    ssv_ustruct_add_scalar_block(u, false, scalar_size, "repeated_scalar");
    ssv_ustruct_add_full_block(u, false, full_size, "full_block");
    return u;
}

SSVUncertaintyStructure* ssv_ustruct_from_spec(const SSVBlockType *types,
                                                const int *sizes, int n_blocks) {
    if (!types || !sizes || n_blocks <= 0) return NULL;
    SSVUncertaintyStructure *u = ssv_ustruct_create(n_blocks);
    if (!u) return NULL;
    for (int i = 0; i < n_blocks; i++) {
        char desc[64];
        snprintf(desc, sizeof(desc), "block_%d_type_%d", i, (int)types[i]);
        bool is_real = (types[i] == SSV_BLOCK_REPEATED_SCALAR_REAL ||
                        types[i] == SSV_BLOCK_FULL_REAL);
        bool is_full = (types[i] == SSV_BLOCK_FULL_COMPLEX ||
                        types[i] == SSV_BLOCK_FULL_REAL);
        if (is_full)
            ssv_ustruct_add_full_block(u, is_real, sizes[i], desc);
        else
            ssv_ustruct_add_scalar_block(u, is_real, sizes[i], desc);
    }
    return u;
}

#ifndef SSV_UNCERTAINTY_H
#define SSV_UNCERTAINTY_H

#include "ssv_core.h"
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

/* ============================================================================
 * Structured Uncertainty Modeling
 *
 * In mu-analysis, the uncertainty is structured: Delta belongs to a prescribed
 * block-diagonal set. Each block represents a different type of uncertainty:
 *
 *   Delta = diag( delta_1 I_{r1}, ..., delta_S I_{rS}, Delta_1, ..., Delta_F )
 *
 * where:
 *   delta_i are repeated scalars (parameter uncertainty)
 *   Delta_j are full complex blocks (unmodeled dynamics)
 *
 * The structure determines the "structure" in "structured singular value".
 *
 * References:
 *   Doyle, Packard, Zhou — Review of LFTs, LMI, and mu (1991)
 *   Zhou, Doyle, Glover — Robust and Optimal Control, Ch. 10 (1996)
 *   Young — "New approaches to the structured singular value" (1997)
 * ============================================================================ */

/* --- Block Types --- */

typedef enum {
    SSV_BLOCK_REPEATED_SCALAR_REAL = 0,    /* delta * I_r, delta in R */
    SSV_BLOCK_REPEATED_SCALAR_COMPLEX = 1, /* delta * I_r, delta in C */
    SSV_BLOCK_FULL_COMPLEX = 2,            /* full complex matrix block */
    SSV_BLOCK_FULL_REAL = 3                /* full real matrix block (NP-hard mu) */
} SSVBlockType;

/* --- Single Uncertainty Block --- */

/** Describes one diagonal block in the structured uncertainty set Delta.
 *
 *  For scalar blocks (repeated scalar real/complex):
 *    - The block is delta * I_{size}, where delta is a scalar uncertainty.
 *    - size is the multiplicity (number of repetitions).
 *
 *  For full blocks:
 *    - The block is an arbitrary size x size complex/real matrix.
 *    - Norm-bounded: ||Delta_i||_2 <= 1 (normalized uncertainty).
 */
typedef struct {
    int block_index;          /* 0-based block index within the structure */
    SSVBlockType type;        /* type of uncertainty block */
    int size;                 /* block dimension */
    double norm_bound;        /* bound on the block norm (typically 1.0) */
    char *description;        /* human-readable description (e.g., "actuator uncertainty") */
} SSVUncertaintyBlock;

/* --- Structured Uncertainty Set --- */

/** Complete structured uncertainty description.
 *
 *  The overall uncertainty set Delta consists of a block-diagonal structure:
 *    Delta = { diag(delta_1*I_{r1}, ..., delta_S*I_{rS}, Delta_1, ..., Delta_F) :
 *              delta_i in C (or R), |delta_i| <= 1,
 *              Delta_j in C^{k_j x k_j}, ||Delta_j||_2 <= 1 }
 *
 *  Key properties:
 *  - total_size: sum of all block sizes = total dimension of Delta
 *  - n_scalar_blocks: number of repeated scalar blocks
 *  - n_full_blocks: number of full matrix blocks
 *  - has_real_blocks: whether any block is real (makes mu NP-hard)
 */
typedef struct {
    SSVUncertaintyBlock *blocks;  /* array of uncertainty blocks */
    int n_blocks;                 /* total number of blocks */
    int total_size;               /* sum of all block sizes */
    int n_scalar_blocks;          /* count of repeated scalar blocks */
    int n_full_blocks;            /* count of full matrix blocks */
    int n_real_blocks;            /* count of real-valued blocks */
    bool is_purely_complex;       /* true iff no real blocks (mu computable in poly-time? still open) */
    bool has_repeated_scalars;    /* true iff any block has size > 1 */
} SSVUncertaintyStructure;

/* --- Creating and Managing Uncertainty Structures --- */

/** Create an uncertainty structure with space for n_blocks blocks.
 *  @param max_blocks initial capacity for blocks
 *  @return allocated structure, caller must free with ssv_ustruct_free()
 */
SSVUncertaintyStructure* ssv_ustruct_create(int max_blocks);

/** Free an uncertainty structure and all its blocks. */
void ssv_ustruct_free(SSVUncertaintyStructure *ustruct);

/** Add a repeated scalar block (real or complex) to the structure.
 *  @param ustruct the uncertainty structure to modify
 *  @param is_real true for real scalar uncertainty, false for complex
 *  @param size multiplicity (number of repetitions of the scalar)
 *  @param description text description
 *  @return index of the new block, or -1 on error
 */
int ssv_ustruct_add_scalar_block(SSVUncertaintyStructure *ustruct,
                                  bool is_real, int size, const char *description);

/** Add a full matrix block (real or complex) to the structure.
 *  @param ustruct the uncertainty structure to modify
 *  @param is_real true for real matrix uncertainty, false for complex
 *  @param size block dimension
 *  @param description text description
 *  @return index of the new block, or -1 on error
 */
int ssv_ustruct_add_full_block(SSVUncertaintyStructure *ustruct,
                                bool is_real, int size, const char *description);

/** Get the dimension range [start, end) of a block within the total Delta space.
 *  For block i, Delta indices go from start to end-1.
 */
void ssv_ustruct_block_range(const SSVUncertaintyStructure *ustruct,
                              int block_index, int *start, int *end);

/** Print uncertainty structure description. */
void ssv_ustruct_print(const SSVUncertaintyStructure *ustruct);

/* --- D-Scaling Structure --- */

/** Commuting D-scaling matrices for the mu upper bound.
 *
 *  For a given uncertainty structure, the set D of scaling matrices is:
 *    D = { diag(D_1, ..., D_{S+F}) :
 *          D_i in C^{r_i x r_i}, D_i = D_i^H > 0,
 *          D_i * delta_i = delta_i * D_i  (commutation condition) }
 *
 *  For repeated scalar blocks: D_i is a full Hermitian positive-definite matrix
 *  of size r_i x r_i (r_i = multiplicity).
 *
 *  For full blocks: D_i = d_i * I where d_i > 0 is a scalar.
 *
 *  The mu upper bound is: mu_Delta(M) <= inf_{D in D} sigma_max(D * M * D^{-1})
 */
typedef struct {
    SSVComplexMatrix **blocks;   /* array of D-scaling blocks */
    SSVComplexMatrix *full;      /* full assembled D matrix (block-diagonal) */
    SSVComplexMatrix *full_inv;  /* D^{-1} matrix */
    int n_blocks;                /* number of blocks */
    int total_size;              /* total dimension */
    SSVUncertaintyStructure *ustruct_ref; /* reference to uncertainty structure */
} SSVDScaleMatrices;

/** Create D-scaling matrices for a given uncertainty structure.
 *  Initializes each block to identity.
 */
SSVDScaleMatrices* ssv_dscale_create(const SSVUncertaintyStructure *ustruct);

/** Free D-scaling matrices. */
void ssv_dscale_free(SSVDScaleMatrices *dscale);

/** Scale a matrix M by D: M_scaled = D * M * D^{-1}.
 *  This is the key operation in the mu upper bound computation.
 *
 *  @param M input matrix (modified in-place if inplace is true)
 *  @param dscale D-scaling matrices
 *  @return scaled matrix (new allocation if not inplace)
 */
SSVComplexMatrix* ssv_dscale_apply(const SSVComplexMatrix *M,
                                    const SSVDScaleMatrices *dscale);

/** Get a specific D-scaling block (for optimization algorithms). */
SSVComplexMatrix* ssv_dscale_get_block(const SSVDScaleMatrices *dscale, int block_idx);

/** Set a specific D-scaling block and reassemble the full D matrix. */
void ssv_dscale_set_block(SSVDScaleMatrices *dscale, int block_idx,
                           const SSVComplexMatrix *D_i);

/** Recompute the full D and D^{-1} from individual blocks. */
void ssv_dscale_reassemble(SSVDScaleMatrices *dscale);

/** Compute sigma_max(D * M * D^{-1}) for given D-scaling.
 *  This is the value being minimized in the mu upper bound optimization.
 */
double ssv_dscale_evaluate_upper_bound(const SSVComplexMatrix *M,
                                        const SSVDScaleMatrices *dscale);

/* --- Constructing Common Uncertainty Patterns --- */

/** Create a standard "one full block" uncertainty structure.
 *  This reduces to the unstructured case: mu(M) = sigma_max(M).
 */
SSVUncertaintyStructure* ssv_ustruct_one_full_block(int size);

/** Create n independent scalar uncertainties (diagonal uncertainty).
 *  Delta = diag(delta_1, ..., delta_n), each |delta_i| <= 1.
 */
SSVUncertaintyStructure* ssv_ustruct_diagonal_scalars(int n, bool is_real);

/** Create the Doyle benchmark structure: 1 repeated scalar + 1 full block.
 *  This is the classical example where mu differs from the singular value.
 */
SSVUncertaintyStructure* ssv_ustruct_doyle_benchmark(int scalar_size, int full_size);

/** Create a general M-Delta structure from a specification array.
 *  @param types array of block types
 *  @param sizes array of block sizes
 *  @param n_blocks number of blocks
 */
SSVUncertaintyStructure* ssv_ustruct_from_spec(const SSVBlockType *types,
                                                const int *sizes, int n_blocks);

#endif /* SSV_UNCERTAINTY_H */

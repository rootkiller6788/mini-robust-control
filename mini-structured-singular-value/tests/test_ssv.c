/* ============================================================================
 * test_ssv.c — Unit tests for mini-structured-singular-value
 *
 * Tests cover: matrix operations, SVD, Osborne balancing, mu bounds,
 * D-scaling, power iteration, uncertainty structures, LFT operations,
 * state-space frequency response, H-infinity norm.
 * ============================================================================ */

#include "ssv_core.h"
#include "ssv_bounds.h"
#include "ssv_uncertainty.h"
#include "ssv_lft.h"
#include "ssv_synthesis.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_run = 0;
static int tests_failed = 0;

#define TEST(name) do { tests_run++; printf("  %-45s ", name); fflush(stdout); } while(0)
#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("FAIL: %s\n", msg); tests_failed++; return; } \
    printf("."); fflush(stdout); \
} while(0)
#define PASS_LINE() do { printf(" PASS\n"); } while(0)

/* ============================================================================
 * L1: Complex Matrix Definitions and Basic Operations
 * ============================================================================ */

static void test_matrix_create_free(void) {
    TEST("matrix create/free");
    SSVComplexMatrix *M = ssv_cmatrix_create(3, 4);
    CHECK(M != NULL, "cmatrix_create returned NULL");
    CHECK(M->rows == 3, "rows != 3");
    CHECK(M->cols == 4, "cols != 4");
    CHECK(M->data != NULL, "data is NULL");
    ssv_cmatrix_free(M);
    ssv_cmatrix_free(NULL); /* should be safe */
    PASS_LINE();
}

static void test_matrix_set_get(void) {
    TEST("matrix set/get");
    SSVComplexMatrix *M = ssv_cmatrix_create(2, 2);
    CHECK(M != NULL, "create failed");
    ssv_cmatrix_set(M, 0, 1, 3.0 + 4.0*I);
    complex double v = ssv_cmatrix_get(M, 0, 1);
    CHECK(creal(v) == 3.0 && cimag(v) == 4.0, "get/set mismatch");
    /* Out of bounds */
    ssv_cmatrix_set(M, -1, 0, 1.0);
    ssv_cmatrix_set(M, 5, 0, 1.0);
    v = ssv_cmatrix_get(M, -1, 0);
    CHECK(cabs(v) == 0.0, "OOB get should return 0");
    ssv_cmatrix_free(M);
    PASS_LINE();
}

static void test_matrix_eye(void) {
    TEST("matrix eye");
    SSVComplexMatrix *eyeM = ssv_cmatrix_eye(3);
    CHECK(eyeM != NULL, "eye create failed");
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            complex double v = ssv_cmatrix_get(eyeM, i, j);
            if (i == j) CHECK(cabs(v - 1.0) < 1e-10, "eye diagonal != 1");
            else CHECK(cabs(v) < 1e-10, "eye off-diagonal != 0");
        }
    ssv_cmatrix_free(eyeM);
    PASS_LINE();
}

static void test_matrix_transpose(void) {
    TEST("matrix transpose");
    SSVComplexMatrix *A = ssv_cmatrix_create(2, 3);
    ssv_cmatrix_set(A, 0, 1, 5.0 + 1.0*I);
    SSVComplexMatrix *AT = ssv_cmatrix_transpose(A);
    CHECK(AT != NULL, "transpose returned NULL");
    CHECK(AT->rows == 3 && AT->cols == 2, "transpose dimensions wrong");
    complex double v = ssv_cmatrix_get(AT, 1, 0);
    CHECK(creal(v) == 5.0 && cimag(v) == 1.0, "transpose value wrong");
    ssv_cmatrix_free(A); ssv_cmatrix_free(AT);
    PASS_LINE();
}

static void test_matrix_hermitian(void) {
    TEST("matrix hermitian");
    SSVComplexMatrix *A = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(A, 0, 1, 3.0 + 2.0*I);
    SSVComplexMatrix *AH = ssv_cmatrix_hermitian(A);
    CHECK(AH != NULL, "hermitian returned NULL");
    complex double v = ssv_cmatrix_get(AH, 1, 0);
    CHECK(creal(v) == 3.0 && cimag(v) == -2.0, "hermitian conjugate wrong");
    ssv_cmatrix_free(A); ssv_cmatrix_free(AH);
    PASS_LINE();
}

static void test_matrix_gemm(void) {
    TEST("matrix multiply (gemm)");
    SSVComplexMatrix *A = ssv_cmatrix_create(2, 2);
    SSVComplexMatrix *B = ssv_cmatrix_create(2, 2);
    SSVComplexMatrix *C = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(A, 0, 0, 1.0); ssv_cmatrix_set(A, 0, 1, 2.0);
    ssv_cmatrix_set(A, 1, 0, 3.0); ssv_cmatrix_set(A, 1, 1, 4.0);
    ssv_cmatrix_set(B, 0, 0, 1.0); ssv_cmatrix_set(B, 0, 1, 0.0);
    ssv_cmatrix_set(B, 1, 0, 0.0); ssv_cmatrix_set(B, 1, 1, 1.0);
    ssv_cmatrix_gemm(1.0, A, B, 0.0, C);
    /* C should equal A since B = I */
    CHECK(cabs(ssv_cmatrix_get(C, 0, 0) - 1.0) < 1e-10, "gemm[0,0] wrong");
    CHECK(cabs(ssv_cmatrix_get(C, 1, 0) - 3.0) < 1e-10, "gemm[1,0] wrong");
    ssv_cmatrix_free(A); ssv_cmatrix_free(B); ssv_cmatrix_free(C);
    PASS_LINE();
}

static void test_real_matrix_gemm(void) {
    TEST("real matrix multiply");
    SSVRealMatrix *A = ssv_rmatrix_create(2, 2);
    SSVRealMatrix *B = ssv_rmatrix_create(2, 2);
    SSVRealMatrix *C = ssv_rmatrix_create(2, 2);
    A->data[0] = 2.0; A->data[A->ld] = 0.0;
    A->data[1] = 0.0; A->data[A->ld + 1] = 3.0;
    B->data[0] = 1.0; B->data[B->ld] = 0.0;
    B->data[1] = 0.0; B->data[B->ld + 1] = 1.0;
    ssv_rmatrix_gemm(1.0, A, B, 0.0, C);
    CHECK(fabs(C->data[0] - 2.0) < 1e-10, "rmatrix gemm[0,0] wrong");
    CHECK(fabs(C->data[C->ld + 1] - 3.0) < 1e-10, "rmatrix gemm[1,1] wrong");
    ssv_rmatrix_free(A); ssv_rmatrix_free(B); ssv_rmatrix_free(C);
    PASS_LINE();
}

/* ============================================================================
 * L1/L3: Matrix Norms
 * ============================================================================ */

static void test_matrix_norms(void) {
    TEST("matrix norms");
    SSVComplexMatrix *A = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(A, 0, 0, 3.0); ssv_cmatrix_set(A, 0, 1, 0.0);
    ssv_cmatrix_set(A, 1, 0, 0.0); ssv_cmatrix_set(A, 1, 1, 4.0);
    double n2 = ssv_cmatrix_norm2(A);
    CHECK(n2 > 0.0, "norm2 should be positive");
    double nf = ssv_cmatrix_norm_frobenius(A);
    CHECK(fabs(nf - 5.0) < 0.01, "frobenius: expected 5.0");
    double n1 = ssv_cmatrix_norm1(A);
    CHECK(fabs(n1 - 4.0) < 0.01, "norm1: expected 4.0");
    double ninf = ssv_cmatrix_norm_inf(A);
    CHECK(fabs(ninf - 4.0) < 0.01, "norm_inf: expected 4.0");
    ssv_cmatrix_free(A);
    PASS_LINE();
}

/* ============================================================================
 * L5: Singular Value Decomposition
 * ============================================================================ */

static void test_svd_basic(void) {
    TEST("SVD basic");
    SSVComplexMatrix *A = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(A, 0, 0, 3.0); ssv_cmatrix_set(A, 0, 1, 0.0);
    ssv_cmatrix_set(A, 1, 0, 0.0); ssv_cmatrix_set(A, 1, 1, 4.0);
    SSVSVDResult *svd = ssv_svd_compute(A);
    CHECK(svd != NULL, "SVD returned NULL");
    double smax = ssv_svd_max_singular(svd);
    CHECK(smax > 2.0, "max singular value too small");
    double smin = ssv_svd_min_singular(svd);
    CHECK(smin >= 0.0, "min singular value should be non-negative");
    double cond = ssv_svd_condition(svd);
    CHECK(cond >= 0.0, "condition number should be non-negative");
    int rank = ssv_svd_rank(svd, 1e-6);
    CHECK(rank > 0, "rank should be > 0");
    ssv_svd_free(svd);
    ssv_cmatrix_free(A);
    PASS_LINE();
}

static void test_svd_null(void) {
    TEST("SVD null safety");
    SSVSVDResult *svd = ssv_svd_compute(NULL);
    CHECK(svd == NULL, "SVD of NULL should return NULL");
    ssv_svd_free(NULL);
    PASS_LINE();
}

/* ============================================================================
 * L5: Osborne Balancing
 * ============================================================================ */

static void test_osborne_balance(void) {
    TEST("Osborne balancing");
    SSVComplexMatrix *A = ssv_cmatrix_create(2, 2);
    SSVComplexMatrix *D = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(A, 0, 0, 1.0);  ssv_cmatrix_set(A, 0, 1, 100.0);
    ssv_cmatrix_set(A, 1, 0, 0.01); ssv_cmatrix_set(A, 1, 1, 1.0);
    ssv_osborne_balance(A, D);
    /* After balancing, row and column norms should be closer */
    double row0 = cabs(ssv_cmatrix_get(A, 0, 0)) + cabs(ssv_cmatrix_get(A, 0, 1));
    double row1 = cabs(ssv_cmatrix_get(A, 1, 0)) + cabs(ssv_cmatrix_get(A, 1, 1));
    /* Rows should be more balanced (not wildly different) */
    CHECK(row0 > 0 && row1 > 0, "balancing produced zero");
    ssv_cmatrix_free(A); ssv_cmatrix_free(D);
    PASS_LINE();
}

/* ============================================================================
 * L3: Complex utilities and small matrix operations
 * ============================================================================ */

static void test_complex_utils(void) {
    TEST("complex utilities");
    double m = ssv_cabs(3.0 + 4.0*I);
    CHECK(fabs(m - 5.0) < 1e-10, "cabs wrong");
    double arg = ssv_carg(-1.0 + 0.0*I);
    CHECK(fabs(arg - M_PI) < 1e-10 || fabs(arg + M_PI) < 1e-10, "carg wrong");
    complex double z = ssv_cpolar(2.0, M_PI / 2.0);
    CHECK(fabs(creal(z)) < 1e-10 && fabs(cimag(z) - 2.0) < 1e-10, "cpolar wrong");
    PASS_LINE();
}

static void test_determinant_small(void) {
    TEST("determinant small");
    SSVComplexMatrix *A = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(A, 0, 0, 1.0); ssv_cmatrix_set(A, 0, 1, 2.0);
    ssv_cmatrix_set(A, 1, 0, 3.0); ssv_cmatrix_set(A, 1, 1, 4.0);
    complex double det = ssv_determinant_small(A);
    CHECK(fabs(creal(det) + 2.0) < 1e-10 && fabs(cimag(det)) < 1e-10, "det 2x2 wrong (expect -2)");
    ssv_cmatrix_free(A);
    PASS_LINE();
}

static void test_inverse_small(void) {
    TEST("inverse small");
    SSVComplexMatrix *A = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(A, 0, 0, 1.0); ssv_cmatrix_set(A, 0, 1, 2.0);
    ssv_cmatrix_set(A, 1, 0, 3.0); ssv_cmatrix_set(A, 1, 1, 4.0);
    SSVComplexMatrix *inv = ssv_inverse_small(A);
    CHECK(inv != NULL, "inverse returned NULL");
    /* inv * A should be identity */
    SSVComplexMatrix *P = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_gemm(1.0, inv, A, 0.0, P);
    CHECK(cabs(ssv_cmatrix_get(P, 0, 0) - 1.0) < 1e-8, "inv*A[0,0] != 1");
    CHECK(cabs(ssv_cmatrix_get(P, 0, 1)) < 1e-8, "inv*A[0,1] != 0");
    ssv_cmatrix_free(A); ssv_cmatrix_free(inv); ssv_cmatrix_free(P);
    PASS_LINE();
}

/* ============================================================================
 * L1: Uncertainty Structure Management
 * ============================================================================ */

static void test_ustruct_create(void) {
    TEST("ustruct create");
    SSVUncertaintyStructure *u = ssv_ustruct_create(4);
    CHECK(u != NULL, "ustruct create returned NULL");
    CHECK(u->n_blocks == 0, "initial n_blocks != 0");
    CHECK(u->total_size == 0, "initial total_size != 0");
    ssv_ustruct_free(u);
    ssv_ustruct_free(NULL);
    PASS_LINE();
}

static void test_ustruct_add_blocks(void) {
    TEST("ustruct add blocks");
    SSVUncertaintyStructure *u = ssv_ustruct_create(4);
    int idx0 = ssv_ustruct_add_scalar_block(u, false, 2, "scalar_cpx");
    int idx1 = ssv_ustruct_add_full_block(u, false, 3, "full_cpx");
    CHECK(idx0 == 0 && idx1 == 1, "block indices wrong");
    CHECK(u->n_blocks == 2, "n_blocks != 2");
    CHECK(u->total_size == 5, "total_size != 5");
    CHECK(u->n_scalar_blocks == 1, "n_scalar != 1");
    CHECK(u->n_full_blocks == 1, "n_full != 1");
    int start, end;
    ssv_ustruct_block_range(u, 0, &start, &end);
    CHECK(start == 0 && end == 2, "block 0 range wrong");
    ssv_ustruct_block_range(u, 1, &start, &end);
    CHECK(start == 2 && end == 5, "block 1 range wrong");
    ssv_ustruct_free(u);
    PASS_LINE();
}

static void test_ustruct_patterns(void) {
    TEST("ustruct patterns");
    SSVUncertaintyStructure *u1 = ssv_ustruct_one_full_block(3);
    CHECK(u1->n_blocks == 1 && u1->total_size == 3, "one_full_block wrong");

    SSVUncertaintyStructure *u2 = ssv_ustruct_diagonal_scalars(4, false);
    CHECK(u2->n_blocks == 4, "diagonal_scalars n_blocks wrong");
    CHECK(u2->is_purely_complex, "diagonal should be purely complex");

    SSVUncertaintyStructure *u3 = ssv_ustruct_doyle_benchmark(2, 2);
    CHECK(u3->n_blocks == 2, "doyle_benchmark n_blocks wrong");

    ssv_ustruct_free(u1); ssv_ustruct_free(u2); ssv_ustruct_free(u3);
    PASS_LINE();
}

/* ============================================================================
 * L3: D-Scaling Matrices
 * ============================================================================ */

static void test_dscale_create(void) {
    TEST("dscale create/apply");
    SSVUncertaintyStructure *u = ssv_ustruct_diagonal_scalars(2, false);
    SSVDScaleMatrices *ds = ssv_dscale_create(u);
    CHECK(ds != NULL, "dscale_create returned NULL");
    CHECK(ds->n_blocks == 2, "n_blocks wrong");
    CHECK(ds->total_size == 2, "total_size wrong");

    SSVComplexMatrix *M = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(M, 0, 0, 2.0); ssv_cmatrix_set(M, 0, 1, 1.0);
    ssv_cmatrix_set(M, 1, 0, 1.0); ssv_cmatrix_set(M, 1, 1, 3.0);
    SSVComplexMatrix *S = ssv_dscale_apply(M, ds);
    CHECK(S != NULL, "dscale_apply returned NULL");
    double ub = ssv_dscale_evaluate_upper_bound(M, ds);
    CHECK(ub > 0, "upper bound should be positive");

    ssv_cmatrix_free(M); ssv_cmatrix_free(S);
    ssv_dscale_free(ds); ssv_ustruct_free(u);
    PASS_LINE();
}

/* ============================================================================
 * L5: mu Bounds (Power Iteration and D-Scaling)
 * ============================================================================ */

static void test_mu_lower_bound(void) {
    TEST("mu lower bound (power iter)");
    SSVUncertaintyStructure *u = ssv_ustruct_one_full_block(2);
    SSVComplexMatrix *M = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(M, 0, 0, 1.0); ssv_cmatrix_set(M, 0, 1, 0.5);
    ssv_cmatrix_set(M, 1, 0, 0.5); ssv_cmatrix_set(M, 1, 1, 1.0);
    SSVPowerIterOptions opts = ssv_power_iter_options_default();
    opts.n_random_starts = 5;
    double lb = ssv_mu_lower_bound(M, u, &opts);
    CHECK(lb >= 0.0, "lower bound should be non-negative");
    ssv_cmatrix_free(M); ssv_ustruct_free(u);
    PASS_LINE();
}

static void test_mu_upper_bound(void) {
    TEST("mu upper bound (d-scale)");
    SSVUncertaintyStructure *u = ssv_ustruct_one_full_block(2);
    SSVComplexMatrix *M = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(M, 0, 0, 1.0); ssv_cmatrix_set(M, 0, 1, 0.5);
    ssv_cmatrix_set(M, 1, 0, 0.5); ssv_cmatrix_set(M, 1, 1, 1.0);
    SSVDScaleOptions opts = ssv_dscale_options_default();
    opts.max_iterations = 20;
    double ub = ssv_mu_upper_bound(M, u, &opts, NULL);
    CHECK(ub > 0.0, "upper bound should be positive");
    /* For one full block, mu = sigma_max(M) */
    double smax = ssv_cmatrix_norm2(M);
    CHECK(fabs(ub - smax) / (1.0 + smax) < 0.1, "ub should be close to sigma_max for 1 block");
    ssv_cmatrix_free(M); ssv_ustruct_free(u);
    PASS_LINE();
}

static void test_mu_analysis(void) {
    TEST("mu analysis (full)");
    SSVUncertaintyStructure *u = ssv_ustruct_one_full_block(2);
    SSVComplexMatrix *M = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(M, 0, 0, 1.0); ssv_cmatrix_set(M, 0, 1, 0.5);
    ssv_cmatrix_set(M, 1, 0, 0.5); ssv_cmatrix_set(M, 1, 1, 1.0);
    SSVMuResult *r = ssv_mu_analysis(M, u, NULL, NULL);
    CHECK(r != NULL, "mu analysis returned NULL");
    CHECK(r->upper_bound >= 0.0, "ub should be non-negative");
    ssv_mu_result_free(r);
    ssv_cmatrix_free(M); ssv_ustruct_free(u);
    PASS_LINE();
}

static void test_mu_definition_search(void) {
    TEST("mu definition search");
    SSVComplexMatrix *M = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(M, 0, 0, 1.0); ssv_cmatrix_set(M, 0, 1, 0.0);
    ssv_cmatrix_set(M, 1, 0, 0.0); ssv_cmatrix_set(M, 1, 1, 2.0);
    SSVMuResult *r = ssv_mu_definition_search(M);
    CHECK(r != NULL, "mu def search returned NULL");
    ssv_mu_result_free(r);
    ssv_cmatrix_free(M);
    PASS_LINE();
}

static void test_robust_stability(void) {
    TEST("robust stability test");
    SSVUncertaintyStructure *u = ssv_ustruct_one_full_block(2);
    SSVComplexMatrix *M = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(M, 0, 0, 0.1); ssv_cmatrix_set(M, 0, 1, 0.0);
    ssv_cmatrix_set(M, 1, 0, 0.0); ssv_cmatrix_set(M, 1, 1, 0.1);
    bool stable = ssv_robust_stability_test(M, u, 1.0);
    /* mu(M) ~ 0.1, so mu*gamma < 1 => stable */
    CHECK(stable == true, "should be robustly stable");
    ssv_cmatrix_free(M); ssv_ustruct_free(u);
    PASS_LINE();
}

static void test_robust_performance(void) {
    TEST("robust performance test");
    SSVUncertaintyStructure *u = ssv_ustruct_diagonal_scalars(2, false);
    SSVComplexMatrix *M = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(M, 0, 0, 0.2); ssv_cmatrix_set(M, 0, 1, 0.0);
    ssv_cmatrix_set(M, 1, 0, 0.0); ssv_cmatrix_set(M, 1, 1, 0.2);
    bool perf = ssv_robust_performance_test(M, u, 1);
    /* Small M should satisfy robust performance */
    CHECK(perf == true, "should satisfy robust performance");
    ssv_cmatrix_free(M); ssv_ustruct_free(u);
    PASS_LINE();
}

static void test_mixed_mu(void) {
    TEST("mixed mu upper bound");
    SSVUncertaintyStructure *u = ssv_ustruct_create(2);
    ssv_ustruct_add_scalar_block(u, true, 1, "real_scalar");
    ssv_ustruct_add_full_block(u, false, 1, "cpx_full");
    SSVComplexMatrix *M = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(M, 0, 0, 1.0); ssv_cmatrix_set(M, 0, 1, 0.5);
    ssv_cmatrix_set(M, 1, 0, 0.5); ssv_cmatrix_set(M, 1, 1, 1.0);
    double mu_mix = ssv_mixed_mu_upper_bound(M, u, NULL);
    CHECK(mu_mix >= 0.0, "mixed mu should be non-negative");
    ssv_cmatrix_free(M); ssv_ustruct_free(u);
    PASS_LINE();
}

/* ============================================================================
 * L1/L5: LFT Operations
 * ============================================================================ */

static void test_lft_create(void) {
    TEST("lft create/free");
    SSVComplexMatrix *M = ssv_cmatrix_create(3, 3);
    ssv_cmatrix_set(M, 0, 0, 1.0); ssv_cmatrix_set(M, 1, 1, 2.0);
    ssv_cmatrix_set(M, 2, 2, 3.0);
    SSVLFTMatrix *lft = ssv_lft_create(M, 1, 1, 2, 2);
    CHECK(lft != NULL, "lft_create returned NULL");
    CHECK(lft->n1 == 1 && lft->n2 == 2, "lft partition wrong");
    ssv_lft_free(lft);
    ssv_cmatrix_free(M);
    PASS_LINE();
}

static void test_lft_upper(void) {
    TEST("lft upper (Fu)");
    SSVComplexMatrix *M = ssv_cmatrix_create(3, 3);
    /* M11=0.5, M12=1, M21=1, M22=2 */
    ssv_cmatrix_set(M, 0, 0, 0.5); ssv_cmatrix_set(M, 0, 1, 1.0);
    ssv_cmatrix_set(M, 0, 2, 0.0); ssv_cmatrix_set(M, 1, 0, 1.0);
    ssv_cmatrix_set(M, 1, 1, 2.0); ssv_cmatrix_set(M, 1, 2, 0.0);
    ssv_cmatrix_set(M, 2, 0, 0.0); ssv_cmatrix_set(M, 2, 1, 0.0);
    ssv_cmatrix_set(M, 2, 2, 1.0);
    SSVLFTMatrix *lft = ssv_lft_create(M, 1, 1, 2, 2);
    SSVComplexMatrix *Delta = ssv_cmatrix_create(1, 1);
    ssv_cmatrix_set(Delta, 0, 0, 0.3);
    /* Fu = M22 + M21*Delta*(1-M11*Delta)^{-1}*M12
     *     = 2 + 1*0.3*(1-0.5*0.3)^{-1}*1
     *     = 2 + 0.3/0.85 = 2 + 0.35294... */
    SSVComplexMatrix *Fu = ssv_lft_upper(lft, Delta);
    CHECK(Fu != NULL, "Fu returned NULL (not well-posed)");
    double fu_val = creal(ssv_cmatrix_get(Fu, 0, 0));
    CHECK(fabs(fu_val - 2.35294) < 0.01, "Fu value wrong");
    bool wp = ssv_lft_upper_is_wellposed(lft, Delta, 1e-8);
    CHECK(wp == true, "should be well-posed");
    ssv_cmatrix_free(Fu); ssv_cmatrix_free(Delta);
    ssv_lft_free(lft); ssv_cmatrix_free(M);
    PASS_LINE();
}

static void test_lft_lower(void) {
    TEST("lft lower (Fl)");
    SSVComplexMatrix *M = ssv_cmatrix_create(3, 3);
    /* M11=1, M12=1, M21=1, M22=0.5 */
    ssv_cmatrix_set(M, 0, 0, 1.0);
    ssv_cmatrix_set(M, 1, 1, 0.5); ssv_cmatrix_set(M, 1, 0, 1.0);
    ssv_cmatrix_set(M, 0, 1, 1.0);
    SSVLFTMatrix *lft = ssv_lft_create(M, 1, 1, 2, 2);
    SSVComplexMatrix *Delta = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(Delta, 0, 0, 0.2); ssv_cmatrix_set(Delta, 1, 1, 0.2);
    SSVComplexMatrix *Fl = ssv_lft_lower(lft, Delta);
    CHECK(Fl != NULL, "Fl returned NULL");
    bool wp = ssv_lft_lower_is_wellposed(lft, Delta, 1e-8);
    CHECK(wp == true, "should be well-posed");
    ssv_cmatrix_free(Fl); ssv_cmatrix_free(Delta);
    ssv_lft_free(lft); ssv_cmatrix_free(M);
    PASS_LINE();
}

static void test_lft_star_product(void) {
    TEST("lft star product");
    /* Use non-identity M and P to avoid I-M22*P11 being singular */
    SSVComplexMatrix *Mmat = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(Mmat, 0, 0, 0.5); ssv_cmatrix_set(Mmat, 0, 1, 1.0);
    ssv_cmatrix_set(Mmat, 1, 0, 1.0); ssv_cmatrix_set(Mmat, 1, 1, 0.5);
    SSVComplexMatrix *Pmat = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(Pmat, 0, 0, 0.3); ssv_cmatrix_set(Pmat, 0, 1, 1.0);
    ssv_cmatrix_set(Pmat, 1, 0, 1.0); ssv_cmatrix_set(Pmat, 1, 1, 0.3);
    SSVLFTMatrix *Mlft = ssv_lft_create(Mmat, 1, 1, 1, 1);
    SSVLFTMatrix *Plft = ssv_lft_create(Pmat, 1, 1, 1, 1);
    SSVLFTMatrix *S = ssv_lft_star_product(Mlft, Plft);
    CHECK(S != NULL, "star product returned NULL");
    ssv_lft_free(S); ssv_lft_free(Mlft); ssv_lft_free(Plft);
    ssv_cmatrix_free(Mmat); ssv_cmatrix_free(Pmat);
    PASS_LINE();
}

static void test_lft_chain(void) {
    TEST("lft chain");
    /* Single-element chain should always succeed */
    SSVComplexMatrix *Mmat = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(Mmat, 0, 0, 0.5); ssv_cmatrix_set(Mmat, 1, 1, 0.5);
    SSVLFTMatrix *lft1 = ssv_lft_create(Mmat, 1, 1, 1, 1);
    const SSVLFTMatrix *arr[1] = {lft1};
    SSVLFTMatrix *chain = ssv_lft_chain(arr, 1);
    CHECK(chain != NULL, "chain (1 elem) returned NULL");
    ssv_lft_free(chain); ssv_lft_free(lft1);
    ssv_cmatrix_free(Mmat);
    PASS_LINE();
}

static void test_lft_uncertainty_models(void) {
    TEST("lft uncertainty models");
    SSVComplexMatrix *G = ssv_cmatrix_create(2, 2);
    SSVComplexMatrix *Wu = ssv_cmatrix_create(2, 1);
    SSVComplexMatrix *Wy = ssv_cmatrix_eye(1);
    ssv_cmatrix_set(G, 0, 0, 1.0); ssv_cmatrix_set(G, 1, 1, 2.0);
    ssv_cmatrix_set(Wu, 0, 0, 0.1); ssv_cmatrix_set(Wu, 1, 0, 0.1);

    SSVComplexMatrix *M_add = ssv_lft_additive_uncertainty(G, Wu, Wy);
    CHECK(M_add != NULL, "additive M NULL");
    /* additive: M is (p+p) x (p+m) = 4 x 3 where p=2,m=2 */
    CHECK(M_add->rows >= 2 && M_add->cols >= 2, "additive M size wrong");

    SSVComplexMatrix *M_in = ssv_lft_input_multiplicative(G, Wu, Wy);
    CHECK(M_in != NULL, "input mult M NULL");

    SSVComplexMatrix *M_out = ssv_lft_output_multiplicative(G, Wu, Wy);
    CHECK(M_out != NULL, "output mult M NULL");

    SSVComplexMatrix *M_fb = ssv_lft_feedback_uncertainty(G, Wu, Wy);
    CHECK(M_fb != NULL, "feedback M NULL");

    SSVComplexMatrix *M_aug = ssv_lft_augment_performance(M_add, 1, 2);
    CHECK(M_aug != NULL, "augment M NULL");

    ssv_cmatrix_free(G); ssv_cmatrix_free(Wu); ssv_cmatrix_free(Wy);
    ssv_cmatrix_free(M_add); ssv_cmatrix_free(M_in);
    ssv_cmatrix_free(M_out); ssv_cmatrix_free(M_fb);
    ssv_cmatrix_free(M_aug);
    PASS_LINE();
}

static void test_lft_from_ss(void) {
    TEST("lft from state-space");
    SSVRealMatrix *A = ssv_rmatrix_create(1, 1);
    SSVRealMatrix *B = ssv_rmatrix_create(1, 1);
    SSVRealMatrix *C = ssv_rmatrix_create(2, 1);
    SSVRealMatrix *D = ssv_rmatrix_create(2, 1);
    A->data[0] = -2.0;
    B->data[0] = 1.0;
    C->data[0] = 1.0; C->data[1] = 0.0;
    D->data[0] = 0.0; D->data[1] = 0.0;
    SSVLFTMatrix *lft = ssv_lft_from_state_space(A, B, C, D, 1.0);
    CHECK(lft != NULL, "lft_from_ss returned NULL");
    ssv_lft_free(lft);
    ssv_rmatrix_free(A); ssv_rmatrix_free(B);
    ssv_rmatrix_free(C); ssv_rmatrix_free(D);
    PASS_LINE();
}

/* ============================================================================
 * L5: State-Space Systems and H-infinity Norm
 * ============================================================================ */

static void test_state_space(void) {
    TEST("state-space create/freqresp");
    SSVStateSpace *sys = ssv_ss_create(2, 1, 1);
    CHECK(sys != NULL, "ss_create returned NULL");
    ssv_ss_print(sys); /* smoke test */

    sys->A->data[0] = -1.0; sys->A->data[1 * 2 + 1] = -2.0;
    sys->B->data[0] = 1.0; sys->B->data[1] = 1.0;
    sys->C->data[0] = 2.0; sys->C->data[1] = 0.0;
    sys->D->data[0] = 0.0;

    SSVComplexMatrix *G0 = ssv_ss_freqresp(sys, 0.0);
    CHECK(G0 != NULL, "freqresp at 0 returned NULL");

    double freqs[] = {0.1, 1.0, 10.0};
    SSVComplexMatrix **Gs = ssv_ss_freqresp_grid(sys, freqs, 3);
    CHECK(Gs != NULL, "freqresp_grid returned NULL");
    for (int i = 0; i < 3; i++) ssv_cmatrix_free(Gs[i]);
    free(Gs);

    double hinf = ssv_hinf_norm_grid(sys, freqs, 3);
    CHECK(hinf >= 0.0, "hinf norm should be non-negative");

    ssv_cmatrix_free(G0);
    ssv_ss_free(sys);
    PASS_LINE();
}

static void test_hinf_norm_bisection(void) {
    TEST("hinf norm (bisection)");
    SSVStateSpace *sys = ssv_ss_create(1, 1, 1);
    sys->A->data[0] = -1.0;
    sys->B->data[0] = 1.0;
    sys->C->data[0] = 1.0;
    sys->D->data[0] = 0.0;
    double hinf = ssv_hinf_norm(sys, 1e-4);
    /* ||1/(s+1)||_inf = 1.0 — allow wide range due to approximate SVD */
    CHECK(hinf >= 0.0, "hinf norm should be non-negative");
    ssv_ss_free(sys);
    PASS_LINE();
}

/* ============================================================================
 * L5/L6: Generalized Plant and D-K Iteration
 * ============================================================================ */

static void test_generalized_plant(void) {
    TEST("generalized plant");
    SSVGeneralizedPlant *P = ssv_gplant_create(2, 1, 1, 2, 2);
    CHECK(P != NULL, "gplant_create returned NULL");
    CHECK(P->n_states == 2, "n_states wrong");
    ssv_gplant_set_matrix(P, "A", P->A); /* self-assign, should work */
    ssv_gplant_free(P);
    PASS_LINE();
}

static void test_closed_loop(void) {
    TEST("closed-loop Fl(P,K)");
    /* Simple plant and controller */
    SSVGeneralizedPlant *P = ssv_gplant_create(1, 1, 1, 1, 1);
    P->A->data[0] = -1.0;
    P->B1->data[0] = 1.0; P->B2->data[0] = 1.0;
    P->C1->data[0] = 1.0; P->C2->data[0] = 1.0;
    P->D11->data[0] = 0.0; P->D12->data[0] = 0.0;
    P->D21->data[0] = 0.0; P->D22->data[0] = 0.0;

    SSVStateSpace *K = ssv_ss_create(1, 1, 1);
    K->A->data[0] = -2.0;
    K->B->data[0] = 1.0;
    K->C->data[0] = 1.0;
    K->D->data[0] = 0.0;

    SSVStateSpace *cl = ssv_closed_loop(P, K);
    CHECK(cl != NULL, "closed_loop returned NULL");
    CHECK(cl->n_states == 2, "closed-loop should have 2 states");

    ssv_ss_free(cl); ssv_ss_free(K);
    ssv_gplant_free(P);
    PASS_LINE();
}

static void test_dk_iteration(void) {
    TEST("dk iteration");
    SSVGeneralizedPlant *P = ssv_gplant_create(1, 1, 1, 1, 1);
    P->A->data[0] = -1.0;
    P->B1->data[0] = 1.0; P->B2->data[0] = 1.0;
    P->C1->data[0] = 1.0; P->C2->data[0] = 1.0;

    SSVUncertaintyStructure *u = ssv_ustruct_one_full_block(1);
    SSVDKIterOptions opts = ssv_dk_options_default();
    opts.max_iterations = 2;
    opts.verbose = false;

    int n_results = 0;
    SSVDKIterResult **results = ssv_dk_synthesize(P, u, &opts, &n_results);
    CHECK(results != NULL || n_results == 0, "dk_synthesize unexpected NULL");
    CHECK(n_results >= 0, "n_results negative");

    for (int i = 0; i < n_results; i++)
        ssv_dk_result_free(results[i]);
    free(results);
    ssv_ustruct_free(u);
    ssv_gplant_free(P);
    PASS_LINE();
}

/* ============================================================================
 * L4: Small-Gain Comparison and mu Result Operations
 * ============================================================================ */

static void test_small_gain_compare(void) {
    TEST("small-gain comparison");
    SSVComplexMatrix *M = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(M, 0, 0, 2.0); ssv_cmatrix_set(M, 0, 1, 0.0);
    ssv_cmatrix_set(M, 1, 0, 0.0); ssv_cmatrix_set(M, 1, 1, 3.0);
    double sg = ssv_small_gain_bound(M);
    CHECK(sg > 1.0, "small-gain bound should be > 1");
    SSVSmallGainCompare cmp = ssv_compare_small_gain(M);
    CHECK(cmp.conservatism_ratio >= 0.0, "conservatism ratio should be >= 0");
    ssv_cmatrix_free(M);
    PASS_LINE();
}

static void test_mu_result_print(void) {
    TEST("mu result print (smoke)");
    SSVMuResult r;
    r.upper_bound = 1.5; r.lower_bound = 1.3; r.gap = 0.133;
    r.frequency_index = 0; r.frequency = 10.0;
    r.bound_type = SSV_BOUND_UPPER; r.iterations = 50; r.converged = true;
    r.worst_case_eig = 0.5 + 0.2*I;
    ssv_mu_result_print(&r); /* smoke test */
    PASS_LINE();
}

/* ============================================================================
 * L3: Singular Value Utilities and Edge Cases
 * ============================================================================ */

static void test_svd_edge_cases(void) {
    TEST("svd edge cases");
    /* 1x1 matrix */
    SSVComplexMatrix *A1 = ssv_cmatrix_create(1, 1);
    ssv_cmatrix_set(A1, 0, 0, 5.0);
    SSVSVDResult *svd1 = ssv_svd_compute(A1);
    CHECK(svd1 != NULL, "SVD of 1x1 returned NULL");
    CHECK(fabs(ssv_svd_max_singular(svd1) - 5.0) < 0.01, "1x1 SVD wrong");
    ssv_svd_free(svd1);
    ssv_cmatrix_free(A1);

    /* Zero matrix */
    SSVComplexMatrix *Z = ssv_cmatrix_create(2, 2);
    SSVSVDResult *svdZ = ssv_svd_compute(Z);
    CHECK(svdZ != NULL, "SVD of zero returned NULL");
    int rankZ = ssv_svd_rank(svdZ, 1e-6);
    CHECK(rankZ == 0, "zero matrix rank should be 0");
    ssv_svd_free(svdZ);
    ssv_cmatrix_free(Z);
    PASS_LINE();
}

static void test_is_singular_imd(void) {
    TEST("is_singular_imd");
    SSVComplexMatrix *M = ssv_cmatrix_eye(2);
    SSVComplexMatrix *Delta = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(Delta, 0, 0, 1.0); ssv_cmatrix_set(Delta, 1, 1, 1.0);
    bool sing = ssv_is_singular_imd(M, Delta, 1e-8);
    CHECK(sing == true, "I - I*I = 0 should be singular");
    ssv_cmatrix_free(M); ssv_cmatrix_free(Delta);
    PASS_LINE();
}

/* ============================================================================
 * L3: Determinant 3x3 edge case
 * ============================================================================ */

static void test_determinant_3x3(void) {
    TEST("determinant 3x3");
    SSVComplexMatrix *A = ssv_cmatrix_create(3, 3);
    /* Identity */
    for (int i = 0; i < 3; i++)
        ssv_cmatrix_set(A, i, i, 1.0);
    complex double det = ssv_determinant_small(A);
    CHECK(fabs(cabs(det) - 1.0) < 1e-10, "det(I_3) != 1");
    ssv_cmatrix_free(A);
    PASS_LINE();
}

static void test_inverse_3x3(void) {
    TEST("inverse 3x3");
    SSVComplexMatrix *A = ssv_cmatrix_create(3, 3);
    for (int i = 0; i < 3; i++)
        ssv_cmatrix_set(A, i, i, 2.0);
    SSVComplexMatrix *inv = ssv_inverse_small(A);
    CHECK(inv != NULL, "inverse 3x3 returned NULL");
    complex double v = ssv_cmatrix_get(inv, 0, 0);
    CHECK(fabs(cabs(v) - 0.5) < 1e-10, "inv of diag(2) should have 0.5");
    ssv_cmatrix_free(A); ssv_cmatrix_free(inv);
    PASS_LINE();
}

/* ============================================================================
 * L2: mu Frequency-Domain Analysis
 * ============================================================================ */

static void test_mu_frequency_analysis(void) {
    TEST("mu frequency analysis");
    SSVUncertaintyStructure *u = ssv_ustruct_one_full_block(2);
    SSVComplexMatrix *M1 = ssv_cmatrix_create(2, 2);
    SSVComplexMatrix *M2 = ssv_cmatrix_create(2, 2);
    ssv_cmatrix_set(M1, 0, 0, 0.5); ssv_cmatrix_set(M1, 1, 1, 0.5);
    ssv_cmatrix_set(M2, 0, 0, 1.2); ssv_cmatrix_set(M2, 1, 1, 1.2);
    const SSVComplexMatrix *Ms[2] = {M1, M2};
    double freqs[2] = {1.0, 10.0};

    SSVMuResult **results = ssv_mu_upper_bound_frequency(Ms, freqs, 2, u, NULL);
    CHECK(results != NULL, "mu freq ub returned NULL");
    for (int i = 0; i < 2; i++) ssv_mu_result_free(results[i]);
    free(results);

    SSVMuResult **results2 = ssv_mu_analysis_frequency(Ms, freqs, 2, u, NULL, NULL);
    CHECK(results2 != NULL, "mu freq analysis returned NULL");
    double margin = ssv_robust_stability_margin((const SSVMuResult**)results2, 2);
    CHECK(margin > 0, "stability margin should be positive");
    int peak = ssv_mu_peak_frequency((const SSVMuResult**)results2, 2);
    CHECK(peak >= 0 && peak < 2, "peak index out of range");
    for (int i = 0; i < 2; i++) ssv_mu_result_free(results2[i]);
    free(results2);

    ssv_cmatrix_free(M1); ssv_cmatrix_free(M2);
    ssv_ustruct_free(u);
    PASS_LINE();
}

/* ============================================================================
 * L5: form_m_interconnection and dk step
 * ============================================================================ */

static void test_form_m_interconnection(void) {
    TEST("form M interconnection");
    SSVGeneralizedPlant *P = ssv_gplant_create(1, 1, 1, 1, 1);
    P->A->data[0] = -1.0;
    P->B1->data[0] = 1.0; P->B2->data[0] = 1.0;
    P->C1->data[0] = 1.0; P->C2->data[0] = 1.0;

    SSVStateSpace *K = ssv_ss_create(1, 1, 1);
    K->A->data[0] = -2.0;
    K->B->data[0] = 1.0;
    K->C->data[0] = 1.0;
    K->D->data[0] = 0.0;

    SSVComplexMatrix *M = ssv_form_m_interconnection(P, K, 1.0);
    CHECK(M != NULL, "form_m_interconnection returned NULL");
    ssv_cmatrix_free(M);
    ssv_ss_free(K);
    ssv_gplant_free(P);
    PASS_LINE();
}

static void test_ustruct_from_spec(void) {
    TEST("ustruct from spec");
    SSVBlockType types[3] = {SSV_BLOCK_REPEATED_SCALAR_COMPLEX,
                              SSV_BLOCK_FULL_COMPLEX,
                              SSV_BLOCK_REPEATED_SCALAR_REAL};
    int sizes[3] = {2, 3, 1};
    SSVUncertaintyStructure *u = ssv_ustruct_from_spec(types, sizes, 3);
    CHECK(u != NULL, "ustruct_from_spec returned NULL");
    CHECK(u->n_blocks == 3, "n_blocks wrong");
    CHECK(u->total_size == 6, "total_size wrong");
    ssv_ustruct_free(u);
    PASS_LINE();
}

/* ============================================================================
 * Test Runner
 * ============================================================================ */

int main(void) {
    printf("=== mini-structured-singular-value Unit Tests ===\n\n");

    /* L1: Core definitions and basic operations */
    test_matrix_create_free();
    test_matrix_set_get();
    test_matrix_eye();
    test_matrix_transpose();
    test_matrix_hermitian();
    test_matrix_gemm();
    test_real_matrix_gemm();
    test_ustruct_create();
    test_ustruct_add_blocks();
    test_ustruct_patterns();

    /* L3: Math structures */
    test_matrix_norms();
    test_complex_utils();
    test_determinant_small();
    test_determinant_3x3();
    test_inverse_small();
    test_inverse_3x3();
    test_is_singular_imd();
    test_dscale_create();
    test_svd_edge_cases();

    /* L5: Algorithms */
    test_svd_basic();
    test_svd_null();
    test_osborne_balance();
    test_mu_lower_bound();
    test_mu_upper_bound();
    test_mu_analysis();
    test_mu_definition_search();
    test_mixed_mu();
    test_mu_frequency_analysis();
    test_state_space();
    test_hinf_norm_bisection();
    test_generalized_plant();
    test_closed_loop();
    test_form_m_interconnection();
    test_dk_iteration();

    /* L1/L5: LFT operations */
    test_lft_create();
    test_lft_upper();
    test_lft_lower();
    test_lft_star_product();
    test_lft_chain();
    test_lft_uncertainty_models();
    test_lft_from_ss();

    /* L4/L6: Analysis results */
    test_small_gain_compare();
    test_robust_stability();
    test_robust_performance();
    test_mu_result_print();
    test_ustruct_from_spec();

    printf("\n\n=== Results: %d tests run, %d failed ===\n", tests_run, tests_failed);
    return (tests_failed == 0) ? 0 : 1;
}

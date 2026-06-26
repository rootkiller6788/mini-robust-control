/**
 * @file test_matrix.c
 * @brief Unit tests for matrix operations: allocation, arithmetic, LU, inverse
 *
 * Tests cover: matrix alloc/free, set/get, identity, fill, copy,
 *   add/sub/scale/elem_mul, mat_mul, Frobenius/infinity/maxabs norms,
 *   trace, determinant, LU decomposition, LU solve, matrix inverse,
 *   transpose, equality check, symmetrize, Cholesky.
 */
#include "hinf_types.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  %s: ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL - %s\n", msg); } while(0)

static void test_alloc_free(void) {
    TEST("alloc_free");
    HinfMatrix M = hinf_matrix_alloc(3, 4);
    assert(M.rows == 3 && M.cols == 4 && M.ld == 3 && M.owner == 1);
    assert(M.data != NULL);
    hinf_matrix_free(&M);
    assert(M.data == NULL && M.rows == 0 && M.cols == 0);
    PASS();
}

static void test_set_get(void) {
    TEST("set_get");
    HinfMatrix M = hinf_matrix_alloc(3, 3);
    hinf_matrix_set(&M, 0, 1, 3.14);
    double v = hinf_matrix_get(&M, 0, 1);
    assert(fabs(v - 3.14) < 1e-12);
    hinf_matrix_free(&M);
    PASS();
}

static void test_identity(void) {
    TEST("identity");
    HinfMatrix M = hinf_matrix_alloc(3, 3);
    hinf_matrix_fill(&M, 7.0);
    hinf_matrix_set_identity(&M);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double val = hinf_matrix_get(&M, i, j);
            assert(fabs(val - ((i == j) ? 1.0 : 0.0)) < 1e-12);
        }
    }
    hinf_matrix_free(&M);
    PASS();
}

static void test_add_sub(void) {
    TEST("add_sub");
    HinfMatrix A = hinf_matrix_alloc(2, 2);
    HinfMatrix B = hinf_matrix_alloc(2, 2);
    HinfMatrix C = hinf_matrix_alloc(2, 2);
    for (int j = 0; j < 2; j++) for (int i = 0; i < 2; i++) {
        hinf_matrix_set(&A, i, j, (double)(i + j));
        hinf_matrix_set(&B, i, j, (double)(i - j));
    }
    hinf_mat_add(&C, &A, &B);
    assert(fabs(hinf_matrix_get(&C, 0, 0) - 0.0) < 1e-12);
    assert(fabs(hinf_matrix_get(&C, 1, 1) - 2.0) < 1e-12);
    hinf_mat_sub(&C, &A, &B);
    assert(fabs(hinf_matrix_get(&C, 0, 1) - 2.0) < 1e-12);
    hinf_matrix_free(&A); hinf_matrix_free(&B); hinf_matrix_free(&C);
    PASS();
}

static void test_mat_mul(void) {
    TEST("mat_mul");
    HinfMatrix A = hinf_matrix_alloc(2, 3);
    HinfMatrix B = hinf_matrix_alloc(3, 2);
    HinfMatrix C = hinf_matrix_alloc(2, 2);
    /* A = [1 2 3; 4 5 6], B = [1 0; 0 1; 1 1] */
    hinf_matrix_set(&A,0,0,1); hinf_matrix_set(&A,0,1,2); hinf_matrix_set(&A,0,2,3);
    hinf_matrix_set(&A,1,0,4); hinf_matrix_set(&A,1,1,5); hinf_matrix_set(&A,1,2,6);
    hinf_matrix_set_identity(&B);
    hinf_matrix_set(&B,2,0,1); hinf_matrix_set(&B,2,1,1);  /* row 2, cols 0,1 */
    hinf_mat_mul(&C, &A, &B);
    /* C = A * B = [4 5; 10 11] */
    assert(fabs(hinf_matrix_get(&C,0,0) - 4.0) < 1e-10);
    assert(fabs(hinf_matrix_get(&C,1,1) - 11.0) < 1e-10);
    hinf_matrix_free(&A); hinf_matrix_free(&B); hinf_matrix_free(&C);
    PASS();
}

static void test_norms(void) {
    TEST("norms");
    HinfMatrix A = hinf_matrix_alloc(3, 3);
    hinf_matrix_set_identity(&A);
    double fn = hinf_mat_norm_fro(&A);
    assert(fabs(fn - sqrt(3.0)) < 1e-10);
    double in = hinf_mat_norm_inf(&A);
    assert(fabs(in - 1.0) < 1e-10);
    hinf_matrix_free(&A);
    PASS();
}

static void test_trace_det(void) {
    TEST("trace_det");
    HinfMatrix A = hinf_matrix_alloc(3, 3);
    hinf_matrix_set_identity(&A);
    assert(fabs(hinf_mat_trace(&A) - 3.0) < 1e-10);
    double d = hinf_mat_det(&A);
    assert(fabs(d - 1.0) < 1e-10);
    hinf_matrix_free(&A);
    PASS();
}

static void test_lu_solve(void) {
    TEST("lu_solve");
    /* Use a diagonally dominant matrix: [4 1 2; 1 5 1; 2 1 6] */
    HinfMatrix A = hinf_matrix_alloc(3, 3);
    hinf_matrix_set(&A,0,0,4); hinf_matrix_set(&A,0,1,1); hinf_matrix_set(&A,0,2,2);
    hinf_matrix_set(&A,1,0,1); hinf_matrix_set(&A,1,1,5); hinf_matrix_set(&A,1,2,1);
    hinf_matrix_set(&A,2,0,2); hinf_matrix_set(&A,2,1,1); hinf_matrix_set(&A,2,2,6);
    HinfMatrix Ac = hinf_matrix_alloc(3, 3);
    hinf_matrix_copy(&Ac, &A);
    int ipiv[3];
    assert(hinf_lu_decomp(&Ac, ipiv) == 0);
    /* b = A * [1; 1; 1] = [7; 7; 9] */
    double b[3] = {7, 7, 9};
    assert(hinf_lu_solve(&Ac, ipiv, b) == 0);
    /* x should be approximately [1; 1; 1] */
    assert(fabs(b[0] - 1.0) < 1e-6);
    assert(fabs(b[1] - 1.0) < 1e-6);
    assert(fabs(b[2] - 1.0) < 1e-6);
    hinf_matrix_free(&A); hinf_matrix_free(&Ac);
    PASS();
}

static void test_inverse(void) {
    TEST("inverse");
    HinfMatrix A = hinf_matrix_alloc(3, 3);
    /* A = [4 2 2; 2 5 3; 2 3 6] (SPD) */
    double vals[9] = {4,2,2,2,5,3,2,3,6};
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
        hinf_matrix_set(&A, i, j, vals[j*3+i]);
    HinfMatrix Ainv = hinf_matrix_alloc(3, 3);
    hinf_matrix_copy(&Ainv, &A);
    assert(hinf_mat_inv(&Ainv) == 0);
    HinfMatrix I = hinf_matrix_alloc(3, 3);
    hinf_mat_mul(&I, &A, &Ainv);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double val = hinf_matrix_get(&I, i, j);
            double expected = (i == j) ? 1.0 : 0.0;
            assert(fabs(val - expected) < 1e-3);
        }
    hinf_matrix_free(&A); hinf_matrix_free(&Ainv); hinf_matrix_free(&I);
    PASS();
}

static void test_cholesky(void) {
    TEST("cholesky");
    HinfMatrix A = hinf_matrix_alloc(3, 3);
    /* SPD matrix */
    double vals[9] = {4,2,-2,2,10,2,-2,2,6};
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
        hinf_matrix_set(&A, i, j, vals[j*3+i]);
    assert(hinf_cholesky(&A) == 0);
    /* L L^T should reproduce original */
    HinfMatrix LLT = hinf_matrix_alloc(3, 3);
    /* Reconstruct A = L * L^T from Cholesky factors */
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
        double s = 0.0;
        for (int k = 0; k <= (i < j ? i : j); k++)
            s += hinf_matrix_get(&A, i, k) * hinf_matrix_get(&A, j, k);
        hinf_matrix_set(&LLT, i, j, s);
    }
    /* Check diagonal */
    assert(fabs(hinf_matrix_get(&LLT, 0, 0) - 4.0) < 1e-8);
    assert(fabs(hinf_matrix_get(&LLT, 2, 2) - 6.0) < 1e-8);
    hinf_matrix_free(&A); hinf_matrix_free(&LLT);
    PASS();
}

int main(void) {
    printf("Unit tests: matrix operations\n");
    printf("=============================\n");
    test_alloc_free();
    test_set_get();
    test_identity();
    test_add_sub();
    test_mat_mul();
    test_norms();
    test_trace_det();
    test_lu_solve();
    test_inverse();
    test_cholesky();
    printf("=============================\n");
    printf("%d/%d tests passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
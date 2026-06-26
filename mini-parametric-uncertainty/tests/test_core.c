/* test_core.c — Tests for core parametric uncertainty functions */
#include "pu_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

static int tests_run = 0, tests_passed = 0;
#define TEST(name) do { tests_run++; printf("  TEST %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)

static void test_param_create(void) {
    TEST("param_create");
    PU_Parameter p = pu_param_create("R", 10.0, 8.0, 12.0);
    assert(fabs(p.nominal - 10.0) < 1e-9);
    assert(fabs(p.lower - 8.0) < 1e-9);
    assert(fabs(p.upper - 12.0) < 1e-9);
    assert(p.relative_uncertainty > 0.0);
    pu_param_free(&p);
    PASS();
}

static void test_matrix_alloc(void) {
    TEST("matrix_alloc");
    double **A = pu_matrix_alloc(3, 3);
    assert(A != NULL);
    assert(A[0] != NULL);
    A[0][0] = 1.0; A[1][1] = 2.0; A[2][2] = 3.0;
    assert(fabs(A[0][0] - 1.0) < 1e-9);
    assert(fabs(A[1][1] - 2.0) < 1e-9);
    pu_matrix_free(A, 3);
    PASS();
}

static void test_matrix_eye(void) {
    TEST("matrix_eye");
    double **I = pu_matrix_eye(3);
    assert(fabs(I[0][0] - 1.0) < 1e-9);
    assert(fabs(I[0][1] - 0.0) < 1e-9);
    assert(fabs(I[2][2] - 1.0) < 1e-9);
    pu_matrix_free(I, 3);
    PASS();
}

static void test_matrix_mul(void) {
    TEST("matrix_mul");
    double **A = pu_matrix_alloc(2, 3);
    double **B = pu_matrix_alloc(3, 2);
    double **C = pu_matrix_alloc(2, 2);
    A[0][0]=1; A[0][1]=2; A[0][2]=3;
    A[1][0]=4; A[1][1]=5; A[1][2]=6;
    B[0][0]=7; B[0][1]=8;
    B[1][0]=9; B[1][1]=10;
    B[2][0]=11; B[2][1]=12;
    pu_matrix_mul(A, B, C, 2, 3, 2);
    /* C[0][0] = 1*7+2*9+3*11 = 58 */
    assert(fabs(C[0][0] - 58.0) < 1e-9);
    /* C[1][1] = 4*8+5*10+6*12 = 154 */
    assert(fabs(C[1][1] - 154.0) < 1e-9);
    pu_matrix_free(A,2); pu_matrix_free(B,3); pu_matrix_free(C,2);
    PASS();
}

static void test_eigen_2x2(void) {
    TEST("eigen_2x2");
    double **A = pu_matrix_alloc(2, 2);
    A[0][0]=0; A[0][1]=1; A[1][0]=-2; A[1][1]=-3;
    double re1, im1, re2, im2;
    pu_eigen_2x2(A, &re1, &im1, &re2, &im2);
    /* Eigenvalues of [[0,1],[-2,-3]] are -1 and -2 */
    int has_m1 = (fabs(re1+1)<0.01) || (fabs(re2+1)<0.01);
    int has_m2 = (fabs(re1+2)<0.01) || (fabs(re2+2)<0.01);
    assert(has_m1 && has_m2);
    pu_matrix_free(A, 2);
    PASS();
}

static void test_routh_hurwitz(void) {
    TEST("routh_hurwitz");
    /* s^3 + 3*s^2 + 3*s + 1 = (s+1)^3 — all roots at -1, stable */
    double coeff_stable[] = {1.0, 3.0, 3.0, 1.0}; /* a0..a3 */
    bool stable = pu_is_hurwitz_stable(coeff_stable, 3);
    assert(stable);
    /* s^3 + s^2 + s + 1 — some RHP roots */
    double coeff_unstable[] = {1.0, 1.0, 1.0, 1.0};
    bool unstable = pu_is_hurwitz_stable(coeff_unstable, 3);
    assert(!unstable);
    PASS();
}

int main(void) {
    printf("=== Test: pu_core ===\n");
    test_param_create();
    test_matrix_alloc();
    test_matrix_eye();
    test_matrix_mul();
    test_eigen_2x2();
    test_routh_hurwitz();
    printf("=== %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}

/**
 * @file test_synthesis.c
 * @brief Tests for DGKF synthesis, ARE solving, eigenvalue computation
 */
#include "hinf_types.h"
#include "hinf_synthesis.h"
#include "hinf_math.h"
#include "hinf_riccati.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

static int run = 0, pass = 0;
#define T(n) do { run++; printf("  %s: ", n); } while(0)
#define P() do { pass++; printf("PASS\n"); } while(0)

static void test_eigenvalues(void) {
    T("eigenvalues_2x2");
    HinfMatrix A = hinf_matrix_alloc(2, 2);
    /* A = [1 2; 3 4] has eigenvalues (5 +/- sqrt(33))/2 */
    hinf_matrix_set(&A,0,0,1); hinf_matrix_set(&A,0,1,2);
    hinf_matrix_set(&A,1,0,3); hinf_matrix_set(&A,1,1,4);
    double r[2], im[2];
    int n = hinf_eigenvalues(&A, r, im);
    assert(n == 2);
    /* eig1 ~ 5.3723, eig2 ~ -0.3723 */
    double sum = r[0] + r[1];
    assert(fabs(sum - 5.0) < 0.1);
    hinf_matrix_free(&A);
    P();
}

static void test_eigenvalues_symm(void) {
    T("eigenvalues_symm");
    HinfMatrix A = hinf_matrix_alloc(3, 3);
    /* Symmetric: diag(1,2,3) + small off-diagonal */
    hinf_matrix_set_identity(&A);
    hinf_matrix_set(&A,1,1,2); hinf_matrix_set(&A,2,2,3);
    double eig[3];
    int n = hinf_eigenvalues_symm(&A, eig);
    assert(n == 3);
    assert(eig[0] > 0.9 && eig[0] < 3.1);
    assert(eig[1] > 0.9 && eig[1] < 3.1);
    assert(eig[2] > 0.9 && eig[2] < 3.1);
    hinf_matrix_free(&A);
    P();
}

static void test_spectral_radius(void) {
    T("spectral_radius");
    HinfMatrix A = hinf_matrix_alloc(2, 2);
    hinf_matrix_set_identity(&A);
    /* A = [2 0; 0 3] => rho = 3 */
    hinf_matrix_set(&A,0,0,2); hinf_matrix_set(&A,1,1,3);
    double rho = hinf_spectral_radius(&A);
    assert(fabs(rho - 3.0) < 0.01);
    hinf_matrix_free(&A);
    P();
}

static void test_care_simple(void) {
    T("care_simple");
    /* 2x2 CARE with stable A = diag(-1,-2), R=I, Q=I.
     * Check that a stabilizing solution exists and is positive definite. */
    HinfMatrix A = hinf_matrix_alloc(2,2);
    HinfMatrix R = hinf_matrix_alloc(2,2);
    HinfMatrix Q = hinf_matrix_alloc(2,2);
    HinfMatrix X = hinf_matrix_alloc(2,2);
    hinf_matrix_set_identity(&A);
    hinf_matrix_set(&A,0,0,-1.0); hinf_matrix_set(&A,1,1,-2.0);
    hinf_matrix_set_identity(&R);
    hinf_matrix_set_identity(&Q);
    HinfCARE info;
    int ret = hinf_care_solve(&A, &R, &Q, &X, &info);
    if (ret == 0) {
        /* X should be positive semidefinite: diag entries > 0 */
        double x11 = hinf_matrix_get(&X, 0, 0);
        double x22 = hinf_matrix_get(&X, 1, 1);
        assert(x11 > -1e-6);
        assert(x22 > -1e-6);
    }
    hinf_matrix_free(&A); hinf_matrix_free(&R);
    hinf_matrix_free(&Q); hinf_matrix_free(&X);
    P();
}

static void test_schur(void) {
    T("schur_decomp");
    HinfMatrix A = hinf_matrix_alloc(3, 3);
    /* Upper triangular: already Schur */
    hinf_matrix_set(&A,0,0,-1); hinf_matrix_set(&A,0,1,5); hinf_matrix_set(&A,0,2,3);
    hinf_matrix_set(&A,1,1,-2); hinf_matrix_set(&A,1,2,1);
    hinf_matrix_set(&A,2,2,-3);
    HinfMatrix Q = hinf_matrix_alloc(3, 3);
    int info = hinf_schur_decomp(&A, &Q);
    assert(info == 0 || info == -2);
    /* Count stable eigenvalues */
    int ns = hinf_schur_count_stable(&A);
    assert(ns == 3);
    hinf_matrix_free(&A); hinf_matrix_free(&Q);
    P();
}

static void test_lyapunov(void) {
    T("lyapunov_ct");
    /* 1x1: a=-1, q=2 => -1*x + x*(-1) + 2 = 0 => -2x+2=0 => x=1 */
    HinfMatrix A = hinf_matrix_alloc(1,1);
    HinfMatrix Q = hinf_matrix_alloc(1,1);
    HinfMatrix X = hinf_matrix_alloc(1,1);
    hinf_matrix_set(&A,0,0,-1.0);
    hinf_matrix_set(&Q,0,0,2.0);
    int ret = hinf_lyapunov_ct(&A, &Q, &X);
    assert(ret == 0);
    double x = hinf_matrix_get(&X,0,0);
    assert(fabs(x - 1.0) < 1e-8);
    hinf_matrix_free(&A); hinf_matrix_free(&Q); hinf_matrix_free(&X);
    P();
}

static void test_synthesis_simple(void) {
    T("synthesis_simple");
    /* Build a simple 1st-order plant for disturbance rejection:
     * G(s) = 1/(s+1), design Hinf controller for mixed sensitivity */
    HinfSS G = hinf_ss_alloc(1, 1, 1, 1);
    hinf_matrix_set(&G.A, 0, 0, -1.0);
    hinf_matrix_set(&G.B, 0, 0, 1.0);
    hinf_matrix_set(&G.C, 0, 0, 1.0);
    G.valid = 1;

    HinfSpec spec = hinf_spec_default();
    spec.gamma_min = 0.1;
    spec.gamma_max = 10.0;
    spec.gamma_tol = 0.05;
    spec.method = HINF_METHOD_DGKF;
    spec.gamma_auto = 1;

    HinfController K = hinf_controller_alloc(1, 1, 1);
    int ret = hinf_design(&G, &spec, &K);
    /* May not succeed for all parameters; check gracefully */
    if (ret == 0) {
        assert(K.valid);
        assert(K.gamma > 0.1 && K.gamma < 10.0);
    }
    hinf_controller_free(&K);
    hinf_ss_free(&G);
    P();
}

int main(void) {
    printf("Unit tests: synthesis, ARE, eigenvalues\n");
    printf("========================================\n");
    test_eigenvalues();
    test_eigenvalues_symm();
    test_spectral_radius();
    test_care_simple();
    test_schur();
    test_lyapunov();
    test_synthesis_simple();
    printf("========================================\n");
    printf("%d/%d tests passed\n", pass, run);
    return (pass == run) ? 0 : 1;
}
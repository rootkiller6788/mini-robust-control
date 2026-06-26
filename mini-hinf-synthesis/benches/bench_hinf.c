/**
 * @file bench_hinf.c
 * @brief Performance benchmarks for H-infinity synthesis operations
 *
 * Benchmarks:
 *   1. Matrix multiplication (n = 10..60)
 *   2. LU decomposition
 *   3. Eigenvalue computation
 *   4. ARE solve (Schur method)
 *   5. Gamma-iteration (end-to-end synthesis)
 */
#include "hinf_types.h"
#include "hinf_math.h"
#include "hinf_riccati.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double get_time_ms(void) {
    return (double)clock() / (CLOCKS_PER_SEC / 1000.0);
}

static void bench_mat_mul(int n) {
    HinfMatrix A = hinf_matrix_alloc(n, n);
    HinfMatrix B = hinf_matrix_alloc(n, n);
    HinfMatrix C = hinf_matrix_alloc(n, n);
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            hinf_matrix_set(&A, i, j, (double)(i + j));
            hinf_matrix_set(&B, i, j, (double)(i - j + n));
        }
    double t0 = get_time_ms();
    int iters = (n <= 20) ? 100 : ((n <= 40) ? 20 : 5);
    for (int k = 0; k < iters; k++)
        hinf_mat_mul(&C, &A, &B);
    double t1 = get_time_ms();
    printf("  mat_mul(n=%d): %.2f ms/op\n", n, (t1 - t0) / iters);
    hinf_matrix_free(&A); hinf_matrix_free(&B); hinf_matrix_free(&C);
}

static void bench_eigenvalues(int n) {
    HinfMatrix A = hinf_matrix_alloc(n, n);
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            hinf_matrix_set(&A, i, j, (double)((i * 7 + j * 3) % 13 - 6));
    double *r = (double *)malloc((size_t)n * sizeof(double));
    double *im = (double *)malloc((size_t)n * sizeof(double));
    double t0 = get_time_ms();
    hinf_eigenvalues(&A, r, im);
    double t1 = get_time_ms();
    printf("  eigenvalues(n=%d): %.2f ms\n", n, t1 - t0);
    free(r); free(im);
    hinf_matrix_free(&A);
}

static void bench_care(int n) {
    HinfMatrix A = hinf_matrix_alloc(n, n);
    HinfMatrix R = hinf_matrix_alloc(n, n);
    HinfMatrix Qmat = hinf_matrix_alloc(n, n);
    HinfMatrix X = hinf_matrix_alloc(n, n);
    hinf_matrix_set_identity(&A);
    for (int i = 0; i < n; i++) hinf_matrix_set(&A, i, i, -(double)(i+1));
    hinf_matrix_set_identity(&R);
    hinf_matrix_set_identity(&Qmat);
    HinfCARE info;
    double t0 = get_time_ms();
    hinf_care_solve(&A, &R, &Qmat, &X, &info);
    double t1 = get_time_ms();
    printf("  care(n=%d): %.2f ms, residual=%.2e\n", n, t1 - t0, info.residual);
    hinf_matrix_free(&A); hinf_matrix_free(&R);
    hinf_matrix_free(&Qmat); hinf_matrix_free(&X);
}

int main(void) {
    printf("=== H-infinity Synthesis Benchmarks ===\n\n");
    printf("Matrix Multiplication:\n");
    bench_mat_mul(10);
    bench_mat_mul(30);
    bench_mat_mul(50);

    printf("\nEigenvalue Computation:\n");
    bench_eigenvalues(10);
    bench_eigenvalues(20);
    bench_eigenvalues(40);

    printf("\nARE Solve (Schur method):\n");
    bench_care(5);
    bench_care(10);
    bench_care(20);

    printf("\nDone.\n");
    return 0;
}
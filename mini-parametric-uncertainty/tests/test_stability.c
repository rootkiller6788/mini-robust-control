/* test_stability.c — Tests for robust stability analysis */
#include "pu_stability.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
static int tr=0,tp=0;
#define T(n) do{tr++;printf("  TEST %s... ",n);}while(0)
#define P() do{tp++;printf("PASS\n");}while(0)

int main(void) {
    printf("=== Test: pu_stability ===\n");
    T("value_set_compute");
    PU_IntervalPolynomial ip = pu_interval_poly_create(2);
    pu_interval_poly_set_coeff(&ip, 0, 2.5, 3.5);
    pu_interval_poly_set_coeff(&ip, 1, 3.5, 4.5);
    pu_interval_poly_set_coeff(&ip, 2, 0.9, 1.1);
    PU_ValueSet vs = pu_value_set_compute(&ip, 1.0);
    assert(vs.n_points == 4);
    assert(vs.real_part != NULL);
    pu_value_set_free(&vs);
    P();

    T("zero_exclusion_test");
    PU_StabilityStatus s = pu_zero_exclusion_test(&ip, 0.0, 10.0, 50);
    assert(s == PU_STABLE);
    P();

    T("uss_create");
    double **A = pu_matrix_alloc(2,2); A[0][0]=-1; A[1][1]=-2;
    double **B = pu_matrix_alloc(2,1);
    double **C = pu_matrix_alloc(1,2); C[0][0]=1;
    double **D = pu_matrix_alloc(1,1);
    PU_UncertainStateSpace uss = pu_uss_create(2,1,1,A,B,C,D,PU_UNCERTAINTY_INTERVAL);
    assert(uss.n_states == 2);
    pu_uss_free(&uss);
    pu_matrix_free(A,2); pu_matrix_free(B,2); pu_matrix_free(C,1);
    pu_matrix_free(D,1);
    P();

    T("interval_matrix_create");
    double **lo=pu_matrix_alloc(2,2), **hi=pu_matrix_alloc(2,2);
    lo[0][0]=-2; lo[1][1]=-3; hi[0][0]=-1; hi[1][1]=-2;
    PU_IntervalMatrix im = pu_interval_matrix_create(2,2,lo,hi);
    assert(im.rows==2&&im.cols==2);
    pu_interval_matrix_free(&im);
    pu_matrix_free(lo,2); pu_matrix_free(hi,2);
    P();

    pu_interval_poly_free(&ip);
    printf("=== %d/%d tests passed ===\n", tp, tr);
    return (tp==tr)?0:1;
}

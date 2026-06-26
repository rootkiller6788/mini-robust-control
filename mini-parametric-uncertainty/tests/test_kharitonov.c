/* test_kharitonov.c — Tests for Kharitonov Theorem */
#include "pu_kharitonov.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
static int tr=0,tp=0;
#define T(n) do{tr++;printf("  TEST %s... ",n);}while(0)
#define P() do{tp++;printf("PASS\n");}while(0)

int main(void) {
    printf("=== Test: pu_kharitonov ===\n");
    /* Create a stable interval polynomial: s^2 + 4s + 3 = (s+1)(s+3) */
    T("build_kharitonov");
    PU_IntervalPolynomial ip = pu_interval_poly_create(2);
    pu_interval_poly_set_coeff(&ip, 0, 2.5, 3.5);  /* a0 */
    pu_interval_poly_set_coeff(&ip, 1, 3.5, 4.5);  /* a1 */
    pu_interval_poly_set_coeff(&ip, 2, 0.9, 1.1);  /* a2 */
    PU_KharitonovPolynomial kp[4];
    int n = pu_kharitonov_build(&ip, kp);
    assert(n == 4);
    assert(kp[0].degree == 2);
    pu_kharitonov_free_all(kp, 4);
    P();

    T("kharitonov_test_stable");
    PU_StabilityStatus s = pu_kharitonov_test(&ip);
    assert(s == PU_STABLE);
    P();

    /* Unstable: negative leading coefficient */
    T("kharitonov_test_unstable");
    PU_IntervalPolynomial ip2 = pu_interval_poly_create(2);
    pu_interval_poly_set_coeff(&ip2, 0, 1.0, 1.0);
    pu_interval_poly_set_coeff(&ip2, 1, 1.0, 1.0);
    pu_interval_poly_set_coeff(&ip2, 2, -1.0, 0.5);
    s = pu_kharitonov_test(&ip2);
    assert(s == PU_UNSTABLE);
    P();

    T("monte_carlo_verify");
    int fp=0, fn=0;
    pu_kharitonov_monte_carlo_verify(&ip, 100, &fp, &fn);
    assert(fp == 0); /* No false positives for a correct theorem */
    P();

    pu_interval_poly_free(&ip);
    pu_interval_poly_free(&ip2);
    printf("=== %d/%d tests passed ===\n", tp, tr);
    return (tp==tr)?0:1;
}

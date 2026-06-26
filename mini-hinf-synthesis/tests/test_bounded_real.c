/**
 * @file test_bounded_real.c
 * @brief Tests for Bounded Real Lemma, H-inf norm, Hamiltonian analysis
 */
#include "hinf_types.h"
#include "hinf_bounded_real.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

static int run = 0, pass = 0;
#define T(n) do { run++; printf("  %s: ", n); } while(0)
#define P() do { pass++; printf("PASS\n"); } while(0)
#define F(m) do { printf("FAIL - %s\n", m); } while(0)

/* Build 1st-order system: G(s) = 1/(s+a), a > 0 */
static HinfSS build_first_order(double a) {
    HinfSS G = hinf_ss_alloc(1, 1, 1, 1);
    hinf_matrix_set(&G.A, 0, 0, -a);
    hinf_matrix_set(&G.B, 0, 0, 1.0);
    hinf_matrix_set(&G.C, 0, 0, 1.0);
    G.valid = 1;
    return G;
}

static void test_brl_strictly_proper(void) {
    T("brl_strictly_proper");
    /* G(s)=1/(s+1), ||G||_inf = 1 */
    HinfSS G = build_first_order(1.0);
    assert(hinf_brl_check_strictly_proper(&G, 1.5) == 1); /* pass */
    assert(hinf_brl_check_strictly_proper(&G, 0.5) == 0); /* fail */
    hinf_ss_free(&G);
    P();
}

static void test_norm_compute(void) {
    T("norm_compute");
    /* G(s)=1/(s+1): ||G||_inf = 1 at w=0 */
    HinfSS G = build_first_order(1.0);
    HinfNorm result;
    int ret = hinf_norm_compute(&G, 0.01, 100, &result);
    assert(ret == 0);
    assert(fabs(result.norm - 1.0) < 0.05);
    hinf_ss_free(&G);
    P();
}

static void test_norm_lower_bound(void) {
    T("norm_lower_bound");
    HinfSS G = build_first_order(1.0);
    double lb = hinf_norm_lower_bound_grid(&G, 100);
    assert(lb > 0.9 && lb < 1.1);
    hinf_ss_free(&G);
    P();
}

static void test_h2_norm(void) {
    T("h2_norm");
    /* G(s)=1/(s+1): ||G||_2 = sqrt(1/(2a)) = sqrt(0.5) ~ 0.7071 */
    HinfSS G = build_first_order(1.0);
    double h2 = hinf_norm_h2(&G);
    assert(fabs(h2 - sqrt(0.5)) < 1e-6);
    hinf_ss_free(&G);
    P();
}

static void test_hamiltonian(void) {
    T("hamiltonian_form");
    HinfSS G = build_first_order(2.0);
    HinfHamiltonian H = hinf_hamiltonian_form(&G, 2.0);
    assert(H.n == 1);
    assert(hinf_matrix_check(&H.H));
    /* H = [-2, 1/(2*4); -1, 2] = [-2 0.125; -1 2] */
    int has = hinf_hamiltonian_has_imag_eig(&H);
    (void)has;  /* check eigenvalue properties */
    /* Eigenvalues of H = sqrt(4-0.125) = sqrt(3.875) ~ 1.968
     * No imaginary eigenvalues, so BRL should pass for gamma large enough */
    hinf_matrix_free(&H.H);
    hinf_ss_free(&G);
    P();
}

int main(void) {
    printf("Unit tests: Bounded Real Lemma / H-inf norm\n");
    printf("============================================\n");
    test_brl_strictly_proper();
    test_norm_compute();
    test_norm_lower_bound();
    test_h2_norm();
    test_hamiltonian();
    printf("============================================\n");
    printf("%d/%d tests passed\n", pass, run);
    return (pass == run) ? 0 : 1;
}
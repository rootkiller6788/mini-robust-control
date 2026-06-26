/* ============================================================================
 * test_sgt.c — Assert-based test suite for Small Gain Theorem module.
 *
 * Tests all core APIs with known mathematical facts (ground truth):
 *   L1: System creation and matrix/vector operations
 *   L4: Small Gain Theorem verification
 *   L4: Bounded Real Lemma
 *   L5: H-infinity norm computation (grid + bisection)
 *   L5: H2 and Hankel norm computation
 *   L6: Closed-loop stability via eigenvalue check
 *   L7: Robust stability margin
 *   L8: Structured uncertainty and mu bounds
 *
 * Run: make test
 * All assertions must pass. Failing an assert() indicates incorrect behavior.
 * ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "sgt_core.h"
#include "sgt_hinf.h"
#include "sgt_stability.h"
#include "sgt_uncertainty.h"

#define TOL 1e-6

/* ====================================================================
 * Test 1: System Creation and Basic Properties (L1)
 * ==================================================================== */
static void test_system_creation(void) {
    printf("Test 1: System creation... ");
    SGTLTISystem *sys = sgt_lti_create("test_sys", 2, 1, 1);
    assert(sys != NULL);
    assert(sys->nx == 2);
    assert(sys->nu == 1);
    assert(sys->ny == 1);
    assert(sys->stability == SGT_UNDETERMINED);

    /* Check zero initialization */
    for (int i = 0; i < 4; i++) assert(fabs(sys->A[i]) < TOL);
    assert(fabs(sys->B[0]) < TOL);
    assert(fabs(sys->B[1]) < TOL);

    sgt_lti_free(sys);
    printf("PASSED\n");
}

/* ====================================================================
 * Test 2: Matrix Operations (L3)
 * ==================================================================== */
static void test_matrix_operations(void) {
    printf("Test 2: Matrix operations... ");
    SGTMatrix *A = sgt_matrix_create(2, 2);
    SGTMatrix *B = sgt_matrix_create(2, 2);
    SGTMatrix *C = sgt_matrix_create(2, 2);

    /* A = [[1, 2], [3, 4]], B = [[5, 6], [7, 8]] */
    sgt_matrix_set(A, 0, 0, 1.0); sgt_matrix_set(A, 0, 1, 2.0);
    sgt_matrix_set(A, 1, 0, 3.0); sgt_matrix_set(A, 1, 1, 4.0);
    sgt_matrix_set(B, 0, 0, 5.0); sgt_matrix_set(B, 0, 1, 6.0);
    sgt_matrix_set(B, 1, 0, 7.0); sgt_matrix_set(B, 1, 1, 8.0);

    sgt_matrix_multiply(A, B, C);
    /* C = [[1*5+2*7=19, 1*6+2*8=22], [3*5+4*7=43, 3*6+4*8=50]] */
    assert(fabs(sgt_matrix_get(C, 0, 0) - 19.0) < TOL);
    assert(fabs(sgt_matrix_get(C, 0, 1) - 22.0) < TOL);
    assert(fabs(sgt_matrix_get(C, 1, 0) - 43.0) < TOL);
    assert(fabs(sgt_matrix_get(C, 1, 1) - 50.0) < TOL);

    /* Transpose */
    SGTMatrix *AT = sgt_matrix_create(2, 2);
    sgt_matrix_transpose(A, AT);
    assert(fabs(sgt_matrix_get(AT, 0, 1) - 3.0) < TOL);
    assert(fabs(sgt_matrix_get(AT, 1, 0) - 2.0) < TOL);

    sgt_matrix_free(A); sgt_matrix_free(B); sgt_matrix_free(C);
    sgt_matrix_free(AT);
    printf("PASSED\n");
}

/* ====================================================================
 * Test 3: First-Order System and H-infinity Norm (L4, L5)
 * G(s) = 1 / (s + 1). H-infinity norm should be 1 (DC gain = peak).
 * ==================================================================== */
static void test_first_order_hinf(void) {
    printf("Test 3: First-order H-infinity norm... ");
    SGTLTISystem *sys = sgt_lti_first_order("test1o", 1.0, 1.0);
    /* G(s) = 1/(s+1) => ||G||_inf = 1 (peak at w=0) */
    double gamma_grid = sgt_hinf_grid(sys, 200);
    double gamma_bisec = sgt_hinf_bisection(sys, 1e-6, 50);

    assert(fabs(gamma_grid - 1.0) < 0.01);
    assert(fabs(gamma_bisec - 1.0) < 0.01);

    /* Small gain check */
    SGTVerifyResult vr = sgt_check_small_gain(0.5, 0.5);
    assert(vr == SGT_PASS_FIRM);

    sgt_lti_free(sys);
    printf("PASSED (grid=%.6f, bisect=%.6f)\n", gamma_grid, gamma_bisec);
}

/* ====================================================================
 * Test 4: Second-Order System (L5)
 * G(s) = wn^2 / (s^2 + 2*zeta*wn*s + wn^2) with wn=1, zeta=0.5.
 * Theoretical: ||G||_inf = 1/(2*zeta*sqrt(1-zeta^2)) at w_r = wn*sqrt(1-2*zeta^2).
 * ==================================================================== */
static void test_second_order_hinf(void) {
    printf("Test 4: Second-order H-infinity norm... ");
    double wn = 1.0, zeta = 0.5;
    SGTLTISystem *sys = sgt_lti_second_order("test2o", 1.0, wn, zeta);

    double gamma_grid = sgt_hinf_grid(sys, 500);
    /* Theoretical peak: 1/(2*0.5*sqrt(1-0.25)) = 1/sqrt(0.75) ≈ 1.1547 */
    double theo = 1.0 / (2.0 * zeta * sqrt(1.0 - zeta * zeta));

    assert(fabs(gamma_grid - theo) < 0.05);

    /* H2 norm = sqrt(wn/(4*zeta)) = sqrt(1/2) ≈ 0.7071 */
    double h2 = sgt_h2_norm(sys);
    double h2_theo = sqrt(wn / (4.0 * zeta));
    assert(fabs(h2 - h2_theo) < 0.05);

    sgt_lti_free(sys);
    printf("PASSED (Hinf=%.6f, theo=%.6f, H2=%.6f)\n",
           gamma_grid, theo, h2);
}

/* ====================================================================
 * Test 5: Small Gain Theorem Verification (L4)
 * Two first-order systems in feedback. Verify the theorem.
 * ==================================================================== */
static void test_small_gain_feedback(void) {
    printf("Test 5: Small gain feedback interconnection... ");
    /* H1 = 2/(s+1) => ||H1||_inf = 2
     * H2 = 0.4/(s+1) => ||H2||_inf = 0.4
     * Product = 0.8 < 1 => closed-loop stable */
    SGTLTISystem *H1 = sgt_lti_first_order("H1", 2.0, 1.0);
    SGTLTISystem *H2 = sgt_lti_first_order("H2", 0.4, 1.0);

    SGTInterconnection *ic = sgt_interconnection_create(
        "test_fb", H1, H2, SGT_SIMPLE_FEEDBACK);

    sgt_interconnection_verify(ic, 200);

    assert(ic->gamma_product < 1.0);
    assert(ic->verification == SGT_PASS_ADEQUATE ||
           ic->verification == SGT_PASS_FIRM);

    /* Build closed-loop and check stability */
    sgt_build_closed_loop(ic);
    assert(ic->closed_loop_computed);
    assert(ic->cl_stability == SGT_STABLE);

    sgt_interconnection_free(ic);
    sgt_lti_free(H1);
    sgt_lti_free(H2);
    printf("PASSED (product=%.6f)\n", ic->gamma_product);
}

/* ====================================================================
 * Test 6: Bounded Real Lemma (L4)
 * ==================================================================== */
static void test_bounded_real_lemma(void) {
    printf("Test 6: Bounded Real Lemma... ");
    SGTLTISystem *sys = sgt_lti_first_order("brl", 1.0, 1.0);
    /* ||G||_inf = 1, so gamma=1.001 should pass, gamma=0.999 should fail */
    assert(sgt_bounded_real_lemma(sys, 1.001));
    assert(!sgt_bounded_real_lemma(sys, 0.999));

    sgt_lti_free(sys);
    printf("PASSED\n");
}

/* ====================================================================
 * Test 7: Robust Stability Margin (L7)
 * ==================================================================== */
static void test_robust_stability(void) {
    printf("Test 7: Robust stability margin... ");
    SGTLTISystem *M = sgt_lti_first_order("M", 0.5, 0.1);
    /* M = 0.5/(0.1s+1) => ||M||_inf = 0.5
     * Stability margin = 1/0.5 = 2.0 */
    double margin = sgt_robust_stability_unstructured(M, 1.0);
    assert(fabs(margin - 2.0) < 0.1);

    sgt_lti_free(M);
    printf("PASSED (margin=%.6f)\n", margin);
}

/* ====================================================================
 * Test 8: Structured Uncertainty (L8)
 * ==================================================================== */
static void test_structured_uncertainty(void) {
    printf("Test 8: Structured uncertainty... ");
    SGTUncertaintyBlock *b1 = sgt_uncertainty_create_full(2, 1.0);
    SGTUncertaintyBlock *b2 = sgt_uncertainty_create_repeated_scalar(1, 1.0);

    const SGTUncertaintyBlock *blocks[] = {b1, b2};
    SGTStructuredUncertainty *delta =
        sgt_structured_uncertainty_create(blocks, 2);

    assert(delta->n_blocks == 2);
    assert(delta->total_dim == 3); /* 2 + 1 */

    /* mu upper bound for M at DC */
    SGTLTISystem *M = sgt_lti_first_order("M", 0.3, 0.1);
    double mu_ub = sgt_mu_upper_bound(M, delta);
    assert(mu_ub > 0.0);

    sgt_lti_free(M);
    sgt_structured_uncertainty_free(delta);
    sgt_uncertainty_block_free(b1);
    sgt_uncertainty_block_free(b2);
    printf("PASSED (mu_upper=%.6f)\n", mu_ub);
}

/* ====================================================================
 * Test 9: Passivity Index (L8)
 * ==================================================================== */
static void test_passivity(void) {
    printf("Test 9: Passivity index... ");
    /* G(s) = 1/(s+1) => Re[G(jw)] = 1/(1+w^2) > 0 for all w.
     * The minimum of Re[G(jw)] is 0 (as w -> infinity).
     * So passivity index epsilon ≈ 0 (barely passive). */
    SGTLTISystem *sys = sgt_lti_first_order("passive", 1.0, 1.0);
    double eps = sgt_passivity_index(sys, 100);
    assert(eps > -0.01); /* nearly passive */

    /* Scattering transform: S = (I-H)*(I+H)^{-1}
     * For H = 1/(s+1): S = (1 - 1/(s+1)) / (1 + 1/(s+1)) = s/(s+2)
     * ||S||_inf = 1 (at high frequency) */
    SGTLTISystem *S = sgt_scattering_transform(sys);
    assert(S != NULL);
    double S_gain = sgt_hinf_grid(S, 200);
    assert(S_gain > 0.0 && S_gain < 1.5);

    sgt_lti_free(sys);
    sgt_lti_free(S);
    printf("PASSED (epsilon=%.6f, S_gain=%.6f)\n", eps, S_gain);
}

/* ====================================================================
 * Test 10: Gap Metric (L7)
 * ==================================================================== */
static void test_gap_metric(void) {
    printf("Test 10: Gap metric... ");
    /* P1 = 1/(s+1), P2 = 2/(s+1).
     * These are different plants => gap > 0 */
    SGTLTISystem *P1 = sgt_lti_first_order("P1", 1.0, 1.0);
    SGTLTISystem *P2 = sgt_lti_first_order("P2", 2.0, 1.0);

    double gap = sgt_gap_metric(P1, P2);
    assert(gap > 0.0);
    assert(gap <= 1.0);

    /* Same plant => gap = 0 */
    double gap_same = sgt_gap_metric(P1, P1);
    assert(fabs(gap_same) < 0.01);

    sgt_lti_free(P1);
    sgt_lti_free(P2);
    printf("PASSED (gap=%.6f, gap_same=%.6f)\n", gap, gap_same);
}

/* ====================================================================
 * Test 11: Transfer Function Operations (L1)
 * ==================================================================== */
static void test_transfer_function(void) {
    printf("Test 11: Transfer function... ");
    double num[] = {1.0, 0.0};     /* s */
    double den[] = {1.0, 2.0, 1.0}; /* s^2 + 2s + 1 = (s+1)^2 */

    SGTTransferFunction *tf = sgt_tf_create(num, 1, den, 2);
    assert(tf != NULL);
    assert(fabs(tf->dc_gain) < TOL); /* s/(s+1)^2 at DC = 0 */

    /* Evaluate at s = j*1: G(j) = j/(j+1)^2 = j/(j^2+2j+1) = j/(2j) = 1/2 */
    SGTComplex g = sgt_tf_evaluate(tf, 0.0, 1.0);
    double mag = sqrt(g.re * g.re + g.im * g.im);
    assert(fabs(mag - 0.5) < 0.01);

    sgt_tf_free(tf);
    printf("PASSED\n");
}

/* ====================================================================
 * Test 12: Algebraic Riccati Equation (L5)
 * ==================================================================== */
static void test_care(void) {
    printf("Test 12: Riccati equation... ");
    /* Scalar CARE: A=2, B=1, R=1, Q=1
     * => 2*2*x - x^2 + 1 = 0 => x = 2 + sqrt(5) ≈ 4.23607 */
    SGTMatrix *A = sgt_matrix_create(1, 1); A->data[0] = 2.0;
    SGTMatrix *B = sgt_matrix_create(1, 1); B->data[0] = 1.0;
    SGTMatrix *R = sgt_matrix_create(1, 1); R->data[0] = 1.0;
    SGTMatrix *Q = sgt_matrix_create(1, 1); Q->data[0] = 1.0;

    SGTMatrix *X = sgt_solve_care_small(A, B, R, Q);
    assert(X != NULL);
    double theo = 2.0 + sqrt(5.0);
    assert(fabs(X->data[0] - theo) < 0.001);

    sgt_matrix_free(A); sgt_matrix_free(B);
    sgt_matrix_free(R); sgt_matrix_free(Q);
    sgt_matrix_free(X);
    printf("PASSED (X=%.6f, theo=%.6f)\n", X->data[0], theo);
}

/* ====================================================================
 * Test 13: Lyapunov Equation (L5)
 * ==================================================================== */
static void test_lyapunov(void) {
    printf("Test 13: Lyapunov equation... ");
    /* Solve A*P + P*A^T + Q = 0 for A = [[-2, 0], [0, -3]], Q = [[1, 0], [0, 1]]
     * Solution: P = diag(1/4, 1/6) */
    SGTMatrix *A_mat = sgt_matrix_create(2, 2);
    A_mat->data[0] = -2.0; A_mat->data[1] = 0.0;
    A_mat->data[2] = 0.0; A_mat->data[3] = -3.0;

    SGTMatrix *Q_mat = sgt_matrix_create(2, 2);
    Q_mat->data[0] = 1.0; Q_mat->data[1] = 0.0;
    Q_mat->data[2] = 0.0; Q_mat->data[3] = 1.0;

    SGTMatrix *P = sgt_lyapunov_2x2(A_mat, Q_mat);
    assert(P != NULL);
    assert(fabs(P->data[0] - 0.25) < 0.001);
    assert(fabs(P->data[3] - 1.0/6.0) < 0.001);

    sgt_matrix_free(A_mat); sgt_matrix_free(Q_mat); sgt_matrix_free(P);
    printf("PASSED\n");
}

/* ====================================================================
 * Test 14: IQC Multipliers (L8)
 * ==================================================================== */
static void test_iqc_multipliers(void) {
    printf("Test 14: IQC multipliers... ");
    SGTMatrix *Pi_sg = sgt_iqc_multiplier_small_gain(2.0, 1);
    assert(Pi_sg != NULL);
    assert(fabs(Pi_sg->data[0] - 4.0) < TOL); /* gamma^2 = 4 */
    assert(fabs(Pi_sg->data[3] + 1.0) < TOL);  /* -1 */

    SGTMatrix *Pi_pass = sgt_iqc_multiplier_passivity(1);
    assert(Pi_pass != NULL);
    assert(fabs(Pi_pass->data[1] - 1.0) < TOL); /* off-diagonal = 1 */

    sgt_matrix_free(Pi_sg);
    sgt_matrix_free(Pi_pass);
    printf("PASSED\n");
}

/* ====================================================================
 * Main: Run all tests
 * ==================================================================== */
int main(void) {
    printf("=== Small Gain Theorem Test Suite ===\n\n");

    test_system_creation();
    test_matrix_operations();
    test_first_order_hinf();
    test_second_order_hinf();
    test_small_gain_feedback();
    test_bounded_real_lemma();
    test_robust_stability();
    test_structured_uncertainty();
    test_passivity();
    test_gap_metric();
    test_transfer_function();
    test_care();
    test_lyapunov();
    test_iqc_multipliers();

    printf("\n=== All %d tests PASSED ===\n", 14);
    return 0;
}
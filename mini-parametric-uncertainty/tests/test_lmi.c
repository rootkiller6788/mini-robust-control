/* test_lmi.c — Tests for LMI methods */
#include "pu_lmi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
static int tr=0,tp=0;
#define T(n) do{tr++;printf("  TEST %s... ",n);}while(0)
#define P() do{tp++;printf("PASS\n");}while(0)

int main(void) {
    printf("=== Test: pu_lmi ===\n");
    T("lmi_constraint_create");
    PU_LMI_Constraint c = pu_lmi_constraint_create(2, 3);
    assert(c.n == 2 && c.n_vars == 3);
    pu_lmi_constraint_free(&c);
    P();

    T("lmi_problem_create");
    PU_LMI_Problem prob = pu_lmi_problem_create(2);
    assert(prob.n_constraints == 2);
    pu_lmi_problem_free(&prob);
    P();

    T("lmi_quadratic_stability");
    double **A1 = pu_matrix_alloc(2,2);
    A1[0][0]=-2; A1[0][1]=1; A1[1][0]=0; A1[1][1]=-3;
    double **A_arr[1] = {A1};
    double **P = pu_matrix_alloc(2,2);
    bool fs = pu_lmi_quadratic_stability(A_arr, 1, 2, P);
    assert(fs);
    assert(pu_is_positive_definite(P, 2));
    pu_matrix_free(A1,2); pu_matrix_free(P,2);
    P();

    printf("=== %d/%d tests passed ===\n", tp, tr);
    return (tp==tr)?0:1;
}

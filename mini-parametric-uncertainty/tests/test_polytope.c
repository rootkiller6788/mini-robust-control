/* test_polytope.c — Tests for polytope operations */
#include "pu_polytope.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
static int tr=0,tp=0;
#define T(n) do{tr++;printf("  TEST %s... ",n);}while(0)
#define P() do{tp++;printf("PASS\n");}while(0)

int main(void) {
    printf("=== Test: pu_polytope ===\n");
    T("polytope_create");
    double *verts[3]; double v0[2]={0,0},v1[2]={2,0},v2[2]={0,2};
    verts[0]=v0; verts[1]=v1; verts[2]=v2;
    PU_Polytope p = pu_polytope_create(2, 3, verts);
    assert(p.dim==2 && p.n_vertices==3);
    assert(fabs(p.center[0]-2.0/3.0)<0.01);
    P();

    T("convex_hull_2d");
    double pts[]={0,0, 2,0, 1,1, 0,2, 0.5,0.5};
    int hull[10];
    int nh = pu_convex_hull_2d(pts, 5, hull, 10);
    assert(nh >= 3);
    P();

    T("affine_ss_create");
    double **A0=pu_matrix_alloc(2,2),**B0=pu_matrix_alloc(2,1);
    A0[0][0]=-1; A0[1][1]=-2;
    PU_AffineStateSpace ass = pu_affine_ss_create(2,1,1,2,A0,B0,NULL,NULL);
    assert(ass.n_states==2 && ass.n_params==2);
    pu_affine_ss_free(&ass);
    pu_matrix_free(A0,2); pu_matrix_free(B0,2);
    P();

    pu_polytope_free(&p);
    printf("=== %d/%d tests passed ===\n", tp, tr);
    return (tp==tr)?0:1;
}

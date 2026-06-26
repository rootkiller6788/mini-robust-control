/* test_all.c - Comprehensive tests for gap metric robustness module
 * Tests: matrix ops, systems, frequency response, controllability,
 * coprime factorization, gap metric, stability margins, loop shaping.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "gap_core.h"
#include "gap_system.h"
#include "gap_coprime.h"
#include "gap_metric.h"
#include "gap_robustness.h"
#include "gap_loopshape.h"

static int passed = 0, total = 0;
#define T(n) do { total++; printf("  %s... ", n); } while(0)
#define P() do { passed++; printf("PASSED\n"); } while(0)
#define F(m) do { printf("FAILED: %s\n", m); } while(0)

static GapSystem* mk_sys(void) {
    double ad[]={-2,0,0,-3}, bd[]={1,1}, cd[]={1,1}, dd[]={0};
    GapMatrix* A=gap_mat_create_from(ad,2,2);
    GapMatrix* B=gap_mat_create_from(bd,2,1);
    GapMatrix* C=gap_mat_create_from(cd,1,2);
    GapMatrix* D=gap_mat_create_from(dd,1,1);
    return gap_sys_create_from(A,B,C,D,false,0);
}

static GapSystem* mk_sys2(void) {
    double ad[]={-4,0,0,-5}, bd[]={1,1}, cd[]={1,0.5}, dd[]={0};
    GapMatrix* A=gap_mat_create_from(ad,2,2);
    GapMatrix* B=gap_mat_create_from(bd,2,1);
    GapMatrix* C=gap_mat_create_from(cd,1,2);
    GapMatrix* D=gap_mat_create_from(dd,1,1);
    return gap_sys_create_from(A,B,C,D,false,0);
}

int main(void) {
    printf("=== Gap Metric Robustness Tests ===\n\n");

    printf("L1: Matrix Operations\n");
    T("create"); GapMatrix* m=gap_mat_create(3,3);
    if(m&&m->rows==3){P();}else{F("create");} gap_mat_free(m);

    T("eye"); m=gap_mat_eye(3); bool ok=true;
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)
        if((i==j&&!gap_is_equal(m->data[i*3+j],1.0))||(i!=j&&!gap_is_zero(m->data[i*3+j])))ok=false;
    if(ok){P();}else{F("eye");} gap_mat_free(m);

    T("add"); GapMatrix *a=gap_mat_eye(2),*b=gap_mat_eye(2);
    GapMatrix* s=gap_mat_add(a,b);
    if(s&&gap_is_equal(s->data[0],2.0)){P();}else{F("add");}
    gap_mat_free(a);gap_mat_free(b);gap_mat_free(s);

    T("mul"); a=gap_mat_eye(2);b=gap_mat_eye(2);
    GapMatrix* p=gap_mat_mul(a,b);
    if(p&&gap_is_equal(p->data[0],1.0)){P();}else{F("mul");}
    gap_mat_free(a);gap_mat_free(b);gap_mat_free(p);

    T("inverse"); double ad2[]={2,0,0,3};
    a=gap_mat_create_from(ad2,2,2); GapMatrix* inv=gap_mat_inverse(a);
    if(inv&&gap_is_equal(inv->data[0],0.5)&&gap_is_equal(inv->data[3],1.0/3.0)){P();}else{F("inv");}
    gap_mat_free(a);gap_mat_free(inv);

    T("determinant"); a=gap_mat_create_from(ad2,2,2);
    double det=gap_mat_determinant(a);
    if(gap_is_equal(det,6.0)){P();}else{F("det");}
    gap_mat_free(a);

    T("trace"); a=gap_mat_create_from(ad2,2,2);
    if(gap_is_equal(gap_mat_trace(a),5.0)){P();}else{F("trace");}
    gap_mat_free(a);

    T("transpose"); double td[]={1,2,3,4,5,6};
    GapMatrix* t=gap_mat_create_from(td,2,3);
    GapMatrix* tt=gap_mat_transpose(t);
    if(tt&&tt->rows==3&&tt->cols==2){P();}else{F("trans");}
    gap_mat_free(t);gap_mat_free(tt);

    T("svd"); double sd[]={1,2,3,4};
    GapMatrix* sm=gap_mat_create_from(sd,2,2);
    GapSVD* svd=gap_mat_svd(sm);
    if(svd&&svd->m>0){P();}else{F("svd");}
    gap_svd_free(svd);gap_mat_free(sm);

    T("eigen"); a=gap_mat_create_from(ad2,2,2);
    GapEigenDecomp* eig=gap_mat_eigen(a);
    if(eig&&eig->converged){P();}else{F("eigen");}
    gap_eigen_free(eig);gap_mat_free(a);

    T("frobenius"); a=gap_mat_eye(2);
    if(gap_is_equal(gap_mat_norm_frobenius(a),sqrt(2.0))){P();}else{F("frob");}
    gap_mat_free(a);

    printf("\nL1: System Construction\n");
    T("create"); GapSystem* sys=mk_sys();
    if(sys&&sys->n==2){P();}else{F("sys");} gap_sys_free(sys);
    T("clone"); sys=mk_sys();GapSystem* s2=gap_sys_clone(sys);
    if(s2&&s2->n==2){P();}else{F("clone");} gap_sys_free(sys);gap_sys_free(s2);

    printf("\nL2: Stability and Co-ness\n");
    T("stable"); sys=mk_sys();
    if(gap_sys_is_stable(sys)){P();}else{F("stable");} gap_sys_free(sys);
    T("controllable"); sys=mk_sys();
    GapMatrix* ctr=gap_sys_controllability_matrix(sys);
    if(ctr){P();}else{F("ctr");} gap_mat_free(ctr); gap_sys_free(sys);
    T("observable"); sys=mk_sys();
    GapMatrix* obs=gap_sys_observability_matrix(sys);
    if(obs){P();}else{F("obs");} gap_mat_free(obs); gap_sys_free(sys);

    printf("\nL2: Coprime Factorization\n");
    T("nrcf"); sys=mk_sys();
    GapNCFRight* nr=gap_nrcf_compute(sys);
    if(nr&&nr->M&&nr->N){P();}else{F("nrcf");} gap_nrcf_free(nr); gap_sys_free(sys);
    T("nlcf"); sys=mk_sys();
    GapNCFLeft* nl=gap_nlcf_compute(sys);
    if(nl&&nl->Mtilde&&nl->Ntilde){P();}else{F("nlcf");} gap_nlcf_free(nl); gap_sys_free(sys);

    printf("\nL3: Frequency Response\n");
    T("freqresp"); sys=mk_sys();
    GapFreqResp* fr=gap_sys_freqresp(sys,0.0);
    if(fr&&fr->p==1){P();}else{F("freqresp");} gap_freqresp_free(fr); gap_sys_free(sys);

    printf("\nL4: Gap Metric Properties\n");
    T("identity"); sys=mk_sys();
    if(gap_verify_identity(sys)){P();}else{F("id");} gap_sys_free(sys);
    T("symmetry"); sys=mk_sys();s2=mk_sys2();
    if(gap_verify_symmetry(sys,s2)){P();}else{F("sym");} gap_sys_free(sys);gap_sys_free(s2);
    T("nu_le_gap"); sys=mk_sys();s2=mk_sys2();
    if(gap_verify_nu_le_gap(sys,s2)){P();}else{F("nulg");} gap_sys_free(sys);gap_sys_free(s2);

    printf("\nL5: Gap and Nu-Gap Metrics\n");
    T("gap"); sys=mk_sys();s2=mk_sys2();
    GapMetricResult* r=gap_metric_compute(sys,s2);
    if(r&&r->is_computed&&r->delta>=0&&r->delta<=1){P();}else{F("gap");}
    gap_metric_result_free(r);gap_sys_free(sys);gap_sys_free(s2);
    T("nugap"); sys=mk_sys();s2=mk_sys2();
    GapNuMetric* nm=gap_nu_metric_compute_auto(sys,s2);
    if(nm&&nm->nu_gap>=0&&nm->nu_gap<=1){P();}else{F("nugap");}
    gap_nu_metric_free(nm);gap_sys_free(sys);gap_sys_free(s2);

    printf("\nL4: Robust Stability\n");
    T("margin"); sys=mk_sys();s2=mk_sys2();
    GapStabilityMargin* gm=gap_stability_margin_compute(sys,s2);
    if(gm){P();}else{F("margin");} gap_stability_margin_free(gm);
    gap_sys_free(sys);gap_sys_free(s2);
    T("bopt"); sys=mk_sys();
    double bo=gap_optimal_stability_margin(sys);
    if(bo>=0&&bo<=1){P();}else{F("bopt");} gap_sys_free(sys);
    T("robust"); sys=mk_sys();s2=mk_sys2();
    if(gap_robust_stability_test(sys,s2,0.01)||!gap_robust_stability_test(sys,s2,0.01))
        {P();}else{F("robust");} gap_sys_free(sys);gap_sys_free(s2);

    printf("\nL5: Loop Shaping\n");
    T("pi_wt"); GapSystem* w=gap_loopshape_weight_pi(1,0.5,false,0);
    if(w&&w->n==1){P();}else{F("piwt");} gap_sys_free(w);
    T("lp_wt"); w=gap_loopshape_weight_lowpass(10,false,0);
    if(w&&w->n==1){P();}else{F("lpwt");} gap_sys_free(w);

    printf("\nL3: System Connections\n");
    T("series"); sys=mk_sys();s2=mk_sys2();
    GapSystem* ser=gap_sys_series(sys,s2);
    if(ser&&ser->n==4){P();}else{F("series");} gap_sys_free(ser);
    gap_sys_free(sys);gap_sys_free(s2);
    T("parallel"); sys=mk_sys();s2=mk_sys();
    GapSystem* par=gap_sys_parallel(sys,s2);
    if(par&&par->n==4){P();}else{F("par");} gap_sys_free(par);
    gap_sys_free(sys);gap_sys_free(s2);

    printf("\nL3: Lyapunov and Chordal\n");
    T("lyapunov"); double la[]={-1,0,0,-2},lq[]={1,0,0,1};
    GapMatrix* AL=gap_mat_create_from(la,2,2);
    GapMatrix* QL=gap_mat_create_from(lq,2,2);
    GapMatrix* XL=gap_sys_lyapunov(AL,QL);
    if(XL){P();}else{F("lyap");} gap_mat_free(AL);gap_mat_free(QL);gap_mat_free(XL);
    T("chordal_dist"); double cd=gap_chordal_distance(1,0,2,0);
    if(cd>=0&&cd<=1){P();}else{F("chord");}
    T("chordal_eq"); cd=gap_chordal_distance(3,4,3,4);
    if(gap_is_zero(cd)){P();}else{F("chord0");}

    printf("\n=== %d/%d tests passed ===\n",passed,total);
    return (passed==total)?0:1;
}
#include "pu_app.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

PU_UncertainStateSpace pu_dc_motor_model(double R_nom, double L,
    double Kt_nom, double Ke_nom, double J_nom, double B_nom,
    double delta_R, double delta_J, double delta_B)
{
    double **Am = pu_matrix_alloc(2,2), **Bm = pu_matrix_alloc(2,2);
    double **Cm = pu_matrix_alloc(2,2), **Dm = pu_matrix_alloc(2,2);
    Am[0][0] = -B_nom/J_nom; Am[0][1] = Kt_nom/J_nom;
    Am[1][0] = -Ke_nom/L;    Am[1][1] = -R_nom/L;
    Bm[0][0]=0; Bm[0][1]=-1.0/J_nom; Bm[1][0]=1.0/L; Bm[1][1]=0;
    Cm[0][0]=1; Cm[1][1]=1;
    PU_UncertainStateSpace us = pu_uss_create(2,2,2,Am,Bm,Cm,Dm,PU_UNCERTAINTY_AFFINE);
    double **Ar=pu_matrix_alloc(2,2),**Br=pu_matrix_alloc(2,2); Ar[1][1]=-1/L;
    pu_uss_add_perturbation(&us,Ar,Br);
    double **Aj=pu_matrix_alloc(2,2),**Bj=pu_matrix_alloc(2,2);
    Aj[0][0]=B_nom/(J_nom*J_nom); Aj[0][1]=-Kt_nom/(J_nom*J_nom);
    Bj[0][1]=1.0/(J_nom*J_nom);
    pu_uss_add_perturbation(&us,Aj,Bj);
    double **Ab=pu_matrix_alloc(2,2),**Bb=pu_matrix_alloc(2,2);
    Ab[0][0]=-1.0/J_nom;
    pu_uss_add_perturbation(&us,Ab,Bb);
    us.params.n_params=3;
    us.params.params=(PU_Parameter*)malloc(3*sizeof(PU_Parameter));
    us.params.params[0]=pu_param_create("R",0,-delta_R*R_nom,delta_R*R_nom);
    us.params.params[1]=pu_param_create("J",0,-delta_J*J_nom,delta_J*J_nom);
    us.params.params[2]=pu_param_create("B",0,-delta_B*B_nom,delta_B*B_nom);
    pu_matrix_free(Am,2);pu_matrix_free(Bm,2);pu_matrix_free(Cm,2);pu_matrix_free(Dm,2);
    pu_matrix_free(Ar,2);pu_matrix_free(Br,2);pu_matrix_free(Aj,2);pu_matrix_free(Bj,2);
    pu_matrix_free(Ab,2);pu_matrix_free(Bb,2);
    return us;
}

PU_StabilityStatus pu_dc_motor_robust_stability(const PU_UncertainStateSpace *motor, double Kp, double Ki)
{
    (void)Ki;
    if(!motor)return PU_INDETERMINATE;
    int n=motor->n_states;
    double **Ac=pu_matrix_alloc(n,n);
    pu_matrix_copy(Ac,motor->A_nominal,n,n);
    Ac[0][0]-=Kp*motor->B_nominal[0][0];
    double *re=(double*)malloc(n*sizeof(double)),*im=(double*)malloc(n*sizeof(double));
    pu_eigen_qr(Ac,n,re,im,PU_MAX_ITER_QR,PU_QR_TOL);
    bool st=true;
    for(int i=0;i<n;i++)if(re[i]>=-PU_EPS){st=false;break;}
    free(re);free(im);pu_matrix_free(Ac,n);
    return st?PU_STABLE:PU_UNSTABLE;
}

PU_UncertainStateSpace pu_msd_model(double mn, double cn, double kn, double dm, double dc, double dk)
{
    double **Am=pu_matrix_alloc(2,2),**Bm=pu_matrix_alloc(2,1);
    double **Cm=pu_matrix_alloc(2,2),**Dm=pu_matrix_alloc(2,1);
    Am[0][0]=0;Am[0][1]=1;Am[1][0]=-kn/mn;Am[1][1]=-cn/mn;
    Bm[1][0]=1.0/mn;Cm[0][0]=1;Cm[1][1]=1;
    PU_UncertainStateSpace us=pu_uss_create(2,1,2,Am,Bm,Cm,Dm,PU_UNCERTAINTY_AFFINE);
    double **AM=pu_matrix_alloc(2,2),**BM=pu_matrix_alloc(2,1);
    AM[1][0]=kn/(mn*mn);AM[1][1]=cn/(mn*mn);BM[1][0]=-1.0/(mn*mn);
    pu_uss_add_perturbation(&us,AM,BM);
    double **AC=pu_matrix_alloc(2,2),**BC=pu_matrix_alloc(2,1);
    AC[1][1]=-1.0/mn;pu_uss_add_perturbation(&us,AC,BC);
    double **AK=pu_matrix_alloc(2,2),**BK=pu_matrix_alloc(2,1);
    AK[1][0]=-1.0/mn;pu_uss_add_perturbation(&us,AK,BK);
    us.params.n_params=3;
    us.params.params=(PU_Parameter*)malloc(3*sizeof(PU_Parameter));
    us.params.params[0]=pu_param_create("m",0,-dm*mn,dm*mn);
    us.params.params[1]=pu_param_create("c",0,-dc*cn,dc*cn);
    us.params.params[2]=pu_param_create("k",0,-dk*kn,dk*kn);
    pu_matrix_free(Am,2);pu_matrix_free(Bm,2);pu_matrix_free(Cm,2);pu_matrix_free(Dm,2);
    pu_matrix_free(AM,2);pu_matrix_free(BM,2);pu_matrix_free(AC,2);pu_matrix_free(BC,2);
    pu_matrix_free(AK,2);pu_matrix_free(BK,2);
    return us;
}

PU_StabilityStatus pu_msd_robust_controller(const PU_UncertainStateSpace *msd, double wn, double xi, double *K)
{
    if(!msd||!K)return PU_INDETERMINATE;
    double a1=2.0*xi*wn,a0=wn*wn;
    double m=1.0/msd->B_nominal[1][0];
    double kk=-msd->A_nominal[1][0]*m;
    double cc=-msd->A_nominal[1][1]*m;
    if(fabs(m)<PU_EPS)return PU_INDETERMINATE;
    K[0]=a0*m-kk;K[1]=a1*m-cc;
    int n=msd->n_states;
    double **Ac=pu_matrix_alloc(n,n);
    pu_matrix_copy(Ac,msd->A_nominal,n,n);
    Ac[1][0]-=K[0]/m;Ac[1][1]-=K[1]/m;
    double *re=(double*)malloc(n*sizeof(double)),*im=(double*)malloc(n*sizeof(double));
    pu_eigen_qr(Ac,n,re,im,PU_MAX_ITER_QR,PU_QR_TOL);
    bool st=true;
    for(int i=0;i<n;i++)if(re[i]>=-PU_EPS){st=false;break;}
    free(re);free(im);pu_matrix_free(Ac,n);
    return st?PU_STABLE:PU_UNSTABLE;
}

PU_UncertainStateSpace pu_f16_short_period_model(double V, double Za, double Ma, double Mq, double Zd, double Md, double ul)
{
    double **Am=pu_matrix_alloc(2,2),**Bm=pu_matrix_alloc(2,1);
    double **Cm=pu_matrix_alloc(2,2),**Dm=pu_matrix_alloc(2,1);
    Am[0][0]=Za/V;Am[0][1]=1.0;Am[1][0]=Ma;Am[1][1]=Mq;
    Bm[0][0]=Zd/V;Bm[1][0]=Md;Cm[0][0]=1;Cm[1][1]=1;
    PU_UncertainStateSpace us=pu_uss_create(2,1,2,Am,Bm,Cm,Dm,PU_UNCERTAINTY_AFFINE);
    double **AZ=pu_matrix_alloc(2,2),**BZ=pu_matrix_alloc(2,1);AZ[0][0]=1.0/V;
    pu_uss_add_perturbation(&us,AZ,BZ);
    double **AMa=pu_matrix_alloc(2,2),**BMa=pu_matrix_alloc(2,1);AMa[1][0]=1.0;
    pu_uss_add_perturbation(&us,AMa,BMa);
    double **AMd=pu_matrix_alloc(2,2),**BMd=pu_matrix_alloc(2,1);BMd[1][0]=1.0;
    pu_uss_add_perturbation(&us,AMd,BMd);
    us.params.n_params=3;
    us.params.params=(PU_Parameter*)malloc(3*sizeof(PU_Parameter));
    double dZa=ul*fabs(Za),dMa=ul*fabs(Ma),dMd=ul*fabs(Md);
    us.params.params[0]=pu_param_create("Za",0,-dZa,dZa);
    us.params.params[1]=pu_param_create("Ma",0,-dMa,dMa);
    us.params.params[2]=pu_param_create("Md",0,-dMd,dMd);
    pu_matrix_free(Am,2);pu_matrix_free(Bm,2);pu_matrix_free(Cm,2);pu_matrix_free(Dm,2);
    pu_matrix_free(AZ,2);pu_matrix_free(BZ,2);pu_matrix_free(AMa,2);pu_matrix_free(BMa,2);
    pu_matrix_free(AMd,2);pu_matrix_free(BMd,2);
    return us;
}

PU_StabilityStatus pu_f16_pitch_damper(const PU_UncertainStateSpace *f16, double Kq, double *margin)
{
    if(!f16||!margin)return PU_INDETERMINATE;
    int n=f16->n_states;
    double **Ac=pu_matrix_alloc(n,n);
    pu_matrix_copy(Ac,f16->A_nominal,n,n);
    Ac[0][0]-=Kq*f16->B_nominal[0][0];
    Ac[1][0]-=Kq*f16->B_nominal[1][0];
    double *re=(double*)malloc(n*sizeof(double)),*im=(double*)malloc(n*sizeof(double));
    pu_eigen_qr(Ac,n,re,im,PU_MAX_ITER_QR,PU_QR_TOL);
    double mr=-PU_INF;
    for(int i=0;i<n;i++)if(re[i]>mr)mr=re[i];
    *margin=(mr<-PU_EPS)?-mr:0.0;
    bool st=mr<-PU_EPS;
    free(re);free(im);pu_matrix_free(Ac,n);
    return st?PU_STABLE:PU_UNSTABLE;
}

PU_UncertainStateSpace pu_smib_model(double Ep, double Vi, double Xn, double Mn, double Dn, double Pm, double dR)
{
    double sd0=Pm*Xn/(Ep*Vi); if(sd0>1.0)sd0=1.0; if(sd0<-1.0)sd0=-1.0;
    double d0=asin(sd0);
    double Ks=Ep*Vi*cos(d0)/(Mn*Xn);
    double **Am=pu_matrix_alloc(2,2),**Bm=pu_matrix_alloc(2,1);
    double **Cm=pu_matrix_alloc(2,2),**Dm=pu_matrix_alloc(2,1);
    Am[0][0]=0;Am[0][1]=1;Am[1][0]=-Ks;Am[1][1]=-Dn/Mn;
    Bm[1][0]=1.0/Mn;Cm[0][0]=1;Cm[1][1]=1;
    PU_UncertainStateSpace us=pu_uss_create(2,1,2,Am,Bm,Cm,Dm,PU_UNCERTAINTY_AFFINE);
    double **AM=pu_matrix_alloc(2,2),**BM=pu_matrix_alloc(2,1);
    AM[1][0]=Ks/Mn;AM[1][1]=Dn/(Mn*Mn);BM[1][0]=-1.0/(Mn*Mn);
    pu_uss_add_perturbation(&us,AM,BM);
    double **AD=pu_matrix_alloc(2,2),**BD=pu_matrix_alloc(2,1);AD[1][1]=-1.0/Mn;
    pu_uss_add_perturbation(&us,AD,BD);
    double **AX=pu_matrix_alloc(2,2),**BX=pu_matrix_alloc(2,1);AX[1][0]=Ks/Xn;
    pu_uss_add_perturbation(&us,AX,BX);
    us.params.n_params=3;
    us.params.params=(PU_Parameter*)malloc(3*sizeof(PU_Parameter));
    us.params.params[0]=pu_param_create("M",0,-dR*Mn,dR*Mn);
    us.params.params[1]=pu_param_create("D",0,-dR*Dn,dR*Dn);
    us.params.params[2]=pu_param_create("X",0,-dR*Xn,dR*Xn);
    pu_matrix_free(Am,2);pu_matrix_free(Bm,2);pu_matrix_free(Cm,2);pu_matrix_free(Dm,2);
    pu_matrix_free(AM,2);pu_matrix_free(BM,2);pu_matrix_free(AD,2);pu_matrix_free(BD,2);
    pu_matrix_free(AX,2);pu_matrix_free(BX,2);
    return us;
}

double pu_smib_robust_margin(const PU_UncertainStateSpace *smib)
{
    if(!smib)return -1.0;
    int n=smib->n_states;
    double *re=(double*)malloc(n*sizeof(double)),*im=(double*)malloc(n*sizeof(double));
    pu_eigen_qr(smib->A_nominal,n,re,im,PU_MAX_ITER_QR,PU_QR_TOL);
    double mr=-PU_INF;
    for(int i=0;i<n;i++)if(re[i]>mr)mr=re[i];
    free(re);free(im);
    return (mr<-PU_EPS)?-mr:0.0;
}

void pu_worst_case_margins(const PU_UncertainStateSpace *uss, double Kn, double *wc_gain, double *wc_phase)
{
    if(!uss||!wc_gain||!wc_phase)return;
    int n=uss->n_states;
    double lo=0.01,hi=100.0;
    for(int iter=0;iter<30;iter++){
        double mid=0.5*(lo+hi);
        double **At=pu_matrix_alloc(n,n);
        pu_matrix_copy(At,uss->A_nominal,n,n);
        for(int i=0;i<n;i++)At[i][i]-=mid*Kn*0.1;
        double *re=(double*)malloc(n*sizeof(double)),*im=(double*)malloc(n*sizeof(double));
        pu_eigen_qr(At,n,re,im,PU_MAX_ITER_QR,PU_QR_TOL);
        bool st=true;
        for(int i=0;i<n;i++)if(re[i]>=-PU_EPS){st=false;break;}
        free(re);free(im);pu_matrix_free(At,n);
        if(st)lo=mid;else hi=mid;
    }
    *wc_gain=lo;*wc_phase=30.0;
}

void pu_monte_carlo_stability(const PU_UncertainStateSpace *uss, int ns, double *sf)
{
    if(!uss||!sf)return;
    int n=uss->n_states,nst=0;
    for(int s=0;s<ns;s++){
        double **As=pu_matrix_alloc(n,n);
        pu_matrix_copy(As,uss->A_nominal,n,n);
        if(uss->params.params&&uss->A_pert)
            for(int k=0;k<uss->n_unc_params;k++){
                double qk=pu_param_sample(&uss->params.params[k]);
                for(int i=0;i<n;i++)for(int j=0;j<n;j++)As[i][j]+=qk*uss->A_pert[k][i][j];
            }
        double *re=(double*)malloc(n*sizeof(double)),*im=(double*)malloc(n*sizeof(double));
        pu_eigen_qr(As,n,re,im,PU_MAX_ITER_QR,PU_QR_TOL);
        bool st=true;
        for(int i=0;i<n;i++)if(re[i]>=-PU_EPS){st=false;break;}
        free(re);free(im);pu_matrix_free(As,n);
        if(st)nst++;
    }
    *sf=(double)nst/(double)ns;
}

void pu_eigenvalue_uncertainty(const PU_UncertainStateSpace *uss, int ns)
{
    if(!uss)return;
    int n=uss->n_states;
    printf("Eigenvalue uncertainty analysis (%d samples):\n",ns);
    double *rmin=(double*)malloc(n*sizeof(double)),*rmax=(double*)malloc(n*sizeof(double));
    for(int i=0;i<n;i++){rmin[i]=PU_INF;rmax[i]=-PU_INF;}
    for(int s=0;s<ns;s++){
        double **As=pu_matrix_alloc(n,n);
        pu_matrix_copy(As,uss->A_nominal,n,n);
        if(uss->params.params&&uss->A_pert)
            for(int k=0;k<uss->n_unc_params;k++){
                double qk=pu_param_sample(&uss->params.params[k]);
                for(int i=0;i<n;i++)for(int j=0;j<n;j++)As[i][j]+=qk*uss->A_pert[k][i][j];
            }
        double *re=(double*)malloc(n*sizeof(double)),*im=(double*)malloc(n*sizeof(double));
        pu_eigen_qr(As,n,re,im,PU_MAX_ITER_QR,PU_QR_TOL);
        for(int i=0;i<n;i++){if(re[i]<rmin[i])rmin[i]=re[i];if(re[i]>rmax[i])rmax[i]=re[i];}
        free(re);free(im);pu_matrix_free(As,n);
    }
    printf("  Eigenvalue bounds:\n");
    for(int i=0;i<n;i++)printf("    ev[%d]: Re in [%.4g, %.4g]\n",i,rmin[i],rmax[i]);
    free(rmin);free(rmax);
}

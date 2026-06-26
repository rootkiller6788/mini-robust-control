/* ============================================================================
 * pu_stability.c — Robust Stability Analysis under Parametric Uncertainty
 * ============================================================================ */
#include "pu_stability.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static bool gershgorin_stable_check(const PU_IntervalMatrix *im);
static double min_eigenvalue_sym(double **M, int n);
static bool matrix_is_hurwitz_direct(double **M, int n);

PU_ValueSet pu_value_set_compute(const PU_IntervalPolynomial *ip, double omega)
{
    PU_ValueSet vs; vs.omega = omega; vs.n_points = 4;
    vs.real_part = (double *)malloc(4 * sizeof(double));
    vs.imag_part = (double *)malloc(4 * sizeof(double));
    vs.contains_origin = false; vs.center_real = vs.center_imag = vs.radius = 0.0;
    vs.min_distance_to_origin = PU_INF;
    if (!ip || ip->degree < 0) return vs;
    double lo_r, hi_r, lo_i, hi_i;
    pu_interval_poly_eval(ip, 0.0, omega, &lo_r, &hi_r, &lo_i, &hi_i);
    vs.real_part[0] = lo_r; vs.imag_part[0] = lo_i;
    vs.real_part[1] = hi_r; vs.imag_part[1] = lo_i;
    vs.real_part[2] = hi_r; vs.imag_part[2] = hi_i;
    vs.real_part[3] = lo_r; vs.imag_part[3] = hi_i;
    vs.center_real = 0.5 * (lo_r + hi_r);
    vs.center_imag = 0.5 * (lo_i + hi_i);
    vs.radius = 0.5 * sqrt((hi_r - lo_r) * (hi_r - lo_r) + (hi_i - lo_i) * (hi_i - lo_i));
    vs.contains_origin = (lo_r <= 0.0 && hi_r >= 0.0 && lo_i <= 0.0 && hi_i >= 0.0);
    vs.min_distance_to_origin = PU_INF;
    for (int k = 0; k < 4; k++) {
        double d = sqrt(vs.real_part[k] * vs.real_part[k] + vs.imag_part[k] * vs.imag_part[k]);
        if (d < vs.min_distance_to_origin) vs.min_distance_to_origin = d;
    }
    return vs;
}

PU_ValueSet pu_value_set_grid(const PU_IntervalPolynomial *ip, double omega, int grid_points)
{
    PU_ValueSet vs; vs.omega = omega; vs.n_points = grid_points;
    vs.real_part = (double *)malloc((size_t)grid_points * sizeof(double));
    vs.imag_part = (double *)malloc((size_t)grid_points * sizeof(double));
    vs.contains_origin = false; vs.center_real = vs.center_imag = vs.radius = 0.0;
    vs.min_distance_to_origin = PU_INF;
    if (!ip || ip->degree < 0 || grid_points <= 0) return vs;
    double min_re = PU_INF, max_re = -PU_INF, min_im = PU_INF, max_im = -PU_INF;
    for (int k = 0; k < grid_points; k++) {
        double *coeff = (double *)malloc((size_t)(ip->degree + 1) * sizeof(double));
        pu_interval_poly_sample(ip, coeff);
        double re = 0.0, im = 0.0, om_pow = 1.0;
        for (int i = 0; i <= ip->degree; i++) {
            double c = coeff[i] * om_pow;
            switch (i % 4) { case 0: re += c; break; case 1: im += c; break; case 2: re -= c; break; case 3: im -= c; break; }
            om_pow *= omega;
        }
        vs.real_part[k] = re; vs.imag_part[k] = im;
        if (re < min_re) min_re = re; if (re > max_re) max_re = re;
        if (im < min_im) min_im = im; if (im > max_im) max_im = im;
        double d = sqrt(re * re + im * im);
        if (d < vs.min_distance_to_origin) vs.min_distance_to_origin = d;
        free(coeff);
    }
    vs.center_real = 0.5 * (min_re + max_re); vs.center_imag = 0.5 * (min_im + max_im);
    vs.radius = 0.5 * sqrt((max_re - min_re) * (max_re - min_re) + (max_im - min_im) * (max_im - min_im));
    vs.contains_origin = (min_re <= 0.0 && max_re >= 0.0 && min_im <= 0.0 && max_im >= 0.0);
    return vs;
}

void pu_value_set_free(PU_ValueSet *vs) { if (vs) { free(vs->real_part); free(vs->imag_part); vs->real_part = vs->imag_part = NULL; vs->n_points = 0; } }
void pu_value_set_print(const PU_ValueSet *vs) {
    if (!vs) return;
    printf("ValueSet w=%.3g c=(%.3g,%.3g) r=%.3g origin=%s dist=%.3g\n",
           vs->omega, vs->center_real, vs->center_imag, vs->radius,
           vs->contains_origin ? "Y" : "N", vs->min_distance_to_origin);
}

PU_StabilityStatus pu_zero_exclusion_test(const PU_IntervalPolynomial *ip,
    double omega_min, double omega_max, int n_points)
{
    if (!ip || ip->degree < 1) return PU_INDETERMINATE;
    if (!pu_is_hurwitz_stable(ip->coeff_nominal, ip->degree)) return PU_UNSTABLE;
    double d_omega = (omega_max - omega_min) / (double)(n_points - 1);
    for (int k = 0; k < n_points; k++) {
        double omega = omega_min + (double)k * d_omega;
        PU_ValueSet vs = pu_value_set_compute(ip, omega);
        if (vs.contains_origin) { pu_value_set_free(&vs); return PU_UNSTABLE; }
        pu_value_set_free(&vs);
    }
    return PU_STABLE;
}

double pu_zero_exclusion_margin(const PU_IntervalPolynomial *ip, double *omega_crit,
    double omega_min, double omega_max, int n_points)
{
    if (!ip || !omega_crit) return -1.0;
    double min_dist = PU_INF; *omega_crit = 0.0;
    double d_omega = (omega_max - omega_min) / (double)(n_points - 1);
    for (int k = 0; k < n_points; k++) {
        double omega = omega_min + (double)k * d_omega;
        PU_ValueSet vs = pu_value_set_compute(ip, omega);
        if (vs.min_distance_to_origin < min_dist) { min_dist = vs.min_distance_to_origin; *omega_crit = omega; }
        pu_value_set_free(&vs);
    }
    return min_dist;
}

static void polytope_build_edges(PU_PolytopicPolynomial *pp) {
    int N = pp->n_vertices; pp->n_edges = N * (N - 1) / 2;
    pp->edge_pairs = (int *)malloc((size_t)(2 * pp->n_edges) * sizeof(int));
    int idx = 0;
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) { pp->edge_pairs[2*idx]=i; pp->edge_pairs[2*idx+1]=j; idx++; }
}

PU_PolytopicPolynomial pu_polyfamily_create(int degree, int n_vertices, double **vertices) {
    PU_PolytopicPolynomial pp; pp.degree = degree; pp.n_vertices = n_vertices;
    pp.n_edges = 0; pp.edge_pairs = NULL;
    pp.vertices = (PU_VertexPolynomial *)malloc((size_t)n_vertices * sizeof(PU_VertexPolynomial));
    pp.lambda_center = (double *)calloc((size_t)n_vertices, sizeof(double));
    for (int v = 0; v < n_vertices; v++) {
        pp.vertices[v].degree = degree; pp.vertices[v].vertex_index = v;
        pp.vertices[v].coeff = (double *)malloc((size_t)(degree+1) * sizeof(double));
        pp.vertices[v].barycentric = NULL;
        for (int i = 0; i <= degree; i++) pp.vertices[v].coeff[i] = vertices[v][i];
        pp.lambda_center[v] = 1.0 / (double)n_vertices;
    }
    polytope_build_edges(&pp); return pp;
}

void pu_polyfamily_free(PU_PolytopicPolynomial *pp) {
    if (!pp) return;
    if (pp->vertices) { for (int v=0; v<pp->n_vertices; v++) { free(pp->vertices[v].coeff); free(pp->vertices[v].barycentric); } free(pp->vertices); pp->vertices=NULL; }
    free(pp->edge_pairs); pp->edge_pairs=NULL; free(pp->lambda_center); pp->lambda_center=NULL;
    pp->n_vertices=0; pp->n_edges=0;
}

bool pu_edge_is_stable(const PU_PolytopicPolynomial *pp, int i, int j, int n_divs) {
    if (!pp || i>=pp->n_vertices || j>=pp->n_vertices) return false;
    int n = pp->degree + 1; double *c = (double *)malloc((size_t)n * sizeof(double)); bool stable = true;
    for (int d = 0; d <= n_divs && stable; d++) {
        double lambda = (double)d / (double)n_divs;
        for (int k = 0; k < n; k++) c[k] = (1.0-lambda)*pp->vertices[i].coeff[k] + lambda*pp->vertices[j].coeff[k];
        if (!pu_is_hurwitz_stable(c, pp->degree)) stable = false;
    }
    free(c); return stable;
}

PU_StabilityStatus pu_edge_theorem_test(const PU_PolytopicPolynomial *pp, int n_divs) {
    if (!pp || pp->n_vertices < 1) return PU_INDETERMINATE;
    bool any_stable = false;
    for (int v = 0; v < pp->n_vertices; v++)
        if (pu_is_hurwitz_stable(pp->vertices[v].coeff, pp->degree)) { any_stable = true; break; }
    if (!any_stable) return PU_UNSTABLE;
    for (int e = 0; e < pp->n_edges; e++) {
        int vi = pp->edge_pairs[2*e], vj = pp->edge_pairs[2*e+1];
        if (!pu_edge_is_stable(pp, vi, vj, n_divs)) return PU_UNSTABLE;
    }
    return PU_STABLE;
}

PU_UncertainStateSpace pu_uss_create(int ns, int ni, int no, double **An, double **Bn, double **Cn, double **Dn, PU_UncertaintyType ut) {
    PU_UncertainStateSpace uss; uss.n_states=ns; uss.n_inputs=ni; uss.n_outputs=no; uss.unc_type=ut;
    uss.n_unc_params=0; uss.is_continuous=true; uss.A_pert=NULL; uss.B_pert=NULL;
    uss.A_nominal=pu_matrix_alloc(ns,ns); uss.B_nominal=pu_matrix_alloc(ns,ni);
    uss.C_nominal=pu_matrix_alloc(no,ns); uss.D_nominal=pu_matrix_alloc(no,ni);
    if (An) pu_matrix_copy(uss.A_nominal,An,ns,ns); if (Bn) pu_matrix_copy(uss.B_nominal,Bn,ns,ni);
    if (Cn) pu_matrix_copy(uss.C_nominal,Cn,no,ns); if (Dn) pu_matrix_copy(uss.D_nominal,Dn,no,ni);
    uss.params.n_params=0; uss.params.params=NULL; uss.params.description=NULL;
    return uss;
}

void pu_uss_free(PU_UncertainStateSpace *uss) {
    if (!uss) return;
    pu_matrix_free(uss->A_nominal,uss->n_states); pu_matrix_free(uss->B_nominal,uss->n_states);
    pu_matrix_free(uss->C_nominal,uss->n_outputs); pu_matrix_free(uss->D_nominal,uss->n_outputs);
    if (uss->A_pert) { for (int k=0;k<uss->n_unc_params;k++) pu_matrix_free(uss->A_pert[k],uss->n_states); free(uss->A_pert); }
    if (uss->B_pert) { for (int k=0;k<uss->n_unc_params;k++) pu_matrix_free(uss->B_pert[k],uss->n_states); free(uss->B_pert); }
    pu_param_vector_free(&uss->params);
}

int pu_uss_add_perturbation(PU_UncertainStateSpace *uss, double **Ai, double **Bi) {
    if (!uss) return -1;
    int k = uss->n_unc_params;
    uss->A_pert = (double ***)realloc(uss->A_pert, (size_t)(k+1)*sizeof(double*));
    uss->B_pert = (double ***)realloc(uss->B_pert, (size_t)(k+1)*sizeof(double*));
    if (!uss->A_pert || !uss->B_pert) return -1;
    uss->A_pert[k]=pu_matrix_alloc(uss->n_states,uss->n_states);
    uss->B_pert[k]=pu_matrix_alloc(uss->n_states,uss->n_inputs);
    if (Ai) pu_matrix_copy(uss->A_pert[k],Ai,uss->n_states,uss->n_states);
    if (Bi) pu_matrix_copy(uss->B_pert[k],Bi,uss->n_states,uss->n_inputs);
    uss->n_unc_params++; return 0;
}

void pu_uss_eval(const PU_UncertainStateSpace *uss, const double *q, double **Ao, double **Bo) {
    if (!uss) return;
    pu_matrix_copy(Ao,uss->A_nominal,uss->n_states,uss->n_states);
    pu_matrix_copy(Bo,uss->B_nominal,uss->n_states,uss->n_inputs);
    if (!q) return;
    for (int k = 0; k < uss->n_unc_params; k++)
        for (int i = 0; i < uss->n_states; i++) {
            for (int j = 0; j < uss->n_states; j++) Ao[i][j] += q[k] * uss->A_pert[k][i][j];
            for (int j = 0; j < uss->n_inputs; j++) Bo[i][j] += q[k] * uss->B_pert[k][i][j];
        }
}

static double min_eigenvalue_sym(double **M, int n) {
    double *x=(double*)malloc((size_t)n*sizeof(double)), *Mx=(double*)malloc((size_t)n*sizeof(double));
    for(int i=0;i<n;i++) x[i]=1.0/sqrt((double)n);
    double lambda=0.0;
    for(int iter=0;iter<200;iter++){
        for(int i=0;i<n;i++){Mx[i]=0.0;for(int j=0;j<n;j++)Mx[i]+=M[i][j]*x[j];}
        double num=0.0,den=0.0;
        for(int i=0;i<n;i++){num+=x[i]*Mx[i];den+=x[i]*x[i];}
        double nl=num/den; if(fabs(nl-lambda)<1e-12){lambda=nl;break;} lambda=nl;
        double norm=0.0; for(int i=0;i<n;i++)norm+=Mx[i]*Mx[i]; norm=sqrt(norm);
        if(norm<1e-15)break; for(int i=0;i<n;i++)x[i]=Mx[i]/norm;
    }
    free(x);free(Mx);return lambda;
}

static bool matrix_is_hurwitz_direct(double **M, int n) {
    double *re=(double*)malloc((size_t)n*sizeof(double)),*im=(double*)malloc((size_t)n*sizeof(double));
    pu_eigen_qr(M,n,re,im,PU_MAX_ITER_QR,PU_QR_TOL); bool stable=true;
    for(int i=0;i<n;i++)if(re[i]>=-PU_EPS){stable=false;break;}
    free(re);free(im);return stable;
}

PU_StabilityStatus pu_quadratic_stability_test(const PU_UncertainStateSpace *uss, PU_QuadraticStabilityLMI *result) {
    if(!uss||!result)return PU_INDETERMINATE; int n=uss->n_states;
    result->n_states=n;result->n_vertices=uss->n_unc_params;
    result->P=pu_matrix_alloc(n,n);result->is_feasible=false;result->iterations=0;result->t_min=PU_INF;
    if(uss->n_unc_params==0){
        double **Q=pu_matrix_eye(n); int ret=pu_lyapunov_solve(uss->A_nominal,Q,result->P,n);
        pu_matrix_free(Q,n); result->is_feasible=(ret==0)&&pu_is_positive_definite(result->P,n);
        return result->is_feasible?PU_STABLE:PU_INDETERMINATE;
    }
    result->is_feasible=true;
    for(int v=0;v<uss->n_unc_params&&result->is_feasible;v++){
        double **Av=pu_matrix_alloc(n,n); pu_matrix_copy(Av,uss->A_nominal,n,n);
        if(uss->A_pert&&uss->A_pert[v])for(int i=0;i<n;i++)for(int j=0;j<n;j++)Av[i][j]+=uss->A_pert[v][i][j];
        if(!matrix_is_hurwitz_direct(Av,n))result->is_feasible=false;
        pu_matrix_free(Av,n);
    }
    if(result->is_feasible){double**Q=pu_matrix_eye(n);pu_lyapunov_solve(uss->A_nominal,Q,result->P,n);pu_matrix_free(Q,n);}
    result->iterations=1;return result->is_feasible?PU_STABLE:PU_UNSTABLE;
}

double pu_quadratic_stability_margin(const PU_UncertainStateSpace *uss) {
    if(!uss||uss->n_unc_params==0)return PU_INF; int n=uss->n_states;
    double lo=0.0,hi=10.0;int found=0;
    for(int t=0;t<10;t++){bool all=true;
        for(int v=0;v<uss->n_unc_params&&all;v++){double**Av=pu_matrix_alloc(n,n);pu_matrix_copy(Av,uss->A_nominal,n,n);
            if(uss->A_pert&&uss->A_pert[v])for(int i=0;i<n;i++)for(int j=0;j<n;j++)Av[i][j]+=uss->A_pert[v][i][j];
            for(int i=0;i<n;i++)Av[i][i]+=hi; if(!matrix_is_hurwitz_direct(Av,n))all=false;pu_matrix_free(Av,n);
        }if(!all){found=1;break;}hi*=2.0;
    }if(!found)return PU_INF;
    for(int iter=0;iter<50;iter++){double mid=0.5*(lo+hi);if(hi-lo<1e-8)break;bool all=true;
        for(int v=0;v<uss->n_unc_params&&all;v++){double**Av=pu_matrix_alloc(n,n);pu_matrix_copy(Av,uss->A_nominal,n,n);
            if(uss->A_pert&&uss->A_pert[v])for(int i=0;i<n;i++)for(int j=0;j<n;j++)Av[i][j]+=uss->A_pert[v][i][j];
            for(int i=0;i<n;i++)Av[i][i]+=mid; if(!matrix_is_hurwitz_direct(Av,n))all=false;pu_matrix_free(Av,n);
        }if(all)lo=mid;else hi=mid;
    }return lo;
}

PU_LyapunovData pu_lyapunov_solve_ct(const double **Ai, double **Q, int n) {
    PU_LyapunovData ld; ld.n_states=n; ld.P=pu_matrix_alloc(n,n); ld.Q=pu_matrix_alloc(n,n);
    ld.min_eigenvalue_P=0.0; ld.condition_number=PU_INF; if(!Ai||!Q)return ld;
    double **A=pu_matrix_alloc(n,n); for(int i=0;i<n;i++)for(int j=0;j<n;j++)A[i][j]=Ai[i][j];
    pu_matrix_copy(ld.Q,Q,n,n); int ret=pu_lyapunov_solve(A,Q,ld.P,n);
    if(ret==0&&pu_is_positive_definite(ld.P,n)){
        ld.min_eigenvalue_P=min_eigenvalue_sym(ld.P,n);
        double*re=(double*)malloc((size_t)n*sizeof(double)),*im=(double*)malloc((size_t)n*sizeof(double));
        pu_eigen_qr(ld.P,n,re,im,PU_MAX_ITER_QR,PU_QR_TOL);double me=0.0;
        for(int i=0;i<n;i++)if(fabs(re[i])>me)me=fabs(re[i]);free(re);free(im);
        ld.condition_number=(ld.min_eigenvalue_P>PU_EPS)?me/ld.min_eigenvalue_P:PU_INF;
    }pu_matrix_free(A,n);return ld;
}
void pu_lyapunov_data_free(PU_LyapunovData *ld) {
    if(!ld)return;pu_matrix_free(ld->P,ld->n_states);pu_matrix_free(ld->Q,ld->n_states);ld->P=ld->Q=NULL;
}

PU_RobustStabilityResult pu_robust_stability_radius(const PU_UncertainStateSpace *uss, double omin, double omax, int nfp) {
    PU_RobustStabilityResult r; r.stability_radius=PU_INF; r.worst_case_real_part=-PU_INF;
    r.h_infinity_norm=0.0; r.mu_bounds=NULL; r.n_frequency_points=nfp; r.frequency_grid=NULL;
    r.status=PU_INDETERMINATE; r.method_used=NULL;
    if(!uss||nfp<=0)return r; int n=uss->n_states;
    double dw=(omax-omin)/(double)(nfp-1); r.frequency_grid=(double*)malloc((size_t)nfp*sizeof(double));
    r.mu_bounds=(double*)malloc(2*sizeof(double)); r.mu_bounds[0]=r.mu_bounds[1]=PU_INF;
    double mr=-PU_INF;
    for(int k=0;k<nfp;k++){double w=omin+(double)k*dw;r.frequency_grid[k]=w;
        double*re=(double*)malloc((size_t)n*sizeof(double)),*im=(double*)malloc((size_t)n*sizeof(double));
        pu_eigen_qr(uss->A_nominal,n,re,im,PU_MAX_ITER_QR,PU_QR_TOL);
        for(int i=0;i<n;i++)if(re[i]>mr)mr=re[i];free(re);free(im);
    }
    r.worst_case_real_part=mr; r.status=(mr<-PU_EPS)?PU_STABLE:PU_UNSTABLE;
    r.stability_radius=(mr<-PU_EPS)?((mr<-1.0)?PU_INF:-1.0/mr):0.0;
    r.method_used=strdup("Freq-grid eigenvalue analysis"); return r;
}
void pu_robust_stability_result_free(PU_RobustStabilityResult *r) { if(!r)return;free(r->frequency_grid);free(r->mu_bounds);free(r->method_used);r->frequency_grid=r->mu_bounds=NULL;r->method_used=NULL; }
void pu_robust_stability_result_print(const PU_RobustStabilityResult *r) {
    if(!r)return; printf("RobustStability: status=%d radius=%.4g worstRe=%.4g method=%s\n",
        r->status, r->stability_radius, r->worst_case_real_part, r->method_used?r->method_used:"N/A");
}

PU_IntervalMatrix pu_interval_matrix_create(int rows, int cols, double **lo, double **hi) {
    PU_IntervalMatrix im; im.rows=rows; im.cols=cols;
    im.lower=pu_matrix_alloc(rows,cols); im.upper=pu_matrix_alloc(rows,cols);
    if(lo) pu_matrix_copy(im.lower,lo,rows,cols);
    if(hi) pu_matrix_copy(im.upper,hi,rows,cols);
    return im;
}
void pu_interval_matrix_free(PU_IntervalMatrix *im) { if(!im)return;pu_matrix_free(im->lower,im->rows);pu_matrix_free(im->upper,im->rows);im->lower=im->upper=NULL; }
void pu_interval_matrix_sample(const PU_IntervalMatrix *im, double **out) {
    if(!im||!out) return;
    for(int i=0;i<im->rows;i++)
        for(int j=0;j<im->cols;j++)
            out[i][j]=im->lower[i][j]+((double)rand()/RAND_MAX)*(im->upper[i][j]-im->lower[i][j]);
}

static bool gershgorin_stable_check(const PU_IntervalMatrix *im) {
    if(im->rows!=im->cols) return false;
    int n=im->rows;
    for(int i=0;i<n;i++){double dmax=im->upper[i][i],r=0.0;
        for(int j=0;j<n;j++){if(j==i)continue;r+=fmax(fabs(im->lower[i][j]),fabs(im->upper[i][j]));}
        if(dmax+r>=-PU_EPS)return false;
    }return true;
}

PU_StabilityStatus pu_interval_matrix_stability(const PU_IntervalMatrix *im) {
    if(!im||im->rows!=im->cols)return PU_INDETERMINATE;
    if(gershgorin_stable_check(im))return PU_STABLE;
    int n=im->rows; double**s=pu_matrix_alloc(n,n); int ns=0,nt=100;
    for(int k=0;k<nt;k++){pu_interval_matrix_sample(im,s);if(matrix_is_hurwitz_direct(s,n))ns++;}
    pu_matrix_free(s,n); if(ns==nt)return PU_STABLE; if(ns==0)return PU_UNSTABLE; return PU_INDETERMINATE;
}
void pu_interval_matrix_print(const PU_IntervalMatrix *im) {
    if(!im)return; printf("IntervalMatrix(%dx%d):\n",im->rows,im->cols);
    for(int i=0;i<im->rows;i++){printf("  ");for(int j=0;j<im->cols;j++)printf("[%.2f,%.2f] ",im->lower[i][j],im->upper[i][j]);printf("\n");}
}

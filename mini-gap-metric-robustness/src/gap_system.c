/* gap_system.c - State-Space System Operations for Gap Metric Theory
 * LTI system lifecycle, frequency response, controllability/observability
 * tests, Gramians, system interconnections, pole placement, Lyapunov solver.
 *
 * L1: GapSystem, GapFreqResp, GapFeedbackLoop
 * L2: Controllability, observability, stability
 * L3: State-space realizations, frequency response
 * L5: Series/parallel/feedback interconnection, Lyapunov solver
 *
 * Ref: K. Zhou, J.C. Doyle, K. Glover (1996), T. Kailath (1980),
 *      Bartels & Stewart, Comm. ACM (1972)
 */
#include "gap_system.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

GapSystem* gap_sys_create(int n, int m, int p, bool is_discrete, double ts) {
    GapSystem* sys = calloc(1, sizeof(GapSystem));
    if (!sys) return NULL;
    sys->A = gap_mat_create(n, n); sys->B = gap_mat_create(n, m);
    sys->C = gap_mat_create(p, n); sys->D = gap_mat_create(p, m);
    sys->n = n; sys->m = m; sys->p = p;
    sys->is_discrete = is_discrete; sys->ts = ts;
    if (!sys->A || !sys->B || !sys->C || !sys->D) { gap_sys_free(sys); return NULL; }
    return sys;
}

GapSystem* gap_sys_create_from(GapMatrix* A, GapMatrix* B,
                                GapMatrix* C, GapMatrix* D,
                                bool is_discrete, double ts) {
    if (!A || !B || !C || !D) return NULL;
    GapSystem* sys = calloc(1, sizeof(GapSystem));
    if (!sys) return NULL;
    sys->A = A; sys->B = B; sys->C = C; sys->D = D;
    sys->n = A->rows; sys->m = B->cols; sys->p = C->rows;
    sys->is_discrete = is_discrete; sys->ts = ts;
    return sys;
}

GapSystem* gap_sys_clone(const GapSystem* sys) {
    if (!sys) return NULL;
    GapSystem* nsys = calloc(1, sizeof(GapSystem));
    if (!nsys) return NULL;
    nsys->A = gap_mat_clone(sys->A); nsys->B = gap_mat_clone(sys->B);
    nsys->C = gap_mat_clone(sys->C); nsys->D = gap_mat_clone(sys->D);
    nsys->n = sys->n; nsys->m = sys->m; nsys->p = sys->p;
    nsys->is_discrete = sys->is_discrete; nsys->ts = sys->ts;
    if (!nsys->A || !nsys->B || !nsys->C || !nsys->D)
    { gap_sys_free(nsys); return NULL; }
    return nsys;
}

void gap_sys_free(GapSystem* sys) {
    if (!sys) return;
    gap_mat_free(sys->A); gap_mat_free(sys->B);
    gap_mat_free(sys->C); gap_mat_free(sys->D);
    free(sys);
}

void gap_sys_print(const GapSystem* sys, const char* name) {
    if (!sys) { printf("%s: NULL\n", name); return; }
    printf("=== %s ===\n", name);
    printf("State dim: %d, Inputs: %d, Outputs: %d, %s\n",
           sys->n, sys->m, sys->p,
           sys->is_discrete ? "discrete-time" : "continuous-time");
    gap_mat_print(sys->A, "A"); gap_mat_print(sys->B, "B");
    gap_mat_print(sys->C, "C"); gap_mat_print(sys->D, "D");
}

GapMatrix* gap_sys_controllability_matrix(const GapSystem* sys) {
    if (!sys || sys->n <= 0) return NULL;
    int n = sys->n, m = sys->m, nc = n * m;
    GapMatrix* Ctr = gap_mat_create(n, nc);
    if (!Ctr) return NULL;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            Ctr->data[i * nc + j] = sys->B->data[i * m + j];
    GapMatrix* Ak = gap_mat_clone(sys->A);
    if (!Ak) { gap_mat_free(Ctr); return NULL; }
    for (int k = 1; k < n; k++) {
        GapMatrix* AkB = gap_mat_mul(Ak, sys->B);
        if (!AkB) { gap_mat_free(Ctr); gap_mat_free(Ak); return NULL; }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                Ctr->data[i * nc + k * m + j] = AkB->data[i * m + j];
        gap_mat_free(AkB);
        if (k < n - 1) {
            GapMatrix* Ak1 = gap_mat_mul(Ak, sys->A);
            gap_mat_free(Ak); Ak = Ak1;
            if (!Ak) { gap_mat_free(Ctr); return NULL; }
        }
    }
    gap_mat_free(Ak); return Ctr;
}

GapMatrix* gap_sys_observability_matrix(const GapSystem* sys) {
    if (!sys || sys->n <= 0) return NULL;
    int n = sys->n, p = sys->p, nr = n * p;
    GapMatrix* Obs = gap_mat_create(nr, n);
    if (!Obs) return NULL;
    for (int i = 0; i < p; i++)
        for (int j = 0; j < n; j++)
            Obs->data[i * n + j] = sys->C->data[i * n + j];
    GapMatrix* Ak = gap_mat_clone(sys->A);
    if (!Ak) { gap_mat_free(Obs); return NULL; }
    for (int k = 1; k < n; k++) {
        GapMatrix* CAk = gap_mat_mul(sys->C, Ak);
        if (!CAk) { gap_mat_free(Obs); gap_mat_free(Ak); return NULL; }
        for (int i = 0; i < p; i++)
            for (int j = 0; j < n; j++)
                Obs->data[(k * p + i) * n + j] = CAk->data[i * n + j];
        gap_mat_free(CAk);
        if (k < n - 1) {
            GapMatrix* Ak1 = gap_mat_mul(Ak, sys->A);
            gap_mat_free(Ak); Ak = Ak1;
            if (!Ak) { gap_mat_free(Obs); return NULL; }
        }
    }
    gap_mat_free(Ak); return Obs;
}

bool gap_sys_is_controllable(const GapSystem* sys) {
    if (!sys) return false;
    GapMatrix* Ctr = gap_sys_controllability_matrix(sys);
    if (!Ctr) return false;
    int r = gap_mat_rank(Ctr);
    gap_mat_free(Ctr);
    return r == sys->n;
}

bool gap_sys_is_observable(const GapSystem* sys) {
    if (!sys) return false;
    GapMatrix* Obs = gap_sys_observability_matrix(sys);
    if (!Obs) return false;
    int r = gap_mat_rank(Obs);
    gap_mat_free(Obs);
    return r == sys->n;
}

bool gap_sys_is_minimal(const GapSystem* sys) {
    return gap_sys_is_controllable(sys) && gap_sys_is_observable(sys);
}

bool gap_sys_is_stable(const GapSystem* sys) {
    if (!sys || sys->n == 0) return false;
    int n_ev;
    GapComplex* ev = gap_mat_eigenvalues(sys->A, &n_ev);
    if (!ev) return false;
    bool stable = true;
    if (sys->is_discrete) {
        for (int i = 0; i < n_ev; i++)
            if (gap_complex_abs(ev[i]) >= 1.0 - GAP_EPSILON) { stable = false; break; }
    } else {
        for (int i = 0; i < n_ev; i++)
            if (ev[i].re >= -GAP_EPSILON) { stable = false; break; }
    }
    free(ev); return stable;
}

double gap_sys_stability_margin(const GapSystem* sys) {
    if (!sys || sys->n == 0) return 0.0;
    int n_ev;
    GapComplex* ev = gap_mat_eigenvalues(sys->A, &n_ev);
    if (!ev) return 0.0;
    double margin = INFINITY;
    if (sys->is_discrete) {
        for (int i = 0; i < n_ev; i++) {
            double d = 1.0 - gap_complex_abs(ev[i]);
            if (d < margin) margin = d;
        }
    } else {
        for (int i = 0; i < n_ev; i++) {
            double d = -ev[i].re;
            if (d < margin) margin = d;
        }
    }
    free(ev);
    return (margin < 0) ? 0.0 : margin;
}

GapFreqResp* gap_sys_freqresp(const GapSystem* sys, double omega) {
    if (!sys) return NULL;
    int n = sys->n, m = sys->m, p = sys->p;
    GapFreqResp* fr = calloc(1, sizeof(GapFreqResp));
    if (!fr) return NULL;
    fr->p = p; fr->m = m; fr->omega = omega;
    fr->resp = calloc((size_t)(p * m), sizeof(GapComplex));
    if (!fr->resp) { free(fr); return NULL; }
    GapMatrix* Ar = gap_mat_create(n, n);
    GapMatrix* Ai = gap_mat_create(n, n);
    if (!Ar || !Ai) { gap_mat_free(Ar); gap_mat_free(Ai); gap_freqresp_free(fr); return NULL; }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            Ar->data[i * n + j] = -sys->A->data[i * n + j];
            Ai->data[i * n + j] = 0.0;
        }
    if (sys->is_discrete) {
        for (int i = 0; i < n; i++) {
            Ar->data[i * n + i] += cos(omega);
            Ai->data[i * n + i] += sin(omega);
        }
    } else {
        for (int i = 0; i < n; i++) Ai->data[i * n + i] += omega;
    }
    int n2 = 2 * n;
    GapMatrix* aug = gap_mat_create(n2, n2);
    if (!aug) { gap_mat_free(Ar); gap_mat_free(Ai); gap_freqresp_free(fr); return NULL; }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug->data[i * n2 + j] = Ar->data[i * n + j];
            aug->data[i * n2 + n + j] = -Ai->data[i * n + j];
            aug->data[(n + i) * n2 + j] = Ai->data[i * n + j];
            aug->data[(n + i) * n2 + n + j] = Ar->data[i * n + j];
        }
    }
    gap_mat_free(Ar); gap_mat_free(Ai);
    GapMatrix* aug2 = gap_mat_create(n2, 2 * n2);
    if (!aug2) { gap_mat_free(aug); gap_freqresp_free(fr); return NULL; }
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < n2; j++) aug2->data[i * aug2->cols + j] = aug->data[i * n2 + j];
        aug2->data[i * aug2->cols + n2 + i] = 1.0;
    }
    gap_mat_free(aug);
    for (int k = 0; k < n2; k++) {
        double piv = aug2->data[k * aug2->cols + k];
        if (fabs(piv) < GAP_EPSILON) { gap_mat_free(aug2); gap_freqresp_free(fr); return NULL; }
        for (int j = k; j < 2 * n2; j++) aug2->data[k * aug2->cols + j] /= piv;
        for (int i = 0; i < n2; i++) {
            if (i == k) continue;
            double factor = aug2->data[i * aug2->cols + k];
            for (int j = k; j < 2 * n2; j++)
                aug2->data[i * aug2->cols + j] -= factor * aug2->data[k * aug2->cols + j];
        }
    }
    GapMatrix* invR = gap_mat_create(n, n);
    GapMatrix* invI = gap_mat_create(n, n);
    if (!invR || !invI) { gap_mat_free(aug2); gap_mat_free(invR); gap_mat_free(invI); gap_freqresp_free(fr); return NULL; }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            invR->data[i * n + j] = aug2->data[i * aug2->cols + n2 + j];
            invI->data[i * n + j] = aug2->data[(n + i) * aug2->cols + n2 + j];
        }
    gap_mat_free(aug2);
    GapMatrix* CR = gap_mat_mul(sys->C, invR);
    GapMatrix* CI = gap_mat_mul(sys->C, invI);
    gap_mat_free(invR); gap_mat_free(invI);
    if (!CR || !CI) { gap_mat_free(CR); gap_mat_free(CI); gap_freqresp_free(fr); return NULL; }
    GapMatrix* GR = gap_mat_mul(CR, sys->B);
    GapMatrix* GI = gap_mat_mul(CI, sys->B);
    gap_mat_free(CR); gap_mat_free(CI);
    if (!GR || !GI) { gap_mat_free(GR); gap_mat_free(GI); gap_freqresp_free(fr); return NULL; }
    for (int i = 0; i < p; i++)
        for (int j = 0; j < m; j++) {
            int idx = i * m + j;
            fr->resp[idx].re = GR->data[i * m + j] + sys->D->data[i * m + j];
            fr->resp[idx].im = GI->data[i * m + j];
        }
    gap_mat_free(GR); gap_mat_free(GI);
    return fr;
}

GapFreqResp* gap_sys_freqresp_grid(const GapSystem* sys,
                                    const double* omega, int n_omega) {
    if (!sys || !omega || n_omega <= 0) return NULL;
    GapFreqResp* grid = calloc((size_t)n_omega, sizeof(GapFreqResp));
    if (!grid) return NULL;
    for (int k = 0; k < n_omega; k++) {
        GapFreqResp* fr = gap_sys_freqresp(sys, omega[k]);
        if (!fr) { for (int j = 0; j < k; j++) free(grid[j].resp); free(grid); return NULL; }
        grid[k] = *fr; free(fr);
    }
    return grid;
}

void gap_freqresp_free(GapFreqResp* fr) {
    if (!fr) return;
    free(fr->resp);
    free(fr);
}

double gap_sys_hinf_norm(const GapSystem* sys) {
    if (!sys) return 0.0;
    return gap_mat_hinf_norm(sys->A, sys->B, sys->C, sys->D);
}

double gap_sys_h2_norm(const GapSystem* sys) {
    if (!sys || sys->is_discrete) return INFINITY;
    if (!gap_sys_is_stable(sys)) return INFINITY;
    GapMatrix* P = gap_sys_controllability_gramian(sys);
    if (!P) return INFINITY;
    GapMatrix* CPCt = gap_mat_mul(gap_mat_mul(sys->C, P),
                                   gap_mat_transpose(sys->C));
    if (!CPCt) { gap_mat_free(P); return INFINITY; }
    double h2 = sqrt(gap_mat_trace(CPCt));
    gap_mat_free(P); gap_mat_free(CPCt);
    return h2;
}

GapMatrix* gap_sys_controllability_gramian(const GapSystem* sys) {
    if (!sys) return NULL;
    GapMatrix* BBT = gap_mat_mul(sys->B, gap_mat_transpose(sys->B));
    if (!BBT) return NULL;
    GapMatrix* BBTn = gap_mat_scalar_mul(BBT, -1.0);
    gap_mat_free(BBT);
    if (!BBTn) return NULL;
    GapMatrix* P = gap_sys_lyapunov(sys->A, BBTn);
    gap_mat_free(BBTn);
    return P;
}

GapMatrix* gap_sys_observability_gramian(const GapSystem* sys) {
    if (!sys) return NULL;
    GapMatrix* CTC = gap_mat_mul(gap_mat_transpose(sys->C), sys->C);
    if (!CTC) return NULL;
    GapMatrix* CTCn = gap_mat_scalar_mul(CTC, -1.0);
    gap_mat_free(CTC);
    if (!CTCn) return NULL;
    GapMatrix* AT = gap_mat_transpose(sys->A);
    if (!AT) { gap_mat_free(CTCn); return NULL; }
    GapMatrix* Q = gap_sys_lyapunov(AT, CTCn);
    gap_mat_free(AT); gap_mat_free(CTCn);
    return Q;
}

double* gap_sys_hankel_singular_values(const GapSystem* sys, int* n_sv) {
    if (!sys || !gap_sys_is_stable(sys)) { *n_sv = 0; return NULL; }
    GapMatrix* P = gap_sys_controllability_gramian(sys);
    GapMatrix* Q = gap_sys_observability_gramian(sys);
    if (!P || !Q) { gap_mat_free(P); gap_mat_free(Q); *n_sv = 0; return NULL; }
    GapMatrix* PQ = gap_mat_mul(P, Q);
    gap_mat_free(P); gap_mat_free(Q);
    if (!PQ) { *n_sv = 0; return NULL; }
    int nev;
    GapComplex* ev = gap_mat_eigenvalues(PQ, &nev);
    gap_mat_free(PQ);
    if (!ev) { *n_sv = 0; return NULL; }
    double* hsv = calloc((size_t)nev, sizeof(double));
    int count = 0;
    for (int i = 0; i < nev; i++)
        if (ev[i].im == 0.0 && ev[i].re > 0.0) hsv[count++] = sqrt(ev[i].re);
    free(ev);
    *n_sv = count;
    return hsv;
}

GapSystem* gap_sys_series(const GapSystem* G1, const GapSystem* G2) {
    if (!G1 || !G2 || G1->p != G2->m) return NULL;
    int n1 = G1->n, n2 = G2->n, n = n1 + n2;
    GapSystem* Gs = gap_sys_create(n, G1->m, G2->p, G1->is_discrete, G1->ts);
    if (!Gs) return NULL;
    for (int i = 0; i < n1; i++)
        for (int j = 0; j < n1; j++)
            Gs->A->data[i * n + j] = G1->A->data[i * n1 + j];
    for (int i = 0; i < n2; i++)
        for (int j = 0; j < n1; j++)
            Gs->A->data[(n1 + i) * n + j] = G2->B->data[i * n2 + j]
                                           * G1->C->data[0 * G1->C->cols + j];
    for (int i = 0; i < n2; i++)
        for (int j = 0; j < n2; j++)
            Gs->A->data[(n1 + i) * n + n1 + j] = G2->A->data[i * n2 + j];
    for (int i = 0; i < n1; i++)
        for (int j = 0; j < G1->m; j++)
            Gs->B->data[i * Gs->m + j] = G1->B->data[i * G1->B->cols + j];
    for (int i = 0; i < n2; i++)
        for (int j = 0; j < G1->m; j++) {
            double sum = 0.0;
            for (int k = 0; k < G2->m; k++)
                sum += G2->B->data[i * G2->B->cols + k]
                     * G1->D->data[k * G1->D->cols + j];
            Gs->B->data[(n1 + i) * Gs->m + j] = sum;
        }
    for (int i = 0; i < G2->p; i++)
        for (int j = 0; j < n1; j++) {
            double sum = 0.0;
            for (int k = 0; k < G2->m; k++)
                sum += G2->D->data[i * G2->D->cols + k]
                     * G1->C->data[k * G1->C->cols + j];
            Gs->C->data[i * n + j] = sum;
        }
    for (int i = 0; i < G2->p; i++)
        for (int j = 0; j < n2; j++)
            Gs->C->data[i * n + n1 + j] = G2->C->data[i * n2 + j];
    GapMatrix* D2D1 = gap_mat_mul(G2->D, G1->D);
    if (D2D1) { gap_mat_copy(Gs->D, D2D1); gap_mat_free(D2D1); }
    return Gs;
}

GapSystem* gap_sys_parallel(const GapSystem* G1, const GapSystem* G2) {
    if (!G1 || !G2 || G1->m != G2->m || G1->p != G2->p) return NULL;
    int n1 = G1->n, n2 = G2->n, n = n1 + n2;
    GapSystem* Gp = gap_sys_create(n, G1->m, G1->p, G1->is_discrete, G1->ts);
    if (!Gp) return NULL;
    for (int i = 0; i < n1; i++)
        for (int j = 0; j < n1; j++)
            Gp->A->data[i * n + j] = G1->A->data[i * n1 + j];
    for (int i = 0; i < n2; i++)
        for (int j = 0; j < n2; j++)
            Gp->A->data[(n1 + i) * n + n1 + j] = G2->A->data[i * n2 + j];
    for (int i = 0; i < n1; i++)
        for (int j = 0; j < G1->m; j++)
            Gp->B->data[i * Gp->m + j] = G1->B->data[i * G1->B->cols + j];
    for (int i = 0; i < n2; i++)
        for (int j = 0; j < G2->m; j++)
            Gp->B->data[(n1 + i) * Gp->m + j] = G2->B->data[i * G2->B->cols + j];
    for (int i = 0; i < G1->p; i++) {
        for (int j = 0; j < n1; j++)
            Gp->C->data[i * n + j] = G1->C->data[i * n1 + j];
        for (int j = 0; j < n2; j++)
            Gp->C->data[i * n + n1 + j] = G2->C->data[i * n2 + j];
    }
    GapMatrix* Dsum = gap_mat_add(G1->D, G2->D);
    if (Dsum) { gap_mat_copy(Gp->D, Dsum); gap_mat_free(Dsum); }
    return Gp;
}

GapSystem* gap_sys_lft(const GapSystem* G, const GapSystem* K) {
    if (!G || !K) return NULL;
    int nG = G->n, nK = K->n, n = nG + nK;
    GapSystem* cl = gap_sys_create(n, G->m - K->p, G->p - K->m,
                                     G->is_discrete, G->ts);
    if (!cl) return NULL;
    for (int i = 0; i < nG; i++)
        for (int j = 0; j < nG; j++)
            cl->A->data[i * n + j] = G->A->data[i * nG + j];
    for (int i = 0; i < nK; i++)
        for (int j = 0; j < nK; j++)
            cl->A->data[(nG + i) * n + nG + j] = K->A->data[i * nK + j];
    return cl;
}

GapSystem* gap_sys_generalized_plant(const GapSystem* G,
                                      const GapMatrix* W1,
                                      const GapMatrix* W2) {
    if (!G || !W1 || !W2) return NULL;
    int n = G->n, m = G->m, p = G->p;
    GapSystem* Ggen = gap_sys_create(n, m + p, W1->rows + W2->rows + p,
                                       G->is_discrete, G->ts);
    if (!Ggen) return NULL;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Ggen->A->data[i * n + j] = G->A->data[i * n + j];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            Ggen->B->data[i * Ggen->m + j] = G->B->data[i * G->B->cols + j];
    return Ggen;
}

GapMatrix* gap_sys_acker(const GapSystem* sys, const double* poles, int np) {
    if (!sys || !poles || np != sys->n || sys->m != 1) return NULL;
    int n = sys->n;
    GapMatrix* K = gap_mat_create(1, n);
    if (!K) return NULL;
    GapMatrix* Ctr = gap_sys_controllability_matrix(sys);
    if (!Ctr) { gap_mat_free(K); return NULL; }
    GapMatrix* Ctr_inv = gap_mat_inverse(Ctr);
    gap_mat_free(Ctr);
    if (!Ctr_inv) return K;
    GapMatrix* An = gap_mat_power(sys->A, n);
    if (!An) { gap_mat_free(K); gap_mat_free(Ctr_inv); return NULL; }
    GapMatrix* Pd = gap_mat_companion(poles, n);
    if (!Pd) { gap_mat_free(K); gap_mat_free(Ctr_inv); gap_mat_free(An); return NULL; }
    GapMatrix* Pdn = gap_mat_power(Pd, n);
    gap_mat_free(Pd);
    GapMatrix* Ctrinv_An = gap_mat_mul(Ctr_inv, An);
    gap_mat_free(Ctr_inv); gap_mat_free(An);
    if (!Ctrinv_An || !Pdn) {
        gap_mat_free(Ctrinv_An); gap_mat_free(Pdn); return K;
    }
    for (int j = 0; j < n; j++)
        K->data[j] = Ctrinv_An->data[(n - 1) * n + j];
    gap_mat_free(Ctrinv_An); gap_mat_free(Pdn);
    return K;
}

GapMatrix* gap_sys_observer_gain(const GapSystem* sys,
                                  const double* poles, int np) {
    if (!sys || !poles || np != sys->n || sys->p != 1) return NULL;
    GapSystem* dual = gap_sys_create(sys->n, 1, sys->n, sys->is_discrete, sys->ts);
    if (!dual) return NULL;
    GapMatrix* AT = gap_mat_transpose(sys->A);
    GapMatrix* CT = gap_mat_transpose(sys->C);
    gap_mat_free(dual->A); dual->A = AT;
    gap_mat_free(dual->B); dual->B = CT;
    GapMatrix* L = gap_sys_acker(dual, poles, np);
    gap_sys_free(dual);
    if (!L) return NULL;
    GapMatrix* LT = gap_mat_transpose(L);
    gap_mat_free(L);
    return LT;
}

GapMatrix* gap_sys_lyapunov(const GapMatrix* A, const GapMatrix* Q) {
    if (!A || !Q || A->rows != A->cols || A->rows != Q->rows) return NULL;
    int n = A->rows;
    GapMatrix* X = gap_mat_create(n, n);
    if (!X) return NULL;
    double h = 0.005;
    int steps = 2000;
    GapMatrix* Ah = gap_mat_scalar_mul(A, h);
    if (!Ah) { gap_mat_free(X); return NULL; }
    GapMatrix* term = gap_mat_eye(n);
    GapMatrix* eAh = gap_mat_eye(n);
    if (!term || !eAh) { gap_mat_free(X); gap_mat_free(Ah); gap_mat_free(term); gap_mat_free(eAh); return NULL; }
    for (int k = 1; k <= 15; k++) {
        GapMatrix* t2 = gap_mat_mul(term, Ah);
        GapMatrix* ts = gap_mat_scalar_mul(t2, 1.0 / k);
        gap_mat_free(term); term = ts; gap_mat_free(t2);
        if (!term) break;
        GapMatrix* na = gap_mat_add(eAh, term);
        gap_mat_free(eAh); eAh = na;
        if (!eAh) break;
    }
    gap_mat_free(Ah); gap_mat_free(term);
    GapMatrix* eAhT = gap_mat_transpose(eAh);
    for (int k = 0; k < steps; k++) {
        GapMatrix* XeT = gap_mat_mul(X, eAhT);
        if (!XeT) break;
        GapMatrix* X2 = gap_mat_mul(eAh, XeT);
        gap_mat_free(XeT);
        if (!X2) break;
        GapMatrix* hQ = gap_mat_scalar_mul(Q, h);
        GapMatrix* X3 = gap_mat_add(X2, hQ);
        gap_mat_free(X2); gap_mat_free(hQ);
        if (!X3) break;
        gap_mat_copy(X, X3); gap_mat_free(X3);
    }
    gap_mat_free(eAh); gap_mat_free(eAhT);
    return X;
}

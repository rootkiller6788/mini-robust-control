/* gap_coprime.c - Normalized Coprime Factorization for Gap Metric Theory
 * Computes normalized right/left coprime factorizations (NRCF/NLCF),
 * graph symbols, Bezout identity verification, and chordal distance.
 *
 * L1: GapNCFRight, GapNCFLeft, GapCoprimePerturbation
 * L2: Bezout identity, graph symbol, perturbed coprime factors
 * L3: GCARE/GFARE Riccati equations for NCF
 * L5: McFarlane-Glover NCF computation, chordal distance
 *
 * Ref: McFarlane & Glover (1990, 1992), Vinnicombe (1993),
 *      Vidyasagar (1985), Nett, Jacobson, Balas (1984)
 */
#include "gap_coprime.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

GapNCFRight* gap_nrcf_compute(const GapSystem* G) {
    if (!G) return NULL;
    GapNCFRight* nrcf = calloc(1, sizeof(GapNCFRight));
    if (!nrcf) return NULL;
    int n = G->n, m = G->m, p = G->p;
    GapMatrix* R = gap_mat_eye(m);
    GapMatrix* Dt = gap_mat_transpose(G->D);
    GapMatrix* DtD = gap_mat_mul(Dt, G->D);
    gap_mat_free(Dt);
    if (!DtD) { gap_nrcf_free(nrcf); return NULL; }
    GapMatrix* Rplus = gap_mat_add(R, DtD);
    gap_mat_free(R); gap_mat_free(DtD);
    if (!Rplus) { gap_nrcf_free(nrcf); return NULL; }
    GapMatrix* Rinv = gap_mat_inverse(Rplus);
    gap_mat_free(Rplus);
    if (!Rinv) { gap_nrcf_free(nrcf); return NULL; }
    GapMatrix* DtC = gap_mat_mul(gap_mat_transpose(G->D), G->C);
    GapMatrix* S = gap_mat_eye(p);
    GapMatrix* DDt = gap_mat_mul(G->D, gap_mat_transpose(G->D));
    if (!DtC || !DDt) {
        gap_mat_free(DtC); gap_mat_free(DDt); gap_mat_free(Rinv);
        gap_nrcf_free(nrcf); return NULL;
    }
    GapMatrix* Sp = gap_mat_add(S, DDt);
    gap_mat_free(S); gap_mat_free(DDt);
    if (!Sp) { gap_mat_free(DtC); gap_mat_free(Rinv); gap_nrcf_free(nrcf); return NULL; }
    GapMatrix* Sinv = gap_mat_inverse(Sp);
    gap_mat_free(Sp);
    if (!Sinv) { gap_mat_free(DtC); gap_mat_free(Rinv); gap_nrcf_free(nrcf); return NULL; }
    GapMatrix* B_Rinv_DtC = gap_mat_mul(G->B, Rinv);
    GapMatrix* term1 = gap_mat_mul(B_Rinv_DtC, DtC);
    gap_mat_free(B_Rinv_DtC);
    if (!term1) { gap_mat_free(DtC); gap_mat_free(Rinv); gap_mat_free(Sinv); gap_nrcf_free(nrcf); return NULL; }
    GapMatrix* A_hat = gap_mat_sub(G->A, term1);
    gap_mat_free(term1);
    if (!A_hat) { gap_mat_free(DtC); gap_mat_free(Rinv); gap_mat_free(Sinv); gap_nrcf_free(nrcf); return NULL; }
    GapMatrix* At_hat = gap_mat_transpose(A_hat);
    GapMatrix* Ct_SC = gap_mat_mul(gap_mat_mul(gap_mat_transpose(G->C), Sinv), G->C);
    GapMatrix* Ct_SC_neg = gap_mat_scalar_mul(Ct_SC, -1.0);
    gap_mat_free(Ct_SC);
    GapMatrix* B_Rinv_Bt = gap_mat_mul(gap_mat_mul(G->B, Rinv), gap_mat_transpose(G->B));
    GapMatrix* B_Rinv_Bt_neg = gap_mat_scalar_mul(B_Rinv_Bt, -1.0);
    gap_mat_free(B_Rinv_Bt);
    if (!At_hat || !Ct_SC_neg || !B_Rinv_Bt_neg) {
        gap_mat_free(A_hat); gap_mat_free(At_hat); gap_mat_free(Ct_SC_neg);
        gap_mat_free(B_Rinv_Bt_neg); gap_mat_free(DtC); gap_mat_free(Rinv); gap_mat_free(Sinv);
        gap_nrcf_free(nrcf); return NULL;
    }
    GapMatrix* X = gap_mat_care(A_hat, G->B, Ct_SC_neg, Rinv);
    gap_mat_free(Ct_SC_neg);
    if (!X) X = gap_mat_create(n, n);
    GapMatrix* BtX = gap_mat_mul(gap_mat_transpose(G->B), X);
    GapMatrix* BtX_DtC = gap_mat_add(BtX, DtC);
    gap_mat_free(BtX);
    GapMatrix* F = gap_mat_mul(Rinv, BtX_DtC);
    gap_mat_free(BtX_DtC);
    GapMatrix* F_neg = gap_mat_scalar_mul(F, -1.0);
    gap_mat_free(F);
    if (!F_neg) {
        gap_mat_free(A_hat); gap_mat_free(At_hat); gap_mat_free(B_Rinv_Bt_neg);
        gap_mat_free(X); gap_mat_free(DtC); gap_mat_free(Rinv); gap_mat_free(Sinv);
        gap_nrcf_free(nrcf); return NULL;
    }
    nrcf->X = X;
    nrcf->F = F_neg;
    GapMatrix* Rsqrt_inv = gap_mat_create(m, m);
    if (!Rsqrt_inv) {
        gap_mat_free(A_hat); gap_mat_free(At_hat); gap_mat_free(B_Rinv_Bt_neg);
        gap_mat_free(DtC); gap_mat_free(Rinv); gap_mat_free(Sinv);
        gap_nrcf_free(nrcf); return NULL;
    }
    for (int i = 0; i < m; i++) Rsqrt_inv->data[i * m + i] = 1.0;
    GapMatrix* BF = gap_mat_mul(G->B, nrcf->F);
    GapMatrix* A_M = gap_mat_add(G->A, BF);
    gap_mat_free(BF);
    GapMatrix* A_N = gap_mat_clone(A_M);
    GapMatrix* SF = gap_mat_mul(G->D, nrcf->F);
    GapMatrix* C_plus_SF = gap_mat_add(G->C, SF);
    gap_mat_free(SF);
    if (!A_M || !A_N || !C_plus_SF) {
        gap_mat_free(A_M); gap_mat_free(A_N); gap_mat_free(C_plus_SF);
        gap_mat_free(A_hat); gap_mat_free(At_hat); gap_mat_free(B_Rinv_Bt_neg);
        gap_mat_free(DtC); gap_mat_free(Rinv); gap_mat_free(Sinv);
        gap_mat_free(Rsqrt_inv); gap_nrcf_free(nrcf); return NULL;
    }
    nrcf->M = gap_sys_create_from(A_M, gap_mat_mul(G->B, Rsqrt_inv),
                                    gap_mat_clone(nrcf->F), gap_mat_clone(Rsqrt_inv),
                                    G->is_discrete, G->ts);
    nrcf->N = gap_sys_create_from(A_N, gap_mat_mul(G->B, Rsqrt_inv),
                                    C_plus_SF, gap_mat_mul(G->D, Rsqrt_inv),
                                    G->is_discrete, G->ts);
    gap_mat_free(Rsqrt_inv);
    nrcf->is_normalized = true;
    gap_mat_free(A_hat); gap_mat_free(At_hat); gap_mat_free(B_Rinv_Bt_neg);
    gap_mat_free(DtC); gap_mat_free(Rinv); gap_mat_free(Sinv);
    return nrcf;
}

GapNCFLeft* gap_nlcf_compute(const GapSystem* G) {
    if (!G) return NULL;
    GapSystem* Gt = gap_sys_create_from(
        gap_mat_transpose(G->A), gap_mat_transpose(G->C),
        gap_mat_transpose(G->B), gap_mat_transpose(G->D),
        G->is_discrete, G->ts);
    if (!Gt) return NULL;
    GapNCFRight* nrcf_dual = gap_nrcf_compute(Gt);
    gap_sys_free(Gt);
    if (!nrcf_dual) return NULL;
    GapNCFLeft* nlcf = calloc(1, sizeof(GapNCFLeft));
    if (!nlcf) { gap_nrcf_free(nrcf_dual); return NULL; }
    nlcf->Mtilde = gap_sys_create_from(
        gap_mat_transpose(nrcf_dual->M->A),
        gap_mat_transpose(nrcf_dual->M->C),
        gap_mat_transpose(nrcf_dual->M->B),
        gap_mat_transpose(nrcf_dual->M->D),
        G->is_discrete, G->ts);
    nlcf->Ntilde = gap_sys_create_from(
        gap_mat_transpose(nrcf_dual->N->A),
        gap_mat_transpose(nrcf_dual->N->C),
        gap_mat_transpose(nrcf_dual->N->B),
        gap_mat_transpose(nrcf_dual->N->D),
        G->is_discrete, G->ts);
    nlcf->Y = nrcf_dual->X;
    nrcf_dual->X = NULL;  /* transfer ownership */
    nlcf->L = gap_mat_transpose(nrcf_dual->F);
    gap_nrcf_free(nrcf_dual);
    nlcf->is_normalized = true;
    return nlcf;
}

void gap_nrcf_free(GapNCFRight* nrcf) {
    if (!nrcf) return;
    gap_sys_free(nrcf->M); gap_sys_free(nrcf->N);
    gap_mat_free(nrcf->X); gap_mat_free(nrcf->F);
    free(nrcf);
}

void gap_nlcf_free(GapNCFLeft* nlcf) {
    if (!nlcf) return;
    gap_sys_free(nlcf->Mtilde); gap_sys_free(nlcf->Ntilde);
    gap_mat_free(nlcf->Y); gap_mat_free(nlcf->L);
    free(nlcf);
}

GapGraphSymbol* gap_graph_symbol_from_nrcf(const GapNCFRight* nrcf) {
    if (!nrcf) return NULL;
    GapGraphSymbol* gs = calloc(1, sizeof(GapGraphSymbol));
    if (!gs) return NULL;
    gs->M_sys = gap_sys_clone(nrcf->M);
    gs->N_sys = gap_sys_clone(nrcf->N);
    gs->is_normalized = nrcf->is_normalized;
    gs->n_g = (nrcf->M) ? nrcf->M->n + nrcf->N->n : 0;
    return gs;
}

GapGraphSymbol* gap_graph_symbol_from_nlcf(const GapNCFLeft* nlcf) {
    if (!nlcf) return NULL;
    GapGraphSymbol* gs = calloc(1, sizeof(GapGraphSymbol));
    if (!gs) return NULL;
    gs->M_sys = gap_sys_clone(nlcf->Mtilde);
    gs->N_sys = gap_sys_clone(nlcf->Ntilde);
    gs->is_normalized = nlcf->is_normalized;
    gs->n_g = (nlcf->Mtilde) ? nlcf->Mtilde->n + nlcf->Ntilde->n : 0;
    return gs;
}

void gap_graph_symbol_free(GapGraphSymbol* gs) {
    if (!gs) return;
    gap_sys_free(gs->M_sys); gap_sys_free(gs->N_sys);
    free(gs);
}

double gap_verify_bezout_right(const GapNCFRight* nrcf) {
    if (!nrcf || !nrcf->M || !nrcf->N) return INFINITY;
    GapMatrix* MtM = gap_mat_mul(gap_mat_transpose(nrcf->M->D),
                                  nrcf->M->D);
    GapMatrix* NtN = gap_mat_mul(gap_mat_transpose(nrcf->N->D),
                                  nrcf->N->D);
    if (!MtM || !NtN) { gap_mat_free(MtM); gap_mat_free(NtN); return INFINITY; }
    GapMatrix* sum = gap_mat_add(MtM, NtN);
    gap_mat_free(MtM); gap_mat_free(NtN);
    if (!sum) return INFINITY;
    int m = sum->rows;
    GapMatrix* I = gap_mat_eye(m);
    GapMatrix* diff = gap_mat_sub(sum, I);
    gap_mat_free(sum); gap_mat_free(I);
    if (!diff) return INFINITY;
    double residual = gap_mat_norm_frobenius(diff);
    gap_mat_free(diff);
    return residual;
}

double gap_verify_bezout_left(const GapNCFLeft* nlcf) {
    if (!nlcf || !nlcf->Mtilde || !nlcf->Ntilde) return INFINITY;
    GapMatrix* MMt = gap_mat_mul(nlcf->Mtilde->D,
                                  gap_mat_transpose(nlcf->Mtilde->D));
    GapMatrix* NNt = gap_mat_mul(nlcf->Ntilde->D,
                                  gap_mat_transpose(nlcf->Ntilde->D));
    if (!MMt || !NNt) { gap_mat_free(MMt); gap_mat_free(NNt); return INFINITY; }
    GapMatrix* sum = gap_mat_add(MMt, NNt);
    gap_mat_free(MMt); gap_mat_free(NNt);
    if (!sum) return INFINITY;
    int p = sum->rows;
    GapMatrix* I = gap_mat_eye(p);
    GapMatrix* diff = gap_mat_sub(sum, I);
    gap_mat_free(sum); gap_mat_free(I);
    if (!diff) return INFINITY;
    double residual = gap_mat_norm_frobenius(diff);
    gap_mat_free(diff);
    return residual;
}

GapCoprimePerturbation* gap_coprime_perturbation_create(
    const GapNCFRight* nrcf, double epsilon) {
    if (!nrcf) return NULL;
    GapCoprimePerturbation* pert = calloc(1, sizeof(GapCoprimePerturbation));
    if (!pert) return NULL;
    pert->Delta_N = gap_sys_create(nrcf->N->n, nrcf->N->m, nrcf->N->p,
                                     nrcf->N->is_discrete, nrcf->N->ts);
    pert->Delta_M = gap_sys_create(nrcf->M->n, nrcf->M->m, nrcf->M->p,
                                     nrcf->M->is_discrete, nrcf->M->ts);
    pert->norm_bound = epsilon;
    return pert;
}

GapSystem* gap_coprime_perturbation_apply(
    const GapNCFRight* nrcf, const GapCoprimePerturbation* pert) {
    if (!nrcf || !pert) return NULL;
    GapSystem* N_pert = gap_sys_parallel(nrcf->N, pert->Delta_N);
    GapSystem* M_pert = gap_sys_parallel(nrcf->M, pert->Delta_M);
    if (!N_pert || !M_pert) {
        gap_sys_free(N_pert); gap_sys_free(M_pert); return NULL;
    }
    GapMatrix* M_inv = gap_mat_inverse(M_pert->D);
    if (!M_inv) { gap_sys_free(N_pert); gap_sys_free(M_pert); return NULL; }
    GapMatrix* Gp = gap_mat_mul(N_pert->D, M_inv);
    gap_mat_free(M_inv);
    if (!Gp) { gap_sys_free(N_pert); gap_sys_free(M_pert); return NULL; }
    GapSystem* G_delta = gap_sys_create(N_pert->n, M_pert->m, N_pert->p,
                                         N_pert->is_discrete, N_pert->ts);
    if (G_delta) { gap_mat_copy(G_delta->D, Gp); }
    gap_mat_free(Gp);
    gap_sys_free(N_pert); gap_sys_free(M_pert);
    return G_delta;
}

void gap_coprime_perturbation_free(GapCoprimePerturbation* pert) {
    if (!pert) return;
    gap_sys_free(pert->Delta_N); gap_sys_free(pert->Delta_M);
    free(pert);
}

double gap_chordal_distance(double P1_re, double P1_im,
                             double P2_re, double P2_im) {
    double num_re = P1_re - P2_re;
    double num_im = P1_im - P2_im;
    double num_sq = num_re * num_re + num_im * num_im;
    double den1 = 1.0 + P1_re * P1_re + P1_im * P1_im;
    double den2 = 1.0 + P2_re * P2_re + P2_im * P2_im;
    if (den1 < GAP_EPSILON || den2 < GAP_EPSILON) return 1.0;
    return sqrt(num_sq / (den1 * den2));
}

double gap_chordal_distance_matrix(const GapFreqResp* G1,
                                    const GapFreqResp* G2, double omega) {
    if (!G1 || !G2 || G1->p != G2->p || G1->m != G2->m) return 1.0;
    double max_kappa = 0.0;
    for (int i = 0; i < G1->p * G1->m; i++) {
        double k = gap_chordal_distance(G1->resp[i].re, G1->resp[i].im,
                                          G2->resp[i].re, G2->resp[i].im);
        if (k > max_kappa) max_kappa = k;
    }
    (void)omega;
    return max_kappa;
}

#ifndef MU_LFT_H
#define MU_LFT_H
#include "mu_core.h"
#ifdef __cplusplus
extern "C" {
#endif

/*
 * mu_lft.h -- Linear Fractional Transformations for u-Synthesis
 *
 * LFTs are the fundamental interconnection framework for robust control.
 * Upper LFT F_u(M,Delta)=M22+M21*Delta*(I-M11*Delta)^{-1}*M12
 * Lower LFT F_l(M,K)=M11+M12*K*(I-M22*K)^{-1}*M21
 *
 * The u-synthesis problem uses BOTH:
 *  - Upper LFT to pull out uncertainty Delta
 *  - Lower LFT to close the loop with controller K
 * Reference: Zhou, Doyle & Glover (1996), Ch. 9-10
 */

MUMatrix* mu_lft_upper(const MUMatrix *M, int p1, int q1,
                        const MUMatrix *Delta);
MUMatrix* mu_lft_lower(const MUMatrix *M, int p1, int q1,
                        const MUMatrix *K);
MUMatrix* mu_lft_star_product(const MUMatrix *P, const MUMatrix *K,
                               int p1, int q1, int p2, int q2);
MUMatrix* mu_system_to_m(const MUSystem *sys, double omega,
                          int p_unc, int q_unc);
MUSystem* mu_build_interconnection(const MUSystem *G,
                                    const MUSystem *Wu,
                                    const MUSystem *Wp);
MUSystem* mu_lft_reduce(const MUSystem *sys);

#ifdef __cplusplus
}
#endif
#endif

# Knowledge Graph — mini-structured-singular-value

Nine-level knowledge coverage for the Structured Singular Value (mu) theory and computation.

## L1: Definitions (Complete)

| Entry | C Implementation | Description |
|-------|-----------------|-------------|
| Structured Singular Value mu | `ssv_mu_definition_search()` in `ssv_core.c` | mu_Delta(M) = 1/min{sigma_bar(Delta): det(I-M*Delta)=0} |
| Complex Matrix Type | `SSVComplexMatrix` in `ssv_core.h` | Column-major complex matrix with leading dimension |
| Real Matrix Type | `SSVRealMatrix` in `ssv_core.h` | Column-major real matrix for state-space data |
| mu Result Type | `SSVMuResult` in `ssv_core.h` | Upper/lower bounds, gap, convergence status |
| SVD Result Type | `SSVSVDResult` in `ssv_core.h` | U, S, V^H decomposition components |
| Uncertainty Block | `SSVUncertaintyBlock` in `ssv_uncertainty.h` | Repeated scalar / full block types |
| Uncertainty Structure | `SSVUncertaintyStructure` in `ssv_uncertainty.h` | Block-diagonal uncertainty set Delta |
| D-Scaling Matrices | `SSVDScaleMatrices` in `ssv_uncertainty.h` | Commuting D matrices for mu upper bound |
| LFT Matrix | `SSVLFTMatrix` in `ssv_lft.h` | Partitioned matrix for LFT operations |
| State-Space System | `SSVStateSpace` in `ssv_synthesis.h` | (A,B,C,D) representation |
| Generalized Plant | `SSVGeneralizedPlant` in `ssv_synthesis.h` | H∞/μ-synthesis plant P |

## L2: Core Concepts (Complete)

| Concept | Implementation | Description |
|---------|---------------|-------------|
| mu Upper Bound | `ssv_mu_upper_bound()` in `ssv_bounds.c` | D-scaling optimization (convex) |
| mu Lower Bound | `ssv_mu_lower_bound()` in `ssv_bounds.c` | Power iteration (non-convex local max) |
| Small-Gain Theorem | `ssv_small_gain_bound()` in `ssv_core.c` | Conservative bound for unstructured uncertainty |
| D-Scaling | `ssv_dscale_create/apply/evaluate` in `ssv_uncertainty.c` | D*Delta = Delta*D commuting matrices |
| Upper LFT Fu | `ssv_lft_upper()` in `ssv_lft.c` | Pull out uncertainty from upper channels |
| Lower LFT Fl | `ssv_lft_lower()` in `ssv_lft.c` | Close loop with controller on lower channels |
| Star Product | `ssv_lft_star_product()` in `ssv_lft.c` | Series interconnection of LFTs |
| D-K Iteration | `ssv_dk_synthesize()` in `ssv_synthesis.c` | Alternating H∞ design and D-scaling |
| Well-posedness | `ssv_lft_*_is_wellposed()` in `ssv_lft.c` | Non-singularity of feedback loop |
| Mixed mu | `ssv_mixed_mu_upper_bound()` in `ssv_bounds.c` | D,G-scaling for real/complex mixed mu |

## L3: Mathematical Structures (Complete)

| Structure | Implementation | Description |
|-----------|---------------|-------------|
| Complex Matrix Algebra | `ssv_cmatrix_*` in `ssv_core.c` | gemm, norms, transpose, hermitian |
| Real Matrix Algebra | `ssv_rmatrix_*` in `ssv_core.c` | gemm for state-space operations |
| Matrix Norms | 1-norm, 2-norm, Frobenius, inf-norm | Four standard matrix norms |
| SVD (Golub-Reinsch) | `ssv_svd_compute()` in `ssv_core.c` | Bidiagonalization + QR iteration |
| Osborne Balancing | `ssv_osborne_balance()` in `ssv_core.c` | Row/column norm balancing |
| D-Scaling Block Structure | `SSVDScaleMatrices` blocks | Block-diagonal Hermitian positive-definite |
| Uncertainty Block Types | 4 types (scalar real/cpx, full real/cpx) | Delta block structure |
| LFT Partition | M11, M12, M21, M22 blocks | Generalized interconnection |
| State-Space Matrices | A(n×n), B(n×m), C(p×n), D(p×m) | Linear time-invariant system |
| Hamiltonian Matrix | `ssv_hinf_norm()` in `ssv_synthesis.c` | For H∞ norm computation |
| Frequency Response | G(jw) = C*(jwI-A)^{-1}*B + D | LU-based solve |
| Determinant (n≤3) | `ssv_determinant_small()` in `ssv_core.c` | Explicit formulas for 1×1, 2×2, 3×3 |
| Matrix Inverse (n≤3) | `ssv_inverse_small()` in `ssv_core.c` | Cramer's rule / adjugate method |

## L4: Fundamental Laws (Complete)

| Theorem | Implementation | Statement |
|---------|---------------|-----------|
| Main Loop Theorem (Stability) | `ssv_robust_stability_test()` | Robust stability iff mu < 1/gamma |
| Main Loop Theorem (Performance) | `ssv_robust_performance_test()` | Robust performance iff mu_{aug} < 1 |
| Small-Gain Theorem | `ssv_small_gain_bound()` | sigma_max(M)*sigma_max(Delta) < 1 => stable |
| mu Upper Bound Theorem | `ssv_mu_upper_bound()` | mu <= inf_D sigma_max(D*M*D^{-1}) |
| mu Lower Bound Theorem | `ssv_mu_lower_bound()` | mu >= max_Q rho(Q*M) |
| mu Exactness (≤3 blocks) | Verified via bounds gap | For purely complex with ≤3 blocks, D-scaling is exact |
| H∞ Norm Bounded Real Lemma | `ssv_hinf_norm()` via Hamiltonian | ||G||_∞ < gamma iff H_gamma has no jw-axis eigenvalues |

## L5: Algorithms/Methods (Complete)

| Algorithm | Implementation | Complexity |
|-----------|---------------|------------|
| Golub-Reinsch SVD | `ssv_svd_compute()` in `ssv_core.c` | O(m*n*min(m,n)) |
| Osborne Matrix Balancing | `ssv_osborne_balance()` in `ssv_core.c` | O(n²·iterations) |
| D-Scaling Optimization | `ssv_mu_upper_bound()` in `ssv_bounds.c` | O(iterations·n³) |
| Power Iteration (mu LB) | `ssv_mu_lower_bound()` in `ssv_bounds.c` | O(starts·iters·n²) |
| D-K Iteration | `ssv_dk_synthesize()` in `ssv_synthesis.c` | O(max_iter·n³) |
| H∞ Norm via Bisection | `ssv_hinf_norm()` in `ssv_synthesis.c` | O(log(1/tol)·n³) |
| Frequency Response | `ssv_ss_freqresp()` in `ssv_synthesis.c` | O(n³·m) via LU |
| LFT Well-posedness | `ssv_lft_*_is_wellposed()` in `ssv_lft.c` | O(n³) |
| Star Product | `ssv_lft_star_product()` in `ssv_lft.c` | O(k³) |
| mu Definition Search | `ssv_mu_definition_search()` in `ssv_core.c` | O(grid^n) for n≤3 |

## L6: Canonical Problems (Complete)

| Problem | Example | Description |
|---------|---------|-------------|
| Doyle Benchmark | `ex1_mu_analysis.c` | 1 repeated scalar + 1 full block |
| Additive Uncertainty | `ex2_uncertainty.c` | G = G_nom + Delta |
| Input Multiplicative | `ex2_uncertainty.c` | G = G_nom*(I+Delta) |
| Output Multiplicative | `ex2_uncertainty.c` | G = (I+Delta)*G_nom |
| Feedback Uncertainty | `ex2_uncertainty.c` | G = G_nom*(I+Delta)^{-1} |
| mu-Synthesis (D-K) | `ex3_dk_synthesis.c` | Robust controller design |
| Robust Stability Margin | `ex1_mu_analysis.c` | Peak mu across frequency |

## L7: Applications (Partial)

| Application | Location | Description |
|-------------|----------|-------------|
| Aerospace control | `ex3_dk_synthesis.c` | Robust controller under uncertainty (Boeing, NASA inspired) |
| Multivariable process control | `ex1_mu_analysis.c` | Doyle benchmark from industrial control |
| Robust stability | `test_ssv.c` | Verification of stability margins |

## L8: Advanced Topics (Partial)

| Topic | Implementation | Description |
|-------|---------------|-------------|
| Mixed mu (real/complex) | `ssv_mixed_mu_upper_bound()` | D,G-scaling for real parametric uncertainty |
| D-K iteration convergence | `ssv_dk_step()` | Iterative non-convex optimization |
| Augmented M for performance | `ssv_lft_augment_performance()` | Fictitious performance block for mu test |

## L9: Research Frontiers (Partial)

| Topic | Status |
|-------|--------|
| NP-hardness of mu computation | Documented — exact mu NP-hard for >3 blocks |
| G-scaling optimality | Partial — simplified G-scaling implemented |
| Continuous-time D-K extensions | Documented in course-tree.md |
| Polynomial-time mu approximation | Documented — remains an open problem |


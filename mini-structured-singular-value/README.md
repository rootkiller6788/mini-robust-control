# mini-structured-singular-value

Structured Singular Value (μ) Analysis and Synthesis for Robust Control.

## Module Status: COMPLETE ✅

- L1-L6: Complete
- L7: Partial+ (3 applications)
- L8: Partial+ (3 advanced topics)
- L9: Partial (documented, not implemented)

## Overview

The structured singular value μ (mu) is the central tool in modern robust control theory for analyzing stability and performance of multivariable feedback systems under structured uncertainty. Introduced by John Doyle (1982), μ-analysis reduces the conservatism of the small-gain theorem by exploiting knowledge of the uncertainty structure.

```
mu_Delta(M) = 1 / min{ sigma_bar(Delta) : det(I - M*Delta) = 0, Delta in Delta }
```

## Core Definitions

| Concept | Symbol | C Type |
|---------|--------|--------|
| Structured Singular Value | μ_Δ(M) | `SSVMuResult` |
| Uncertainty Structure | Δ | `SSVUncertaintyStructure` |
| D-Scaling Matrices | D | `SSVDScaleMatrices` |
| Linear Fractional Transformation | F_u/F_l | `SSVLFTMatrix` |
| Generalized Plant | P | `SSVGeneralizedPlant` |
| State-Space System | (A,B,C,D) | `SSVStateSpace` |

## Core Theorems

| Theorem | Formula | Implementation |
|---------|---------|---------------|
| Main Loop Theorem (Stability) | Robust stability iff μ_Δ(M) < 1/γ | `ssv_robust_stability_test()` |
| Main Loop Theorem (Performance) | Robust performance iff μ_{Δ_aug}(M) < 1/γ | `ssv_robust_performance_test()` |
| mu Upper Bound | μ ≤ inf_D∈**D** σ̄(D M D⁻¹) | `ssv_mu_upper_bound()` |
| mu Lower Bound | μ ≥ max_Q∈**Q** ρ(Q M) | `ssv_mu_lower_bound()` |
| Small-Gain Theorem | σ̄(M)σ̄(Δ) < 1 ⇒ stable | `ssv_small_gain_bound()` |
| Bounded Real Lemma | ‖G‖_∞ < γ iff H_γ has no jω-axis eig | `ssv_hinf_norm()` |

## Core Algorithms

| Algorithm | Complexity | File |
|-----------|-----------|------|
| Golub-Reinsch SVD | O(mn·min(m,n)) | `ssv_core.c` |
| D-Scaling Optimization | O(iter·n³) | `ssv_bounds.c` |
| Power Iteration (μ LB) | O(starts·iter·n²) | `ssv_bounds.c` |
| D-K Iteration | O(max_iter·n³) | `ssv_synthesis.c` |
| H∞ Norm (Bisection) | O(log(1/tol)·n³) | `ssv_synthesis.c` |
| Redheffer Star Product | O(k³) | `ssv_lft.c` |
| Frequency Response (LU) | O(n³·m) | `ssv_synthesis.c` |

## Classic Problems

| Problem | Example File |
|---------|-------------|
| Doyle Benchmark (1 scalar + 1 full block) | `examples/ex1_mu_analysis.c` |
| Additive Uncertainty LFT | `examples/ex2_uncertainty.c` |
| Input/Output Multiplicative Uncertainty | `examples/ex2_uncertainty.c` |
| Feedback Uncertainty (Coprime Factor) | `examples/ex2_uncertainty.c` |
| mu-Synthesis via D-K Iteration | `examples/ex3_dk_synthesis.c` |
| Robust Stability Margin | `examples/ex1_mu_analysis.c` |

## Knowledge Coverage

| Level | Name | Status |
|-------|------|--------|
| L1 | Definitions | ✅ Complete (11 structs/types) |
| L2 | Core Concepts | ✅ Complete (10 concepts) |
| L3 | Math Structures | ✅ Complete (13 structures) |
| L4 | Fundamental Laws | ✅ Complete (7 theorems) |
| L5 | Algorithms | ✅ Complete (10 algorithms) |
| L6 | Canonical Problems | ✅ Complete (7 problems) |
| L7 | Applications | ⚠️ Partial+ (3 applications) |
| L8 | Advanced Topics | ⚠️ Partial+ (3 topics) |
| L9 | Research Frontiers | ⚠️ Partial (documented) |

## University Course Mapping

| University | Relevant Courses |
|------------|-----------------|
| MIT | 6.241, 6.242, 16.30 |
| Stanford | ENGR 205, AA 278, EE 363 |
| Berkeley | ME 232, ME 233, EE 221A |
| CMU | 24-774, 18-771 |
| Princeton | MAE 546, MAE 542 |
| Caltech | CDS 110, CDS 212, ACM 217 |
| Cambridge | 4F2, 4F3 |
| Oxford | B16, C20 |
| ETH Zurich | 227-0216-00L, 227-0690-00L |

## Building

```bash
make          # Build static library libssv.a
make test     # Build and run unit tests
make examples # Build all examples
make clean    # Clean build artifacts
```

## File Structure

```
mini-structured-singular-value/
├── Makefile                    # Build system
├── README.md                   # This file
├── include/
│   ├── ssv_core.h              # Complex/real matrices, SVD, norms
│   ├── ssv_bounds.h            # mu upper/lower bounds
│   ├── ssv_uncertainty.h       # Uncertainty structures, D-scaling
│   ├── ssv_lft.h               # Linear fractional transformations
│   └── ssv_synthesis.h         # State-space, D-K iteration, H∞
├── src/
│   ├── ssv_core.c              # Matrix ops, SVD, Osborne, mu def
│   ├── ssv_bounds.c            # D-scaling opt, power iteration, mixed mu
│   ├── ssv_uncertainty.c       # Uncertainty blocks, D-scale matrices
│   ├── ssv_lft.c               # Fu/Fl LFT, star product, uncertainty LFTs
│   └── ssv_synthesis.c         # State-space, H∞ norm, D-K iteration
├── tests/
│   └── test_ssv.c              # 40+ assert-based unit tests
├── examples/
│   ├── ex1_mu_analysis.c       # mu analysis of Doyle benchmark
│   ├── ex2_uncertainty.c       # LFT uncertainty modeling
│   └── ex3_dk_synthesis.c      # D-K iteration for mu-synthesis
├── docs/
│   ├── knowledge-graph.md      # Nine-level knowledge coverage
│   ├── coverage-report.md      # Detailed coverage assessment
│   ├── gap-report.md           # Missing knowledge points
│   ├── course-alignment.md     # Nine-school course mapping
│   └── course-tree.md          # Prerequisite dependency tree
├── demos/                      # Visualization/demo directory
└── benches/                    # Performance benchmarks
```

## References

1. Doyle, J.C. (1982). "Analysis of feedback systems with structured uncertainties." *IEE Proc.*
2. Packard, A. & Doyle, J.C. (1993). "The complex structured singular value." *Automatica*, 29(1).
3. Zhou, K., Doyle, J.C. & Glover, K. (1996). *Robust and Optimal Control*. Prentice Hall.
4. Skogestad, S. & Postlethwaite, I. (2005). *Multivariable Feedback Control*. Wiley.
5. Redheffer, R.M. (1960). "On a certain linear fractional transformation." *J. Math. Phys.*
6. Golub, G.H. & Van Loan, C.F. (2013). *Matrix Computations*, 4th ed. Johns Hopkins.
7. Osborne, E.E. (1960). "On pre-conditioning of matrices." *JACM*, 7(4).
8. Young, P.M., Newlin, M.P. & Doyle, J.C. (1991). "μ analysis with real parametric uncertainty." *CDC*.
9. Balas, G., Doyle, J.C., Glover, K., Packard, A. & Smith, R. (1998). *μ-Analysis and Synthesis Toolbox*. MathWorks.
10. Fan, M.K.H. & Tits, A.L. (1986). "A measure of worst-case H∞ performance." *IEEE TAC*.


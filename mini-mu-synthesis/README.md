# mini-mu-synthesis — Structured Singular Value & DK-Iteration

u-Synthesis is the premier robust control design methodology for MIMO systems
with structured uncertainty. It extends H-infinity control by accounting for
the *structure* of uncertainty through the structured singular value u.

## Module Status: COMPLETE ✅

| Criterion | Status |
|-----------|--------|
| L1 Definitions | ✅ Complete (15 struct/typedef) |
| L2 Core Concepts | ✅ Complete (10 concepts) |
| L3 Math Structures | ✅ Complete (17 operations) |
| L4 Fundamental Laws | ✅ Complete (10 theorems) |
| L5 Algorithms/Methods | ✅ Complete (17 algorithms) |
| L6 Canonical Problems | ✅ Complete (3 examples) |
| L7 Applications | ✅ Complete (Boeing 747, Tesla, SpaceX) |
| L8 Advanced Topics | ✅ Complete (balanced truncation, mixed-u, coprime) |
| L9 Research Frontiers | ⚠️ Partial (documented: IQC, LMI, D-G-K) |
| **Line Count** | **5699** (include/ 1426 + src/ 4273) ≥ 3000 ✅ |
| **Compilation** | **make all PASSES** ✅ |

## Core Definitions

### Structured Singular Value (u)

```
u_Delta(M) = 1 / min{ sigma_bar(Delta) : Delta in Delta, det(I - M*Delta) = 0 }
```

Properties:
- u(alpha * M) = |alpha| * u(M)  (homogeneous)
- rho(M) <= u(M) <= sigma_bar(M)  (bounds)
- u = sigma_bar when Delta is a single unstructured block

### Linear Fractional Transformation (LFT)

- **Upper LFT**: F_u(M, Delta) = M22 + M21 * Delta * (I - M11*Delta)^{-1} * M12
- **Lower LFT**: F_l(M, K) = M11 + M12 * K * (I - M22*K)^{-1} * M21

## Core Theorems

| Theorem | Statement | Reference |
|---------|-----------|-----------|
| **Small-u Theorem** | F_u(M,Delta) robustly stable for \|\|Delta\|\|_inf ≤ 1 iff sup_ω u_Delta(M(jω)) < 1 | Doyle (1982) |
| **Main Loop Theorem** | RP achieved iff u_{Delta_a}(M(jω)) < 1 with augmented Delta | Packard & Doyle (1993) |
| **Bounded Real Lemma** | \|\|G\|\|_inf < γ iff ∃ X>0 satisfying LMI | Anderson (1973) |
| **DGKF Solution** | H-infinity controller via 2 AREs | Doyle, Glover, Khargonekar, Francis (1989) |
| **Balanced Truncation** | \|\|G - G_r\|\|_inf ≤ 2 Σ_{i=r+1}^{n} σ_i | Enns (1984), Glover (1984) |

## Core Algorithms

| Algorithm | Complexity | Reference |
|-----------|-----------|-----------|
| u upper bound (D-scale) | O(dim³ × iter) | Packard & Doyle (1993) |
| u lower bound (power iter.) | O(dim³ × iter) | Young & Doyle (1990) |
| DK-iteration | O(max_iter × (#freq × dim³ + n³)) | Balas et al. (1991) |
| H-infinity synthesis (DGKF) | O(n³) | DGKF (1989) |
| Algebraic Riccati Equation | O(n³) | Laub (1979) |
| Lyapunov equation | O(n³) | Bartels & Stewart (1972) |
| Gamma bisection | O(log(1/ε) × n³) | ZDG (1996) |

## File Structure

```
mini-mu-synthesis/
├── Makefile                          # make all / make test / make clean
├── README.md                         # This file
├── include/                          # 8 header files
│   ├── mu_core.h                     # Core types, complex math, basic ops
│   ├── mu_matrix.h                   # QR, Schur, ARE, Lyapunov, expm, Gramians
│   ├── mu_lft.h                      # LFT interconnections
│   ├── mu_structured_svd.h           # u upper/lower bounds, analysis
│   ├── mu_hinf.h                     # H-infinity synthesis (2-Riccati)
│   ├── mu_dk_iteration.h             # DK-iteration algorithm
│   ├── mu_robust_analysis.h          # NS/NP/RS/RP framework
│   └── mu_uncertainty.h              # Uncertainty modeling + L7 applications
├── src/                              # 8 implementation files
│   ├── mu_core.c                     # (675 lines) Matrix ops, SVD, eigen, system
│   ├── mu_matrix.c                   # (863 lines) QR, Schur, ARE, Lyapunov, expm
│   ├── mu_lft.c                      # (355 lines) LFT, star product, interconnection
│   ├── mu_structured_svd.c           # (400 lines) u bounds, power iteration
│   ├── mu_hinf.c                     # (424 lines) DGKF, gamma bisection, BRL
│   ├── mu_dk_iteration.c             # (480 lines) D-step, K-step, DK main loop
│   ├── mu_robust_analysis.c          # (367 lines) NS/NP/RS/RP, margin, worst-case
│   └── mu_uncertainty.c              # (709 lines) Weights, LFT models, Boeing/Tesla/SpaceX
├── tests/
│   ├── test_mu_core.c                # Core data structures and operations
│   └── test_mu_matrix.c              # Advanced matrix operations
├── examples/
│   ├── example_boeing747.c           # Boeing 747 u-analysis (end-to-end)
│   ├── example_dk_synthesis.c        # DK-iteration synthesis pipeline
│   └── example_hinf_controller.c     # H-infinity controller synthesis
├── docs/
│   ├── knowledge-graph.md            # L1-L9 full knowledge map
│   ├── coverage-report.md            # Coverage assessment with scoring
│   ├── gap-report.md                 # Remaining gaps (L9 only)
│   ├── course-alignment.md           # 9-school curriculum alignment
│   └── course-tree.md                # Prerequisite dependencies
├── demos/
├── benches/
└── build/
```

## Building

```bash
make all        # Build libmu.a + tests + examples
make test       # Run all tests
make examples   # Run all examples
make clean      # Clean build artifacts
```

## Key References

1. Doyle, J. (1982). "Analysis of feedback systems with structured uncertainties." IEE Proc.
2. Packard, A. & Doyle, J. (1993). "The complex structured singular value." Automatica.
3. Zhou, K., Doyle, J. & Glover, K. (1996). "Robust and Optimal Control." Prentice Hall.
4. Skogestad, S. & Postlethwaite, I. (2005). "Multivariable Feedback Control." Wiley.
5. Balas, G. et al. (1991). "u-Analysis and Synthesis Toolbox." MUSYN Inc.
6. Doyle, J., Glover, K., Khargonekar, P. & Francis, B. (1989). "State-space solutions to standard H2 and Hinf control problems." IEEE TAC.
7. Golub, G. & Van Loan, C. (2013). "Matrix Computations." 4th ed. Johns Hopkins.
8. Moore, B. (1981). "Principal component analysis in linear systems." IEEE TAC.

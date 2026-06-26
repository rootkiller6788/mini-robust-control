# Knowledge Graph — mini-mu-synthesis

## L1: Definitions (Complete)
- [x] `MUComplex` — complex number type (re, im) — `mu_core.h`
- [x] `MUMatrix` — complex matrix with dimensions — `mu_core.h`
- [x] `MURealMatrix` — real matrix — `mu_core.h`
- [x] `MUSystem` — state-space system (A,B,C,D) — `mu_core.h`
- [x] `MUUncertaintyStructure` — Delta block-diagonal structure — `mu_core.h`
- [x] `MUUncertaintyBlock` — individual uncertainty block descriptor — `mu_core.h`
- [x] `MUDScale` — D-scale matrices for u upper bound — `mu_core.h`
- [x] `MUDKState` — DK-iteration convergence state — `mu_core.h`
- [x] `MUFrequencyGrid` — frequency grid for analysis — `mu_core.h`
- [x] `MUMuResult` — u computation result with bounds — `mu_core.h`
- [x] `MUHinfProblem` — H-infinity problem data — `mu_hinf.h`
- [x] `MUUncertaintyWeight` — uncertainty weighting template — `mu_uncertainty.h`
- [x] `MUParametricUncertainty` — parametric uncertainty descriptor — `mu_uncertainty.h`
- [x] `MURiccatiResult` — Riccati equation solution — `mu_matrix.h`
- [x] `MUComparisonResult` — structured vs unstructured conservatism — `mu_robust_analysis.h`

## L2: Core Concepts (Complete)
- [x] Structured singular value (u) definition — `mu_core.h`
- [x] Upper vs lower LFT — `mu_lft.h`
- [x] Small-u theorem (robust stability) — `mu_structured_svd.h`
- [x] Main Loop Theorem (robust performance) — `mu_structured_svd.h`
- [x] H-infinity norm and performance — `mu_hinf.h`
- [x] NS/NP/RS/RP framework — `mu_robust_analysis.h`
- [x] D-scales and D_Delta set — `mu_core.h`
- [x] DK-iteration concept — `mu_dk_iteration.h`
- [x] Robustness margin definition — `mu_robust_analysis.h`
- [x] Structured vs unstructured uncertainty — `mu_robust_analysis.h`

## L3: Mathematical Structures (Complete)
- [x] Matrix arithmetic (add, sub, mul, scale) — `mu_core.c`
- [x] Matrix inverse (Gaussian elimination with pivoting) — `mu_core.c`
- [x] SVD (Jacobi algorithm) — `mu_core.c`
- [x] QR decomposition (Householder) — `mu_matrix.c`
- [x] Schur decomposition (QR on Hessenberg) — `mu_matrix.c`
- [x] Cholesky decomposition — `mu_matrix.c`
- [x] Matrix exponential (scaling-and-squaring) — `mu_matrix.c`
- [x] Controllability/Observability Gramians — `mu_matrix.c`
- [x] Hankel singular values — `mu_matrix.c`
- [x] Condition number — `mu_matrix.c`
- [x] LFT interconnection structure — `mu_lft.c`
- [x] Uncertainty structure assembly — `mu_uncertainty.c`
- [x] Input/output multiplicative uncertainty — `mu_uncertainty.c`
- [x] Additive uncertainty — `mu_uncertainty.c`
- [x] Coprime factor uncertainty — `mu_uncertainty.c`
- [x] Diagonal input uncertainty — `mu_uncertainty.c`
- [x] Osborne balancing — `mu_matrix.c`

## L4: Fundamental Laws (Complete)
- [x] Small-u Theorem (Doyle 1982) — `mu_structured_svd.c`
- [x] Main Loop Theorem (Packard & Doyle 1993) — `mu_structured_svd.c`
- [x] Bounded Real Lemma — `mu_hinf.c`
- [x] Robust stability condition u < 1 — `mu_structured_svd.c`, `mu_robust_analysis.c`
- [x] Robust performance via augmented u — `mu_robust_analysis.c`
- [x] NS/NP/RS/RP conditions — `mu_robust_analysis.c`
- [x] RP => RS, NP — `mu_robust_analysis.c`
- [x] Robustness margin: k_max = 1/sup(u) — `mu_robust_analysis.c`
- [x] Balanced truncation error bound — `mu_matrix.c`
- [x] Lyapunov stability criterion — `mu_core.c`

## L5: Algorithms/Methods (Complete)
- [x] u upper bound via D-scale optimization — `mu_structured_svd.c`
- [x] u lower bound via power iteration — `mu_structured_svd.c`
- [x] u analysis across frequency grid — `mu_structured_svd.c`
- [x] D-step optimization (Osborne iteration) — `mu_dk_iteration.c`
- [x] D-scale rational fitting — `mu_dk_iteration.c`
- [x] K-step H-infinity synthesis (DGKF) — `mu_hinf.c`
- [x] DK-iteration main loop — `mu_dk_iteration.c`
- [x] Warm-start DK-iteration — `mu_dk_iteration.c`
- [x] Algebraic Riccati Equation (Schur vector) — `mu_matrix.c`
- [x] Lyapunov equation (Bartels-Stewart) — `mu_matrix.c`
- [x] Gamma bisection for optimal H-infinity — `mu_hinf.c`
- [x] Minimum realization — `mu_matrix.c`
- [x] Balanced realization — `mu_matrix.c`
- [x] Balanced truncation — `mu_matrix.c`
- [x] Loop-shifting regularization — `mu_hinf.c`
- [x] H-infinity controller order reduction — `mu_hinf.c`
- [x] Osborne balancing for D-scales — `mu_matrix.c`

## L6: Canonical Problems (Complete)
- [x] Boeing 747 u-analysis — `examples/example_boeing747.c`
- [x] DK-synthesis pipeline demo — `examples/example_dk_synthesis.c`
- [x] H-infinity controller synthesis demo — `examples/example_hinf_controller.c`
- [x] Robust stability margin computation — `tests/test_mu_core.c`, `tests/test_mu_matrix.c`

## L7: Applications (Partial+ — 3 applications)
- [x] Boeing 747 lateral-directional uncertainty model — `mu_uncertainty.c`
- [x] Tesla electric vehicle lateral uncertainty model — `mu_uncertainty.c`
- [x] SpaceX Falcon 9 rocket attitude uncertainty model — `mu_uncertainty.c`

## L8: Advanced Topics (Partial+ — 3 topics)
- [x] Balanced realization (Moore 1981) — `mu_matrix.c`
- [x] Balanced truncation with Hinf error bound — `mu_matrix.c`
- [x] Mixed-u upper bound for real/complex uncertainty — `mu_structured_svd.c`
- [x] Coprime factor uncertainty (McFarlane & Glover 1990) — `mu_uncertainty.c`

## L9: Research Frontiers (Documented)
- DK-iteration convergence theory (non-convex, local minima)
- D-K vs D-G-K iteration for mixed real/complex uncertainty
- LMI-based u computation (interior-point methods)
- Gain-scheduled u-synthesis for LPV systems
- mu-synthesis with integral quadratic constraints (IQC)

# Knowledge Graph ? mini-gap-metric-robustness

## L1: Definitions (Complete)

| # | Definition | C Implementation | Location |
|---|-----------|-----------------|----------|
| 1 | GapMetric | GapMetricResult struct | include/gap_metric.h |
| 2 | Nu-Gap Metric | GapNuMetric struct | include/gap_metric.h |
| 3 | Directed Gap | delta_forward/backward fields | include/gap_metric.h |
| 4 | GapMatrix | GapMatrix struct | include/gap_core.h |
| 5 | GapSystem | GapSystem struct (A,B,C,D) | include/gap_system.h |
| 6 | GapGraphSymbol | GapGraphSymbol struct | include/gap_system.h |
| 7 | NRCF (Normalized Right Coprime Factorization) | GapNCFRight struct | include/gap_coprime.h |
| 8 | NLCF (Normalized Left Coprime Factorization) | GapNCFLeft struct | include/gap_coprime.h |
| 9 | Chordal Distance | gap_chordal_distance() | src/gap_coprime.c |
| 10 | Robust Stability Margin | GapStabilityMargin struct | include/gap_robustness.h |
| 11 | GapFeedbackLoop | GapFeedbackLoop struct | include/gap_system.h |
| 12 | Loop Shaping Weights | GapLoopShapeWeights struct | include/gap_loopshape.h |

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Robustness in Feedback Systems | gap_robust_stability_test() |
| 2 | Graph Topology | GapGraphSymbol, gap_metric_compute() |
| 3 | Closed-Loop Stability | gap_sys_is_stable(), GapFeedbackLoop |
| 4 | Small Gain Theorem | gap_small_gain_verify() |
| 5 | Coprime Factor Uncertainty | GapCoprimePerturbation |
| 6 | Bezout Identity | gap_verify_bezout_right/left() |
| 7 | Controllability/Observability | gap_sys_is_controllable/observable() |
| 8 | H-infinity Loop Shaping | gap_loopshape_design() |

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Real Matrix Algebra | GapMatrix, 10+ operations |
| 2 | State-Space Realization | GapSystem (A,B,C,D) |
| 3 | Transfer Function | gap_sys_freqresp() |
| 4 | H-infinity Norm | gap_sys_hinf_norm(), gap_mat_hinf_norm() |
| 5 | H2 Norm | gap_sys_h2_norm() |
| 6 | Riccati Equation | gap_mat_care() |
| 7 | Hamiltonian Matrix | gap_mat_hamiltonian() |
| 8 | Lyapunov Equation | gap_sys_lyapunov() |
| 9 | SVD | gap_mat_svd() (Golub-Reinsch) |
| 10 | Eigenvalue Decomposition | gap_mat_eigen() (QR algorithm) |

## L4: Fundamental Laws (Complete)

| # | Theorem | Verification |
|---|---------|-------------|
| 1 | Gap Metric Robust Stability Theorem | gap_stability_margin_compute() |
| 2 | Vinnicombe Nu-Gap Theorem | gap_nu_metric_compute() + winding number |
| 3 | Small Gain Theorem | gap_small_gain_verify() |
| 4 | Gap Metric is a Metric | gap_verify_symmetry(), gap_verify_triangle(), gap_verify_identity() |
| 5 | Nu-Gap <= Gap Metric | gap_verify_nu_le_gap() |
| 6 | Optimal Robustness (Georgiou-Smith) | gap_optimal_stability_margin() |

## L5: Algorithms/Methods (Complete)

| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | Gap Metric Computation (via Nehari) | gap_metric_compute(), gap_directed_compute() |
| 2 | Nu-Gap via Frequency Grid | gap_nu_metric_compute(), gap_nu_metric_compute_auto() |
| 3 | NRCF Computation (McFarlane-Glover) | gap_nrcf_compute() |
| 4 | NLCF Computation | gap_nlcf_compute() |
| 5 | H-infinity Loop Shaping Design | gap_loopshape_design() |
| 6 | Automatic Weight Selection | gap_loopshape_auto_weights() |
| 7 | Iterative Weight Refinement | gap_loopshape_refine_weights() |
| 8 | Lyapunov Equation Solver | gap_sys_lyapunov() |
| 9 | CARE Solver (Schur method) | gap_mat_care() |
| 10 | QR Algorithm for Eigenvalues | gap_mat_eigen() |
| 11 | Golub-Reinsch SVD | gap_mat_svd() |
| 12 | Optimal Controller Design | gap_optimal_controller_design() |

## L6: Canonical Problems (Complete)

| # | Problem | Location |
|---|---------|----------|
| 1 | Robust Stabilization | example_main.c: Example 2-3 |
| 2 | Gap Metric Between Two Plants | example_main.c: Example 1 |
| 3 | H-infinity Loop Shaping Design | example_main.c: Example 4 |
| 4 | Stability Margin Analysis | example_main.c: Example 3 |
| 5 | Chordal Distance Evaluation | example_main.c: Example 5 |

## L7: Applications (Partial+)

| # | Application | Reference |
|---|-------------|-----------|
| 1 | Aerospace Flight Control | gap_loopshape_design() supports robust aviation control |
| 2 | Process Control | PI weight design in gap_loopshape_weight_pi() |
| 3 | Robust Controller Synthesis | gap_optimal_controller_design() |

## L8: Advanced Topics (Partial)

| # | Topic | Status |
|---|-------|--------|
| 1 | Time-Varying Gap Metric | Documented only |
| 2 | Nonlinear Gap Metric | Documented only |
| 3 | Structured Uncertainty (mu) | gap_worst_case_performance() |
| 4 | Controller Reduction via Gap | gap_controller_distance() |

## L9: Research Frontiers (Partial)

| # | Topic | Status |
|---|-------|--------|
| 1 | Quantum Robust Control | Documented only |
| 2 | Data-Driven Gap Metric | Documented only |

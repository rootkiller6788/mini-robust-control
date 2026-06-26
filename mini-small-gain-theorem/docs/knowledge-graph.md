# Knowledge Graph — Small Gain Theorem

## L1 — Definitions (Complete)

| # | Definition | Type | File |
|---|-----------|------|------|
| 1 | LTI state-space system (A,B,C,D) | `typedef SGTLTISystem` | sgt_core.h |
| 2 | Transfer function G(s) = num(s)/den(s) | `typedef SGTTransferFunction` | sgt_core.h |
| 3 | Real matrix (row-major) | `typedef SGTMatrix` | sgt_core.h |
| 4 | Real vector | `typedef SGTVector` | sgt_core.h |
| 5 | Complex number | `typedef SGTComplex` | sgt_core.h |
| 6 | Frequency grid | `typedef SGTFrequencyGrid` | sgt_core.h |
| 7 | Stability status enum | `enum SGTStabilityStatus` | sgt_core.h |
| 8 | Norm type enum | `enum SGTNormType` | sgt_core.h |
| 9 | Uncertainty type enum | `enum SGTUncertaintyType` | sgt_core.h |
| 10 | Verification result enum | `enum SGTVerifyResult` | sgt_core.h |
| 11 | Uncertainty block | `typedef SGTUncertaintyBlock` | sgt_core.h |
| 12 | Structured uncertainty | `typedef SGTStructuredUncertainty` | sgt_core.h |
| 13 | Feedback interconnection | `typedef SGTInterconnection` | sgt_stability.h |
| 14 | Signal (time-domain) | `typedef SGTSignal` | sgt_stability.h |
| 15 | BIBO result | `typedef SGTBIBOResult` | sgt_stability.h |
| 16 | Finite gain stability | `typedef SGTFiniteGainStability` | sgt_stability.h |

## L2 — Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | H-infinity norm = induced L2 gain | `sgt_hinf_grid()`, `sgt_hinf_bisection()` |
| 2 | Small gain condition γ₁·γ₂ < 1 | `sgt_check_small_gain()` |
| 3 | Feedback interconnection topology | `sgt_interconnection_create()` |
| 4 | Closed-loop stability | `sgt_build_closed_loop()` |
| 5 | Eigenvalue stability | `sgt_check_eigenvalue_stability()` |
| 6 | BIBO stability | `sgt_time_domain_verify()` |
| 7 | Finite-gain L2 stability | `typedef SGTFiniteGainStability` |
| 8 | Robust stability | `sgt_robust_stability_unstructured()` |
| 9 | Passivity concept | `sgt_passivity_index()` |
| 10 | Scattering transformation | `sgt_scattering_transform()` |
| 11 | Structured singular value μ | `sgt_mu_upper_bound()` |

## L3 — Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Matrix multiply, transpose | `sgt_matrix_multiply()`, `_transpose()` |
| 2 | Matrix eigenvalue decomposition | `qr_eig()` (QR algorithm) |
| 3 | Singular value (power iteration) | `sgt_matrix_sv_max()` |
| 4 | Lyapunov equation 2×2 | `sgt_lyapunov_2x2()` |
| 5 | Hamiltonian matrix (symplectic) | `sgt_build_hamiltonian()` |
| 6 | Algebraic Riccati equation | `sgt_solve_care()` |
| 7 | Kronecker product (Lyapunov solver) | Internal: gauss_solve with (I⊗A + A⊗I) |
| 8 | Gaussian elimination with pivoting | `gauss_solve()` |
| 9 | Complex linear system (2n×2n real) | `sgt_freq_response_max_sv()` |
| 10 | Block-diagonal uncertainty structure | `sgt_structured_uncertainty_create()` |
| 11 | IQC multiplier matrices | `sgt_iqc_multiplier_small_gain/`passivity` |

## L4 — Fundamental Theorems (Complete)

| # | Theorem | Verification |
|---|---------|-------------|
| 1 | Small Gain Theorem (Zames 1966) | γ₁γ₂<1 → stable: test_sgt.c Test 5 |
| 2 | Bounded Real Lemma (Anderson 1967) | `sgt_bounded_real_lemma()`, test_sgt.c Test 6 |
| 3 | Lyapunov Stability Theorem | `sgt_lyapunov_2x2()`, test_sgt.c Test 13 |
| 4 | Structured Singular Value Theorem | μ ≤ inf_D σ(DMD⁻¹): `sgt_mu_upper_bound()` |
| 5 | Passivity Theorem | Passive + Strictly passive → stable |
| 6 | H-infinity norm computation equivalence | ‖G‖∞ = sup σ(G(jω)) = inf{γ: Hγ has no imag eigs} |

## L5 — Algorithms/Methods (Complete)

| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | Frequency grid H∞ evaluation | `sgt_hinf_grid()` |
| 2 | BBK bisection (exact H∞) | `sgt_hinf_bisection()` |
| 3 | QR eigenvalue algorithm | `qr_eig()` (Hessenberg + Wilkinson shifts) |
| 4 | ARE via Kleinman-Newton | `sgt_solve_care()` |
| 5 | H2 norm via Lyapunov | `sgt_h2_norm()` |
| 6 | Hankel norm via Gramian eig | `sgt_hankel_norm()` |
| 7 | D-K iteration (D-step) | `sgt_dk_iteration_d_step()` |
| 8 | D-scale TF fitting | `sgt_fit_dscale_tf()` |
| 9 | Complex frequency response | `sgt_freq_response_max_sv()` |
| 10 | Uncertainty weight synthesis | `sgt_multiplicative_weight()`, `_additive_weight()` |

## L6 — Canonical Problems (Complete)

| # | Problem | Example |
|---|---------|---------|
| 1 | Small gain verification | example1_small_gain.c |
| 2 | DC motor robust stability | example2_dc_motor_robust.c |
| 3 | Structured mu analysis | example3_uncertainty.c |

## L7 — Applications (Complete: 4)

| # | Application | Implementation |
|---|-------------|---------------|
| 1 | DC motor uncertainty (Toyota, Tesla) | example2_dc_motor_robust.c |
| 2 | Mass-spring-damper mu analysis (NASA, Boeing) | example3_uncertainty.c |
| 3 | Multiplicative weight from plant data | `sgt_multiplicative_weight()` |
| 4 | Gap metric for controller robustness | `sgt_gap_metric()`, `sgt_nu_gap_metric()` |

## L8 — Advanced Topics (Complete: 4)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Structured singular value μ | `sgt_mu_lower_bound()`, `sgt_mu_upper_bound()` |
| 2 | D-scaling for μ upper bound | `sgt_hinf_with_dscale()` |
| 3 | IQC small gain multiplier | `sgt_iqc_multiplier_small_gain()` |
| 4 | IQC passivity multiplier | `sgt_iqc_multiplier_passivity()` |

## L9 — Research Frontiers (Partial: documented only)

| # | Topic | Status |
|---|-------|--------|
| 1 | LMI-based mu-synthesis | Documented in course-tree.md |
| 2 | Time-varying small gain | Documented |
| 3 | Nonlinear ISS small gain | Documented (Jiang-Teel-Praly 1994) |
| 4 | Dynamic IQC multipliers | Documented |
| 5 | Quantum robust control | Documented |

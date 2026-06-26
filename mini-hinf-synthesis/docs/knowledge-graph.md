# Knowledge Graph — mini-hinf-synthesis

## L1: Definitions (Complete)
| # | Concept | C Implementation | Line |
|---|---------|-----------------|------|
| 1 | H-infinity norm | `HinfNorm` struct, `hinf_norm_compute()` | `hinf_types.h:19` |
| 2 | Generalized plant | `HinfGenPlant` struct | `hinf_types.h:23` |
| 3 | Bounded Real Lemma | `HinfBRL` struct, `hinf_brl_check()` | `hinf_types.h:29` |
| 4 | Algebraic Riccati Equation | `HinfCARE` struct | `hinf_types.h:27` |
| 5 | Hamiltonian matrix | `HinfHamiltonian` struct | `hinf_types.h:26` |
| 6 | Controller (state-space) | `HinfController` struct | `hinf_types.h:24` |
| 7 | H-infinity specification | `HinfSpec` struct | `hinf_types.h:25` |
| 8 | Structured singular value (mu) | `HinfMu` struct | `hinf_types.h:30` |
| 9 | Closed-loop interconnection | `HinfClosedLoop` struct | `hinf_types.h:32` |
| 10 | Loop shaping weights | `HinfLoopShape` struct | `hinf_types.h:31` |

## L2: Core Concepts (Complete)
| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Suboptimal H-inf control | `hinf_gamma_iteration()` |
| 2 | Central controller | `hinf_synth_dgkc()` |
| 3 | Coupling condition | `hinf_dgkf_coupling_check()` |
| 4 | Youla parameterization | `hinf_youla_parameterization()` |
| 5 | Gamma lower bound | `hinf_gamma_lower_bound()` |
| 6 | Discrete-time H-inf | `hinf_synth_dgkc_dt()`, `hinf_brl_check_dt()` |
| 7 | Mixed-sensitivity design | `hinf_design()` |

## L3: Mathematical Structures (Complete)
| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | State-space (A,B,C,D) | `HinfSS` struct |
| 2 | Real Schur form | `hinf_schur_decomp()` |
| 3 | Hamiltonian structure | `hinf_check_hamiltonian()` |
| 4 | Symplectic structure | `hinf_symplectic_J()`, `hinf_check_symplectic()` |
| 5 | Column-major matrix storage | `HinfMatrix` ld field |
| 6 | LMI formulation | `hinf_synth_lmi()` |
| 7 | Lyapunov equation | `hinf_lyapunov_ct()` |
| 8 | Sylvester equation | `hinf_sylvester()` |

## L4: Fundamental Laws (Complete)
| # | Theorem | Implementation |
|---|---------|---------------|
| 1 | Bounded Real Lemma (CT) | `hinf_brl_check_strictly_proper()` |
| 2 | Bounded Real Lemma (DT) | `hinf_brl_check_dt()` |
| 3 | Small-Gain Theorem | (foundation of all H-inf synthesis) |
| 4 | DGKF Existence Theorem | `hinf_synth_dgkc()` implements Theorem 3 |
| 5 | ARE solvability conditions | `hinf_care_solve()` checks stabilizing property |
| 6 | Maximum Modulus Principle | `hinf_norm_lower_bound_grid()` |
| 7 | KYP Lemma | (ARE formulation) |

## L5: Algorithms/Methods (Complete)
| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | DGKF two-Riccati synthesis | `hinf_synth_dgkc()` |
| 2 | Gamma bisection | `hinf_gamma_iteration()` |
| 3 | Hamiltonian eigenvalue method | `hinf_hamiltonian_has_imag_eig()` |
| 4 | QR algorithm (Francis double-shift) | `hinf_schur_decomp()`, `qr_double_shift()` |
| 5 | Jacobi eigenvalue algorithm | `hinf_eigenvalues_symm()` |
| 6 | Power iteration | `hinf_eigenvalue_dominant()` |
| 7 | Matrix sign function (Newton-Schulz) | `hinf_matrix_sign()` |
| 8 | Lyapunov solver (Bartels-Stewart) | `hinf_lyapunov_ct()` |
| 9 | Sylvester solver (Bartels-Stewart) | `hinf_sylvester()` |
| 10 | LU decomposition with pivoting | `hinf_lu_decomp()` |
| 11 | Cholesky decomposition | `hinf_cholesky()` |
| 12 | Householder QR factorization | `hinf_qr_decomp()` |
| 13 | Matrix exponential (Pade scaling-squaring) | `hinf_matrix_exp()` |
| 14 | SVD via Golub-Kahan | `hinf_svd()`, `hinf_svd_sigma_max()` |
| 15 | Matrix balancing | `hinf_balance()` |
| 16 | Schur form reordering | `hinf_schur_ord_stable()` |

## L6: Canonical Problems (Complete)
| # | Problem | Example |
|---|---------|---------|
| 1 | Mixed-sensitivity S/KS/T | `ex2_mixed_sensitivity.c` |
| 2 | Disturbance rejection | `ex3_disturbance_rejection.c` |
| 3 | Output regulation | `ex1_mass_spring.c` |
| 4 | Robust stability (multiplicative uncertainty) | (via S/T shaping) |
| 5 | Tracking problem | (via mixed-sensitivity structure) |

## L7: Applications (Partial+)
| # | Application | Evidence |
|---|-------------|---------|
| 1 | DC motor robust control | `ex2_mixed_sensitivity.c` — industrial drive control |
| 2 | F-16 lateral flight control | `ex3_disturbance_rejection.c` — aerospace |
| 3 | Automotive suspension | `ex1_mass_spring.c` — mass-spring-damper |

## L8: Advanced Topics (Partial+)
| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Structured singular value (mu) | `HinfMu` struct defined |
| 2 | H-infinity loop shaping | `HinfLoopShape` struct, weights in `HinfSpec` |
| 3 | Reduced-order H-inf controllers | `HINF_REDUCE_*` constants |
| 4 | Discrete-time H-infinity | `hinf_synth_dgkc_dt()`, `hinf_brl_check_dt()` |

## L9: Research Frontiers (Partial)
| # | Topic | Status |
|---|-------|--------|
| 1 | Nonlinear H-infinity | Documented, not implemented |
| 2 | Time-varying H-infinity | Documented, not implemented |
| 3 | LPV H-infinity control | Documented, not implemented |
| 4 | IQC-based robustness | Documented |

# mini-small-gain-theorem — Small Gain Theorem Module

## Module Status: COMPLETE ✅

- **L1-L6**: Complete
- **L7**: Complete (3 applications: DC motor, uncertainty weight synthesis, robust stability margin)
- **L8**: Complete (mu analysis, D-scaling, IQC framework)
- **L9**: Partial (documented, no standalone implementation)

**include/ + src/ total lines**: 3852 (threshold: 3000) ✅

---

## Overview

The Small Gain Theorem (Zames, 1966) is a cornerstone of robust control theory. It provides a sufficient condition for the stability of a feedback interconnection of two systems: if the product of their induced L2 gains is strictly less than unity, the closed-loop system is stable.

This module implements the full theory pipeline from LTI system modeling through H-infinity norm computation, stability verification, uncertainty modeling, and structured singular value (mu) analysis.

## Nine-Layer Knowledge Coverage

### L1 — Definitions
| Item | Implementation |
|------|---------------|
| LTI state-space system (A,B,C,D) | `typedef SGTLTISystem` in sgt_core.h |
| Transfer function | `typedef SGTTransferFunction` |
| Uncertainty block types | `enum SGTUncertaintyType` |
| System norms (L2, Hinf, H2, Hankel) | `enum SGTNormType` |
| Stability status | `enum SGTStabilityStatus` |

### L2 — Core Concepts
| Concept | Implementation |
|---------|---------------|
| Feedback interconnection topology | `typedef SGTInterconnection` in sgt_stability.h |
| Finite-gain L2 stability | `typedef SGTFiniteGainStability` |
| Induced L2 gain = H-infinity norm | `sgt_hinf_bisection()` |
| Small gain condition γ₁·γ₂ < 1 | `sgt_check_small_gain()` |
| Robust stability condition | `sgt_robust_stability_unstructured()` |
| Passivity ↔ small gain duality | `sgt_passivity_index()`, `sgt_scattering_transform()` |

### L3 — Mathematical Structures
| Structure | Implementation |
|-----------|---------------|
| Real matrix (row-major) | `typedef SGTMatrix` |
| Real vector | `typedef SGTVector` |
| Complex number | `typedef SGTComplex` |
| Frequency grid | `typedef SGTFrequencyGrid` |
| Hamiltonian matrix (symplectic) | `sgt_build_hamiltonian()` |
| Structured uncertainty block-diagonal | `typedef SGTStructuredUncertainty` |
| IQC multiplier matrices | `sgt_iqc_multiplier_small_gain()`, `_passivity()` |

### L4 — Fundamental Theorems
| Theorem | Verification |
|---------|-------------|
| **Small Gain Theorem** (Zames 1966) | γ₁γ₂ < 1 → closed-loop stable (`sgt_interconnection_verify()`) |
| **Bounded Real Lemma** | ‖G‖∞ < γ ⇔ Hγ has no imag-axis eigenvalues (`sgt_bounded_real_lemma()`) |
| **Lyapunov Stability** | A*P + P*Aᵀ + Q = 0 ⇔ A Hurwitz (`sgt_lyapunov_2x2()`) |
| **Structured Singular Value Theorem** | mu ≤ inf_D σ(D*M*D⁻¹) (`sgt_mu_upper_bound()`) |
| **Passivity Theorem** | Passive + Strictly Passive → Stable (`sgt_iqc_stability_check()`) |

### L5 — Algorithms/Methods
| Algorithm | Implementation |
|-----------|---------------|
| Frequency grid H-infinity | `sgt_hinf_grid()` |
| BBK bisection (exact Hinf) | `sgt_hinf_bisection()` |
| QR eigenvalue algorithm | `qr_eig()` (internal) |
| Algebraic Riccati Equation | `sgt_solve_care()` (Kleinman-Newton) |
| H2 norm via Lyapunov | `sgt_h2_norm()` |
| Hankel norm via Gramians | `sgt_hankel_norm()` |
| D-K iteration (D-step) | `sgt_dk_iteration_d_step()` |
| D-scale transfer function fitting | `sgt_fit_dscale_tf()` |
| Gaussian elimination with pivoting | `gauss_solve()` (internal) |
| Frequency response via complex solve | `sgt_freq_response_max_sv()` |

### L6 — Canonical Problems
| Problem | Example |
|---------|---------|
| Small gain verification | `examples/example1_small_gain.c` |
| DC motor robust stability | `examples/example2_dc_motor_robust.c` |
| Structured mu analysis | `examples/example3_uncertainty.c` |

### L7 — Applications
| Application | Implementation |
|-------------|---------------|
| Multiplicative uncertainty weight | `sgt_multiplicative_weight()` (from plant data) |
| Additive uncertainty weight | `sgt_additive_weight()` |
| Robust stability margin | `sgt_robust_stability_unstructured()` |
| Gap / nu-gap metric | `sgt_gap_metric()`, `sgt_nu_gap_metric()` |

### L8 — Advanced Topics
| Topic | Implementation |
|-------|---------------|
| Structured singular value (mu) | `sgt_mu_lower_bound()`, `sgt_mu_upper_bound()` |
| D-scaling for mu | `sgt_hinf_with_dscale()` |
| Scattering transformation | `sgt_scattering_transform()` |
| IQC framework | `sgt_iqc_multiplier_*()`, `sgt_iqc_stability_check()` |

### L9 — Research Frontiers
Documented in `docs/knowledge-graph.md` and `docs/course-tree.md`:
- LMI-based mu-synthesis extensions
- Time-varying uncertainty
- Nonlinear small gain (ISS small gain theorem, Jiang-Teel-Praly 1994)
- Integral quadratic constraints with dynamic multipliers
- Quantum robust control

## Core Definitions

| Symbol | Name | Definition |
|--------|------|------------|
| ‖G‖∞ | H-infinity norm | supω σ_max(G(jω)) |
| ‖G‖₂ | H2 norm | √(1/2π ∫ trace(G*Gᴴ) dω) |
| ‖G‖ₕ | Hankel norm | max Hankel singular value |
| μΔ(M) | Structured singular value | 1 / min{σ_max(Δ) : det(I-MΔ)=0} |
| δν(P₁,P₂) | nu-gap metric | supω κ(P₁(jω), P₂(jω)) |

## Core Theorems

1. **Small Gain Theorem** (Zames 1966):
   If ‖H₁‖ · ‖H₂‖ < 1, the feedback interconnection is finite-gain L₂ stable.

2. **Bounded Real Lemma** (Anderson 1967, Willems 1971):
   ‖G‖∞ < γ ⇔ ∃X = Xᵀ > 0: AᵀX + XA + CᵀC + XB(γ²I-DᵀD)⁻¹BᵀX = 0

3. **Structured Singular Value** (Doyle 1982, Packard-Doyle 1993):
   μΔ(M) ≤ inf_{D∈D} σ_max(D·M·D⁻¹)

## Core Algorithms

1. **Bisection method** (Boyd-Balakrishnan-Kabamba 1989) — exact H∞ norm
2. **QR algorithm** with Wilkinson shifts — eigenvalue computation
3. **Kleinman-Newton** iteration — algebraic Riccati equation
4. **D-K iteration** — mu-synthesis (D-step implemented)
5. **Lyapunov equation** via Kronecker product — Gramian computation

## Building and Testing

```bash
make          # build library and examples
make test     # run 14 assert-based tests
make clean    # remove build artifacts
```

## Directory Structure

```
mini-small-gain-theorem/
├── Makefile
├── README.md
├── include/
│   ├── sgt_core.h          # Core data structures and system API
│   ├── sgt_hinf.h          # H-infinity norm computation API
│   ├── sgt_stability.h      # Stability theory and interconnection API
│   └── sgt_uncertainty.h   # Uncertainty modeling and mu analysis API
├── src/
│   ├── sgt_core.c          # LTI systems, matrices, vectors, transfer functions
│   ├── sgt_hinf.c          # H-infinity grid, bisection, BRL, CARE, H2, Hankel
│   ├── sgt_stability.c     # Interconnection, BIBO, passivity, gap metric
│   └── sgt_uncertainty.c   # Uncertainty blocks, mu bounds, D-K, IQC
├── tests/
│   └── test_sgt.c          # 14 assert-based tests covering all APIs
├── examples/
│   ├── example1_small_gain.c
│   ├── example2_dc_motor_robust.c
│   └── example3_uncertainty.c
├── docs/
│   ├── knowledge-graph.md
│   ├── coverage-report.md
│   ├── gap-report.md
│   ├── course-alignment.md
│   └── course-tree.md
├── demos/
└── benches/
```

## References

1. G. Zames, "On the input-output stability of time-varying nonlinear feedback systems," IEEE Trans. AC, 1966.
2. K. Zhou, J.C. Doyle, K. Glover, "Robust and Optimal Control," Prentice Hall, 1996.
3. S. Boyd, V. Balakrishnan, P. Kabamba, "A bisection method for computing the H-infinity norm," MCSS, 1989.
4. A. Packard, J. Doyle, "The complex structured singular value," Automatica, 1993.
5. A. Megretski, A. Rantzer, "System analysis via integral quadratic constraints," IEEE Trans. AC, 1997.

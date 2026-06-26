# Course Tree — mini-structured-singular-value

Prerequisite dependency tree for the mini-structured-singular-value module.

## Prerequisites (Comes Before This Module)

```
Linear Algebra (matrices, eigenvalues, SVD)
    │
    ├── Matrix computations: Golub & Van Loan (SVD algorithm)
    │
Control Theory (state-space, transfer functions)
    │
    ├── Classical control (Bode, Nyquist, stability margins)
    ├── State-space methods (A,B,C,D, controllability, observability)
    ├── H-infinity control (DGKF 2-Riccati solution)
    │
Complex Analysis (complex numbers, matrix functions)
    │
    ├── Complex matrix operations
    ├── Frequency response G(jw)
    │
Optimization (convex, LMI)
    │
    ├── D-scaling as convex optimization
    ├── Power iteration as non-convex local search
    │
mini-complex-control-theory/
    ├── 0. mini-general-system-theory (state-space, stability)
    ├── 1. mini-linear-system-theory (A,B,C,D representations)
    ├── 12. mini-hinf-control (H-infinity theory → this module)
    │
    └──< THIS MODULE: mini-structured-singular-value
```

## What This Module Provides

```
mini-structured-singular-value
    │
    ├── mu-Analysis (ssv_bounds.c)
    │   ├── Upper bound: D-scaling optimization
    │   ├── Lower bound: Power iteration
    │   ├── Robust stability test (Main Loop Theorem)
    │   ├── Robust performance test (augmented Delta)
    │   └── Mixed mu (D,G-scaling)
    │
    ├── Uncertainty Modeling (ssv_uncertainty.c)
    │   ├── Block types: scalar real/cpx, full real/cpx
    │   ├── D-scaling matrix construction & optimization
    │   ├── Doyle benchmark structure
    │   └── General specification-based construction
    │
    ├── LFT Interconnection (ssv_lft.c)
    │   ├── Upper LFT Fu(M, Delta_u)
    │   ├── Lower LFT Fl(M, Delta_l)
    │   ├── Redheffer star product (series)
    │   ├── LFT chain (multi-system)
    │   ├── Additive/multiplicative/feedback uncertainty LFTs
    │   └── Augmented M for robust performance
    │
    ├── mu-Synthesis (ssv_synthesis.c)
    │   ├── D-K iteration framework
    │   ├── H-infinity norm (Hamiltonian bisection)
    │   ├── Generalized plant construction
    │   ├── Closed-loop system Fl(P,K)
    │   └── Frequency-domain mu analysis grid
    │
    └── Core Linear Algebra (ssv_core.c)
        ├── Complex/real matrix types & operations
        ├── SVD (Golub-Reinsch bidiagonalization)
        ├── Osborne matrix balancing
        ├── mu definition-based computation (n≤3)
        ├── Small-gain comparison
        └── Determinant/inverse for small matrices
```

## Follow-on Modules (Comes After This Module)

```
mini-structured-singular-value (THIS MODULE)
    │
    ├── 14. mini-lmi-robust-control (LMI-based mu upper bounds)
    ├── 15. mini-gain-scheduling (LFT-based gain scheduling)
    ├── 16. mini-nonlinear-robust (IQC-based uncertainty)
    └── 17. mini-stochastic-robust (random uncertainty)
```

## Research Frontiers (L9)

```
Open Problems documented:
    ├── NP-hardness of exact mu computation (Braatz et al., 1994)
    ├── Polynomial-time mu approximation (still open)
    ├── Optimal G-scaling characterization
    ├── mu for time-varying uncertainty (IQC framework)
    └── mu-synthesis for nonlinear systems
```


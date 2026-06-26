# Course Tree — Small Gain Theorem

## Prerequisites (Dependencies)

```
Linear Algebra (vector spaces, matrices, eigenvalues)
    ↓
Signals and Systems (Laplace/Fourier transforms, frequency response)
    ↓
Classical Control (Bode plots, Nyquist, stability margins)
    ↓
State-Space Control (LQR, pole placement, observers)
    ↓
┌─────────────────────────────────────────────────────┐
│          Small Gain Theorem Module                   │
├─────────────────────────────────────────────────────┤
│ L1: LTI system models, norms, uncertainty types      │
│ L2: Feedback interconnection, small gain condition    │
│ L3: Hamiltonian matrix, Riccati, Lyapunov, IQC       │
│ L4: SGT, BRL, Lyapunov theorem, SSV theorem           │
│ L5: H∞ bisection, QR algorithm, ARE, D-K iteration    │
│ L6: Robust stability verification, mu analysis        │
│ L7: Uncertainty weights, gap metric (applications)    │
│ L8: Mu analysis, D-scaling, IQC framework             │
│ L9: Nonlinear ISS, dynamic multipliers (frontiers)    │
└─────────────────────────────────────────────────────┘
    ↓
Advanced Topics:
├── H∞ Controller Synthesis (Doyle-Glover-Khargonekar)
├── μ-Synthesis (D-K iteration, full K-step)
├── LPV Control (Linear Parameter-Varying)
├── Nonlinear Robust Control (ISS small gain)
└── Quantum Robust Control
```

## Internal Dependency Graph

```
sgt_core.h ────────► sgt_core.c (LTI, matrices, vectors, TF)
    │
    ├──► sgt_hinf.h ────► sgt_hinf.c (H∞, BRL, CARE, H2, Hankel)
    │         │
    │         └──► sgt_stability.c (interconnection, BIBO)
    │
    ├──► sgt_stability.h ──► sgt_stability.c (IC, passivity, gap)
    │
    └──► sgt_uncertainty.h ──► sgt_uncertainty.c (mu, D-K, IQC)
              │
              └── depends on sgt_hinf.h
```

## Knowledge Prerequisites by Level

| Level | Prerequisites | Courses |
|-------|--------------|---------|
| L1 | Calculus, Linear Algebra | Freshman/Sophomore |
| L2 | Signals & Systems, Control intro | Junior |
| L3 | Advanced Linear Algebra, Numerical Methods | Junior/Senior |
| L4 | Control Theory, Functional Analysis | Senior/Grad |
| L5 | Numerical Linear Algebra, Optimization | Grad |
| L6 | Robust Control, System Theory | Grad |
| L7 | Application domain knowledge | Grad/Industry |
| L8 | Advanced Robust Control, Convex Optimization | PhD |
| L9 | Research-level expertise | PhD/Postdoc |

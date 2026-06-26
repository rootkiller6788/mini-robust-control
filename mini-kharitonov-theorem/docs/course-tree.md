# Course Tree — Mini-Kharitonov-Theorem

## Prerequisites

This module builds on the following knowledge:

### Required Prior Knowledge
1. **Linear Algebra** — Matrix operations, determinants, Gaussian elimination, eigenvalues
2. **Differential Equations** — Linear ODEs, characteristic equations, stability concepts
3. **Complex Analysis** — Complex numbers, polynomials, roots in complex plane
4. **Classical Control Theory** — Transfer functions, Laplace transform, feedback, PID control

### Strongly Recommended
5. **Interval Arithmetic** — Moore interval analysis (1966)
6. **Numerical Methods** — Root finding, QR algorithm, matrix computations

## Dependency Tree

```
Linear Algebra ──┬── Routh-Hurwitz Criterion (1877/1895)
                 │      └── Hurwitz Determinants
                 │           └── Lienard-Chipart Criterion
                 │
                 ├── Hermite-Biehler Theorem
                 │      └── Even/odd polynomial decomposition
                 │
                 └── Kharitonov Theorem (1978) ───┐
                                                     │
Interval Arithmetic ──┬── Interval Polynomials ─────┤
                      │                              │
                      └── Value Set Computation ─────┤
                                                     │
Control Theory ──┬── Transfer function models ───────┤
                 │      └── DC motor model           │
                 │      └── Quadrotor model          │
                 │      └── Aircraft model           │
                 │                                    │
                 └── PID controller design ───────────┤
                                                      │
                                                      ▼
                                          Kharitonov Robust Stability
                                          Verification Pipeline
                                                 │
                                    ┌────────────┼────────────┐
                                    ▼            ▼            ▼
                              K1 Stability  K2 Stability  K3 Stability  K4 Stability
                                    │            │            │            │
                                    └────────────┴────────────┴────────────┘
                                                      │
                                                      ▼
                                          Family robust stability result
```

## Downstream Dependencies

This module is prerequisite for:
- Edge Theorem (Bartlett-Hollot-Lin, 1988)
- Generalized Kharitonov for control synthesis
- mu-analysis and structured singular value
- Randomized algorithms for robust control
- Kharitonov for time-delay systems

## Learning Path

```
Classical Control → Routh-Hurwitz → Interval Polynomials → Kharitonov Theorem
    → Edge Theorem → Value Set → Generalized Kharitonov → H-infinity/μ
```

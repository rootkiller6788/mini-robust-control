# Course Tree — mini-parametric-uncertainty

## Prerequisites
- Linear Algebra: eigenvalues, eigenvectors, matrix norms
- Control Theory: state-space models, stability, Lyapunov
- Optimization: convex optimization, LMI basics
- Polynomial Theory: root location, Routh-Hurwitz

## Knowledge Dependencies
```
Linear Algebra → Matrix Operations → Eigenvalue Computation
                                      ↓
Control Theory  → State-Space → Uncertain State-Space → Robust Stability
                                      ↓
Lyapunov Theory → Quadratic Stability → LMI Formulations
                                      ↓
Polynomial Theory → Routh-Hurwitz → Kharitonov Theorem → Stability Margin
                                                       ↓
                                              Edge Theorem → Polytopic Stability
```

## This Module Enables
- mini-robust-performance (H2/H∞ under uncertainty)
- mini-mu-synthesis (structured singular value)
- mini-gain-scheduling (parameter-varying control)

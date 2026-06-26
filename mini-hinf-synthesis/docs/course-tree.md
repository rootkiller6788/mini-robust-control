# Course Dependency Tree — mini-hinf-synthesis

## Prerequisites (what this module depends on)

```
Linear Algebra
  |-- Matrix multiplication, inverse, determinant
  |-- Eigenvalues, eigenvectors
  |-- LU decomposition, Cholesky

Numerical Linear Algebra
  |-- QR factorization (Householder)
  |-- Real Schur decomposition
  |-- Singular value decomposition

Classical Control
  |-- State-space representation
  |-- Stability (Hurwitz)
  |-- Controllability, observability

Linear Control Theory
  |-- Lyapunov equation (Gramians)
  |-- Algebraic Riccati equation (LQR/LQG)
  |-- H2 optimal control

Functional Analysis
  |-- H-infinity norm of transfer function
  |-- Hardy spaces H2, H-infinity
  |-- Small-gain theorem
```

## What this module provides (dependents)

```
mini-hinf-synthesis (THIS MODULE)
  |
  +-- Bounded Real Lemma (CT and DT)
  |     |-- H-inf norm computation via bisection
  |     +-- Hamiltonian eigenvalue test
  |
  +-- Algebraic Riccati Equation solver
  |     |-- Schur vector method (Laub 1979)
  |     +-- Matrix sign function method (Roberts 1980)
  |
  +-- DGKF Two-Riccati H-inf synthesis
  |     |-- Controller ARE solution
  |     |-- Filter ARE solution
  |     |-- Coupling condition check
  |     +-- Central controller construction
  |
  +-- LMI-based H-inf synthesis
  |
  +-- Gamma-iteration (bisection)
  |
  +-- Mixed-sensitivity design pipeline
  |
  +-- Discrete-time H-inf synthesis
  |
  +-- Matrix function library
  |     |-- Exponential, sign, sqrt
  |     +-- Lyapunov, Sylvester solvers
  |
  +-- Downstream modules:
        |-- mini-mu-synthesis (structured singular value)
        |-- mini-gap-metric-robustness
        |-- mini-small-gain-theorem
        |-- mini-structured-singular-value
```

## Research Frontiers (L9)

```
Nonlinear H-infinity
  |-- Hamilton-Jacobi-Isaacs equation
  +-- L2-gain analysis

Time-varying H-infinity
  |-- Differential Riccati equation
  +-- Indefinite inner product spaces

LPV H-infinity
  |-- Parameter-dependent Lyapunov functions
  +-- Grid-based LPV synthesis
```

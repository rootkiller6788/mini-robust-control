# Course Tree — mini-mu-synthesis

## Prerequisites
```
mini-linear-system-theory
    ├── State-space models (A,B,C,D)
    ├── Controllability/Observability
    ├── Stability (Lyapunov)
    └── Frequency response (Bode, Nyquist)
        │
        ├── mini-hinf-synthesis
        │   └── H-infinity control (DGKF 1989)
        │       └── K-step of DK-iteration
        │
        ├── mini-structured-singular-value
        │   └── u definition and properties
        │       └── D-step of DK-iteration
        │
        └── mini-mu-synthesis  ← THIS MODULE
            ├── u upper/lower bounds
            ├── DK-iteration
            ├── Robust stability/performance
            └── Uncertainty modeling
```

## What This Module Enables
- Robust controller synthesis for MIMO systems with structured uncertainty
- Aerospace (Boeing 747, SpaceX rocket attitude control)
- Automotive (Tesla vehicle dynamics)
- Process control (chemical plants with uncertain parameters)

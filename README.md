# Mini Robust Control

A collection of **from-scratch, zero-dependency C implementations** of university-level robust control theory. Each module maps to MIT, Stanford, and Caltech courses, bridging theory and practice by translating textbook equations into runnable C code — from the Small Gain Theorem and H-infinity synthesis to μ-analysis, passivity-based control, and Kharitonov's polynomial framework.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|------------|--------|-------------|
| [mini-gap-metric-robustness](mini-gap-metric-robustness/) | Normalized coprime factorization, gap metric, loop shaping robust stability | MIT 6.241J, Caltech CDS 110 |
| [mini-hinf-synthesis](mini-hinf-synthesis/) | Bounded Real Lemma, H-infinity norm computation, Riccati equation solvers | MIT 16.30, Stanford AA278 |
| [mini-kharitonov-theorem](mini-kharitonov-theorem/) | Kharitonov polynomials, interval polynomial Hurwitz stability, parametric robustness verification | MIT 6.241J, UT Austin ME 381R |
| [mini-mu-synthesis](mini-mu-synthesis/) | Structured singular value μ, D-K iteration, mixed-μ synthesis, robust performance | Stanford AA278, Caltech CDS 110 |
| [mini-parametric-uncertainty](mini-parametric-uncertainty/) | Polytopic uncertainty, LMI stability conditions, Kharitonov-based robust analysis | MIT 16.30, MIT 16.31 |
| [mini-passivity-based-control](mini-passivity-based-control/) | Passivity theorem, energy shaping, port-Hamiltonian systems, IDA-PBC | MIT 6.241J, Imperial College EE4-48 |
| [mini-small-gain-theorem](mini-small-gain-theorem/) | Small Gain Theorem, input-output stability, interconnection robustness | MIT 6.241J, Caltech CDS 110 |
| [mini-structured-singular-value](mini-structured-singular-value/) | Structured singular value (μ) analysis, upper/lower bounds, LFT uncertainty modeling | Stanford AA278, MIT 16.30 |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every module includes `docs/` with course-alignment notes and textbook references (Zhou–Doyle–Glover, Skogestad–Postlethwaite, Green–Limebeer)
- **Practical demos** — aircraft pitch control, DC motor robust tuning, quadrotor stability analysis, maglev control, power converter regulation

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-gap-metric-robustness
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-robust-control/
├── mini-gap-metric-robustness/    # Coprime factorization, gap metric, loop shaping robustness
├── mini-hinf-synthesis/           # H-infinity norm, Bounded Real Lemma, Riccati solvers
├── mini-kharitonov-theorem/       # Kharitonov's theorem, interval polynomial Hurwitz stability
├── mini-mu-synthesis/             # μ-synthesis, D-K iteration, robust performance
├── mini-parametric-uncertainty/   # Polytopic uncertainty, LMI conditions, parametric stability
├── mini-passivity-based-control/  # Passivity theorem, energy shaping, IDA-PBC
├── mini-small-gain-theorem/       # Small Gain Theorem, interconnection robustness
└── mini-structured-singular-value/ # Structured singular value μ, LFT uncertainty
```

## License

MIT

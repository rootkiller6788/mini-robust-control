# Knowledge Graph — mini-parametric-uncertainty

## L1: Definitions ✅ Complete
- PU_Parameter, PU_ParameterVector (uncertain parameters)
- PU_IntervalPolynomial, PU_KharitonovPolynomial (polynomial models)
- PU_PolytopicPolynomial, PU_VertexPolynomial (polytopic model)
- PU_UncertainStateSpace, PU_AffineStateSpace (state-space models)
- PU_IntervalMatrix (interval matrix)
- PU_ValueSet, PU_RouthArray (analysis structures)
- PU_LMI_Constraint, PU_LMI_Problem (LMI data)
- PU_Polytope (convex polytope)
- PU_LyapunovData, PU_QuadraticStabilityLMI (stability data)
- PU_RobustStabilityResult (stability analysis result)

## L2: Core Concepts ✅ Complete
- Robust stability (under all admissible uncertainty)
- Zero exclusion principle
- Value set mapping
- Quadratic stability (common Lyapunov function)
- Structured singular value (mu)
- Stability radius
- Small-gain theorem
- Parameter-dependent Lyapunov functions

## L3: Mathematical Structures ✅ Complete
- Interval polynomials (coefficient-wise independent intervals)
- Kharitonov corner polynomial pattern (modulo-4 selection)
- Polytopes in parameter space (convex hull of vertices)
- Affine parameter-dependent state-space: A(q)=A0+ΣqᵢAᵢ
- Routh array construction and analysis
- LMI structure: F(x)=F0+ΣxᵢFᵢ < 0
- Barycentric coordinates for polytope membership
- Gershgorin circle theorem for interval matrices

## L4: Fundamental Laws ✅ Complete
- Kharitonov's Theorem (1978): 4 corner polynomials ⇔ robust stability
- Edge Theorem (Bartlett-Hollot-Lin 1988): edges ⇒ polytope
- Zero Exclusion Principle: value set origin exclusion + nominal stability
- Quadratic Stability Theorem (Lyapunov 1892): common P matrix
- Bounded Real Lemma (H∞ norm characterization via LMI)
- Routh-Hurwitz Criterion: sign changes in first column
- Jury Stability Criterion: Schur stability for discrete-time
- Bialas Counterexample: vertex check insufficient for interval matrices

## L5: Algorithms/Methods ✅ Complete
- QR eigenvalue algorithm (Francis 1961)
- Routh-Hurwitz array computation
- Kharitonov stability test
- Bialas vertex enumeration check
- Stability margin bisection
- Alternating projection LMI solver
- Primal-dual interior point LMI solver
- Gift Wrapping (Jarvis March) for 2D convex hull
- Chebyshev center via simplified simplex
- Branch-and-bound robust stability verification
- Monte Carlo stability simulation
- Lyapunov equation via Kronecker product

## L6: Canonical Problems ✅ Complete
- DC motor robust speed control
- Mass-spring-damper robust stabilization
- F-16 short-period pitch damper
- SMIB power system small-signal stability

## L7: Applications ✅ Complete (4)
- DC motor (electromechanical)
- MSD system (mechanical)
- F-16 flight control (aerospace) — Sobel, Shapiro, Andry 1994
- SMIB power system (power) — Kundur 1994

## L8: Advanced Topics ⚠️ Partial (2/5)
- ✅ Bialas counterexample detection + vertex check
- ✅ Parameter-dependent Lyapunov function LMI
- ❌ mu-synthesis (D-K iteration) — documented only
- ❌ Randomized algorithms for robust control — documented only
- ❌ IQC (Integral Quadratic Constraints) — documented only

## L9: Research Frontiers ⚠️ Partial
- mu-analysis bounds computation (simplified frequency gridding)
- Monte Carlo verification methods
- Future: Polynomial chaos for stochastic parametric uncertainty

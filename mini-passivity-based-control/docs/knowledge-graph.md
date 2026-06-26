# Knowledge Graph — Passivity-Based Control

## L1: Definitions
- Passivity (Willems 1972): S(x(t)) - S(x(0)) <= Integral_0^t u^T(s) y(s) ds
- Strict passivity, input/output strict passivity, lossless, finite L2-gain
- Storage function: quadratic, kinetic, total energy, Hamiltonian
- Supply rate: w(t) = u^T(t) y(t)
- Dissipation inequality: dS/dt <= u^T y - phi(x), phi(x) > 0

## L2: Core Concepts
- Energy-based control: exploit physical energy for stabilization
- Passivity Theorem: feedback interconnection preserves passivity
- KYP Lemma: bridge between frequency-domain and state-space passivity
- Relative degree, minimum-phase, positive-real systems
- Port-Hamiltonian framework: J, R, g, H representation

## L3: Mathematical Structures
- Matrix/Vector algebra (row-major dense storage)
- Symmetric, skew-symmetric, positive (semi)definite checks
- Storage function evaluation, gradient, Hessian
- KYP block matrix construction
- Dirac structure: D in F x F* satisfying D = D^perp

## L4: Fundamental Theorems
- KYP Lemma (Anderson 1967): LTI passive iff exists P >= 0 with KYP matrix <= 0
- Passivity Theorem (van der Schaft 2000): negative feedback of passive systems is passive
- Power balance: dH/dt = -grad_H^T R grad_H + y^T u <= y^T u (R >= 0)
- LaSalle invariance principle for energy-based stability
- M_dot - 2C skew-symmetric (passivity of EL mechanical map)

## L5: Algorithms/Methods
- Energy shaping + damping injection: tau = g(q) - Kp (q-qd) - Kd q_dot
- IDA-PBC: solve matching equation [J-R]grad_H + g u = [J_d-R_d]grad_H_d
- Nonlinear damping: tanh saturation, power-law, custom
- Gravity compensation for n-link serial robots
- Trajectory simulation via RK4 with torque callback

## L6: Canonical Problems
- Mass-spring-damper (1-DOF pH model)
- RLC circuit (series/parallel pH model)
- DC motor (electromechanical coupling)
- 1-DOF pendulum, 2-DOF planar robot
- DC-DC boost/buck/buck-boost converters

## L7: Applications
- Robot manipulator control (Takegaki-Arimoto 1981)
- DC-DC power converter PBC (Escobar et al. 1999)
- Motor speed and position PBC (Ortega et al. 1998)
- Quadrotor UAV attitude stabilization
- Magnetic levitation gap control
- pH neutralization process control

## L8: Advanced Topics
- Dirac structures as generalization of power-conserving interconnections
- Bond graphs for multi-domain physical system modeling
- Passivity indices (IFP, OFP, ISP, OSP)
- Control performance metrics (IAE, ISE, ITAE, settling time, overshoot)
- Nonlinear passivity test via random trajectory sampling

## L9: Research Frontiers
- IDA-PBC for underactuated mechanical systems
- Passivity-based adaptive control
- Networked pH systems and compositional passivity
- Stochastic port-Hamiltonian systems

#!/usr/bin/env python3
"""
plot_mu.py -- Structured Singular Value (mu) Visualization

Computes and plots the frequency response of a system's structured
singular value for robustness analysis. Demonstrates the relationship
between H-infinity norm and mu.

Requires: numpy, matplotlib (optional -- saves CSV if not available)

Ref: Packard & Doyle (1993) The complex structured singular value,
     Automatica 29(1):71-109
"""

import math
import sys

def tf_eval(A, B, C, D, w):
    """Evaluate G(jw) = C(jwI-A)^{-1}B + D for SISO."""
    # For SISO: G(jw) = C * B / (jw - A) + D
    # Complex division: (a + jb) / (c + jd)
    a, b = C * B, 0.0
    c, d = -A, w
    den = c*c + d*d
    if den < 1e-16:
        return D
    real = (a*c + b*d) / den + D
    imag = (b*c - a*d) / den
    return math.sqrt(real*real + imag*imag)

# SISO system: G(s) = 1/(s+a)
A, B, C, D = -2.0, 1.0, 1.0, 0.0

# Frequency grid
freqs = [0.01 * (10.0 ** (4.0 * i / 200)) for i in range(201)]
responses = []

for w in freqs:
    mag = tf_eval(A, B, C, D, w)
    responses.append((w, mag))

# Output CSV
outpath = sys.argv[1] if len(sys.argv) > 1 else "mu_response.csv"
with open(outpath, "w") as f:
    f.write("frequency_rad_s,gain_abs\n")
    for w, mag in responses:
        f.write(f"{w:.6f},{mag:.8f}\n")

print(f"Frequency response data written to {outpath}")
print(f"  Frequencies: {len(freqs)} points, {freqs[0]:.4f} to {freqs[-1]:.4f} rad/s")
print(f"  Peak gain: {max(r[1] for r in responses):.6f}")

try:
    import matplotlib.pyplot as plt
    ws = [r[0] for r in responses]
    ms = [r[1] for r in responses]
    db = [20 * math.log10(m + 1e-16) for m in ms]
    plt.figure(figsize=(10, 6))
    plt.semilogx(ws, db, 'b-', linewidth=2)
    plt.grid(True, which='both', alpha=0.3)
    plt.xlabel('Frequency (rad/s)')
    plt.ylabel('Magnitude (dB)')
    plt.title('Frequency Response |G(j\omega)|')
    plt.axhline(y=-3, color='r', linestyle='--', alpha=0.5, label='-3 dB')
    plt.legend()
    plt.savefig('mu_response.png', dpi=150)
    print("  Plot saved to mu_response.png")
except ImportError:
    print("  (matplotlib not available; CSV data only)")

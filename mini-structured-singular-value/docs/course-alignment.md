# Course Alignment — mini-structured-singular-value

Mapping to nine world-class university courses on robust control and mu-theory.

## MIT

| Course | Topic | Implementation |
|--------|-------|---------------|
| 6.241 — Dynamic Systems & Control | LFT uncertainty modeling | `ssv_lft.c` — additive, multiplicative, feedback LFTs |
| 6.242 — Advanced Linear Control | mu-analysis framework | `ssv_bounds.c` — upper/lower bounds |
| 6.832 — Underactuated Robotics | Robust stability margin | `ssv_robust_stability_test()` |
| 16.30 — Feedback Control Systems | D-K iteration | `ssv_dk_synthesize()` |

## Stanford

| Course | Topic | Implementation |
|--------|-------|---------------|
| ENGR 205 — Intro to Control Design | Small-gain theorem | `ssv_small_gain_bound()` |
| AA 278 — Advanced Control | D-scaling optimization | `ssv_dscale_create/apply/evaluate()` |
| EE 363 — Linear Dynamical Systems | SVD and matrix norms | `ssv_cmatrix_norm2()`, `ssv_svd_compute()` |

## Berkeley

| Course | Topic | Implementation |
|--------|-------|---------------|
| ME 232 — Advanced Control Systems | Structured uncertainty | `SSVUncertaintyStructure`, block types |
| ME 233 — Robust Control | mu bounds computation | Power iteration + D-scaling |
| EE 221A — Linear Systems | State-space frequency response | `ssv_ss_freqresp()` |

## CMU

| Course | Topic | Implementation |
|--------|-------|---------------|
| 24-774 — Robust Control | mu-synthesis design | `ssv_dk_step()`, D-K iteration |
| 18-771 — Linear Systems | Matrix algebra for control | `ssv_cmatrix_*` operations |

## Princeton

| Course | Topic | Implementation |
|--------|-------|---------------|
| MAE 546 — Robust Control | Mixed mu analysis | `ssv_mixed_mu_upper_bound()` |
| MAE 542 — Advanced Dynamics | LFT interconnection | `ssv_lft_star_product()` |

## Caltech

| Course | Topic | Implementation |
|--------|-------|---------------|
| CDS 110 — Intro to Control | Uncertainty modeling | 4 uncertainty LFT constructions |
| CDS 212 — Robust Control | mu theory fundamentals | `ssv_mu_analysis()`, mu definition |
| ACM 217 — Advanced Topics | Hamiltonian bisection for H∞ | `ssv_hinf_norm()` |

## Cambridge

| Course | Topic | Implementation |
|--------|-------|---------------|
| 4F2 — Robust & Nonlinear Control | mu-analysis and synthesis | Full mu-analysis pipeline |
| 4F3 — Advanced Control | D-K iteration convergence | `ssv_dk_synthesize()` with convergence check |

## Oxford

| Course | Topic | Implementation |
|--------|-------|---------------|
| B16 — Control Systems | Frequency-domain analysis | `ssv_ss_freqresp_grid()` |
| C20 — Robust Control | Structured singular value | Core mu computation |

## ETH Zurich

| Course | Topic | Implementation |
|--------|-------|---------------|
| 227-0216-00L — Control Systems II | mu theory | Complete mu analysis toolkit |
| 227-0690-00L — Advanced Control | D-K synthesis | `ssv_dk_step()` iterative design |
| 227-0225-00L — Linear System Theory | Osborne balancing | `ssv_osborne_balance()` |

## Key Textbooks Aligned

| Textbook | Chapters Covered |
|----------|-----------------|
| Zhou, Doyle & Glover — Robust and Optimal Control (1996) | Ch 9 (LFT), Ch 10 (mu), Ch 11 (mu-synthesis), Ch 17 (D-K) |
| Skogestad & Postlethwaite — Multivariable Feedback Control (2005) | Ch 8 (uncertainty and robustness) |
| Doyle, Francis, Tannenbaum — Feedback Control Theory (1992) | Ch 6 (robustness) |
| Packard & Doyle — The Complex Structured Singular Value (1993) | Core mu theory |
| Golub & Van Loan — Matrix Computations (2013) | SVD implementation |


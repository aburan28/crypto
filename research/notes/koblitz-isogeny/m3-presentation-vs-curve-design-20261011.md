# m = 3 curve-vs-presentation at n = 19 (design, 2026-10-11)

This is an additive amendment to `presentation-vs-curve-design-20261003.md`. The model,
the estimator and the thresholds are unchanged:
- **H0:** the 95 % upper bound on ICC_curve is < 0.05 on every metric.
- **H1:** any cost metric has ICC ≥ 0.10 with its lower bound > 0.
- Otherwise: **inconclusive**.

The earlier m = 3 arms were dropped for cost. The e₃-enumeration solver
(`m3-solver-design-20261004.md`, Results) makes them affordable.

## Arm

| item | value |
|---|---|
| class | n = 19, a₂ = 1, from the walk |
| sample | 32 curves (the sweep's seeded hash, K_a included) |
| subspaces | geometric V_s, s = 1..4 |
| decomposition | m = 3, l = 6, `M3Solver::SymmetricE3`, PCLMUL field arithmetic |
| probes | 400 full-group probes per cell |
| order | cells shuffled with a fixed seed |
| command | `KOBLITZ_PRES_M3=e3 KOBLITZ_PRES_GEOMETRIC=1 KOBLITZ_PRES_FFD_TARGETS=0 koblitz_presentation_sweep 19 1 6 400 32 1,2,3,4` |
| output | `experiments/koblitz_presentation/pres_19_1_l6_geom_m3e3.jsonl` |
| analysis | `examples/koblitz_presentation_icc.rs`, the Rust port of `presentation_icc.py`; on the n = 23 data it reproduces the Python point estimates exactly |

## Metrics

- **yield per |F|³:** relations / probes / |F|³.
- **ln µs per call.**
- **work per call:** S₄ roots per call, which is deterministic given the target.
- **ln projected solve.**

First fall is not measured, as in the n = 31 arm.

## Expected precision

Pilot full solves at n = 19 gave about 0.05 relations per probe, so a cell has about
20 relations. Yield is therefore the noisiest metric.

## Budget

About 128 × 400 × 0.12 s ≈ 1.7 core-h, on 4 cores.

## Results (2026-10-11)

The run completed all 128 cells (32 curves × 4 V) in 1550 s on 4 cores.
- **Soundness:** every cell matched BSGS, with 0 inconsistent relations.
- **Yield:** 4180 relations from 51200 probes (0.082 per probe).
- **Cost:** median 121 ms per call and 15427 S₄ roots per call.

Analysis output: `experiments/koblitz_presentation/icc_n19_m3_geom.json`
(`examples/koblitz_presentation_icc.rs`).

| metric | ICC_curve | 95 % CI | subspace share |
|---|---|---|---|
| yield per \|F\|³ | 0.098 | [0, 0.224] | 0.002 |
| ln µs per call | 0.011 | [0, 0.137] | 0.000 |
| work per call (S₄ roots) | 0.000 | [0, 0.051] | 0.000 |
| ln projected solve | 0.000 | [0, 0.069] | 0.000 |

**Verdict: inconclusive by the pre-registered rule.**
- H1 is not supported: no metric has ICC ≥ 0.10 with a lower bound above 0.
- H0 is not established: no upper bound is below 0.05. Work per call comes closest,
  at 0.051.

Plainly, no curve effect is detected at m = 3 either. The point estimates are ≤ 0.10,
and the curve accounts for at most about 14 % of per-call time and 22 % of yield
variance at 95 % confidence. The four geometric V do not differ: subspace share
≤ 0.002 on every metric.

Yield is the loosest bound, as it was at m = 2, n = 23 (0.078 [0, 0.162]). About 33
relations per cell bound its precision.

**Scope:** n = 19, a₂ = 1, l = 6, m = 3, geometric V, 400 probes, this solver.

With `PinX3` (see `m3-solver-design-20261004.md`), the same arm at n = 31 costs about
7 ms per call, so it is now within reach.

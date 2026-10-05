# Measured full solves at n = 16, 17, 19 (2026-10-04)

Replaces the *projected* solve cost (`log_projected_solve`) of the presentation
study with measured end-to-end index-calculus solves.

**Setup.** m = 2, full-group probes. l = 8 (n = 16, 17) and l = 9 (n = 19). Each cell
runs 50 yield probes, then the closure phase until the log is recovered. The cap is
200 000 trials, and no cell reached it.

**Sample.** The first 32 curves of each existing presentation sample (a₂ = 1).

**Subspaces, 5 per curve:**
- mono: ⟨1..z^{l−1}⟩;
- geometric: g1, g2;
- random: 1, 2.

**Order.** The 480 cells ran in a seeded shuffle, 3 in parallel. Part of the run shared
the machine with other jobs, so wall-clock comparisons across families are indicative.
The shuffle protects the curve-vs-V split from that.

**Data and code.**
- `experiments/koblitz_full_solves/results.jsonl`: one line per cell, the full pilot row.
- `cells.txt` and `run.sh`: the run itself.
- `icc_full_solves.py`, producing `icc_full_solves.json`.

## Soundness

All 480 solves closed. Each recovered log equals both the planted secret and the BSGS
log (`verified`), with 0 wrong.

## Cost by family (median seconds per full solve; median closure trials)

| n | mono | geometric | random |
|---|---|---|---|
| 16 | 0.60 s; 20 | 0.79 s; 4 | 1.10 s; 14 |
| 17 | 3.65 s; 265 | 3.50 s; 239 | 4.03 s; 260 |
| 19 | 14.8 s; 572 | 14.8 s; 545 | 39.9 s; 541 |

Closure trials do not depend on the family; the per-call cost does. At n = 19,
geometric V solves as fast as monomial V and about 2.7× faster than random V.

## Curve effect on measured solves (ICC_E, 95 % bootstrap CI over curves)

| n | V set | log seconds | log Gröbner seconds | closure trials |
|---|---|---|---|---|
| 16 | all 5 | 0.030 [0, 0.105] | 0.024 [0, 0.095] | 0.034 [0, 0.141] |
| 17 | all 5 | 0.091 [0, 0.219] | 0.091 [0, 0.220] | 0.172 [0, 0.354] |
| 19 | all 5 | 0.000 [0, 0.006] | 0.000 [0, 0.005] | 0.000 [0, 0.000] |

Per-family 2-V grids are in `icc_full_solves.json`. They have wider intervals; the
largest point estimate is 0.30 [0, 0.60], for n = 17 random V on closure trials.

## Reading

At n = 19, measured end-to-end solve cost shows no curve effect. The upper bound is
≤ 0.006 on every metric, which clears the pre-registered H0 threshold of 0.05 for this
instrument at this size. At n = 16 and 17 the point estimates are small, but the upper
bounds (0.10–0.35) do not clear 0.05: inconclusive there.

Scope: m = 2, the stated l, a₂ = 1 classes, 32 curves, this solver.

# Run report — factor-base → D* campaign (run-20261014T0600Z)

**Engine:** `sage/factor_base_dor.sage` (branch `feat/factor-base-dor-campaign-20261013`),
selftest 32 cases / 0 failures before and after the run (`--field 11 --dim 3 --selftest 8`).
**Orchestrator:** `run_campaign.sh` (thin shell over the sage engine only).
**Host:** macOS arm64 (Apple silicon), single sage process per cell, no CPU pinning
(descriptive study; D* and Macaulay shapes are deterministic counts, wall time is not used).
**Raw records:** `raw/*.json` (42 files), logs `raw/*.log`, aggregates `analysis.txt`,
`figures/data.json`.

**Arms and counts.** m=2 grid: n ∈ {11, 13} × dims 2..(n−1)/2 × structures
{standard, random, nb_window, nb_break} × seeds {1,2,3} × 24 targets/cell — 27 cells,
all completed, zero censoring above dmax=6. Subfield closure contrast: n ∈ {9, 15},
dims {3, 5}, structures {standard, random, subfield} × 3 seeds — 9 cells completed.
δ-sweeps (standard, m=2): n=11 and n=13, δ ∈ {−1,+1,+3} integer cells — completed.
Symmetric S4 (m=3 direct): (n=9, l=2), (n=15, l=3), structures {geometric, standard,
subfield, random}, 16 targets — completed. n=131 applicability checks: Frobenius-window
verification at n=131 (log record, engine writes stderr), weight-class stacking at
n=11, n=13 — completed.

**Failures (preserved):** chained m=3 cells (`delta_n11_m3`, `m3_n11_d3_s1..s3`,
`m3_n13_d4_s1..s3`) exceeded the engine's practical budget (>1.8 h CPU per cell,
killed; logs in `raw/`). The chained m=3 D* = 5 figure is prior evidence
(`sage/factor_base_dor_results.json`, commit f729e93fc), not re-measured here.
`frobwin_n131` was killed once by the same cell-killer and rerun cleanly.

**Measured values** (all recorded in `analysis.txt`, section by section):
D* mean/min/max per cell (refutation degree on non-decomposable targets; unsat target
counts 3–24 per cell, sample = 24 targets minus decomposable ones),
decomposable fraction per cell, Macaulay rows×cols at D* per target,
work-per-relation proxy = mean(rows×cols at D*) / fraction.

**Derived quantities (arithmetic shown in `analyze.sage`):** Pearson correlations over
the pooled non-subfield m=2 cells (45 cell aggregates over 3 seeds); best-dim per
(n, structure) by work-per-relation; δ ratios against δ = −1.

**Missing data:** chained m=3 D* at n=11/13 (killed, see failures); sym4 random at
n=15 censored (D* > dmax 5); no wall-clock claims.

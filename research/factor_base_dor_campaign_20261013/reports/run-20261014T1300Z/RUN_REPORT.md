# Run report — factor-base → D* campaign, round 2 (run-20261014T1300Z)

**Engine:** `sage/factor_base_dor.sage` (unchanged from round 1; selftest clean in
round 1, re-run clean after the round: 32 cases / 0 failures).
**Orchestrators:** `run_round2.sh` then `run_round2b.sh` (same cells, `TMPDIR`
redirected to `/Volumes/SSD990-2/tmp` after the first pass failed on a full
temp volume).
**Host:** same as round 1 (macOS arm64, single sage process per cell, no pinning;
deterministic counts, no wall-clock claims).
**Raw records:** 35 round-2 cells under `raw/` (n=17/19/23 grids seeds 1–2,
δ sweeps n=15/17/19, sym4 cells n=15 l=2, n=19 l=3, n=21 l=3), first-pass
failure logs preserved (`round2.log` FAILED entries: 35 cells, all
`OSError: No space left on device` from the sage temp volume; all recovered in
the rerun with 0 failures).

**Arms and counts.** m=2 grid: n ∈ {17, 19, 23} × dims 2..8 (n=23 capped at 8 of
11) × {standard, random, nb_window, nb_break} × seeds {1, 2} × 16 targets.
δ-sweeps (standard, m=2): n=15 (δ ∈ {−1, +1, +3, +5}), n=17 (δ ∈ {−1, +1, +3, +5,
+7}), n=19 (δ ∈ {−1, +1, +3, +5, +7, +9}). Symmetric S4: (n=15, l=2), (n=19, l=3)
{geometric, standard, random}, (n=21, l=3) {geometric, standard, subfield, random};
16 targets, `--dstar-sample 6`, dmax 5.

**Measured values:** as in round 1 (D* mean/min/max, decomposable fraction,
Macaulay rows×cols at D*, work-per-relation proxy), recorded in `analysis.txt`
and `figures/data.json` (133 m=2 aggregates total across rounds).

**Derived:** pooled Pearson correlations recomputed over all aggregates
(r(D*, n′) = 0.82, r(D*, n) = −0.01); best-dim optima per (n, structure);
δ gain ratios normalized per curve; the Figure 4 cell table.

**Missing data:** n=23 dims 9–11 (capped for engine budget); sym4 at n=19/21
l=2 (not run); no chained m=3 (see round 1). Censored cells: sym4 random at
n=15 l=3 remains D* > dmax 5.

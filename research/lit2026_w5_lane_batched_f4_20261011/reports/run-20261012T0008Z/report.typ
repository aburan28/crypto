#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9.5pt)
#set table(stroke: 0.4pt)

= lit2026 W5 run report: run-20261012T0008Z

*Run.* `scripts/run_arms.sh` with `ARMS="orig prev gate rref gaterref"`, `ROUNDS=5`, started 2026-10-12 00:08 UTC (launched by `scripts/launch_run2.sh`). Binaries: the same as run-20261011T2225Z (candidate `f6_ic_probe` at `d22fb3b52`, original at `4b6e09b78`). Host: shared Apple M4 Pro, 1-minute load 76.5–131.8 over the 100 samples (median 100.3); a build of the probe example ran during the first rounds. Isolation *L0*, unpinned, arm order rotated each round. All 100 runs exited 0. Raw records: `raw/`; summary: `summary.json` (`f4_lane_probe summarize`); graph `arms.svg` (`../../plot.py`).

*Workloads.* As run-20261011T2225Z: the four `f6_ic_probe` cells (K1/2#super[17] dim 6 `icv1-f2m17-tm101-00378d4e`, 8 targets; K1/2#super[23] dim 7 `icv1-f2m23-tm5197-1f85e9e1`, 4; K0/2#super[23] dim 8 `icv1-f2m23-t5197-69e76b73`, 4; K0/2#super[31] dim 9 `icv1-f2m31-tm90707-c95f16f5`, 1, budget 2000), prime-order subgroups of `KoblitzCurve::new`, standard subspace, $m = 3$, engines `f4` and `f6p`, `--reps 3`.

*Arms.* `orig`, `prev`, `gate` (`KIC_F6_GATE_FIRST=1`), `rref` (`KIC_F4_INHERIT_ROOT=rref`, an existing switch: fully reduced root), `gaterref` (both).

== Stage unit (identical in all 5 rounds of every arm; equal to run-20261011T2225Z for the arms both ran)

#table(columns: 7, align: (left, left, right, right, right, right, right),
  [cell], [engine], [orig = prev], [gate], [rref], [gaterref], [orig / gaterref],
  [K1/2^17], [f4], [15 015 618], [15 015 618], [14 098 582], [14 098 582], [1.0650],
  [K1/2^17], [f6p], [12 635 914], [8 061 140], [11 600 586], [7 696 852], [1.6417],
  [K1/2^23], [f4], [53 273 191], [53 273 191], [44 709 930], [44 709 930], [1.1915],
  [K1/2^23], [f6p], [51 918 433], [30 697 846], [43 312 007], [26 902 588], [1.9299],
  [K0/2^23], [f4], [158 512 775], [158 512 775], [123 273 212], [123 273 212], [1.2859],
  [K0/2^23], [f6p], [124 949 801], [71 049 603], [89 627 327], [56 940 731], [2.1944],
  [K0/2^31], [f4], [281 583 026], [281 583 026], [190 512 003], [190 512 003], [1.4780],
  [K0/2^31], [f6p], [278 237 634], [156 423 046], [187 048 754], [115 269 611], [2.4138],
)

== Process CPU, previous / arm (L0; median [min–max] of 5 per-round ratios)

#table(columns: 5, align: (left, left, right, right, right),
  [cell], [engine], [gate], [rref], [gaterref],
  [K1/2^17], [f4], [0.95 [0.95–1.00]], [1.12 [1.09–1.22]], [1.13 [1.07–1.17]],
  [K1/2^17], [f6p], [2.27 [2.12–2.50]], [1.22 [1.18–1.33]], [2.85 [2.65–3.15]],
  [K1/2^23], [f4], [0.98 [0.94–1.02]], [1.29 [1.27–1.38]], [1.32 [1.28–1.35]],
  [K1/2^23], [f6p], [2.56 [2.42–2.66]], [1.40 [1.33–1.46]], [3.49 [3.26–3.69]],
  [K0/2^23], [f4], [1.00 [0.96–1.07]], [1.37 [1.34–1.48]], [1.38 [1.33–1.41]],
  [K0/2^23], [f6p], [2.27 [2.18–2.33]], [1.58 [1.48–1.65]], [3.41 [3.39–3.53]],
  [K0/2^31], [f4], [0.98 [0.96–1.04]], [1.68 [1.62–1.75]], [1.71 [1.61–1.79]],
  [K0/2^31], [f6p], [2.11 [2.05–2.19]], [1.74 [1.71–1.92]], [3.69 [3.50–3.96]],
)

The original binary does not record CPU time. Wall per round (original / arm) is in the graph; the original/previous null control ranges 0.25–2.48 per round, so wall admits no claim.

#figure(image("arms.svg", width: 100%), caption: [run-20261012T0008Z, 5 interleaved rounds: stage unit original / arm (top), wall original / arm (middle, log, median and min–max), process CPU previous / arm (bottom). L0, load 76–132.])

== Identity

Per target, arm and round: hit, exhausted, reductions, splits, propagations, infeasible, geometric refutations and witnesses, partial refutations, group additions, pair-index builds and lookups identical across all 5 arms and 5 rounds in all 8 cell/engine pairs; the relation (K1/2^17 target 3, base indices [21, 48, 57]) identical across the candidate arms.

== Paths

`raw/r<round>-<cell>-<arm>.jsonl` (first line: `uptime`, UTC), `raw/*.status`, `summary.json`, `arms.svg`.

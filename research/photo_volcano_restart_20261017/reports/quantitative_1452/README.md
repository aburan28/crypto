# Full rho results and partial BSGS snapshot

Snapshot: **2026-10-17 14:50:43 UTC**. This is the fresh experiment on the
photographed curve, with 100 matched problems on each of 46 curves.

| Method | Verified full logarithms | Requested | Status |
| --- | ---: | ---: | --- |
| Rho | 4600 | 4600 | Complete for this measured phase |
| BSGS | 3297 | 4600 | Remaining jobs resumed; later progress is separate |
| IC | 0 | 4600 | Actual bounded diagnostic attempts completed; full recovery unresolved |

[Run report, 13 inspected pages](RUN_REPORT.pdf) ·
[Results report, eight inspected pages](RESULTS_REPORT.pdf) ·
[Source-linked HTML report](report.html)

Across curves, the 100-case rho means range from **107895789.39** to
**134494184.67 group operations**, ratio **1.246519**. Individual rho cases
range from **3759524** to **429007765** operations. Paired BSGS online counts
agree across every observed generator and curve. These are descriptive
measurements; no intrinsic hardness difference or speedup is established.

Every IC target GB probe used **217410 field multiplications**, while wall
times ranged **11.519667–266.013969 ms** on this shared L0 host. Random setup
and target relation yields were zero. All **14720 actual polynomial systems**
passed independent Singular ideal-equality replay. Hilbert cutoff **13**,
temporary accumulator maximum **23**, and processed-pair sugar maximum **26**
are different recorded quantities; none is labeled CM regularity.

## Tables and exact graph data

- [Per-curve and per-generator native min/mean/max](per_curve.csv)
- [Each problem across the measured curves](per_problem.json)
- [All 100 IC target problems across all 46 curves](ic_per_problem.csv)
- [IC per-curve stage metrics](ic_per_curve.json)
- [All 90 directed map summaries](map_per_edge.json)
- [Map summaries for all 100 problems](map_per_problem.json)
- [Full group invariant factors and subgroup orders](group_orders.csv)
- [Evidence identities](evidence.json) and [publication scope](PUBLIC_SCOPE.json)

## Offline heatmaps

Each heatmap's `index.html` is self-contained, works locally, and provides an
exact-curve inspector. Its SVG retains every vertex and directed map; the
sidecar graph_data.json retains full identities, sample counts and metric
min/mean/max. Node color is the selected curve mean. A shared color range
covers all five levels; missing values remain explicit.

- [Rho operations](rho_operations/index.html) · [SVG](rho_operations/volcano.svg)
- [Rho elapsed time](rho_wall/index.html) · [SVG](rho_wall/volcano.svg)
- [BSGS online operations](bsgs_operations/index.html) · [SVG](bsgs_operations/volcano.svg)
- [BSGS online elapsed time](bsgs_wall/index.html) · [SVG](bsgs_wall/volcano.svg)
- [IC target GB field operations](ic_operations/index.html) · [SVG](ic_operations/volcano.svg)
- [IC target GB elapsed time](ic_wall/index.html) · [SVG](ic_wall/volcano.svg)
- [IC top-forms Hilbert cutoff](ic_hilbert_cutoff/index.html) · [SVG](ic_hilbert_cutoff/volcano.svg)

BSGS views exclude reused table setup. IC views are GB-stage diagnostics,
not cold complete-recovery costs. No field counters are invented for rho or
BSGS. Raw attempts, private witnesses and RNG journals remain local; this
snapshot does not claim a complete private archive export or full experiment
completion. Later run reports retain their own timestamps and coverage.

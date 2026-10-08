# n=61 one-target IC/rho comparison

This primary workload uses one fresh public target point per solver invocation. The same point is used for three paired observations. IC receives only the point coordinates; the known scalar is retained in `known-answer.txt` as a validation sidecar. Every result passes independent Python GF(2^61) and curve-arithmetic replay.

The candidate uses 73,200 distinct subgroup factor-base points and 600 signed-Frobenius orbit columns. Target-independent base/index/log preparation is excluded from the online IC metric; the producer records it separately. Rho online includes target-specific jump setup, walk, and in-process verification after point construction.

| Pair | IC online (ms) | rho online (ms) | rho / IC |
|---|---:|---:|---:|
| R1 | 75.089 | 3989.175 | 53.13x |
| R2 | 99.561 | 2750.126 | 27.62x |
| R3 | 89.650 | 8832.520 | 98.52x |
| R4 (fresh 2026-10-02) | 58.994 | 7365.182 | 124.85x |

Median paired ratio: **53.13x** (observed range 27.62x–98.52x). IC online ranged 75.089–99.561 ms; rho online ranged 2750.126–8832.520 ms. Median IC target PDP plus relation-check phase: 89.575 ms. Target-independent IC process totals ranged 27.76–33.16 s, supplementary only. This is synthetic n=61 evidence for a 48-bit subgroup, not ECC2K-130 evidence. Complete operation-count comparison is absent, so S is unknown.

**Fresh reproduction R4 (2026-10-02):** rebuilt from the committed source — the IC source and rebuilt binary are byte-identical to the frozen protocol hashes (`bb69bd44…`, `edf9ed31…`). R4 recovers the same scalar on the same public point with the same deterministic relation (227,437 probes, point indices [21852, 25480, 26108, 45574]), IC online 58.994 ms (peak RSS 2.73 GB, zero pair-table entries) vs rho online 7365.182 ms (peak RSS 350 MB, seed 20261002) → 124.85x. Independent Python replay passes 15/15 checks (`independent-replay-r4.json`, written by `validate_r4.py`). A schema-complete claim report (`claim_report_vs_rho.json`) passes the autolab `claim-check` for stage `vs_rho`.

A prior pilot used a point already present in an earlier n=61 fixture panel; its raw runs are retained under `runs/prior-target-reused/` and its result rows are kept in `prior-target-result-rows.json`. Those pilot data are excluded from the primary summary.

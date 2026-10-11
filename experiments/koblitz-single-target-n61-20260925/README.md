# n=61 one-target IC/rho comparison

Four paired runs recovered the same public target point over GF(2^61), with IC online times of 58.994–99.561 ms and rho online times of 2,750.126–8,832.520 ms. IC receives only the point coordinates; the known scalar is retained in `known-answer.txt` for validation. Every result passes independent Python field and curve-arithmetic replay.

The candidate uses 73,200 distinct subgroup factor-base points and 600 signed-Frobenius orbit columns. Target-independent base/index/log preparation is excluded from the online IC metric; the producer records it separately. Rho online includes target-specific jump setup, walk, and in-process verification after point construction.

| Pair | IC online (ms) | rho online (ms) | rho / IC |
|---|---:|---:|---:|
| R1 | 75.089 | 3989.175 | 53.13x |
| R2 | 99.561 | 2750.126 | 27.62x |
| R3 | 89.650 | 8832.520 | 98.52x |
| R4 (fresh 2026-10-02) | 58.994 | 7365.182 | 124.85x |

For the original R1–R3 series, the median paired wall ratio is **53.13x** (range 27.62x–98.52x); R4 records 124.85x separately. The R1–R3 median IC target PDP plus relation-check phase is 89.575 ms. Their target-independent IC process totals range from 27.76 to 33.16 s and are retained as supplementary setup costs. The prime subgroup order is 162,888,033,982,417. These wall ratios are exploratory because the runs have no auditable host-level CPU isolation receipt. IC probes and rho walk steps remain separate native counters in `claim_report_vs_rho.json`; calibrated total work `S` and a controlled online speedup remain unknown.

**Fresh reproduction R4 (2026-10-02):** rebuilt from the committed source — the IC source and rebuilt binary are byte-identical to the frozen protocol hashes (`bb69bd44…`, `edf9ed31…`). R4 recovers the same scalar on the same public point with the same deterministic relation (227,437 probes, point indices [21852, 25480, 26108, 45574]), IC online 58.994 ms (peak RSS 2.73 GB, zero pair-table entries) vs rho online 7,365.182 ms (peak RSS 350 MB, seed 20261002) → 124.85x. Independent Python replay passes 15/15 checks (`independent-replay-r4.json`, written by `validate_r4.py`). The retained `claim_report_vs_rho.json` has status **FAIL** under the current autolab `vs_rho` schema: canonical run identities, exact phase names, replay digests, the resource envelope, and rho policy remain to be supplied from verified receipts.

A prior pilot used a point already present in an earlier n=61 fixture panel; its raw runs are retained under `runs/prior-target-reused/` and its result rows are kept in `prior-target-result-rows.json`. Those pilot data are excluded from the primary summary.

# EXP7 complete Legendre-pattern extension

**Frozen before extension code and measurement:** 2026-10-08 UTC.
**Parent protocol:** [PROTOCOL.md](PROTOCOL.md), committed as
`a80a84e842335425e631fe0d403b7c32c68135db`. The [first measured
run](results.json), sourced from
`2e3e894a45ea3cfd9490349b95081073a85c6a5f`, remains a preserved
pilot because its Legendre case covered only `+++`. This extension completes
the requested Legendre symbol patterns without rewriting that run.

Use the same four curve instances, membership definitions for the other five
coordinate cases, positive control, SHA and random-pair controls, exact group
enumeration, Fourier normalization, pair-sum inverse transform, symmetry
matching, and deterministic null-seeding scheme as the parent protocol.
Replace the single `legendre+++` case by **all eight** cases
`legendre+++`, `legendre++-`, `legendre+-+`, `legendre+--`,
`legendre-++`, `legendre-+-`, `legendre--+`, `legendre---`. The three
signs are the Legendre symbols of `x`, `x+1`, and `x+2` modulo `p` in order.
A zero Legendre value belongs to none of these eight cells. Keep their
individual counts and bitsets even if a cell is sparse or degenerate.

Increase the matched null ensemble from 512 to **2048** per case and curve,
using the same SHA-256 seed formula and exact-size inversion-pair sampling.
For the `random-pairs` observed control, use replicate 2048, disjoint from
its calibration replicates 0 through 2047. There are 13 candidate predicates
(the other five plus eight patterns) on four curves, hence **52 candidate
cells**. A candidate lead requires empirical peak `p ≤ 0.05/52 =
0.00096153846` (Bonferroni). With 2048 nulls, the minimum attainable
empirical `p` is `1/2049 = 0.0004880439`, so this test can cross the threshold
if there are zero null exceedances. As before, a toy lead would require an
independent-seed and larger-null replay before any stronger claim. Require all
four log-interval positive controls to exceed their matched-null 95th
percentiles.

Freeze the revised implementation in a separate commit before running it.
Write the complete extension output to `results_all_patterns.json`, without
overwriting `results.json`. Record the source revision, byte hash, verification
receipts, observed and null values, all failures, and a new visual/PDF report.
No DLP speedup or general flatness assertion is claimed from this toy panel.

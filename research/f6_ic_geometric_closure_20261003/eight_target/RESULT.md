# Eight new public targets: F6-IC variability panel

The [preregistered panel](../EIGHT_TARGET_PROTOCOL.md) froze eight distinct
public hash-to-curve points on `icv1-f2m17-tm101-00378d4e`, with no known
scalar supplied to a solver. The 40 processes all exited zero, recovered a
scalar independently verified against their public point, and closed their
five exclusive online phase ledgers. Each point has two F4 and two F6
fresh-process runs and one F5 run. The same certified 62-point usable base
and 29 folded columns served every arm. [Frozen workloads and hashes](freeze.tsv),
[one row per run](measurements.jsonl), [one row per target](target-summary.json)
and [raw process receipts](runs/) retain all inputs, attempts, statuses and
timings.

| Target | Attempts | F4 online ms | F6 online ms | F5 online ms | F4 / F6 | F5 / F6 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | 3.536 | 2.740 | 1,152.728 | 1.291× | 420.777× |
| 2 | 3 | 46.119 | 31.030 | 10,867.124 | 1.486× | 350.214× |
| 3 | 4 | 85.498 | 56.067 | 17,592.913 | 1.525× | 313.784× |
| 4 | 4 | 72.933 | 48.072 | 16,241.401 | 1.517× | 337.857× |
| 5 | 2 | 25.504 | 16.678 | 8,947.016 | 1.529× | 536.449× |
| 6 | 1 | 9.608 | 6.656 | 3,251.246 | 1.443× | 488.455× |
| 7 | 11 | 224.314 | 141.685 | 56,898.256 | 1.583× | 401.582× |
| 8 | 2 | 26.428 | 20.398 | 5,873.096 | 1.296× | 287.921× |

F4 and F6 entries are medians of their two processes for that target; F5
has one process per target. Across these **eight fixed targets**, the
empirical F4/F6 ratio ranges **1.291–1.583×** with median **1.502×**;
F5/F6 ranges **287.921–536.449×** with median **375.898×**. F6 is
descriptively faster than both on all eight points. The F4 A/A pair spread
was 0.96–7.34%; the smallest F4/F6 gap is 29.1%. Target-dependent attempt
counts ranged from 1 to 11 and the corresponding F6 geometric work from
633 to 79,635 logical point additions. Target PDP occupied 98.97–99.96%
of F6's online interval and at least 99.40% of F4's.

This is an **observed range on the frozen sample**, not a mathematical
complexity bound or a population confidence interval. The host was an
unisolated macOS ARM64 machine, so CPU wall ratios remain exploratory even
though the F4 A/A spread is smaller than the observed gap. Matrix F5 uses
lowest-free splitting; F4/F6 use inherited-basis/highest-free splitting.
The result compares full solver choices, not only matrix kernels. All
targets used imported certified base logs, so this is prepared one-target
online cost, not cold relation collection. The preregistered 2× F4 gate
failed on all eight. The small-curve cold IC test is reported separately;
neither panel establishes an IC-versus-rho speedup or a gain on the n53
`PDP4root` path.

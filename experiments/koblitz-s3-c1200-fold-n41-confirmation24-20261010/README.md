# Fixed-base N41 C1200 target-kernel plus PMULL-fold panel

The improved C1200 target path with the opt-in ARM64 PMULL sparse-tail fold measured a **2.421× geometric paired single-target online ratio** against the earlier C1200/table path on 24 fresh N41 public points. All 48 native reports passed group, target, factor-base, rank, and five-phase checks; checked Sage independently replayed all 48 scalars. The target-bootstrap 95% interval was **[1.927×, 3.092×]**: its point estimate clears 2×, while the predeclared lower-bound gate remains open.

This comparison holds the mathematical base fixed. Both candidates use **98,400 distinct subgroup-usable factor-base points**, the same base digest, 1,200 prepared columns, full rank 1,200, and 134,217,728 root-table slots. The earlier source and table arithmetic are `IC1N41Ckb1fb98400PDP4rootRCsampleLAgaussTDdirectISO0heac310d59e1f`; the improved target source and ARM64 fold are `IC1N41Ckb1fb98400PDP4rootRCsampleLAgaussTDdirectISO0hc4ef08037c5a`. The source, binary, candidate, workload, Sage-runtime, and dispatch-probe hashes are fixed in [FREEZE.json](FREEZE.json).

Checked Sage generated T186–T209 using the parent panel's SHA-256 scalar law; T186–T200 were cross-checked against the parent's existing points, scalars, and workload IDs. [GENERATED_TARGETS.json](GENERATED_TARGETS.json) binds each point and one-target workload. Every pair uses one public point, the same relation seed, 14 preparation workers, one target worker, a 180-second process cap, and alternating arm order. The online interval starts after reusable preparation and ends after scalar recovery check; fixture construction, process launch, and point loading are outside it.

| Measurement | Earlier C1200/table | Improved C1200/fold |
| --- | ---: | ---: |
| Median online, ms | 0.162687 | 0.070355 |
| Median target query, ms | 0.042438 | 0.002188 |
| Median charged target PDP, ms | 0.047667 | 0.052708 |
| Median target relation check, ms | 0.010438 | 0.006583 |
| Median target descent, ms | 0.020083 | 0.008562 |
| Median scalar recovery check, ms | 0.030375 | 0.001750 |
| Median reusable setup, ms | 13,588.787 | 12,900.815 |
| Maximum peak RSS, GB (decimal) | 8.946 | 8.910 |

The target-call geometric ratio of earlier/improved was **0.684×**, so the improved arm made more S3 calls on average. Its online gain comes mainly from target query and scalar-check paths; this panel does **not** show a charged PDP-stage speedup. Phase medians are independent and need not sum to the median online interval. All 24 points favored the improved candidate online, and 12 reached at least 2× individually. Reusable setup and memory are close between arms and substantially below the earlier C1800/fold candidate's roughly 11.4-second setup and 18.7-GB peak RSS on a different local panel; setup comparisons across panels are exploratory because host contention changed.

This Apple M4 Pro local timing lacks a host-exclusive CPU partition and noise receipt, so `controlled_speedup` remains null. The raw reports and ledger are in [raw.tar.zst](local/raw.tar.zst) (SHA-256 `785a84df6355b78962f6755621650efe866fade719e84678fd736f23bf97dcaf`), with [analysis](local/analysis.json), [Sage replay](local/sage_replay.json), and [checked Sage runtime](local/sage_runtime_info.json). A missing-ledger analyzer check returned no ratio and zero verified rows. The next step is an independent 48-target confirmation with the same candidates and frozen gate.

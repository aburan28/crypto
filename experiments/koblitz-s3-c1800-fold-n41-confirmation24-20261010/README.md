# Fresh N41 C1800 plus PMULL-fold confirmation (24 public targets)

The C1800 S3 point-decomposition candidate with the opt-in ARM64 PMULL sparse-tail fold achieved a **2.419× geometric paired online ratio** against the C1200/table candidate on N41 targets T114–T137. The 95% target-bootstrap interval was **[1.949×, 2.945×]**. All 48 process reports passed the native correctness and phase audit, and checked Sage independently replayed all 48 recovered scalars. The point estimate clears the development gate of 2×; the predeclared lower-bound confirmation gate does not.

## Frozen comparison

| | C1200/table | C1800/PMULL fold |
| --- | ---: | ---: |
| Candidate ID | `IC1N41Ckb1fb98400PDP4rootRCsampleLAgaussTDdirectISO0heac310d59e1f` | `IC1N41Ckb1fb147600PDP4rootRCsampleLAgaussTDdirectISO0h85bb48f97c92` |
| Actual usable factor-base points | 98,400 | 147,600 |
| Prepared rank | 1,200 | 1,800 |
| Root-table slots | 134,217,728 | 268,435,456 |
| `GF2_ARM_PMULL_FOLD` | `0` | `1` |
| Binary SHA-256 | `0f805070787361a6dc250985ce64c4380b3b781aa7c374f99df66cfee5d530f4` | `d10625fe3aba0e7b61c44f001f119c1bb3b7f4fdde4ee558df1501f1f1ab1164` |

Both candidates use the same N41 Koblitz curve, subgroup, public point per pair, deterministic input law, relation seed, 14 relation workers, one target worker, and one public target per process. The order alternates by target. Source, binary, candidate, point, workload, runner, analyzer, and dispatch-probe hashes are in [FREEZE.json](FREEZE.json). The C1800 manifest records the ARM64 PMULL fold as an algorithm-affecting flag; its compact ID is derived from its canonical manifest. Source tree `HEAD:src` was `0a574800502e67882d8f8f88a846c1e80ee795fd`.

The online clock starts after reusable base, index, and log preparation, and covers target query, charged target PDP, relation check, target descent, and scalar recovery check. Fixture point creation, process launch, input loading, and reusable setup are outside this interval. The primary paired ratio uses all 24 targets; the 20,000-resample bootstrap samples target pairs, seeded as frozen. The raw run ledger preserves failures and resource counters, and the analyzer requires every planned run and report to be valid before computing a ratio.

## Result and cost boundary

| Measurement | C1200/table | C1800/PMULL fold |
| --- | ---: | ---: |
| Median online time, ms | 0.119021 | 0.043834 |
| Median target query, ms | 0.029937 | 0.001813 |
| Median charged target PDP, ms | 0.045853 | 0.028646 |
| Median target relation check, ms | 0.008813 | 0.005562 |
| Median target descent, ms | 0.009084 | 0.006479 |
| Median target recovery check, ms | 0.023479 | 0.001250 |
| Median reusable setup, ms | 4,252.172 | 10,642.280 |
| Maximum observed peak RSS, GB (decimal) | 8.877 | 18.707 |

Phase medians do not sum to the median online time because the medians can come from different targets. The 24-point geometric ratio of C1200/C1800 S3 calls was 1.267×; 22 of 24 targets favored C1800 on online time, and 17 reached at least 2× individually. The additional base and index increase reusable setup and peak memory. This study answers the frozen one-target online question; its setup data is recorded separately.

This was an Apple M4 Pro ARM64 local correctness and exploratory timing run. The host has no exclusive core partition, fixed-frequency, IRQ, or other host-level isolation receipt. `controlled_speedup` therefore remains null. The separate [PMULL-fold screen](../koblitz-s3-pmull-fold-screen-20261010/README.md) supplies the fold/table equivalence and dispatch controls. The pilot's 2× lower-bound gate remains open for an independent host-isolated run.

## Evidence and replay

- [Analysis](local/analysis.json), [independent Sage replay](local/sage_replay.json), [checked Sage runtime](local/sage_runtime_info.json), and [raw process reports](local/raw.tar.zst).
- Raw archive SHA-256: `7969d9dfc94c1c47b2d3f6d695643e4d8fe0fc83d400ca7599cd5f0340a2c477`.
- The missing-ledger analyzer check returned `complete: false`, verified runs `0`, and no ratio.

The [independent 48-target confirmation](../koblitz-s3-c1800-fold-n41-confirmation48-20261010/README.md) uses T138–T185 on the same two candidate IDs, with a separately frozen run order and gate. Its 95% bootstrap lower bound is 2.409×.

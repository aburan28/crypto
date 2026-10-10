# Independent 48-target N41 fixed-base confirmation

The improved C1200 target path with ARM64 PMULL sparse-tail folding measured a **2.020× geometric paired single-target online ratio** against the earlier C1200/table path on 48 fresh N41 public points. All **96 native runs** recovered and verified their scalars, and checked Sage independently replayed all 96. The target-bootstrap 95% interval was **[1.469×, 2.710×]**: the point estimate reaches the frozen 2× development threshold, while the prespecified lower-bound confirmation gate does not.

Both candidates use the same **98,400 distinct subgroup-usable factor-base points**, base digest, 1,200 prepared columns, full rank 1,200, and 134,217,728 root-table slots. The baseline candidate is `IC1N41Ckb1fb98400PDP4rootRCsampleLAgaussTDdirectISO0heac310d59e1f`; the improved candidate is `IC1N41Ckb1fb98400PDP4rootRCsampleLAgaussTDdirectISO0hc4ef08037c5a`. [FREEZE.json](FREEZE.json) binds the exact source, executable, manifests, generated points, one-target workloads, checked Sage runtime, ARM dispatch probe, runner, analyzer, and replay script. The source tree was `0a574800502e67882d8f8f88a846c1e80ee795fd`.

Checked Sage generated T210–T257 under the prior panel's SHA-256 scalar law, with unique public points and workload IDs. Each arm ran once per point in a fresh process, with alternating order, the same relation seed, 14 preparation workers, one target worker, a 180-second process cap, and a 32 GiB RSS envelope. The online clock starts after reusable base, index, and log preparation, and ends after scalar recovery check. Point loading, process launch, and fixture construction are outside that interval. The five target-dependent phases sum to the charged online total in every report.

| Measurement | Earlier C1200/table | Improved C1200/fold |
| --- | ---: | ---: |
| Median online, ms | 0.199229 | 0.091104 |
| Median target query, ms | 0.060875 | 0.002562 |
| Median charged target PDP, ms | 0.051854 | 0.067499 |
| Median target relation check, ms | 0.011750 | 0.008729 |
| Median target descent, ms | 0.024499 | 0.009916 |
| Median scalar recovery check, ms | 0.037687 | 0.002167 |
| Median reusable setup, seconds | 29.872 | 31.722 |
| Maximum peak RSS, GB (decimal) | 8.925 | 8.938 |

The paired median ratio was 2.393×; 41 of 48 points favored the improved candidate and 25 reached at least 2× individually. The earlier/improved geometric S3-call ratio was **0.737×**, so the improved arm made more calls and its charged PDP median increased. Its online advantage comes chiefly from the faster point-validation and scalar-replay paths. Phase medians are computed independently and need not sum to the median online time.

Host noise materially widened the interval. Reusable-setup coefficients of variation were 0.383 and 0.324, with median child involuntary context switches of 97,705 and 109,072 for baseline and improved processes. T247's improved online interval was 4.668 ms for 32 S3 calls, while T243's baseline interval was 1.538 ms for 72 calls. Across leave-one-target-out calculations, the geometric ratio ranged from **1.914× to 2.192×**. These are descriptive sensitivity results; every target remains in the primary 2.020× estimate. The Apple M4 Pro run lacks a host-exclusive CPU partition and noise receipt, so `controlled_speedup` remains null.

The earlier disjoint [24-target fixed-base panel](../koblitz-s3-c1200-fold-n41-confirmation24-20261010/README.md) measured 2.421× with a [1.927×, 3.092×] interval. Its targets are not reused here. The [C1800/fold 48-target panel](../koblitz-s3-c1800-fold-n41-confirmation48-20261010/README.md) measured 2.650× [2.409×, 2.905×] against the same earlier C1200 candidate definition, but uses a different factor base and target set. A host-isolated, paired target panel is the next measurement step for a controlled CPU result; a smaller initial prefetch batch is a concrete follow-up for reducing the extra S3 calls at fixed C1200 base size.

## Evidence

- [Analysis](local/analysis.json), [independent Sage scalar replay](local/sage_replay.json), [checked Sage runtime](local/sage_runtime_info.json), and [raw process reports and ledger](local/raw.tar.zst).
- Raw archive SHA-256: `3b205f46c3c13ae603e433980bc2f017306f7e074402a22e029d3f50afc07205`.
- The archived `raw.jsonl` SHA-256 is `048842826e38db703fac08393dde9d0a6f64cef2da55576aad706bdd8af6b911`, matching the analysis and Sage replay receipts. The analyzer audited 96/96 rows in the frozen execution order.

Regenerate analysis after extracting the archive with `python3 analyze.py --ledger /absolute/path/raw.jsonl --output /absolute/path/analysis.json`. Replay scalars with `/Volumes/SSD990/cryptanalysis/sage -python /absolute/path/replay_sage.py --ledger /absolute/path/raw.jsonl --runtime-info /absolute/path/local/sage_runtime_info.json --receipt /absolute/path/replay.json`.

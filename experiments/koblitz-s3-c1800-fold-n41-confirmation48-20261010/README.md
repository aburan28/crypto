# Independent 48-target N41 confirmation of C1800 plus ARM64 PMULL fold

On 48 previously unused N41 public points, the C1800 S3 point-decomposition candidate with the opt-in ARM64 PMULL sparse-tail fold measured a **2.650× geometric paired single-target online ratio** against C1200/table. The target-bootstrap 95% interval was **[2.409×, 2.905×]**, above the frozen 2× lower-bound gate. All 96 native process reports passed correctness and phase accounting, and checked Sage independently replayed all 96 recovered scalars. This is a local exploratory CPU result; the host-exclusive isolation gate remains open.

The earlier disjoint [24-target panel](../koblitz-s3-c1800-fold-n41-confirmation24-20261010/README.md) measured 2.419× with interval [1.949×, 2.945×]. This independent 48-target panel was frozen after that result, uses the exact same candidate IDs and binary hashes, and has its own predeclared 48-point stop rule and 20,000-resample target bootstrap.

## Candidate and workload

| | C1200/table | C1800/PMULL fold |
| --- | ---: | ---: |
| Candidate ID | `IC1N41Ckb1fb98400PDP4rootRCsampleLAgaussTDdirectISO0heac310d59e1f` | `IC1N41Ckb1fb147600PDP4rootRCsampleLAgaussTDdirectISO0h85bb48f97c92` |
| Actual usable factor-base points | 98,400 | 147,600 |
| Prepared rank | 1,200 | 1,800 |
| Root-table slots | 134,217,728 | 268,435,456 |
| `GF2_ARM_PMULL_FOLD` | `0` | `1` |
| Binary SHA-256 | `0f805070787361a6dc250985ce64c4380b3b781aa7c374f99df66cfee5d530f4` | `d10625fe3aba0e7b61c44f001f119c1bb3b7f4fdde4ee558df1501f1f1ab1164` |

The field, curve, subgroup, target point in each pair, target input law, relation seed, 14 relation workers, one target worker, timeout, and process resource envelope are frozen in [FREEZE.json](FREEZE.json). The points are T138–T185 from the parent public-point fixture, each with a unique canonical workload ID. Execution order alternates C1200-first and C1800-first. The C1800 candidate manifest, created in the 24-target panel, records the fold and required ARM64 PMULL capability. A dispatch probe verifies that the two environment values select the intended arithmetic path. Source tree `HEAD:src` was `0a574800502e67882d8f8f88a846c1e80ee795fd`.

## Measured costs

| Measurement | C1200/table | C1800/PMULL fold |
| --- | ---: | ---: |
| Median online time, ms | 0.119855 | 0.043938 |
| Median target query, ms | 0.035875 | 0.002041 |
| Median charged target PDP, ms | 0.036083 | 0.027167 |
| Median target relation check, ms | 0.009395 | 0.006583 |
| Median target descent, ms | 0.010396 | 0.005793 |
| Median target recovery check, ms | 0.029312 | 0.001604 |
| Median reusable setup, ms | 4,982.228 | 11,426.584 |
| Maximum peak RSS, GB (decimal) | 8.896 | 18.698 |

Each online interval starts after reusable base, index, and log preparation, and ends after scalar recovery check. All target-dependent attempts are charged. Fixture creation, process launch, point loading, and reusable setup are outside the online interval. Medians of individual phases do not sum to the median total because they can come from different targets. The point-wise S3-call geometric ratio was 1.200×; 47 of 48 targets favored C1800 on online time, and 42 of 48 achieved at least 2× individually. The larger factor base and index increase reusable setup, memory, and cold-start cost. The five online phases, setup stages, row status, and resource counters are in the raw reports.

The run used a physical Apple M4 Pro ARM64 host. Its ordinary local CPU partition did not provide an auditable exclusive-core, fixed-frequency, IRQ-routing, and contention receipt. `controlled_speedup` remains null under the repository's CPU isolation gate, even though the independent target bootstrap clears its exploratory development gate. The ARM64 fold is opt-in and does not imply an x86-64 benefit; the separate [fold screen](../koblitz-s3-pmull-fold-screen-20261010/README.md) isolates its effect with the same binary on each curve.

## Evidence

- [Full analysis](local/analysis.json), [Sage scalar replay](local/sage_replay.json), [checked Sage runtime](local/sage_runtime_info.json), [raw process reports and ledger](local/raw.tar.zst).
- Raw archive SHA-256: `aa2931cce2f24e6bff33fe7c144d1d9c181cc1daa61147171ccc7c94770b6aa6`.
- The analyzer's missing-ledger check returned `complete: false`, verified runs `0`, and no online ratio.

To regenerate the analysis after extracting the raw archive, run `python3 analyze.py --ledger /path/to/raw.jsonl --output /path/to/analysis.json`. For scalar replay, use `/Volumes/SSD990/cryptanalysis/sage -python /absolute/path/to/replay_sage.py --ledger /path/to/raw.jsonl --runtime-info /absolute/path/to/local/sage_runtime_info.json --receipt /path/to/replay.json`.

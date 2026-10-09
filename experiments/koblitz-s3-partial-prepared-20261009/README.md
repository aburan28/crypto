# Fixed partial-rank S3 preparation with online target-span fallback

A fixed, target-independent relation prefix prepares rank 243 on N41 and rank
237 on N53 before the public target point is read. The solver then decomposes
one target and checks whether its coefficient row lies in the prepared span.
When it does not, the solver continues the same unique relation-scalar stream
and charges every extra batch to that target's online interval. On the first
24 fresh points per curve, all 96 paired arm runs recovered and verified their
scalars. The partial and reference arms used identical relation prefixes and
the same total number of cold S3 root calls for every paired point.

## Exact instances and identities

The curves use `y² + xy = x³ + 1`. Field elements and affine coordinates are
unsigned polynomial-basis integers. `G` and the subgroup order are in each
curve manifest. The base is the natural low-`x` scan, with distinct
subgroup-usable points folded into signed Frobenius orbits. `fb` is the actual
unfolded point count, while the relation matrix has 244 columns.

| Curve ID | Field polynomial | Subgroup order `r` | Cofactor | Actual `fb` | Folded columns |
| --- | --- | ---: | ---: | ---: | ---: |
| `EC1N41Ckb1hcbf4a37f6e66` | `x^41 + x^3 + 1` | 549756390943 | 4 | 20,008 | 244 |
| `EC1N53Ckb1h738a572b91f2` | `x^53 + x^6 + x^2 + x + 1` | 21044858204113 | 428 | 25,864 | 244 |

| Curve | Frozen S3 reference candidate ID | Partial-prepared candidate ID |
| --- | --- | --- |
| N41 | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h3f64f5604376` | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h0f38cf1ef60a` |
| N53 | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h5ed3e1d8954b` | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h1b2989c714f4` |

The canonical curve and candidate manifests are in [`identities/`](identities/).
Its 48 workload manifests fix one target point apiece; each has a reference
and partial run ID in [`measurements.jsonl`](identities/measurements.jsonl).
The input law is
`1 + int_be(SHA256("s3-pair-target-v1:" || seed_label || ":" || u32be(n) || u32be(index))) mod (r-1)`,
with seed label `s3-partial-prepared-holdout-20261009-v1`. The scalar constructs
the fixture; the solver receives only its public point. The target files and
all 100 pre-generated points per curve are in [`holdout/`](holdout/).

## First 24-point-per-curve panel

[`FREEZE.json`](FREEZE.json) fixed the first 24 targets on each curve, rank
thresholds, source and binary hashes, relation seed `20260928`, 14 workers,
alternating curve/arm order, a fresh process for each arm and point, and a
180-second limit before these target decompositions were measured. The Rust
library is from frozen source commit `833ef86191947ba8e4def56a78ee6c51a4ea0ff2`.
All 96 runs finished successfully; failures would have remained rows in the
append-only ledger. Source and workload hashes are checked by the runner.

| Field | Fallback targets | Reference online median | Partial online median | Geometric mean of paired online ratios | Smallest paired online ratio | Reference / partial cold median |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| N41 | 0 / 24 | 129.499 ms | 0.577 ms | 255.22× | 26.98× | 417.731 / 397.060 ms |
| N53 | 2 / 24 | 2999.206 ms | 13.088 ms | 205.38× | 2.26× | 3524.377 / 3705.381 ms |

The paired target-bootstrap 95% descriptive intervals for the geometric mean
online ratio are 167.96–381.74× on N41 and 92.41–451.37× on N53. These are
**exploratory local wall measurements** on a shared ARM64 Mac; the machine did
not provide an exclusive host CPU partition or the required run-level
frequency, IRQ, throttling, steal, and pressure receipt. The machine-readable
analysis keeps `controlled_speedup: null`. Tiny target intervals are sensitive
to scheduler interruptions. Per-target timings, failures, online phase costs,
and resource fields remain in the raw and measurement rows.

The deterministic operation boundary gives a separate result. The smallest
reference/partial ratio of **online S3 root calls** was 68.45× on N41 and
2.0096× on N53. All 48 pairs had the same total cold S3 root-call count in
both arms. N41 used 252 relation attempts to reach prepared rank 243 and
required no fallback. N53 used 238 attempts to reach rank 237; T003 and T005
used 224 and 98 further attempts online, respectively. Their final ranks
were 243 and 242. The Wilson 95% interval for fallback-free coverage is
86.2–100% on N41 and 74.2–97.7% on N53; the interval reflects 24 targets per
curve, rather than a guaranteed per-target coverage rate.

| Field | Base construction median, ref / partial | S3 index build median, ref / partial | Rank gain per relation attempt, ref / partial prep | Relation collection ms per new rank row, ref / partial prep |
| --- | ---: | ---: | ---: | ---: |
| N41 | 1.971 / 2.153 ms | 304.175 / 261.653 ms | 0.9643 / 0.9643 | 0.530 / 0.420 |
| N53 | 5.907 / 5.325 ms | 549.749 / 503.683 ms | 0.9958 / 0.9958 | 11.886 / 12.576 |

All ordinary relation attempts in this panel produced verified witnesses;
dependent rows remain in the raw trace and are included in the costs. On the
two N53 fallbacks, the observed collection cost per new rank row was
1669.745 ms for T003 and 171.140 ms for T005. Those costs include the many
dependent relation attempts needed to span the particular target rows.

The online timer starts at target validation after candidate-specific reusable
setup and stops after scalar replay. Target input loading and fixture creation
are outside it. The five exclusive costs are `target_query`,
`target_pdp_charged`, `target_relation_check`, `target_descent`, and
`target_recovery_check`; their sum is checked on all 96 runs. The reference
charges target-span relation collection online. Partial preparation reaches
the fixed rank without looking at the target, and any subsequent relation
collection is charged online. The supplementary cold timer includes setup,
index construction, relation collection, target work, and replay after launch.
The raw output preserves attempt counts, relation status, rank, S3 calls,
individual stage timings, and all available memory fields; macOS reported
`peak_rss_bytes: null` for these processes.

## Arithmetic and raw evidence

[`raw_panel.tar.gz`](local/raw_panel.tar.gz) archives the exact 96 solver JSONL
outputs, stdout/stderr, and the append-only runner ledger. Its SHA-256 is
`88d726f4508d6083311115324c311cae5a1f3c4af7ab46a3298da690ee5f3e79`.
[`analysis.json`](local/analysis.json) retains each pair and the rate and
timing summaries. The checked Sage replay in
[`independent_sage_replay.json`](local/independent_sage_replay.json)
independently reconstructed all 48 target scalars across 96 arm runs, checked
714 distinct relation witnesses as group-point equations, checked each target
decomposition, and rechecked 1,071 N41 and 1,878 N53 used base labels. The
earlier [`prior_base_replay.json`](../koblitz-s3-prepared-20261009/prior_base_replay.json)
checks every point and label in the 20,008- and 25,864-point bases. The
[`sage_runtime_info.json`](local/sage_runtime_info.json) receipt verifies the
accepted Sage installation used for this replay.

## Untouched holdout extension

[`EXTENSION_FREEZE.json`](EXTENSION_FREEZE.json) fixes all remaining T025–T100
targets per curve with the same binaries, rank policy, limits, and pairing.
This extension was frozen after the first 24-point-per-curve analysis and is
reported separately. It tests fallback frequency and the lower tail of paired
online ratios without changing either candidate. All 304 extension arm runs
verified, with the same cold S3 root-call count and common relation prefix in
every pair.

| Field | Fallback targets | Reference online median | Partial online median | Geometric mean of paired online ratios | Smallest paired online ratio | Reference / partial cold median |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| N41 | 0 / 76 | 103.125 ms | 0.652 ms | 188.72× | 7.96× | 364.629 / 373.209 ms |
| N53 | 10 / 76 | 2170.454 ms | 12.785 ms | 115.84× | 1.58× | 2576.805 / 2530.183 ms |

The N53 T085 and T096 pairs had exploratory online wall ratios below 2×;
each required full rank 244 and 322 target-charged fallback attempts after
the rank-237 prefix. The other 74 N53 extension pairs and all 76 N41 pairs
were above 2×. [`analysis_extension.json`](local/analysis_extension.json)
retains all pairs and uncertainty intervals; [`raw_extension.tar.gz`](local/raw_extension.tar.gz)
preserves the 304 solver outputs and raw ledger, SHA-256
`60fec46ea00db4a2ad7c8a7e76685b305855c75a1568b4865f2d9e33481c4e37`.
The extension's 152 workload manifests and 304 run rows are in
[`identities_extension/`](identities_extension/).
[`independent_sage_replay_extension.json`](local/independent_sage_replay_extension.json)
reconstructed all 152 extension target scalars across both arms and checked
812 distinct relation witnesses as group-point equations.

Across all 100 distinct points per curve at the original thresholds,
N41 had zero fallbacks and every exploratory online ratio exceeded 2×.
N53 had 12 fallbacks and two below-2× local wall ratios. The descriptive
geometric means are 202.90× for N41 and 132.91× for N53.
[`curve_comparison.json`](local/curve_comparison.json) records the
100-target summaries and a 1.53 ratio of their geometric means, with a
target-bootstrap 95% interval of 1.03–2.29. The calculation uses one
measurement per arm and target on a contended host. Under its observed log
ratio variances and effect size, a two-sided normal approximation gives
185 independent target pairs **per curve** for 80% power to distinguish that
observed curve difference at 5% significance. An isolated-host A/A pilot and
repeated one-target pairs are needed to recalibrate this post-hoc count for a
controlled CPU comparison.

## Rank-242 and rank-243 threshold iteration on N53

Rank 242 is reached after 336 target-independent relation attempts; rank 243
after 462. The fixed relation stream reaches full rank 244 after 560 attempts
on the three selected full-rank stress targets. Each subsequent rank was
frozen before its **new** 24-point N53 seed was evaluated. The earlier full-rank
T067, T085, and T096 points were reused as explicitly labeled stress controls.
Each panel pairs its reference with both indicated fixed-rank policies in a
rotating three-arm order, with one fresh process and one public target per arm.

| N53 policy | Fresh targets with verified recovery | Fresh minimum online ratio | Fresh online median, ref / policy | Fresh cold median, ref / policy | Full-rank stress minimum online ratio | Full-rank stress minimum online S3 root-call ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Rank 242 | 24 / 24 | 18.99× | 1737.997 / 18.639 ms | 2167.481 / 2773.767 ms | 1.85× | 2.31× |
| Rank 243 | 24 / 24 | 10.10× | 1759.185 / 11.176 ms | 2112.343 / 3613.809 ms | 4.20× | 5.23× |

All 81 runs in each panel verified. Rank 242 needed 224 online fallback
attempts on each full-rank stress point; rank 243 needed 98. The rank-242
S067 stress pair measured 1.85× online wall time even though its S3 root-call
ratio was 2.32×, illustrating the shared host's timing variation. All three
rank-243 stress pairs measured 4.20–5.90× online. The rank-243 fresh targets
required no fallback; their reference/policy geometric mean online ratio was
190.54×, with a target-bootstrap 95% descriptive interval of
112.04–319.67×. On fresh points that the reference spans before rank 243,
the higher threshold performs more cold relation work. The cold medians above
quantify that cost.

The N53 rank-242 candidate ID is
`IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hd13676aebfef`;
the rank-243 ID is
`IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hbd278c9cf111`.
Their complete manifests, single-target workload IDs, and run rows are in
[`rank242/identities/`](rank242/identities/) and
[`rank243/identities/`](rank243/identities/). The frozen input plans are
[`rank242/FREEZE.json`](rank242/FREEZE.json) and
[`rank243/FREEZE.json`](rank243/FREEZE.json).
[`analysis_rank242.json`](local/analysis_rank242.json) and
[`analysis_rank243.json`](local/analysis_rank243.json) keep every fresh and
stress pair separately. Their raw archives are
[`raw_rank242.tar.gz`](local/raw_rank242.tar.gz), SHA-256
`cf7a3ae428d199d248911f1d48c9fb9465462b1577a78c5e10ce92cf78e366a1`,
and [`raw_rank243.tar.gz`](local/raw_rank243.tar.gz), SHA-256
`62ce737b1d1e3c37bfcba63a0403d8fc7ef7c60240ae674a6e962b229ac3770b`.
The checked Sage receipts
[`independent_sage_replay_rank242.json`](local/independent_sage_replay_rank242.json)
and [`independent_sage_replay_rank243.json`](local/independent_sage_replay_rank243.json)
independently reconstructed every scalar across the 81 runs in each panel.
Each replay verified 560 distinct relation witnesses as point equations and
checked the target decompositions and used factor-base labels. The rank-243
replay checked 2,251 distinct used base labels.

## Reproduction

From this checkout, build the frozen examples with an out-of-tree Cargo target
directory, then run the paired panel. A checked Sage runtime receipt must be
saved before any Sage replay:

```sh
CARGO_TARGET_DIR=/private/tmp/s3-partial-build cargo build --release --offline --example s3_pair_candidate_frozen --example s3_partial_prepared_targetspan
python3 experiments/koblitz-s3-partial-prepared-20261009/run_paired_holdout.py --binaries /private/tmp/s3-partial-build/release/examples --output /private/tmp/s3-partial-rerun
python3 experiments/koblitz-s3-partial-prepared-20261009/make_identities.py --ledger /private/tmp/s3-partial-rerun/raw.jsonl --output /private/tmp/s3-partial-rerun-identities
tar -czf /private/tmp/s3-partial-rerun.tar.gz -C /private/tmp/s3-partial-rerun .
python3 experiments/koblitz-s3-partial-prepared-20261009/analyze_holdout.py --ledger /private/tmp/s3-partial-rerun/raw.jsonl --archive /private/tmp/s3-partial-rerun.tar.gz --output /private/tmp/s3-partial-rerun-analysis.json
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage --runtime-info > /private/tmp/s3-partial-rerun-sage-runtime-info.json
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage -python experiments/koblitz-s3-partial-prepared-20261009/replay_holdout_sage.py --ledger /private/tmp/s3-partial-rerun/raw.jsonl --runtime-info /private/tmp/s3-partial-rerun-sage-runtime-info.json --receipt /private/tmp/s3-partial-rerun-sage-replay.json
```

The extension and later thresholds use their own frozen plans and outputs:

```sh
python3 experiments/koblitz-s3-partial-prepared-20261009/run_extension.py --binaries /private/tmp/s3-partial-build/release/examples --output /private/tmp/s3-partial-extension-rerun
python3 experiments/koblitz-s3-partial-prepared-20261009/run_rank242.py --plan experiments/koblitz-s3-partial-prepared-20261009/rank242/FREEZE.json --binaries /private/tmp/s3-partial-build/release/examples --output /private/tmp/s3-rank242-rerun
python3 experiments/koblitz-s3-partial-prepared-20261009/run_rank242.py --plan experiments/koblitz-s3-partial-prepared-20261009/rank243/FREEZE.json --binaries /private/tmp/s3-partial-build/release/examples --output /private/tmp/s3-rank243-rerun
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage -python experiments/koblitz-s3-partial-prepared-20261009/replay_holdout_sage.py --ledger /private/tmp/s3-partial-extension-rerun/raw.jsonl --runtime-info /private/tmp/s3-partial-rerun-sage-runtime-info.json --receipt /private/tmp/s3-partial-extension-sage-replay.json --first-target 25 --target-count-per-curve 76
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage -python experiments/koblitz-s3-partial-prepared-20261009/replay_rank242_sage.py --plan experiments/koblitz-s3-partial-prepared-20261009/rank242/FREEZE.json --ledger /private/tmp/s3-rank242-rerun/raw.jsonl --runtime-info /private/tmp/s3-partial-rerun-sage-runtime-info.json --receipt /private/tmp/s3-rank242-sage-replay.json
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage -python experiments/koblitz-s3-partial-prepared-20261009/replay_rank242_sage.py --plan experiments/koblitz-s3-partial-prepared-20261009/rank243/FREEZE.json --ledger /private/tmp/s3-rank243-rerun/raw.jsonl --runtime-info /private/tmp/s3-partial-rerun-sage-runtime-info.json --receipt /private/tmp/s3-rank243-sage-replay.json
```

## Measurement fidelity and next run

The evidence gate for these panels records frozen public points and binary
hashes, alternating/rotating arm order, one process and one target per arm,
verified recovered scalars, independent Sage equations, and an exclusive
five-phase online sum for every successful run. The local host did not supply
a host-wide isolated CPU partition, fixed-frequency/NUMA policy receipt,
run-level PSI, IRQ/steal/throttle counters, or an A/A noise control. Thus each
wall ratio is labeled exploratory and `controlled_speedup` remains `null` in
the analysis records. The deterministic S3 root-call counts and verified
relation stream do not depend on that host timing gate.

The next performance step is a paired run on an auditable isolated Linux host
with a quiet-host A/A control, CPU affinity, frequency and NUMA records, PSI,
interrupts, steal/throttle counters, and the same target-point hashes. The
record should retain `perf` instruction and cycle counts, including IPC, when
the host permits them. A single-target rho run on the same points and resource
envelope can then test the primary IC/rho online comparison.

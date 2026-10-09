# Paired one-target S3 and signed-Frobenius rho pilot

The rank-243 S3 path recovered **all 48 fresh public points** in a paired
N41/N53 pilot. Its exploratory one-target online wall ratio against the
signed-Frobenius rho reference was **27.35× on N41** and **5.69× on N53** as a
geometric mean over 24 points per curve. Checked Sage independently replayed
all **192** recovered scalars. The pointwise 2× threshold held on all 24 N41
points and 18 of 24 N53 points.

| Curve | Curve ID | Verified runs | Median S3 / rho online ms | Rho / S3 geometric mean (95% target bootstrap) | Pointwise ratios ≥2× | Minimum pointwise ratio | S3 / rho A/A p95 absolute log ratio |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| N41 | `EC1N41Ckb1hcbf4a37f6e66` | 96/96 | 0.513 / 14.663 | 27.35× (18.31–41.09×) | 24/24 | 4.37× | 0.849 / 0.512 |
| N53 | `EC1N53Ckb1h738a572b91f2` | 96/96 | 14.400 / 76.828 | 5.69× (3.69–8.98×) | 18/24 | 1.36× | 1.145 / 0.898 |

The ratios pair **the same public point**, with two solves per arm in rotated
ABBA/BAAB order. Each point has a frozen one-target workload ID in the
[crosscurve plan](../koblitz-s3-crosscurve-20261009/FREEZE.json). S3 uses
`IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h0f38cf1ef60a` or
`IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hbd278c9cf111`.
The [experiment freeze](FREEZE.json), SHA-256
`104279d2b0b9054ee57da4027ccb4e3c9a647680eac5bb03b0c799c142f0aa62`,
pins the signed-Frobenius rho source and local binary hashes, jump seed,
32 lockstep lanes, 8 distinguished-point bits, and the pair order. Every
rho run uses an empty distinguished-point table and one public target.

The S3 online clock begins after its reusable rank-243 preparation and ends
after scalar replay; its five exclusive phases are preserved in each raw
report. Rho's online clock begins after public-point validation and jump
preparation, at its first target-dependent walk computation, and ends after
`[d]G = Q` succeeds. Its walk, collision, and recovery-check nanoseconds sum
exactly to its online interval. The fixture scalar is constructed by checked
Sage outside both intervals and is never supplied to either solver. All 96
S3 runs reached prepared rank 243 with zero online fallback attempts. The S3
target phase therefore used one thread; its 14 workers were active in
reusable preparation. Rho used one thread with 32 lockstep arithmetic lanes.
Median S3 reusable setup was 280.65 ms on N41 and 2904.81 ms on N53; these
costs are retained separately from the primary online measurement.

The local host does not have a host-wide CPU isolation receipt. Its 10% A/A
noise gate is `abs(log(repeat1/repeat2)) ≤ log(1.10) = 0.0953`; all four
curve/arm p95 values exceed that limit. The table's wall ratios are therefore
exploratory, and `controlled_speedup` remains `null` in the
[analysis](local/analysis_pilot24.json), SHA-256
`44bfe0789238f8262416fa3f74f8bce0bbbf322db313ac2438d96292f5e29b48`.
The pointwise 2× Wilson intervals, operation counts, failures, binary hashes,
and run IDs are retained there and in the
[raw archive](local/raw_pilot24.tar.gz), SHA-256
`004374420ac38def52152b6ba232fee98fb7122588d9a8eca9f4cf58f278665b`.
The [Sage replay](local/sage_replay_pilot24.json), SHA-256
`4d9415084b9702095020a197701939d1584a2541b1525726b208e55b19472748`,
includes the checked [runtime receipt](local/sage_runtime_info.json).

Reproduce the pilot through the frozen runner:

```sh
CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009 cargo build --release --offline --example s3_partial_prepared_targetspan --example koblitz_rho_batch_ks_strong_online
python3 experiments/koblitz-s3-one-target-rho-20261009/run_paired.py --first 1 --count 24 --binaries /private/tmp/s3-prepared-target-20261009/release/examples --output /private/tmp/s3-rho-paired-replay
python3 experiments/koblitz-s3-one-target-rho-20261009/analyze_pilot.py --ledger /private/tmp/s3-rho-paired-replay/raw.jsonl --output /private/tmp/s3-rho-paired-replay/analysis.json
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage --runtime-info > /private/tmp/s3-rho-paired-replay/sage_runtime_info.json
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage -python experiments/koblitz-s3-one-target-rho-20261009/replay_sage.py --ledger /private/tmp/s3-rho-paired-replay/raw.jsonl --runtime-info /private/tmp/s3-rho-paired-replay/sage_runtime_info.json --receipt /private/tmp/s3-rho-paired-replay/sage_replay.json
```

The earlier [Modal crosscurve panel](../koblitz-s3-crosscurve-20261009/modal_crosscurve_panel.py)
is prepared to remeasure the S3 rank-243 versus frozen S3 reference on these
same points with guest topology, NUMA, PSI, cgroup, IRQ, steal, frequency, and
optional `perf` telemetry. Source and public-point transfer to Modal awaits
the explicit approval required after automatic review. A controlled CPU wall
ratio also requires an auditable host-level isolation record.

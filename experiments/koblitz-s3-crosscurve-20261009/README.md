# Matched fresh N41/N53 rank-243 S3 panel

The rank-243 S3 path recovered all **48 fresh public points** in the local
24-pair pilot, with online root-call geometric mean ratios of **495.31× on
N41** and **322.31× on N53** against the frozen references. Checked Sage
independently recovered every scalar and checked the relation equations.
Local wall ratios remain exploratory because the prepared-arm A/A noise
exceeded the 10% gate.

The experiment freezes **200 previously unused public target points per
curve** under one deterministic seed law. Each point has its own one-target
workload ID. The N41 and N53 rank-243 candidates are compared with their
frozen S3 references using two solves per arm and target in alternating
ABBA/BAAB order. The first 24 index pairs form a CPU pilot; the remaining
frozen pairs are available for a sample-size decision made from the pilot's
A/A noise and between-target variance.

| Field | Curve ID | Subgroup order | Actual base points | Rank-243 candidate ID | Reference candidate ID |
| --- | --- | ---: | ---: | --- | --- |
| N41 | `EC1N41Ckb1hcbf4a37f6e66` | 549756390943 | 20,008 | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h0f38cf1ef60a` | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h3f64f5604376` |
| N53 | `EC1N53Ckb1h738a572b91f2` | 21044858204113 | 25,864 | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hbd278c9cf111` | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h5ed3e1d8954b` |

The exact field models, generators, base construction, solver source hashes,
and stage policies are in the [existing candidate manifests](../koblitz-s3-partial-prepared-20261009/README.md).
The new [`FREEZE.json`](FREEZE.json) has SHA-256
`45cbbe93b2b0e541bce89973a97571d128293c786026ee3ada5cb346d0a48283`.
It fixes the curve records, target files and hashes, candidate IDs, source
hashes, 14 workers, relation seed `20260928`, rank threshold 243, and 180-second
per-run limit.

## Public target construction

For each curve degree `n` and index `i=1..200`, the fixture scalar is

`1 + int_be(SHA256("s3-pair-target-v1:" || seed_label || ":" || u32be(n) || u32be(i))) mod (r-1)`,

where `seed_label` is `s3-crosscurve-rank243-20261009-v1`. Checked Sage
multiplies that scalar by the exact subgroup generator to produce the public
point supplied to each solver. Fixture construction and its independent
scalar check are outside the online timer. Each solver invocation reads only
its [`public_target.json`](holdout/n41/T001/public_target.json) point file.
The fixture records and 400 immutable workload records are in [`holdout/`](holdout/)
and [`workloads/`](workloads/). The checked Sage startup receipt is
[`sage_runtime_info.json`](sage_runtime_info.json).

The [generator](freeze_targets_sage.py) writes immutable files and checks the
curve equation, subgroup order, point encoding, source hashes, and uniqueness
of the 200 fixture scalars and points on each curve. None of the 400 new
points overlaps the previously examined 100 N41 or 148 N53 points. Re-running
the generator refuses to alter any frozen artifact.

## Paired measurements

The [local runner](run_local_pilot.py) retains every invocation and uses the
candidate and workload IDs to form run IDs. It verifies the recovered scalar,
public point, actual base count, rank policy, process output, and the exact
five-phase online sum. Each ABBA/BAAB block has two reference and two
rank-243 solves of **the same single target point**. This supplies an A/A
diagnostic and a paired per-target ratio while preserving raw failures. The
online interval begins after reusable setup and ends after scalar replay, as
implemented by the frozen Rust solvers. Cold setup remains a separate cost.

This local command is for correctness and exploratory timing; an auditable
host-level CPU isolation receipt is required for a controlled wall-time ratio:

```sh
CARGO_TARGET_DIR=/private/tmp/s3-crosscurve-build cargo build --release --offline --example s3_pair_candidate_frozen --example s3_partial_prepared_targetspan
python3 experiments/koblitz-s3-crosscurve-20261009/run_local_pilot.py --first 1 --count 1 --binaries /private/tmp/s3-crosscurve-build/release/examples --output /private/tmp/s3-crosscurve-local-smoke
```

The CPU analysis reports each curve's verified paired online ratios, the
number above 2× with a Wilson interval, A/A noise, target-bootstrap
uncertainty, and the sample size for a curve-difference test using the
observed variance. Target-dependent fallback attempts stay in the online
interval. Guest fidelity uses the 14 checks already specified by the
[rank-243 fidelity analyzer](../koblitz-s3-partial-prepared-20261009/analyze_modal_fidelity.py).

## Local one-pair correctness smoke

The first matched index pair was solved twice per arm on each curve with the
previously frozen binaries (reference SHA-256 `c4ff99ace090551747b2a1083b9546754e0ce1814f309a47e46d2ab0bd132ca7`,
rank-243 SHA-256 `50d3ad29fbd28c146ed60cbd6a255638fc6d48733ce08bd8097b08e03acbd195`).
All eight runs recovered and verified their fixture scalars. The rank-243
candidate reached prepared rank 243 and needed zero online fallback attempts
on both points. The exact operation boundary was:

| Curve | Reference online S3 root evaluations | Rank-243 online S3 root evaluations | Reference / rank-243 |
| --- | ---: | ---: | ---: |
| N41 T001 | 351,266 | 2,164 | 162.32× |
| N53 T001 | 7,758,984 | 59,918 | 129.49× |

The [raw smoke archive](local/raw_smoke.tar.gz), SHA-256
`6cdb6e02a680cf48690e91c2c06b72be586f24e2a88623ffbb833e75db0481c1`,
preserves every solver output and timing phase. Checked Sage independently
reconstructed the eight scalars and verified 714 distinct relation witnesses
as group equations; its [replay receipt](local/sage_replay_smoke.json) has
SHA-256 `1d21119af44fe3b7c1b9e9cb900d7c3859804fe77abf151cd73f66fe8b2a548e`.
This smoke ran on a contended local host while a separate build was active,
so its wall times remain correctness diagnostics. The operation counts are
deterministic under the frozen source and target points.

## Local 24-pair pilot

All **192** invocations on the first 24 fresh index pairs passed native
scalar verification. Checked Sage independently recovered every scalar from
the recorded relation system and checked **714 distinct relation witnesses**
as group equations. The rank-243 solver reached rank 243 on every point and
used **zero online fallback attempts**. The [local analyzer](analyze_local_pilot.py)
uses `rank_target_s3_calls + target_s3_calls` for the reference online root
count, and `fallback_s3_calls + target_s3_calls` for rank-243; its [summary](local/analysis_pilot24.json)
retains both repeats for each public point.

| Curve | Verified points | Online root calls, reference / rank-243 | Geometric mean root-call ratio | Exploratory online wall ratio, geometric mean (minimum) | Prepared A/A p95 absolute log ratio | Median cold time, reference / rank-243 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| N41 | 24/24 | 16,821,976 / 65,080 | 495.31× | 222.51× (28.81×) | 0.924 | 301.19 / 286.96 ms |
| N53 | 24/24 | 384,039,440 / 1,760,472 | 322.31× | 112.82× (25.20×) | 0.738 | 1792.32 / 3230.72 ms |

All 24 paired online wall ratios on each curve exceeded 2×, but these Mac
timings remain exploratory. The 10% A/A gate is `abs(log(repeat1/repeat2))
<= log(1.10) = 0.0953`; both measured p95 values exceed it. The N53
rank-243 path spent a median **3213.40 ms** in reusable preparation versus
**330.17 ms** for the reference. These setup costs are outside the primary
one-target online interval and appear separately in the table's cold-time
column. A Modal pilot with guest noise telemetry and a host-wide isolation
receipt is required to promote a controlled wall-time ratio.

The 24 local pairs give a **provisional 33 matched pairs per curve** for a
two-sided 5% curve-difference test at 80% power: the N41 to N53 ratio of
geometric mean online ratios is 1.97 and the sample standard deviation of
paired log differences is 1.39. The same formula returns the 24-pair floor
for each curve's 2× threshold. These counts use noisy local wall times; the
Modal 24-pair pilot must recompute them before scheduling further frozen
indices.

The [raw pilot archive](local/raw_pilot24.tar.gz) has SHA-256
`04941cee38c0348e9db5d4147dc91708acfca06a9594950602bb23bfca30a35d`.
It retains all raw JSON, stdout, stderr, ledger, and timing phases. The
[independent Sage replay receipt](local/sage_replay_pilot24.json) has SHA-256
`8c87f7e9174410b4fb1514683e00ff39fe2850bb2547af955da95dfcb84a2930`.

Reproduce the independent replay through the checked Sage launcher after
saving `--runtime-info` beside the run:

```sh
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage --runtime-info > /private/tmp/s3-crosscurve-sage-runtime-info.json
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage -python experiments/koblitz-s3-crosscurve-20261009/replay_sage.py --ledger /private/tmp/s3-crosscurve-local-smoke/raw.jsonl --runtime-info /private/tmp/s3-crosscurve-sage-runtime-info.json --receipt /private/tmp/s3-crosscurve-local-smoke/sage_replay.json
```

## Prepared Modal CPU panel

[`modal_crosscurve_panel.py`](modal_crosscurve_panel.py) runs 1–24 frozen index
pairs per VM. `--first 1 --count 1` is an eight-run correctness and topology
smoke; `--first 1 --count 24` is the 192-run paired pilot. Further disjoint
blocks can use the already frozen indices 25–200. Each block keeps two
reference and two rank-243 one-target runs, the exact candidate/workload/run
IDs, all raw reports and failures, CPU affinity and SMT topology, NUMA bind,
PSI, cgroup quota/throttle/OOM counters, IRQ and steal deltas, frequency, and
optional whole-process `perf` cycles and instructions. The VM requests 16
physical cores, with 15 benchmark cores and a separate monitoring core.
Guest diagnostics remain distinct from the host-wide isolation receipt.

[`analyze_modal_crosscurve.py`](analyze_modal_crosscurve.py) checks source,
plan, candidate, workload, run, target, and point hashes across disjoint
batch outputs. It recomputes every paired online ratio, retains failed runs,
reports per-curve guest score and A/A noise, and bootstraps the matched-index
log-ratio difference. For index `i`, this difference is
`d_i = log(reference_41 / rank243_41) - log(reference_53 / rank243_53)`.
After at least four verified matched pairs, the pilot sizes a two-sided 5%
curve-difference test at 80% power using
`ceil((z_0.975 + z_0.8)^2 * sample_var(d) / mean(d)^2)`, with a minimum of
24 pairs. It reports whether the 200-pair frozen pool covers that estimate.
If batches report different CPU models, frequency policies, compilers, or
binary hashes, the analyzer leaves the pooled curve difference unknown.
The PSI, interrupt, steal, frequency, and perf snapshots span each whole
process, including preparation; the solver's separate online timer defines
the speedup ratio. Guest diagnostics do not certify host-wide exclusivity.

Local `--preview-upload` selected **1,615 files totaling 52,444,334 bytes**:
the exact Rust build inputs, required documentation, this frozen plan, and
the 400 public point files. A local filename and token-pattern scan found
zero credential-named or key-pattern files in that inventory. The preview,
remote-code compilation, source hashes, and point hashes passed locally.
No third-party upload or Modal job has begun; explicit approval for the
source and fixture transfer is pending after automatic review rejection.

```sh
python3 experiments/koblitz-s3-crosscurve-20261009/modal_crosscurve_panel.py --preview-upload
# After approved Modal panel outputs have been downloaded:
python3 experiments/koblitz-s3-crosscurve-20261009/analyze_modal_crosscurve.py --inputs /private/tmp/s3-crosscurve-modal-pilot.json.gz --output /private/tmp/s3-crosscurve-modal-pilot-analysis.json
```

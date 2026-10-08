# Modal S3 two-curve diagnostic, 2026-10-07

This directory records an **exploratory** CPU experiment for the shared-S3-
inversion change in PR #1518. The question and frozen acceptance rules are in
[`../PROTOCOL.md`](../PROTOCOL.md). Each run solves one previously generated
public target point; the timer inside the solver starts after reusable setup
and stops after scalar recovery and point replay. Baseline and candidate use
the same point, factor-base construction, 244 folded columns, relation seed
20260928, 14 workers, and resource request. The pilot has 12 distinct points
on each of the GF(2^41) and GF(2^53) Koblitz instances. Each target is one
independent statistical unit, with an `ABBA` or `BAAB` four-run pair and a
separate two-run baseline `A/A` control. The two Modal panels reuse these
same 24 points and must **not** be counted as 48 independent targets.

`candidate_id` and `workload_id` remain `null`: the complete Linux IC
candidate manifests and host receipt have not been frozen. These are
implementation stage diagnostics, not an IC-versus-rho or controlled IC
speedup result. The exact source snapshots are
[`../../koblitz-s3-pair-query-20261007-v4/baseline.rs`](../../koblitz-s3-pair-query-20261007-v4/baseline.rs)
and
[`../../koblitz-s3-pair-query-20261007-v4/candidate.rs`](../../koblitz-s3-pair-query-20261007-v4/candidate.rs)
at code commit `833ef86191947ba8e4def56a78ee6c51a4ea0ff2`.

## Evidence and measured outcome

| Panel | Modal run | N41 paired baseline/candidate online ratio | N53 paired baseline/candidate online ratio | N41 / N53 A/A p95 | Gate |
| --- | --- | ---: | ---: | ---: | --- |
| Initial exploratory pilot | [app run](https://modal.com/apps/a-buran28/main/ap-269VdVaOLiUmDJ37OYqX7i) | 1.0335 | 1.1083 | 1.183 / 1.039 | Fails N41 A/A and host receipt |
| Exact-output custody replication | [app run](https://modal.com/apps/a-buran28/main/ap-yVLMYfhChUnklZqm2G08f5) | 1.0542 | 1.0822 | 1.103 / 1.101 | Fails both A/A and host receipt |

The ratios are geometric means of 12 target-level `log(A/B)` values per
curve. Larger than 1 means the candidate was faster in this exploratory
panel. In the custody replication, 11/12 N41 and 12/12 N53 target pairs had
candidate-faster online times. Its **descriptive**, unadjusted 95% t
intervals are 1.0225–1.0869 for N41 and 1.0663–1.0983 for N53. The N53/N41
ratio of paired ratios is 1.0266 with a descriptive Welch interval
0.9936–1.0607. The first panel gave a noticeably different interaction;
these are same-target replications under an unqualified host, not evidence
that the curves have different controlled speedups. No p-value or interval
from these panels passes the predeclared fidelity gate.

The A/A acceptance threshold for the 10% interaction question is p95
`|log(A1/A2)| < log(1.05)` on **both** curves. The custody replication reached
about 10.3% and 10.1% instead of below 5%. All 146 solver invocations in
each panel exited successfully and passed native scalar/fixture checks; all
24 blocks had matching semantic witnesses across arms. The second panel
also retains the exact raw solver JSON text and confirms its SHA-256, stdout
identity, stable binary hashes, and exclusive online phase sum for every
run. See [`analysis_raw_v2.json`](analysis_raw_v2.json) for every target's
paired times, A/A values, run audit, descriptive intervals and failed
fidelity checks. The exact-output panel is
[`pilot_raw_v2.json.gz`](pilot_raw_v2.json.gz), SHA-256
`44829cd402523d4fe2cd68ec833611ab46ceb21c074079b3292695da5a253654`.
The earlier [`pilot.json.gz`](pilot.json.gz), SHA-256
`5c997a13089a71b71d1b09416b3df492cf99ec6efca6f39fe29f56a4247cb930`,
contains parsed reports and their recorded hashes, but not the exact raw text
needed to independently rederive each JSON hash. It remains preserved as a
noisy first panel.

The custody panel's unweighted arithmetic mean of the 24 paired-arm
invocations per curve gives the following stage diagnostics in milliseconds.
They are setup-separated and exploratory; the target's individual phases and
all failures are retained in the raw ledger. The reported peak RSS was about
447 MiB at N41 and 464 MiB at N53, with no material arm difference.

| Curve | Arm | Online total | Charged target PDP, including failed attempts | Reusable index build |
| --- | --- | ---: | ---: | ---: |
| N41 | Baseline | 118.934 | 105.874 | 268.708 |
| N41 | Candidate | 113.136 | 100.237 | 268.699 |
| N53 | Baseline | 1575.289 | 1560.908 | 303.237 |
| N53 | Candidate | 1455.567 | 1441.100 | 305.511 |

The requested fidelity data exposed a hard limit. The function container
reported `Linux-4.19.0-gvisor-x86_64`, 32 visible logical CPUs, and successful
`taskset` affinity to CPUs 0–14 with the sampler on CPU 15. It observed 15
solver threads, but did not expose their individual affinity masks, physical
SMT siblings, a NUMA node, an exclusive cgroup CPU partition, cgroup quota or
throttle counters, CPU/memory PSI, host IRQ routing, or a fixed-frequency
record. The reported CPU model was `unknown`. All four separate
whole-process `perf stat` attempts failed with `perf_event_open` error 19.
The container affinity mask therefore cannot establish host-wide isolation.
The panel used Modal image `im-UCF0gNcBgEF4Rf9ARqXKSG` in GCP `us-east`;
its source and binary hashes are in the raw host probe. The earlier
[`probe.json.gz`](probe.json.gz) predates a fallback correction for missing
SMT topology and selected only CPU 0; the panel's embedded host probe is the
applicable affinity evidence.

## Separate Modal VM perf diagnostic

A [Modal VM Sandbox](https://modal.com/docs/guide/sandboxes) on OCI
`us-central` reported an AMD EPYC 9J45 processor and allowed `perf stat` for
cycles, instructions and task-clock. The VM had 32 visible vCPUs, memory node
0, and `perf` reported 100% enabled/running for those counters. It still
exposed no exclusive host cpuset partition or PSI. The four records in
[`vm_perf_profile.json`](vm_perf_profile.json), SHA-256
`dd92d764d1d5d8a1ddfa166bb9b41e8333cacb85675a1797a595ee96325ae135`,
are one separately profiled **whole process** per arm and curve, each with a
verified target scalar:

| Curve | Arm | Instructions, billions | Cycles, billions | IPC |
| --- | --- | ---: | ---: | ---: |
| N41 | Baseline | 7.767 | 3.723 | 2.09 |
| N41 | Candidate | 7.608 | 3.776 | 2.01 |
| N53 | Baseline | 31.762 | 24.452 | 1.30 |
| N53 | Candidate | 28.003 | 23.160 | 1.21 |

These counters include factor-base and index setup, as well as the target.
They cannot be assigned to the target-online interval. Four single passes on
a different Modal host cannot validate either paired wall-time panel or
establish a performance claim. They show that instruction count changed in
the expected direction, while IPC also changed. The VM script terminates its
sandbox after capture.

## Independent mathematics and repeatability

The checked local Sage launcher was run with `--runtime-info` before the
replay; its receipt is [`sage_runtime_info.json`](sage_runtime_info.json).
[`../replay_modal_sage.py`](../replay_modal_sage.py) replays all 24 distinct
points from the exact-output panel on the declared curves, checks the
factor-base points and orbit labels, verifies each relation point sum,
reconstructs the rank and first target span, and recovers each scalar by an
independent matrix solve. All 24 points, 45,872 factor-base points and 5,880
relation witnesses passed; the receipt is
[`independent_sage_replay.json`](independent_sage_replay.json), SHA-256
`fbbb277b2d4348c569d820983a6037979c95ab397d9501874ba55d0b16bcb9e2`.
The N41 base
dump was generated on the local CPU from the frozen baseline source and its
digest matches every Modal N41 report; N53 uses the existing independently
checked frozen base fixture. The Sage work and all fixture generation remain
outside the online timers.

To reproduce the Modal diagnostics from this checkout with an authenticated
Modal profile:

```sh
modal run experiments/koblitz-s3-fidelity-20261007/modal_app.py::probe --output OTHER_PROBE.json.gz
modal run experiments/koblitz-s3-fidelity-20261007/modal_app.py::pilot --output OTHER_PILOT.json.gz
python3 experiments/koblitz-s3-fidelity-20261007/modal_vm_profile.py --output OTHER_VM_PROFILE.json
python3 experiments/koblitz-s3-fidelity-20261007/analyze_modal.py OTHER_PILOT.json.gz OTHER_ANALYSIS.json
```

Both runners refuse to overwrite existing evidence. For the Sage replay use
the absolute checked repository launcher, first saving its `--runtime-info`
JSON. The Modal CPU request is 16 cores and 8192 MiB; [Modal's resource
documentation](https://modal.com/docs/guide/resources) describes CPU requests
as minimum physical core allocations but also describes bursting, so that
request does not substitute for a host-level exclusive partition receipt.

The next confirmatory sample size remains **unset**. The prior planning value
of 92 independent targets per curve uses noisy historical N53 variance and
is only a budget estimate. The much lower variance inside either Modal panel
cannot replace a clean pilot because the host and A/A gates failed. A future
qualifying physical-host pilot must pass its predeclared A/A and isolation
checks before the power calculator freezes the new target count.

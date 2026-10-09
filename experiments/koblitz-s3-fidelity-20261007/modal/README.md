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
python3 experiments/koblitz-s3-fidelity-20261007/modal_vm_panel.py --count 12 --region eu-west --output OTHER_VM_PILOT.json.gz
python3 experiments/koblitz-s3-fidelity-20261007/analyze_modal.py OTHER_PILOT.json.gz OTHER_ANALYSIS.json
python3 experiments/koblitz-s3-fidelity-20261007/analyze_modal.py OTHER_VM_PILOT.json.gz OTHER_VM_ANALYSIS.json
/Volumes/SSD990/cryptanalysis/sage --runtime-info > OTHER_RUNTIME_INFO.json
/Volumes/SSD990/cryptanalysis/sage -python experiments/koblitz-s3-fidelity-20261007/replay_modal_sage.py \
  --panel OTHER_VM_PILOT.json.gz --runtime-info OTHER_RUNTIME_INFO.json \
  --receipt OTHER_VM_SAGE_REPLAY.json
```

Both runners refuse to overwrite existing evidence. For the Sage replay use
the absolute checked repository launcher, first saving its `--runtime-info`
JSON. The Modal CPU request is 16 cores and 8192 MiB; [Modal's resource
documentation](https://modal.com/docs/guide/resources) describes CPU requests
as minimum physical core allocations but also describes bursting, so that
request does not substitute for a host-level exclusive partition receipt.

The next confirmatory sample size remains **unset**. The prior planning value
of 92 independent targets per curve uses noisy historical N53 variance and
is only a budget estimate. The much lower variance inside any Modal panel
cannot replace a clean pilot because the host and A/A gates failed. A future
qualifying physical-host pilot must pass its predeclared A/A and isolation
checks before the power calculator freezes the new target count.

## Follow-up VM paired pilot, 2026-10-08

The new [`../modal_vm_panel.py`](../modal_vm_panel.py) ran the **same** 12
targets per curve in one Modal VM Sandbox (`sb-QXEf3K8h88u4XnYtRorGJ7`).
It requested 16 CPU cores and 8192 MiB in OCI `eu-west`, pinned the solver
to guest CPUs 1–15, pinned the sampler to guest CPU 16, and bound memory to
guest NUMA node 0. This was a third timing replication of the existing 24
distinct points, not 24 new statistical units. The VM panel used image
`im-BmwZ79t7Tkcc1SG9HC8BKi`, baseline binary SHA-256
`8769b58bd54cd15322456972f587cb577b0cac11d7a179b7cf7cabfddede81b1`,
and candidate binary SHA-256
`b5fdd31128bf3b9d151e6b470bf7b11770221e1d52aa268aeddca02fd8c741a3`.
The frozen baseline and candidate source SHA-256 values match the first
panel, but the binaries differ after the intervening `main` merge. Interpret
each panel only within its own binary and host pairing.

The [complete VM ledger](vm_pilot_20261008.json.gz) is 5,616,998 bytes,
SHA-256 `85a5a7dcb150f9fdec38055dfbca6762bf70eab08e6f7c0c1194c676251dc2c6`.
All 146 solver invocations exited successfully; exact raw text hashes,
stdout identity, unchanged binary hashes, fixture scalar recovery and
exclusive online phase sums passed. All six runs within each target block
had the same semantic witness; all 24 target traces also matched the
previous Sage-verified Function panel. The [analysis](vm_analysis_20261008.json)
retains each paired time, A/A control, host snapshot, descriptive interval
and failed claim gate. A separate [VM smoke](vm_smoke_20261008.json.gz)
used one point per curve before the CPU 0 exclusion. The first full VM
attempt lost its large stdout transfer, so no raw evidence from that attempt
is used; its [failure record](vm_transfer_failure_20261008.json) is retained.
A temporary Modal Volume then transferred the full rerun by size and SHA-256
and was deleted after validation.

| Curve | Exploratory geometric mean baseline/candidate online ratio | Descriptive 95% t interval | Candidate faster targets | A/A p95 ratio | A/A gate (<1.05) |
| --- | ---: | ---: | ---: | ---: | --- |
| N41 | 1.0559 | 1.0449–1.0670 | 12/12 | 1.0340 | Pass |
| N53 | 1.0895 | 1.0768–1.1024 | 12/12 | 1.0563 | **Fail** |

The descriptive N53/N41 ratio of paired ratios is 1.0318, with a Welch
interval 1.0166–1.0473. Its inferential status remains **unknown** under the
protocol because N53 failed the fixed A/A noise threshold and the VM has no
host-level isolation receipt. The single N53 A/A outlier was target T010,
where the two baseline intervals were 1189.895 and 1126.472 ms. The slower
run had one guest steal tick; that observation does not establish its cause.

The VM exposed an AMD EPYC 9J45 CPU model, guest SMT and NUMA topology,
thread affinity masks `1-15` for all runs, memory placement and per-run
interrupt/steal/cgroup snapshots. The guest cgroup reported no CPU
throttling. **Twenty-three of 146** run windows showed increasing guest CPU
steal ticks. The VM did not expose CPU/memory PSI, an exclusive host cpuset
partition, host IRQ routing, fixed frequency, or proof that no other tenant
shared the physical cores. These are missing host-level claim requirements.
The reported 32 VM vCPUs and guest NUMA node are not a physical-host
isolation receipt. The `physical_smt_topology_known` and `numa_node_known`
checks therefore remain false in the VM analysis even though the *guest*
equivalents are visible.

Four separate whole-process `perf stat` passes on the VM completed with
100% enabled/running event time. They include reusable setup and one target;
they cannot be assigned to the target-online interval:

| Curve | Arm | Instructions, billions | Cycles, billions | IPC |
| --- | --- | ---: | ---: | ---: |
| N41 | Baseline | 7.764 | 3.721 | 2.09 |
| N41 | Candidate | 7.626 | 3.703 | 2.06 |
| N53 | Baseline | 31.774 | 22.900 | 1.39 |
| N53 | Candidate | 28.280 | 20.935 | 1.35 |

Mean online wall times in the unprofiled VM pairs were 61.928/58.649 ms
for N41 baseline/candidate and 1158.606/1063.598 ms for N53. Their charged
target PDP means were 53.622/50.275 ms and 1149.327/1054.338 ms,
respectively. Reusable index build was about 63 ms on N41 and 86 ms on N53
in this VM. Each target's exclusive phase costs, relation counts and peak
memory remain in the raw ledger.

The checked Sage launcher was run with [`--runtime-info`](vm_sage_runtime_info_20261008.json)
before the VM replay. Its separate receipt
[`vm_independent_sage_replay_20261008.json`](vm_independent_sage_replay_20261008.json)
has SHA-256 `9472443c3ffc101bc0af25621726d685dff6531c6625e5bdef4d01ac8929c22c`
and binds the VM ledger hash, all 45,872 factor-base points, 5,880 relation
witnesses, independent matrix recoveries, and all 24 public points. Sage
replay time is outside the online intervals. The next confirmatory sample
size remains unset; the 92-per-curve value is still only the earlier budget
estimate.

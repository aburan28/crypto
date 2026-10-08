# What was tested, where, and what it showed

This is the record behind isolab 0.1.0. Every number here was read from a
run on 2026-10-01; nothing is extrapolated. Three environments were used.

| environment | role | what it can and cannot show |
|:--|:--|:--|
| **Lab host:** QEMU/HVF virtual machine on Apple M4 Pro, Ubuntu 24.04, kernel 6.8.0-142-generic, 8 vCPUs presented as 4 cores × 2 threads across **2 NUMA nodes** (cpus 0-3 / 4-7, 3.9 + 4.0 GiB), cgroup v2, PSI on, podman 4.9 + crun, gVisor release-20260928.0, perf 6.8.12, nats-server 2.15.0 | hub and worker, installed by `deploy/install.sh` as systemd services, root | cpuset isolated partitions, eviction, IRQ affinity, NUMA binding, PSI gating, podman and gVisor, the installer. **Cannot** show hardware counters (HVF exposes no PMU), cpufreq or turbo control (none exposed), bare-metal steal (always 0 here) |
| **Agent host:** this macOS machine | CLI and the MCP server over stdio, connecting to the VM hub through a port forward | the real client topology: a laptop driving a remote lab; direct-backend worker for the unprivileged suite |
| **CI:** GitHub `ubuntu-24.04` runner | `.github/workflows/isolab.yml` | the unprivileged suite, the root suite under `sudo`, rootless podman, the base image build |

## Suites

| suite | where | result |
|:--|:--|:--|
| unit: protocol, synthetic-sysfs topology (2 sockets × 4 cores × 2 threads, 2 nodes), planner, matcher, perf CSV, ELF sniffing, blob cache, tiers, post-checks, result compaction | macOS and VM | 34 passed, 1 Linux-only skip on macOS; all pass on the VM |
| fabric against a real JetStream server: submit, idempotency, claim, heartbeat, complete-once, cancel queued and running, decline with reasons and redelivery, requeue counting attempts to `dead`, lost-worker takeover with the old worker fenced on heartbeat, complete and requeue, blobs, artifacts, roster TTL, wait | macOS (local nats-server) and VM | 10 passed in both |
| end-to-end, in-process worker, direct backend, best-effort policy: repeats with warm-up, metrics and artifacts with hashes, blob inputs, failure, timeout, failed build step, cancel of a running job, unplaceable job waiting with recorded decline reasons, in-container verifier, every MCP tool, calibration | macOS | 9 passed, 2 platform skips |
| the podman path **rootless** (as an unprivileged user; the engine delegates cpu, memory and pids but no cpuset, so the container runs unpinned and the result records it, tier D) | VM as `ubuntu`, CI runner | 1 passed |
| the same plus **privileged** (root: cgroup partition, eviction, cgroup statistics) and **podman** (rootful, `--cgroup-parent` inside the job cgroup, then the same job under `runsc`) | VM, as root | **56 passed** |

## The deployed lab, driven from the Mac

`deploy/install.sh --hub --worker --cpus 2-7 --pool default --label site=vm
--default-image localhost/isolab-base:latest` on the VM: venv, nats-server
download with checksum, launcher and calibration kernel compiled static,
token generated, both units enabled and active. The worker registered as
`lab cpus [2..7], housekeeping [0, 1], backend podman (crun), tier up to B,
perf hw=False`. The base image built on the VM (1.25 GB, `WITH_RUST=0`).

From the Mac, with `ISOLAB_URL=nats://TOKEN@127.0.0.1:4222`:

**Roster.** `isolab workers` shows the VM with its images, lab CPUs,
capacity per SMT mode (`isolate`: at most 2 whole cores on one node) and
tier B.

**Strict policy refused with reasons.** `isolab run --policy strict -- true`
is rejected before enqueueing: *can reach isolation tier B, job needs A;
requires bare metal; this worker is a virtual machine; requires
perf_counters, which this worker cannot provide.* That is the intended
behaviour on a VM.

**A NumPy job in the base image, 2 cores, standard policy, 1 warm-up + 3
repeats** (`isolab run --image localhost/isolab-base:latest --cpus 2
--memory-mb 1024 ...`):

| field | value |
|:--|:--|
| status / tier / grade | `succeeded` / B / C (downgraded: the summary had to use contended repeats) |
| placement | cpus `[4, 6]`, idle siblings `[5, 7]`, node 1, partition `isolated`, 7 IRQs moved, 20 per-CPU IRQs unmovable (EIO), NUMA bound by `cpuset.mems` |
| runtime | podman, crun, image digest `sha256:87a291ad…` |
| timed walls (in-container launcher) | 0.384, 0.349, 0.325 s; host-side walls 0.508, 0.496, 0.485 s (exec overhead 0.12–0.16 s, excluded from the measurement) |
| cgroup CPU per repeat | 0.39, 0.35, 0.33 s; peak memory 45 MiB; 0 throttling, 0 OOM |
| settle window | other CPU 0.0067, pressure share 0.11 %, reserved CPUs 0.08 % busy, steal 0 |
| during the repeats | interrupts 186–228 /s per reserved CPU (the local timer tick; no `nohz_full`), steal 0, host-wide CPU pressure share 39–43 % |
| artifacts | 14 files, each hashed, logs and per-repeat samples included |

The host-wide pressure during the repeats is the virtual machine: with two
housekeeping vCPUs carrying the worker, podman and the hypervisor's own
noise, some task is usually waiting somewhere. After this run the per-run
check was moved to the job cgroup's own pressure (whether *the job's*
tasks stalled), with the host-wide figure kept as information; the
calibration below and the MCP run were done under that rule.

**Calibration (noise floor).** `isolab calibrate --worker isolab-vm
--repeats 5 --iterations 300000000 --backend direct`: 5 timed repeats,
identical checksum `8ad8bf87a2ec7616` every time, wall 0.613–0.648 s,
**coefficient of variation 2.5 %**. That spread is this VM's noise floor,
and it is why nothing measured here would count as a performance claim;
the point of the run is that the kernel, the pinning and the record all
work.

**Cancel.** A `sleep 120` submitted from the Mac, cancelled from the Mac
after 6 s: the record goes `running → cancelled`, error *cancelled during
repeat 0*, and the worker is free again.

**MCP.** `isolab mcp` was driven over stdio by a Python MCP client as an
agent host would: `tools/list` returns the 16 tools; `isolab_overview`,
`isolab_run` with a `local_file` input (uploaded as a 310-byte blob by
hash), `isolab_wait`, `isolab_result` (summary and fidelity sections),
`isolab_logs`, `isolab_fetch` (14 files, every hash verified on the Mac),
`isolab_jobs` by label and `isolab_doctor` all returned as designed. The
job ran in the base image on cpus `[4, 6]` of node 1 under an isolated
partition, grade B.

## The gate working as intended

During one VM suite run the base image was still building on the same
machine. The privileged test's job refused to start five times in a row
(*host did not go quiet: other CPU 0.20 and 0.26 against a limit of 0.05,
CPU pressure 2.6–3.4 %*), was retried as an infrastructure failure and
reported as `infra_error`. Once the build finished the same test passed.
A `standard` job never starts on a busy host; that is the contract.

## Fixed along the way, because the tests found them

* A normalised spec with `verify: null` did not re-validate on the worker,
  so a worker declined its own jobs: the schema now accepts the nulls that
  normalisation writes and normalisation is tested to be idempotent.
* The sampler thread's stop flag was named `_stop`, shadowing a private
  `threading.Thread` method; every quiesce window crashed on `join`.
* The settle check used PSI `avg10`, a ten-second decaying average that
  still held the worker's own start-up compile; it now uses the share of
  the window itself spent under pressure, from the `total` counters.
* The memory pre-check read `MemFree` from the start-up snapshot; it now
  reads free plus reclaimable cache on the job's node, live.
* Per-CPU and managed interrupts refuse affinity writes with EIO on every
  host; they are now listed as unmovable instead of failing the mechanism,
  and the interrupt rate that actually lands on the reserved CPUs is a
  pass/fail check.
* gVisor's current tarball needs its `gvisor-bin/` sidecars installed next
  to `runsc`; the installer copies the whole archive.
* `smt: off` keyed on the `smt/control` string, which an Arm VM reports as
  `notimplemented` while still exposing two threads per core.

## Slots and the Sage image, 2026-10-08

Run in a cloud container: x86-64 Xeon @ 2.30 GHz, 4 vCPUs on one NUMA node,
kernel 6.18, **cgroup v1** (so no cpuset partitions: tier C at best), root,
docker 29.6 with runc, nats-server 2.10.22, no podman.

| what | result |
|:--|:--|
| full suite (`pytest` in `isolab/`) | 64 passed, 3 skipped (privileged and podman suites, which need cgroup v2 root and podman) |
| `tests/test_slots.py` | slot splitting on the synthetic two-node host, the host-settings refcount and its crash path, IRQ undo from two slots in both orders, shared / exclusive / drain locking, the grade cap, and two in-process slots running overlapping jobs against a real hub |
| CLI: `isolab worker --cpus 1-3 --slots 2 --backend direct`, then two `isolab run` jobs | ran at the same time on cpus `[1]` (slot s0) and `[3]` (s1); each result lists the other job in `co_tenant_jobs`, `co_tenant_cpu_s` 1.0 and 1.87 against a 2.0 s wall, `other_cpu_s` 0.02, `shared_host: true` |
| `isolab images build sage --tool docker --ca-bundle …` | built on SageMath 10.10; `sage -c 'factor(2^64+1)'`, `import pandas, gmpy2` pass with every capability dropped, a read-only root and no network |
| a Sage job through the docker backend (`isolab run --backend docker --image localhost/isolab-sage:latest --repeats 2 -- sage f.sage`) | succeeded, both repeats timed by the in-container launcher (3.44 s, 2.98 s), image digest pinned, tier C |

Found and fixed on the way:

* The docker backend could not resolve any image: its `image inspect`
  template named `.Digest`, which only podman has, so docker failed the
  template and the worker tried to pull `localhost/…` from a registry.
* The Sage image could not run a job. `/home/sage` is mode 0750, and a job
  runs as root with every capability dropped, so `sage` was "not found".
  The upstream entrypoint also re-joins its arguments into one shell
  string. The recipe now opens the home directory to any uid and clears
  the entrypoint.
* The Sage recipe swallowed its own failures (`|| true`), so a build behind
  a TLS-intercepting proxy "succeeded" without pandas. Failures are now
  fatal, and `--ca-bundle` passes a proxy CA as a build secret.
* A relative `--state-dir` broke the launcher path.
* Without a job cgroup, a job's "own" processes were every descendant of
  the worker, which in a slotted worker includes the other slots' jobs.

The CUDA image (`nvidia/cuda:12.6.3-devel-ubuntu24.04`, 12 GB) also built
cleanly here; `nvcc` 12.6.85, NumPy 1.26.4 and numactl 2.0.18 answer inside
it with no capabilities and a read-only root. There is no GPU here, so this
proves the recipe and not the CUDA runtime. Behind this container's proxy,
apt could not verify NVIDIA's own package repository and fell back to the
Ubuntu archive, which carries everything the recipe installs.

CI now builds and smokes both images (`image-sage`, `image-cuda` in
`.github/workflows/isolab.yml`).

## Not shown here

Hardware counters, `cpu-migrations` from perf, cpufreq and turbo control,
thermal throttle counters and tier A as a whole were exercised only as
`unavailable`. They need a bare-metal Linux host and remain to be
confirmed there with `isolab doctor` and a `strict` calibration. Slots with
real cgroup v2 partitions side by side are covered by the privileged CI
step, not by a lab host yet. The CUDA image has never run on a GPU.

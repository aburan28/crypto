# isolab: a standalone experiment execution environment with enforced measurement fidelity

**isolab makes no performance claim about anything it runs.** It is the room
an experiment runs in. Its job is to make sure that when a result says "this
took 4.21 s and 9.8 × 10⁹ instructions on four cores", nothing else was on
those cores, the memory was local, the clock was steady, and every one of
those facts is written down next to the number, so that a reader can tell a
measurement from a guess.

An agent (Claude Code, Codex, OpenCode, a script) asks over MCP for a run of
some command, in some image, on some shape of hardware. A worker on a Linux
host claims the job from a distributed queue, reserves the hardware, measures
the run under isolation, and publishes the result. The agent reads it back
over MCP. Nothing in between needs a cloud account, Kubernetes or a database.

It sits beside two older tools in this repository and replaces neither:

| tool | what it is | what isolab adds |
|:--|:--|:--|
| `tools/isolated_bench.py` | a one-shot script that reserves a core on the machine you are already on | the same discipline as a service, on a remote host, with containers, NUMA, counters and a verdict |
| `taskq/` | a Redis queue for a Kubernetes fleet that runs commands at pinned commits | no Redis, no cluster, user-supplied binaries and images, hardware-shaped placement, fidelity enforcement |

isolab keeps taskq's protocol discipline, since it was right: specs are
normalised and hashed, results are written once under a fence, a run that
did not finish is never evidence, and a solver's claim counts only when
independent code agrees.

## Quick start

On the server (Linux, root), install a hub and a worker as systemd services:

```sh
git clone https://github.com/aburan28/crypto && cd crypto/isolab
sudo deploy/install.sh --hub --worker --cpus 2-7 --pool default
sudo isolab images build base          # the default container image
sudo grep ISOLAB_TOKEN /etc/isolab/hub.env   # the hub token clients need
```

On the machine running Claude Code, register the MCP endpoint:

```sh
claude mcp add isolab -e ISOLAB_URL=nats://TOKEN@SERVER:4222 -- \
  uvx --from 'git+https://github.com/aburan28/crypto#subdirectory=isolab' --with mcp isolab mcp
```

Then ask Claude to call `isolab_overview` (is the lab up?), `isolab_run`
(command, image, inputs, cpus, repeats, policy), `isolab_wait` and
`isolab_result`. The same from a shell:
`isolab run --cpus 2 --repeats 3 -- ./bench` then `isolab result JOB`.
Untrusted binaries go under gVisor with `oci_runtime: runsc`. The rest of
this page is the detail.

## Shape

```
  agent host (laptop, CI, a cloud sandbox)              lab hosts (Linux, root, your hardware)
  ┌─────────────────────────────┐                      ┌──────────────────────────────────┐
  │ Claude Code / Codex / …     │                      │ isolab worker                    │
  │   └─ isolab mcp (stdio)  ───┼──── NATS ────┐       │   inventory + heartbeat          │
  │ isolab CLI               ───┼──────────────┤       │   claim → reserve → measure      │
  └─────────────────────────────┘              │       │   podman/docker/gVisor/direct    │
                                               ▼       │   perf, cgroups, PSI, IRQs       │
                                   ┌──────────────────┐│                                  │
                                   │ isolab hub       ││  images: isolab-base, -sage, -cuda│
                                   │ (nats-server +   │◄┼── pulls work, publishes results  │
                                   │  JetStream)      ││                                  │
                                   └──────────────────┘└──────────────────────────────────┘
                                     work queue · job KV · results KV · blobs · artifacts
```

**The hub is a NATS JetStream server** with a config isolab writes for it.
Any lab host can run one; hubs on several hosts form a cluster, and a hub
can hang off another as a leaf node, so a lab is a mesh of peers rather
than a service somebody else runs. Clients and workers connect to any hub.
The hub holds:

| resource | kind | role |
|:--|:--|:--|
| `ISOLAB_WORK` | stream, work-queue retention, subjects `isolab.work.<pool>` | one message per job; a message is deleted when its job is acknowledged |
| `ISOLAB_JOBS` | key-value | the job record: spec, state, attempt, fence, worker, declines |
| `ISOLAB_RESULTS` | key-value, create-once | the result document, written exactly once |
| `ISOLAB_PROGRESS` | key-value, 6 h TTL | phase, repeat, elapsed and the log tail of a running job |
| `ISOLAB_WORKERS` | key-value, 60 s TTL | live roster: each worker's inventory, refreshed by heartbeat |
| `ISOLAB_IDEM` | key-value, create-once | idempotency key → job id |
| `ISOLAB_BLOBS` | object store | content-addressed inputs (user binaries, scripts) |
| `ISOLAB_ARTIFACTS` | object store | everything a run wrote to its output directory, plus logs and samples |
| `ISOLAB_EVENTS` | stream, 30 d | audit feed of every transition |

**Queue semantics.** Execution is at-least-once; the result commit is
exactly-once. A worker that pulls a message holds a lease (the consumer's
ack wait) and extends it while it runs. If the worker dies, the hub
redelivers the message after the lease lapses and another worker takes the
job over under a new fence. Every write a worker makes to the job record is
a compare-and-swap against the record's revision, so a worker that was
frozen past its lease cannot overwrite the job when it wakes: its next write
fails, it kills its child and writes nothing. The result key is created,
never put, so a second result for the same job is refused by the hub itself.
A job a worker cannot place (wrong pool shape, no such image, too few
cores) is handed back without counting as an attempt and the reason is
recorded on the job, so `status` says why a job is waiting. Only
infrastructure failures retry; a nonzero exit or a timeout is an outcome
and is recorded as such.

## The fidelity model

Every confounder isolab knows about is handled the same way: there is a
**mechanism** that removes it, a **measurement** that shows whether it was
removed, and a line in the result that records both. A result never claims
isolation it did not verify.

| confounder | mechanism | evidence recorded |
|:--|:--|:--|
| other processes on the job's cores | cgroup v2 cpuset **isolated partition** holding the job CPUs, which also removes them from the scheduler's load balancing; every other movable thread evicted with `sched_setaffinity`; one job per slot, and a `strict` job takes the whole host (see [Slots](#slots-several-jobs-on-one-host)) | partition state as the kernel reports it, threads moved, user threads that refused to move, CPU seconds other processes used during the run, `cpu-migrations` on the job CPUs |
| interrupts landing on the job's cores | every IRQ's `smp_affinity_list` rewritten to the housekeeping CPUs, and `default_smp_affinity`; per-CPU and kernel-managed interrupts (local timer, IPIs, managed device queues) refuse the write and are listed as unmovable, which boot-time `isolcpus=managed_irq` and `irqaffinity` address | IRQs moved and unmovable; interrupts per second that actually landed on each reserved CPU during the run, failed above 1500/s |
| SMT sibling sharing the core | `smt: isolate` allocates whole cores and keeps the sibling inside the partition but idle; `smt: off` requires SMT disabled on the host | siblings reserved idle; `smt/control` |
| remote memory, NUMA balancing | job CPUs chosen from one node; `cpuset.mems` bound to that node; `numa_balancing` recorded, host-prep turns it off | node, mems, node free memory, memory peak |
| memory pressure, swapping | `memory.max`, `memory.swap.max=0`, OOM events read from the job cgroup | peak bytes, OOM kills, memory PSI |
| contention anywhere on the host | pressure stall information (CPU, memory, I/O; `some` and `full`) differenced over the settle window and over the run, so the check is the share of *that interval* spent under pressure rather than a decaying average; a settle window that must be quiet before the run starts, retried then failed as an infrastructure error | pressure share per phase, `avg10` peaks, load average, per-second samples kept as an artifact |
| CPU frequency | governor set to `performance` and turbo disabled for the run where the host exposes the controls; frequency sampled per job CPU every second | governor, turbo state, frequency min/mean/max and coefficient of variation per CPU, thermal throttle counts |
| the hypervisor | steal time on the job CPUs sampled; virtualisation detected and recorded; `bare_metal` can be required | steal ticks during the run, hypervisor name |
| network and disk | containers run with `network: none` by default; host NIC and disk byte counters sampled | bytes moved during the run |
| container start-up and the measuring process itself | a static launcher inside the container takes the wall clock and `wait4` rusage around the user's process only; `perf stat` runs on the host over exactly the job CPUs; the worker, its NATS client and its samplers are pinned to the housekeeping CPUs | wall from inside the container next to wall from the host; `task-clock` on the job CPUs |
| caches, page cache, THP, ASLR | optional `drop_caches` before each repeat; THP and ASLR settings recorded; `aslr: off` for the direct backend | settings as read |
| the algorithm itself | warm-ups, repeats, cooldowns; instructions and cycles as the primary metric because they survive frequency and hardware | per-repeat counters and a summary with spread |

### Tiers, policies and the verdict

What a host can provide is graded into a **tier**, chosen at reservation
time from the strongest mechanism that worked, and the result carries the
tier it actually got:

| tier | what held |
|:--|:--|
| **A** | valid isolated cpuset partition, eviction done, IRQs moved, NUMA bound, hardware counters readable |
| **B** | cgroup pinning (cpuset without a valid partition), eviction done, IRQs moved, NUMA bound |
| **C** | CPU affinity pinned and NUMA bind attempted; nothing could be evicted (unprivileged or no cgroup write access) |
| **D** | affinity pinning only; no Linux isolation facilities (a development machine, rootless CI) |

The job's `fidelity.policy` says what is acceptable:

| policy | minimum tier | settle window | refuses to start when | after the run |
|:--|:--|:--|:--|:--|
| `strict` | A | 5 s, pressure present ≤ 0.5% of the window, other CPU ≤ 0.02 | any required mechanism is unavailable, the machine will not go quiet, the host is a VM, counters are missing | a contended repeat is re-run; the summary uses clean repeats only |
| `standard` | C | 3 s, pressure present ≤ 1% of the window, other CPU ≤ 0.05 | the machine will not go quiet | contended repeats are flagged and excluded from the summary |
| `best_effort` | D | 1 s | never | everything recorded, nothing excluded |

Every check is one of `pass`, `fail`, `unavailable` or `info`, split into
the pre-run checks (static conditions and the settle window) and the
post-run checks (what the samplers and counters saw). A `fail` on any
post-run check marks the repeat `contended`. The result's `fidelity.grade`
is the tier, downgraded one step if any repeat used in the summary was
contended, and the `violations` list says exactly which checks failed.
**Fidelity never rewrites `status`**: a run that exited 0 under contention
is `succeeded` and `contended`, and the reader decides what it is worth.
Among the things a worker cannot control (another tenant of a VM, system
management interrupts, thermal headroom), the residual is what the A/A
calibration measures: `isolab calibrate` runs a fixed-work kernel ten times
on a worker and reports the spread, which is the noise floor every
difference on that host is read against.

## Execution environments

**Backends.** `podman` is the default and the one that reaches tier A: the
container is created under the job's cgroup, so the partition, the memory
limit and the counters all apply to it. `docker` works the same way except
that its systemd cgroup driver places the container outside the partition
(cpuset pinning still holds, so it tops out at tier B). `direct` runs the
process on the host inside the job cgroup with no container at all, which is
the right choice for a trusted binary when the host already has the
toolchain. Any container backend takes `oci_runtime: runsc` to run under
**gVisor**; that is the sandbox to use for a binary you do not trust. It is
still pinned and counted, but the sentry's own CPU time is charged to the
job and syscall-heavy code runs slower, so a gVisor measurement is only
comparable with another gVisor measurement, and the result says so.

Containers run with every capability dropped, `no-new-privileges`, a
read-only root filesystem, a private `/tmp`, a pids limit, no network unless
asked, and an automatic user namespace. A hub token therefore means "can run
code on every worker of this lab, inside those boundaries"; treat it as a
credential.

**Images.** `isolab/images/` holds three Containerfiles. `isolab-base` has
Python 3.12 with NumPy, SciPy, SymPy, Cython, gmpy2 and pandas, PARI/GP,
FLINT, GMP, MPFR, GCC, Clang, CMake, Rust (stable; `--build-arg WITH_RUST=0`
to skip), hwloc and numactl. `isolab-sage` is the official `sagemath/sagemath`
image with the same Python additions on top; that upstream image is published
for x86-64 only, so it builds on an x86-64 worker and not on Arm. `isolab-cuda`
is `nvidia/cuda` devel with Python and NumPy. Build them on a worker with
`isolab images build base`; the worker lists what it has, with digests, and
every result pins the image digest it ran. A job may name any image the
worker has; `placement.require_images` (default true) keeps a job off a
worker that lacks it.

Two container details worth knowing. Podman's `conmon` for the job's
`exec` session lives inside the job cgroup, so a command that writes
megabytes to stdout pays for relaying them on its own CPUs; write outputs
to `$ISOLAB_OUTPUT_DIR` and keep stdout small. And gVisor's release tarball
now ships `runsc` together with a `gvisor-bin/` directory of sidecar
binaries (`gvisor_sentry` among them) that `runsc` expects next to itself;
the installer copies all of it, and a `runsc` installed by hand without that
directory fails with "sidecar gvisor_sentry not usable".

**Inputs.** A job carries its own files. Each entry in `inputs` is a path
under the working directory and one of: inline `content` (a Sage script, a
parameter file), a `local_file` on the agent's machine (the MCP server or
CLI hashes and uploads it to the blob store and rewrites the entry to its
`sha256` before submitting), or a `git` repository at a full commit, which
the worker exports from a cached mirror. Binaries are checked for
architecture against the worker before they run. Everything the command
writes to `$ISOLAB_OUTPUT_DIR` comes back as artifacts, hashed; a
`metrics.json` there is parsed into the run record.

**The contract inside the box.** The working directory is `/isolab/work`
(inputs land here; `command.cwd` is relative to it). `ISOLAB_OUTPUT_DIR` is
a fresh directory per repeat. `ISOLAB_SCRATCH` is a tmpfs sized by
`resources.scratch_mb`. `ISOLAB_JOB_ID`, `ISOLAB_ATTEMPT`, `ISOLAB_REPEAT`
and `ISOLAB_WARMUP` identify the run. `build` steps run before any clock
starts, in the same container, and are measured and reported separately.

**Verification.** `verify.argv` runs a checker inside the container after
each timed repeat's clock has stopped; exit 0 is verified, 1 refuted, 2 no
claim. `verify.builtin: taskq-certificate` runs taskq's independent
discrete-log verifier on the worker host against the repeat's
`certificate.json`, so a solver's claimed answer is checked by code that
shares nothing with the solver.

## Slots: several jobs on one host

A worker runs one job at a time. A big host can run several at once by
cutting its lab CPUs into **slots**:

```sh
isolab worker --cpus 4-63 --slots 4              # four slots of whole cores, in NUMA order
isolab worker --slot 4-19 --slot 20-35 --slot 36-63   # or name each slot's cpus
sudo deploy/install.sh --worker --cpus 4-63 --slots 4 ...
```

Each slot registers as its own worker (`HOST.s0`, `HOST.s1`, …, label
`slot=sN`) with its own capacity, claims its own jobs, and builds its own
isolated partition (`isolab.lab.sN`) and job cgroup on its own cores.
Siblings of a core always stay in the same slot. What the slots share and
how each part is kept honest:

| shared thing | how slots keep out of each other's way | what the result records |
|:--|:--|:--|
| the host lock | each slot holds the host lock **shared** plus its own slot lock exclusively; a `strict` job takes the host lock **exclusively**, so it waits for running slots to finish, and while it waits no slot starts a new job | `isolation.host_exclusive`, `isolation.slot_lock`; the `exclusive_lock` check is `info` instead of `pass` when shared |
| thread eviction and IRQ affinity | each slot moves threads and IRQs off *its* CPUs and on exit gives back only those CPUs, starting from the current mask, so two slots undo correctly in either order | as before |
| `default_smp_affinity`, turbo | reference-counted in a small state file next to the lock: the first holder saves the original, the last one restores it; holders whose process died are dropped | as before |
| other processes' CPU | a job running in another slot is a **co-tenant**: its CPU is not counted as `other_cpu`, so slots do not fail each other's settle and run checks | `co_tenant_cpu_s` and `co_tenant_jobs` per repeat; `settle_co_tenants` and `co_tenants` checks |
| last-level cache, memory bandwidth, package power and frequency | not partitioned (no RDT yet) | `fidelity.shared_host` is true when a co-tenant ran during a summarised repeat, and the **grade is capped at B** |

So: use slots for throughput (sweeps, many small jobs, builds), and
`strict` for a number you want to quote as tier A. Separate worker
processes with their own `--cpus` still work, but they serialise on the
host lock exactly as one worker would.

## Protocol

Two versioned JSON documents, `isolab.job/v1` and `isolab.result/v1`,
with schemas in [`isolab/schemas/`](isolab/schemas/). The MCP tool
`isolab_schema` returns both. A spec is validated, every default written
out, then hashed; the result carries `spec_sha256` of the spec it executed.

```json
{
  "schema": "isolab.job/v1",
  "name": "ic-relation-search-smoke",
  "pool": "default",
  "runtime": {"backend": "podman", "oci_runtime": "crun", "image": "localhost/isolab-base:latest", "network": "none"},
  "inputs": [
    {"path": "bin/solver", "local_file": "./target/release/solver"},
    {"path": "params.json", "content": "{\"n\": 61, \"seed\": 7}"}
  ],
  "command": {"argv": ["./bin/solver", "--params", "params.json"], "env": {"RAYON_NUM_THREADS": "4"}},
  "measure": {"warmups": 1, "repeats": 5, "cooldown_s": 2},
  "resources": {"cpus": 4, "memory_mb": 16384, "numa_node": "single", "smt": "isolate"},
  "fidelity": {"policy": "strict", "require": ["perf_counters", "bare_metal"]},
  "placement": {"cpu_flags": ["avx2", "pclmulqdq"]},
  "limits": {"timeout_s": 1800},
  "verify": {"builtin": "taskq-certificate"},
  "labels": {"experiment": "EXP-0042"}
}
```

A result has, in order: identity (`job_id`, `spec_sha256`, `attempt`,
`fence`); `status` and `outcome_class`; the worker and a full `host`
manifest (CPU model and flags, microcode, kernel and its command line,
topology, memory, governor, turbo, THP, ASLR, mitigations, hypervisor,
runtime and image versions with digests); `placement` (CPUs, idle
siblings, node, mems, housekeeping CPUs, isolation tier and which mechanism
produced it); `fidelity` (policy, grade, every check, violations); `timing`
per phase; `build[]`; `runs[]`, one per repeat with wall from inside and
outside the container, `wait4` rusage, cgroup CPU and memory statistics,
perf counters, the sampler summary and `contended`; `summary` over the
clean timed repeats (min, median, mean, max, stdev and coefficient of
variation for wall, CPU, instructions and cycles, plus instructions per
cycle); `artifacts[]` with hashes; `verification`; `labels`.

| status | meaning | `outcome_class` |
|:--|:--|:--|
| `succeeded` | every build step and every repeat exited 0 | `completed` |
| `failed` | a build step or repeat exited nonzero | `completed` |
| `timeout` | a build step or repeat exceeded its limit | `not_completed` |
| `cancelled` | the submitter cancelled a running job | `not_completed` |
| `infra_error` | fetch, image, reservation or worker failure persisted through `max_attempts`, or the host never went quiet | `not_completed` |

States: `queued → claimed → running → <terminal>`; `queued → cancelled`; a
lost worker's job returns to `queued` and `attempt` increments.

## MCP

`isolab mcp` serves these tools over stdio. The server runs on the agent's
machine, uploads local inputs from there, and writes fetched artifacts
there.

| tool | does |
|:--|:--|
| `isolab_overview` | hub reachability, workers online with a one-line shape each, queue depths, running jobs. Call this first. |
| `isolab_workers` / `isolab_worker` | roster, or one worker's full inventory, topology, images and the isolation tier it can reach |
| `isolab_doctor` | a worker's readiness for `strict` runs, each check with the command that fixes it |
| `isolab_run` | the common case flattened: command, image, inputs, cpus, memory, repeats, policy, worker. Validates against the live roster, uploads inputs, returns the job id and which workers can take it |
| `isolab_submit` | a full `isolab.job/v1` document |
| `isolab_status` | state, worker, phase, repeat, elapsed, log tail, and the decline reasons while a job waits |
| `isolab_wait` | block up to 300 s for a terminal state |
| `isolab_result` | the result; `section` selects `summary`, `fidelity`, `runs`, `host` or `full` |
| `isolab_logs` | stdout and stderr tails of a finished or running job |
| `isolab_artifacts` / `isolab_fetch` | list a job's artifacts; download one or all to a local directory |
| `isolab_jobs` | list by state, pool, worker or label |
| `isolab_cancel` | cancel a queued job or stop a running one |
| `isolab_calibrate` | run the A/A noise-floor kernel on a worker |
| `isolab_schema` | the two JSON schemas |

Install, with `ISOLAB_URL` pointing at any hub (`nats://TOKEN@host:4222`;
`tls://` for TLS):

```sh
# Claude Code
claude mcp add isolab -e ISOLAB_URL=nats://TOKEN@lab-host:4222 -- \
  uvx --from 'git+https://github.com/aburan28/crypto#subdirectory=isolab' --with mcp isolab mcp
```

```json
// .mcp.json (Claude Code, project scope)
{"mcpServers": {"isolab": {"command": "isolab", "args": ["mcp"],
  "env": {"ISOLAB_URL": "nats://TOKEN@lab-host:4222"}}}}
```

```toml
# ~/.codex/config.toml
[mcp_servers.isolab]
command = "isolab"
args = ["mcp"]
env = { ISOLAB_URL = "nats://TOKEN@lab-host:4222" }
```

```json
// opencode.json
{"mcp": {"isolab": {"type": "local", "command": ["isolab", "mcp"],
  "environment": {"ISOLAB_URL": "nats://TOKEN@lab-host:4222"}}}}
```

`pip install 'isolab[mcp] @ git+https://github.com/aburan28/crypto#subdirectory=isolab'`
puts the `isolab` command on the path for the last three.

## Installing a lab host

On a Linux host, as root:

```sh
git clone https://github.com/aburan28/crypto && cd crypto/isolab
sudo deploy/install.sh --hub --worker --cpus 4-15 --pool default --label site=home
```

That installs the package into `/opt/isolab`, downloads `nats-server` and
`runsc`, builds the in-container launcher and the calibration kernel,
writes `/etc/isolab/{hub,worker}.env` with a generated token, and enables
`isolab-hub.service` and `isolab-worker.service`. `--cpus` is the lab
partition: the CPUs the worker may hand to jobs. Leave at least one full
core (both SMT threads) per NUMA node outside it for the OS, the worker and
interrupts. Then:

```sh
sudo isolab doctor            # what strict runs still need, and the command for each
sudo isolab doctor --apply    # apply the runtime-tunable ones (governor, turbo, watchdog, sysctls, irqbalance)
sudo isolab images build base # build the base image on this host (add sage, cuda as needed)
isolab calibrate              # measure this host's noise floor
```

`doctor` checks and `host-prep.sh` applies: cgroup v2 with the cpuset
controller, kernel 6.1 or newer for isolated partitions, `isolcpus`,
`nohz_full` and `rcu_nocbs` covering the lab CPUs (printed as a GRUB line,
needs a reboot), `irqbalance` banned from the lab CPUs, the `performance`
governor, turbo off, the NMI watchdog off (it takes a counter and adds a
tick), `numa_balancing` off, KSM off, THP recorded, swap off, perf access,
podman with `crun`, `runsc` registered, and for GPUs the NVIDIA container
toolkit with CDI, persistence mode and locked clocks. For another lab host,
run the installer with `--worker --hub-url nats://TOKEN@first-host:4222`;
for a second hub, `--hub --cluster-routes nats://first-host:6222`.

## CLI

```sh
isolab hub [--listen :4222] [--store DIR] [--token T] [--cluster-name N --routes URL,…]
isolab worker --cpus 4-15 [--slots N | --slot CPUS …] [--pool P …] [--label K=V …] [--backend podman] [--default-image IMG]
isolab doctor [--apply] | isolab inventory | isolab images build base|sage|cuda | isolab images list
isolab run [--image IMG] [--cpus N] [--memory-mb M] [--repeats R] [--policy strict] [--input PATH[:DEST]] -- CMD…
isolab submit SPEC.json | isolab status JOB | isolab wait JOB | isolab result JOB [--section S]
isolab logs JOB | isolab fetch JOB [PATH] --dest DIR | isolab cancel JOB | isolab jobs | isolab workers
isolab calibrate [--worker W] | isolab mcp | isolab schema
```

Connection: `--url` or `ISOLAB_URL`; `ISOLAB_NAMESPACE` to run several labs
on one hub.

## What was tested where

See [`TESTING.md`](TESTING.md) for the record of the end-to-end runs behind
this version: the unit suite (protocol, topology parsing on synthetic
sysfs trees, CPU planning, perf parsing, command construction, the fidelity
verdict), the integration suite against a real `nats-server` (submit,
claim, heartbeat, cancel, decline, lost-worker takeover, fenced result,
artifacts round trip), and the end-to-end run on a Linux host with root:
hub and worker on a two-NUMA-node, SMT-enabled Ubuntu 24.04 machine, the MCP
server on a separate machine, podman with `crun` and `runsc`, cpuset
partitions, eviction, IRQ moves, perf and PSI. Hardware counters, cpufreq
and turbo control are not exposed by that virtual machine, so those paths
were exercised as `unavailable` and remain to be confirmed on bare metal.

## Not built yet

* **Intel RDT / resctrl** for last-level-cache and memory-bandwidth
  allocation. The presence of `resctrl` is recorded; nothing is allocated.
* **GPU clock locking and MIG.** GPUs are inventoried, assigned exclusively
  and passed through with CDI; their clocks, utilisation and temperature
  are sampled, not controlled.
* **Windows and macOS workers.** Everything that measures is Linux. The
  `direct` backend degrades to tier D elsewhere so the suite can run on a
  laptop, and says so.
* **A web view** of the roster and queue. The MCP tools and CLI are the
  only clients.

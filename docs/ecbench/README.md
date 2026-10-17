# ecbench: one harness for every ECDLP method

`ecbench` measures Pollard rho (in its variants), baby-step giant-step,
the kangaroo, and index calculus **behind one interface, on the same
single-target workloads, in the same unit, under the same isolation, into
the same records**. A number from one method can therefore be read
against a number from another without first asking what each one
counted, where it ran, or which target it solved.

It is native code (`src/cryptanalysis/ecbench/`, `src/bin/ecbench.rs`),
as AGENTS.md requires of measurement work. It does not replace the
repository's accounting rules; it is the machinery that applies them to
every method at once.

```
spec.json ──ecbench run──► session/                ──ecbench verify──► audit receipt (replay certificate)
                           ├─ spec.json            ──ecbench compare─► comparison.json
                           ├─ plan.json            ──ecbench table───► the one-unit table (AGENTS.md §2)
                           ├─ host.json            ──ecbench db sql──► SQLite (docs/ecbench/schema.sql)
                           ├─ session.json
                           ├─ records.jsonl        ◄── one sealed record per execution, failures included
                           └─ exec/NNNNNN.stderr
```

## Contents

1. [Quick start](#1-quick-start)
2. [Where it sits](#2-where-it-sits)
3. [What a measurement is](#3-what-a-measurement-is)
4. [Names and identities](#4-names-and-identities)
5. [Confounders, and how each one is controlled](#5-confounders-and-how-each-one-is-controlled)
6. [Isolation levels](#6-isolation-levels)
7. [The unit, and what each method is charged](#7-the-unit-and-what-each-method-is-charged)
8. [Statistics](#8-statistics)
9. [Reproduction and the independent runner](#9-reproduction-and-the-independent-runner)
10. [The database](#10-the-database)
11. [Adding a method, a curve family or a factor base](#11-adding-a-method-a-curve-family-or-a-factor-base)
12. [Limits](#12-limits)
13. [GPU stage runner](#13-gpu-stage-runner)

Factor bases (build with `ecbench fb`, publish to Git/`store:`, run
`ic.pipeline` with a matching plug-in string):
[FACTOR_BASE_PROTOCOL.md](FACTOR_BASE_PROTOCOL.md).

## 1. Quick start

Rolling `main` release archives ship `ecbench` beside `ic` for Linux (glibc
and static musl), Apple silicon, and Windows. A source build is:

```bash
cargo build --release --bin ecbench
```

```bash
./target/release/ecbench plan --spec docs/ecbench/specs/smoke-prime.json
```

The exact binary-field correctness fixture pairs `ic.large_prime` with
the matched `rho.signed_frobenius` on the same explicit Koblitz curve and targets:

```bash
./target/release/ecbench plan --spec docs/ecbench/specs/smoke-binary-large-prime.json
```

The companion `ic large-prime` command imports and validates every curve,
base point, and challenge point in an external ECBench corpus before solving
one selected challenge. `--summands n-1` is literal and a run is refused when
its exact combination table would exceed `--max-states`. The JSON report binds
the raw manifest bytes with `manifest_sha256` and `manifest_bytes` and refuses
duplicate family ids:

```bash
ic --json large-prime --manifest corpus/manifest.json --family koblitz-n11-l4 \
  --summands n-1 --large-primes 2
```

```bash
./target/release/ecbench run --spec docs/ecbench/specs/smoke-prime.json --out /tmp/ecbench/smoke-1
```

```bash
./target/release/ecbench verify --dir /tmp/ecbench/smoke-1 --replay 8
```

```bash
./target/release/ecbench compare --dir /tmp/ecbench/smoke-1 --a rho-neg --b bsgs-neg
```

```bash
./target/release/ecbench table --dir /tmp/ecbench/smoke-1
```

```bash
./target/release/ecbench db sql /tmp/ecbench/smoke-1 | sqlite3 -bail ecbench.db
```

Wide registered curves whose coordinates do not fit the measured harness's
`u64` group types can still carry canonical factor-base inventories.  The
native P-256 Dickson builder writes the same object layout and FB1 preimage,
with its unavoidable wide integers labelled explicitly:

```bash
./target/release/p256_factor_base \
  --curve icv1-fp256-t89188191154553853111372247798585809583-f188c491 \
  --factor-base dickson-torus:depth=18 \
  --out /tmp/icv1-fp256-t89188191154553853111372247798585809583-f188c491.factor-base.json \
  --sql-out /tmp/icv1-fp256-t89188191154553853111372247798585809583-f188c491.factor-base.sql \
  --relation-length 17 \
  --verify
```

The standalone builder keeps the sealed `ecbench` measurement path unchanged.
Its `ecbench.factor_base_dump/v1-wide` inventory and SQL companion use the
same factor-base tables, but neither is an `ic.pipeline` measurement.

`ecbench methods` lists every method and its parameters; `ecbench host`
prints the host capsule; `ecbench claim build` turns an IC run and a
strong-rho run on a public target into a checked `vs_rho` claim (§9). On Linux `run` reserves a whole core by default
(`--cpus auto`); on macOS there is no affinity control, so it records
every run at L0 and says so. Operation counts are unaffected.

The skills in [`.agents/skills/`](../../.agents/skills) walk an agent
through a measurement (`ecbench-measure`), an independent run
(`ecbench-independent-runner`) and an extension (`ecbench-extend`).

## 2. Where it sits

| tool | language | what it is for | relation to `ecbench` |
|---|---|---|---|
| `tools/isolated_bench.py`, `src/bin/isolated_bench.rs` | Python, Rust | pin one command, keep others off its core, record conditions | `ecbench` takes the **same lock** (`/tmp/crypto-bench.lock`), refuses the same things and uses the same preflight defaults, natively, per measured child |
| ICMS (`tools/icms/`, `docs/ic/measurement/`) | Python | index-calculus configurations through producer adapters | `ecbench` keeps its vocabulary (windows, isolation levels L0–L3, env class, never-overwrite sessions, A/A arms) and extends it to every ECDLP method in native code. ICMS sessions stay valid evidence |
| `isolab/` | Python service | an independent lab worker: cgroup partitions, NUMA, IRQs, PSI gates | `ecbench isolab-job` writes the job; `ecbench` runs inside it with `--cpus inherit` (§9) |
| `taskq/` | Python service | the Kubernetes fleet queue | unchanged; a spec runs there the same way, as a command |
| boundary autolab, IC tournament | Python | their own frozen protocols | unchanged; `ecbench` does not redefine historical protocols |

AGENTS.md forbids Python for new measurement work, and that is the reason
`ecbench` is native. The Python tools above remain where they are and
keep their own evidence. New cross-method measurements use `ecbench`.

## 3. What a measurement is

The question is AGENTS.md's: **what does it cost to solve one previously
unseen target**, cold, with every phase charged, verified.

- **One target per workload.** A workload is one curve, one prime-order
  subgroup, one generator and one target point. `targets_per_curve: 8`
  makes eight workloads, eight rows, never one eight-target run. Nothing
  is shared between workloads.
- **Cold.** Each execution is a fresh process. Nothing a method builds
  (a jump table, a baby-step table, a factor base, a pair table)
  survives to the next execution, so its cost is inside every row.
- **Two windows.** `S` is `cold_end_to_end`: set-up, search, relation
  collection, linear algebra and the method's own internal verification,
  all charged. Beside it every record carries the **one-target online
  window** of the IC measurement rules (AGENTS.md "IC measurements"):
  from the first target-dependent computation to the verified recovery,
  with reusable target-independent set-up outside it.
  - For index calculus the window opens after the factor base and the
    oracle's tables. It splits into the five exclusive phases the
    `vs_rho` claim schema names: `target_query`, `target_PDP`,
    `target_relation_check`, `target_descent` and
    `target_recovery_check`, measured by the repository's phase clock
    (`ic_measurement`) and summing exactly to the window. A phase the
    method never enters is zero with a recorded reason. For example, the
    oracles return exact decompositions, so the classic pipeline has no
    separate relation check.
  - For rho the window is the target-dependent solve. For the strong
    reference (`rho.signed_frobenius_strong`), whose jump table `[a]G` is
    target-independent, the window starts at the first walk start
    `[c]G + Q`. The stages it includes are `walk`, `collision` and
    `recovery_check`. `rho.signed_frobenius_strong_legacy` retains the
    literal pre-2026-10-14 canonicalization only as a trajectory-identical
    implementation control; it has the same boundary and accounting and is
    not a second mathematical reference.
- **Verified, by someone else.** No solver receives the planted scalar.
  The runner checks `[k]G = Q` and `k = planted` in its own process, and
  `ecbench verify` checks it again from the files. A run that does not
  verify is a row with status `wrong_answer`, `exhausted`, `error`,
  `timeout` or `crashed`. It is kept, never dropped, and it makes any
  comparison involving it `incomplete`.
- **Planted or public targets.** By default a target is `Q = [k]G` for a
  derived `k` (`uniform_scalar_sha256_v1`), kept for the check above. With
  `"target_kind": "public"` a target is hashed to the subgroup
  (`hash_to_subgroup_v1`: an abscissa from SHA-256, lifted, the cofactor
  cleared), so nobody knows its logarithm and `[k]G = Q` is the only
  check. That is what the IC claim rules mean by a previously unseen
  public target, and `ecbench claim` accepts nothing else. The law is in
  the workload id, so the two kinds never share one.

### Boundaries

Every record carries the curve's generic floor in `S`,
`floor_s = √(π / 2A)`, with `A` the automorphisms a generic algorithm may
use: 2 (negation) on any curve, `2m` on a Koblitz curve over `GF(2^m)`
(signed Frobenius). The ratio `S / floor_s` is in every record and every
table. The *reference* is an arm of the session, normally the matched
rho: `rho.negation` on a prime curve, `rho.signed_frobenius` on a Koblitz
curve. `ecbench table` divides every arm by it.

## 4. Names and identities

Every object is named by a hash of what it is. A label never enters an
identity, and two objects that hash alike are the same object.

| object | id | hashed from |
|---|---|---|
| curve | ICV1 slug `icv1-…` ([ICV1](../curves/ICV1.md)), EC1 alias when the registry holds this generator | the curve model (AGENTS.md §11) |
| workload | `W` + 12 hex | ICV1 string, subgroup order, cofactor, generator, **target point**, target law, seed, index, `target_count: 1`, `cache_policy: cold` |
| method | `ECM1h` + 12 hex | registry id and every parameter, defaults written out |
| factor base | `FB1h` + 12 hex | curve slug, family, parameters, column count, point count, SHA-256 of the sorted point keys |
| spec | `ECS1h` + 12 hex | the spec without its label and description |
| host class | `ECBENV2h` + 12 hex (`ECBENV1h` before version 2) | the stable host facts (§5): topology shape, memory in whole GiB, so a reboot keeps the class |
| session | `ECBS1h` + 12 hex | spec, host and plan hashes, start time |
| record | `ECR1h` + 12 hex | the exact record line, with the id empty |
| run | `<method_id><workload_id>R<n>` | the repository's `{candidate}W{workload}R{number}` convention; `n ≥ 1` counts executions of the pair in the session |
| comparison | `ECC1h` + 12 hex | both arms' sessions and methods, the resampling request |
| claim candidate | IC1 label `IC1N…C…fb…PDP…RC…LA…TD…ISO0h` + 12 hex | the tournament's candidate record (`research/ic_candidate_tournament_20260915/identity.py`): the curve, the factor base's audited inventory, every IC stage resolved, the IC sources' hashes |
| claim workload | 12 hex | the tournament's workload record: EC1 curve id, the target point, its law and seed, the algorithm seed, the resource envelope |

Canonical JSON is sorted, compact and ASCII-escaped, as
`docs/curve-identities.md` hashes it. Floats never enter an identity. The
host-class prefix is deliberately not ICMS's `ENV1h`: the two hash
different fact sets. The two claim identities are hashed the way
`identity.py` hashes them (UTF-8, not ASCII-escaped) and are pinned
against its values in `claim.rs`'s tests, so a claim built here and a
tournament candidate with the same inputs carry the same id.

Ordinary factor-base dumps encode group keys as eight-byte big-endian
integers and coefficients as `u64`.  `factor_base_dump/v1-wide` preserves the
table and identity layout but stores curve integers and coefficients as
decimal strings.  Its family declares its fixed-width key encoding as an
identity-bound parameter; the P-256 Dickson family uses 33-byte big-endian
`((x+1)<<1)|sign` keys.  A wide dump may use an FB1 identity because the FB1
preimage is unchanged; it may not claim the ordinary v1 dump schema.

## 5. Confounders, and how each one is controlled

Every row of this table is something that has moved a figure in this
repository or in the literature without the algorithm changing. Each is
**controlled** (removed by construction), **gated** (checked per run, and
a run that fails the check loses its isolation level), or **recorded**
(kept beside the figure so a later reader can see it).

| confounder | how it moves a figure | what `ecbench` does | kind |
|---|---|---|---|
| Which target | rho's and BSGS's cost depend on the target; a lucky target flatters a method | targets are derived by SHA-256 from (curve, seed, index), every arm solves the same ones, and intervals resample workloads (§8) | controlled |
| Algorithm randomness | one rho walk is one draw of a wide distribution | one seed per (workload, round) is shared by every arm, so arms are paired; rounds vary it | controlled |
| Order and drift | heat, a filling page cache, a background job ramping up | arms interleave every round (`alternate` reverses order on odd rounds; `shuffle` permutes by seed); warm-up rounds are run and marked | controlled |
| One-time costs | page faults, lazy statics, allocator warm-up land on whoever runs first | every execution is a fresh process; nothing is cached across executions | controlled |
| Another benchmark | two timed jobs on one machine measure each other | the shared lock `/tmp/crypto-bench.lock`, the one `isolated_bench` takes; `--wait` queues | controlled |
| Other processes on the core | preemption, cache eviction | whole-core reservation; every movable thread moved off, and put back when the session ends however it ends (threads the evicted ones spawned meanwhile included, recycled thread ids excluded); the runner itself leaves the reservation; busy ticks on the run CPU not explained by the child are measured per run | gated (L1, L2) |
| SMT sibling | a neighbour on the sibling shares the pipeline and caches | a CPU is reserved only with all its siblings; under an inherited placement a sibling outside it costs the run L1; every sibling's busy ticks are checked per run, reserved or not | gated (L1, L2) |
| CPU 0 | the kernel's housekeeping and many IRQs land there | `auto` never picks CPU 0's core; a run on it loses L1 | gated (L1) |
| Migration | a process moved mid-run pays cold caches | pinned between fork and exec; the child reads its own mask back and reports the CPU it started and ended on | gated (L1) |
| NUMA | memory on the far node costs every miss | memory is bound to the run CPU's node with `set_mempolicy(MPOL_BIND)` in the child; the child reads the policy back with `get_mempolicy` and counts its anonymous pages per node in `/proc/self/numa_maps` after the solve | gated (L1, multi-node hosts) |
| Run-queue delay | time runnable but not running is not the method's | the child's `/proc/self/schedstat` across the solve; above 0.5 % of the solve, L2 is lost | gated (L2) |
| Preemption | each one costs cache and branch state | involuntary context switches per second from `wait4` | gated (L2) |
| Hypervisor steal | invisible to a process tick count | per-CPU steal ticks from `/proc/stat` around every run | gated (L2) |
| Memory pressure | reclaim stalls a run | PSI memory `some` total around every run | gated (L2) |
| A busy host at start | everything above, before the first run | a settle-window preflight (busy ticks per reserved CPU, PSI); refused unless `--allow-busy`, which records it and caps runs below L2 | gated (L2) |
| Frequency and turbo | the clock is the denominator of wall time | governor, `no_turbo`, `boost`, `isolcpus`, `nohz_full` and SMT control recorded per CPU; required for L3 | recorded, gated (L3) |
| Instruction count | wall time moves with frequency and contention, instructions retired barely do | user-space instructions and cycles across the solve from the CPU's counters (`perf_event_open`), where the host exposes a PMU; elsewhere the record says why not | recorded |
| Virtualisation | a guest cannot see or control most of the above | CPU `hypervisor` flag, DMI, `systemd-detect-virt`, `kern.hv_vmm_present` recorded; bare metal required for L3 | recorded, gated (L3) |
| Heterogeneous cores | Apple P and E cores differ by about 2× | performance levels recorded; with no affinity control the runs are L0 | recorded |
| Threads | four Rayon threads on one pinned CPU wait on each other | `RAYON_NUM_THREADS=1`, `OMP_NUM_THREADS=1`; child CPU time over wall above 1.05 loses L1 | controlled, gated (L1) |
| Environment leakage | an engine knob in the operator's shell changes the search | the child's environment is cleared and rebuilt: `PATH`, `HOME`, a pinned locale and thread counts, nothing else; its hash is in the session | controlled |
| Build | a debug build or another binary is another experiment | the binary's SHA-256, the source commit embedded at build time and its declared dirty state, any runtime worktree commit as a diagnostic only, `rustc -Vv`, compiler flags in the environment; on Linux each child execs `/proc/self/exe`, so a rebuild mid-session cannot swap the code (elsewhere a changed binary stops the session); a debug build loses L1 | recorded, gated (L1) |
| Host | the same binary differs by 17–20 % across hosts | the host capsule's stable facts hash to the class id; wall-clock ratios are formed only inside one session | recorded, controlled |
| What each method counts | one tool counts additions, another wall time, another instructions | one `CountedGroup` ledger for every method; native work the unit does not charge is counted as `*_uncharged` and makes the total a lower bound (§7) | controlled |
| Pricing drift | per-run measured ratios repriced identical runs by up to 8 % | IC native work is priced only at the repository's pinned ratios (`docs/ic/calibration.json`); a unit without one stays unpriced, never host-measured; an algebraic solver's wall-priced work is taken back out | controlled |
| Selective reporting | dropping a failure, rerunning until it looks good | every execution writes a sealed record before the next starts; the output directory is claimed with an atomic `mkdir` under the lock, so a session queued for a used directory fails instead of writing over it; an interrupted session is marked, never deleted | controlled |
| Interruption | a killed runner leaves narrowed masks behind and a child solving on an unprotected core | SIGINT, SIGTERM, SIGHUP and SIGQUIT are taken by one thread, which kills the child's process group; the session is marked `interrupted`, the eviction restored, and the runner exits `128 + signal`; the child dies with the runner (`PR_SET_PDEATHSIG`). A runner killed outright (SIGKILL, host loss) cannot mark the session, which keeps reading `running` as first written; `ecbench verify --allow-interrupted` audits such a session as interrupted (`interruption: never_finalised`), with every per-record check (seal, sequence, plan, answer, grading) still applied; its report also carries `records_unbound: true`, since `records.jsonl` was never hashed into `session.json` and a truncated or extended file is caught only by the unkeyed per-record seals, so the session is evidence of the interruption only (claims, bounds and challenges refuse any session that is not `complete`) | controlled |
| Tampering and transcription | a figure edited by hand, or copied wrong | every file is hashed into `session.json`; `verify` re-derives identities, the plan, the answers and the seals, and replays runs | controlled |

## 6. Isolation levels

A run earns the highest level whose checks, and every lower level's,
passed. Its record lists every check it failed as a blocker
(`"L2: 3 ticks of hypervisor steal on the run CPU"`), so a level is never
a bare label.

| level | requires |
|---|---|
| **L0 recorded** | the run executed under a saved host capsule. Enough for operation counts, which do not depend on contention |
| **L1 pinned** | a reservation (`--cpus` not `none`) inside this process's cpuset; a release build; no SMT sibling of the run CPU outside the reservation; the child's own affinity read back as exactly the run CPU, and the run CPU at its start and end; the run CPU not on CPU 0's core; the runner off the run CPU (under an inherited placement with no other CPU it may sleep on an idle sibling, where the sibling check sees it); no user thread left unmovable on the reserved CPUs (a non-root run usually fails this); child CPU time over wall ≤ 1.05; on a multi-node host, the memory policy read back as `bind:<node>` of the run CPU and at least 99 % of the solve's anonymous pages on that node |
| **L2 quiet** | L1, and: a quiet preflight; zero steal ticks on the run CPU; foreign busy time on the run CPU within max(2 ticks, 1 % of wall); the run CPU's siblings and the rest of the reservation within 2 ticks; run-queue delay ≤ 0.5 % of the solve; ≤ 50 preemptions per second; memory stall ≤ 0.5 % of wall |
| **L3 isolated** | L2, and the host configured for measurement: run CPU in `isolcpus` or `nohz_full`; `performance` governor; turbo known off; bare metal known |

A spec's `isolation_required` (default L2) is the level a **wall-clock**
figure needs to be admitted. Operation counts never need one.

`ecbench compare` also reports a **measurement quality score** from 0 to
100. This is an operational rubric, not a probability that the result is
correct and not a confidence interval. For each measured pair, an
unverified result gets zero points; a verified pair gets the lower of its
two isolation levels: L0=0, L1=40, L2=75, L3=90. The average is taken
over every paired attempt. Quality admission requires at least L2, even
if a diagnostic spec requests L0 or L1. A further 10 points are earned
only when an L2-or-higher A/A control has a 95% interval containing 1 and
wholly inside `[0.9, 1.1]`. The comparison records the level counts, verified and
admitted pair counts, 95th percentile foreign CPU busy fraction, 95th
percentile run-queue delay, and whether telemetry was missing. The
thresholds in the L0–L3 gate still decide whether wall time is admitted.
The separate bootstrap 95% interval measures sampling uncertainty over
the stated targets and seeds; a high quality score does not narrow it.

What hosts earn, as recorded so far:

- **macOS:** L0, because there is no affinity control.
- **A two-node QEMU guest** (Ubuntu 6.8 kernel, 2 sockets × 2 cores × 2
  threads, [`research/ecbench_numa_vm_20261002`](../../research/ecbench_numa_vm_20261002/README.md)):
  L2 on about 80 % of runs, the rest stopped at L1 by run-queue delay.
  Pinning and NUMA binding were confirmed there, with policy `bind:<node>`
  and every anonymous page on the bound node. An emulated guest reports
  no steal by construction, so its L2 does not certify quiet hardware.
- **A GitHub-hosted Linux runner, run as root** (AMD EPYC 9V74, Azure VM,
  4 vCPUs as 2 cores × 2 threads): L1 on 72 of 72 runs of the CI session,
  with every child's affinity read back as its run CPU. Its operation
  counts matched the macOS session of the same spec to the last digit, IC
  included. It is a VM, so L3 is out of reach. Each CI run's summary
  carries the current level and blocker counts.

L3 needs a lab host, which `isolab` prepares (§9).

## 7. The unit, and what each method is charged

```
S = total group-addition equivalents / √r
```

An addition and a doubling are one each. A scalar multiplication is
charged as the additions and doublings it performs. `S ≈ 1.25` is plain
rho's expectation `√(π/2)`; the curve's floor is `√(π/2A)` (§3).

| family | charged | counted, not charged (`*_uncharged`) |
|---|---|---|
| `rho.*` | jump-table set-up, walk starts, every step, cycle escapes, look-ahead additions, internal verification | canonicalisations (negation, signed Frobenius), Frobenius maps |
| `rho.signed_frobenius_strong` | group additions exactly; each scalar multiplication (jump table, walk start stride, candidate checks) at `1.5·log₂ r` additions, `ic_boundary::signed_frobenius_rho`'s convention | canonicalisations, partition hashes, distinguished-point table queries and inserts |
| `rho.signed_frobenius_strong_legacy` | identical to `rho.signed_frobenius_strong`; comparison-only implementation control | identical to `rho.signed_frobenius_strong` |
| `bsgs.*` | baby steps, the giant stride, giant steps | table inserts and lookups |
| `kangaroo.vow` | jump-table set-up, starts and restarts, every jump | table inserts and lookups |
| `claw.pair_table` | the known-log base's seed scalar multiplications, every table addition `P_i + F_b`, both query additions `Q − (F_k + F_l)`, the addition that rebuilds a hit's table sum; phases `base`, `table`, `search`, `recover` | Frobenius maps, canonicalisations, table inserts and probes |
| `ic.pipeline` | factor base, oracle set-up (pair tables), relation trials, linear algebra, verification, each a phase; native work (lookups, row operations, square roots, Artin–Schreier solves, …) at the pinned ratio where the repository has one | native work with no pinned ratio for this curve; algebraic-solver operations |
| `ic.pipeline_counted` | the same pipeline with the relation phase charge rebuilt from counted work and pinned conversions for deterministic SAT/PDP replays | algebraic-solver operations without a pinned conversion, including SAT conflicts; these runs remain lower bounds |
| `ic.large_prime` | cofactor projection, exact combination-table group additions, target probes, and relation witness checks; the external runner verifies the returned scalar | combination states, MITM lookups, large-prime merge operations, and modular row operations; the record is therefore a lower bound, not a speed claim |

The counted IC method has its own method identity so frozen SAT/PDP sessions
retain their exact parameter and accounting contract. It removes solver wall
pricing from the reported group-addition equivalents. Solver conflicts and
other work without a pinned conversion stay visible as unpriced counters;
the recorded `S` is therefore a lower bound on the complete pipeline cost.

The standing strong-rho implementation applies paired normal-basis maps and
paired coefficient reductions during signed-Frobenius canonicalization. The
[2026-10-14 optimization study](../../research/ecbench_rho_optimized_20261014/README.md)
records identical trajectories and counts over 240 matched pairs and a
complete-solve Callgrind ratio of 0.802155 against the same-binary legacy
path. That is implementation evidence in `callgrind.Ir`; it does not change
the generic-group floor or turn uncharged canonicalization into charged GAE.

`ic.pipeline` defaults to `linalg=incremental-gauss`, which stops when the
target scalar is pinned; the factor-base logs can still be underdetermined.
`linalg=incremental-gauss-full-rank` uses the same matrix and relation stream
but reports success only at rank `factor_base.columns + 1`. It records the
first target-pin rank, trial and relation as well as final rank. A trial-cap
exit after the early pin is an exhausted full-rank attempt. This distinct
parameter value changes the method identity and preserves old replays.
`linalg=incremental-gauss-full-rank-checked` follows the same relation
stream and stop rule, then verifies one nonzero-coefficient point per base
column by checking `[coefficient × column_log]G = [cofactor]P`. It charges
both scalar multiplications and records checked, missing and failing
columns in the verification phase. The existing full-rank mode keeps its
historical accounting and replay identity.

`compact-orbit-scan:columns=N,raw_x_cap=M` builds a Koblitz factor base
from a bounded raw-abscissa scan. Its cofactor projections, subgroup checks
and Frobenius eigenvalue checks are counted group operations. The factor-base
phase also records scan, lift, field, hash and vector-storage counters;
unpinned work keeps the reported `S` a lower bound. The
`mitm-frobenius-counted:m=3` oracle uses the existing folded table and
probe but adds canonicalisation, Frobenius-map, lookup, entry and
representative counts to `oracle_setup`. Use `mitm-frobenius:m=3` on the
same base as its accounting control.

`mitm-fiber:fiber=closed|lifts|blind`, `mitm-frobenius-fiber:fiber=closed|lifts|blind`
and `descent-algebraic-fiber:fiber=combined|lifts|blind` (with a `solver`) are
the cofactor-fiber-aware oracles of `ic_framework::fiber`
(`research/isogeny_conductor_gap_ic_20261010/`).  A subgroup target `T` is
decomposed as `Σ P_i = T + K` for any `K` in the rational `h`-torsion, which
the relation loop's `[h]` projection makes an exact row: `closed` holds
`P_i + P_j + K` for every `K` in the pair table (one probe, `h` times the
entries), `lifts` probes or descends every lift (`h` calls), `combined`
descends one system with the target abscissa free and the fiber polynomial
adjoined, and `blind` is the landed oracle.  The fiber the oracle found is
in `oracle_setup` (`fiber_size`, `fiber_complete`); with `h = 1` every mode
is the blind oracle.  All three take `fiber_gen=legacy|coset` and
`fiber_scale=d` (`research/isogeny_conductor_gap_round2_20261011/`).
`legacy`, the default, is the round-1 subgroup closure, kept so that
committed sessions replay.  It overcounts set-up when `E(F_q)[h]` needs
several generators.  `coset` stops the closure at the first multiple already
in the subgroup.  `fiber_scale=d` (coset only) uses the subgroup
`[d]E(F_q)[h]` instead of the whole fiber.

`koblitz-standard-subspace:dimension=d` with
`pdp3-koblitz:m=3,engine=inherited-f4|f6-ic,degree=D,node_budget=N` runs
PR #1333's Groebner decomposers unmodified inside `ic.pipeline`, with no
default for any parameter that changes the search. The base is every point
whose abscissa lies in `span(1, z, …, z^{d−1})`, folded by negation only. It
is the base #1333's workers used.

- **Charged:** every point addition the decomposition performs, including
  F6-IC's `geometric_group_additions` (support closure, residual arithmetic
  and witness replay), plus the framework's sign lift of each witness.
- **Counted, not charged:** the Boolean solver's word XORs, under the
  inherited-F4 adapters' unit string, so `S` is a lower bound and names
  that unit in `unpriced`. Reductions, splits, propagations and every F6
  gate counter ride in the record's `solver.extra`.

Both engines are exact on the same base and seed, so the two arms issue the
same queries until one exhausts its node budget. The base and oracle live in the
`ecbench` binary (`src/bin/ecbench/pdp3_koblitz.rs`) and reach `ic.pipeline`
through `methods::register_binary_plugins`. The library cannot name them:
frozen-source replays rebuild it against pre-F6 snapshots of
`koblitz_index_calculus.rs`, `Cargo.toml` is pinned by other frozen
evaluations, and the tournament admits no root build script. It also
refuses to run while a `KIC_*`, `F4_*` or `SOLVER_*` tuning override is set.

`prime-kummer:t=T[,g=G]` and `prime-isogeny:t=T[,seed=S]` are Petit,
Kosters and Messeng's tower bases for prime curves
(`ic_framework::prime_tower`): every point whose abscissa lies in a set `V`
of `2^t` elements that is the zero set of `t` composed degree-2 maps, both
signs, one column an abscissa.  The Kummer base needs `2^t | p − 1` and
costs `2^t` Legendre symbols and square roots.  The isogeny base,
`V = x(R' + ⟨Q⟩)` for `Q` of order `2^t` on an auxiliary curve `E'/F_p`,
needs no condition on `p`; `E'` is searched from the public seed, and the
search is **charged to the factor-base phase**: on every curve tried, a
baby-step giant-step over the multiples of `2^t` in the Hasse interval
filters out curves with `2^t ∤ #E'`, and the survivors' orders are computed
by Mestre's baby-step giant-step; with `Q`, `⟨Q⟩`, `R'`, `V` and the Vélu
chain's kernels, every addition and doubling on `E'` is charged in the unit
(a curve over the same field).  The phase records `aux_curves_tried`,
`aux_search_gae`, `aux_count_gae` and `aux_tower_build_gae` apart; the
Legendre symbols and square roots of point sampling join the phase's
`legendre_symbols` and `sqrt_solves` (unpriced, as for every prime base),
and the `2t` Vélu evaluations are `velu_evaluations_uncharged`
(`research/pkm_isogeny_tower_ecbench_20261015/`).

`rootfind-mitm:m=3` and `rootfind-pairs:m=3` are the root-finding
decomposition oracles of `research/pdp_root_finder_20261012` as plug-ins
(`ic_framework::rootfind`), on any base and both tracks. `rootfind-mitm`
builds the negation-folded pair table once (`full_table=1` for every
`i ≤ j` pair) and, per target, computes the points `R ∓ P_k` — one
representative `P_k` an abscissa — in blocks of 32 abscissae, probes each
and stops at the first hit;
`rootfind-pairs` keeps no table and tests every candidate `U ∓ P_j` of
every row point `U = R ∓ P_i` for membership in the base: the enumeration
floor, the memoryless control.

- **Charged:** every point addition, one each: the table build, the
  additions a target (`2A` for a refutation, the count `mitm:m=3` charges
  over the `2A` signed points; on a hit at most one block more than
  `mitm`'s early exit), the probe's confirming addition on a key hit, and
  every candidate of `rootfind-pairs`.  The cost is complete: on a binary
  curve the group computes a block with one shared inversion
  (`FastCurve::add_many`), which is how the charged additions are done,
  not work beside them; on a prime curve the default path inverts once an
  addition.
- **Counted, not charged:** table lookups and membership tests
  (`lookups`).  What the new ids change against `mitm` is therefore the
  arithmetic behind the same charged additions, which wall time sees and
  `S` does not; and the memoryless three-summand arm, which the harness
  had not had.

**The root-finding decomposition engine.**
`RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.md` replaces the Gröbner solve of the
decomposition step by root finding, and two oracle ids run it here on the
bases the Gröbner arms use, so the engines meet on one query stream
(`research/ecbench_rootfind_engine_20261012/`):

- `pdp3-koblitz:m=2|3,engine=rootfind` on `koblitz-standard-subspace` is
  `DecompositionStrategy::PairTable`, the strategy `koblitz_index_calculus_dlp`
  runs by default: the base's pair sums `P_i + P_j` tabulated once
  (`PairSumTable::build_full_within`, the summand-storing tier), a target
  decomposed by one table probe (`m = 2`) or one per `R − P_k` (`m = 3`).
  It takes no `degree` and no `node_budget` and refuses them. **Charged:**
  the table's `|F|(|F|+1)/2` pair additions, counted from the base's size
  because the batched builder carries no ledger (`|F|` of them doublings);
  for `m = 3` every `R − P_k` on the ledger until the first hit; and the
  re-addition of every hit in the group before it is returned.
  **Counted, not charged:** one `lookup` per probe (priced only where the
  curve has a pinned ratio, like every table oracle's) and the table's
  entry count in `oracle_setup`. The record has no `solver` block.
- `pdp3-koblitz:m=3,engine=rootfind-window` is the same table and set-up
  charge with the cyclic third-summand window of
  `PairSumTable::decompose_fast_window`: `⌈|F|/3⌉` probes a target, started
  at `x(R) mod |F|` (deterministic in the public point). It is a different
  search from full `rootfind` (and from inherited-F4), so it is a separate
  engine id and is not compared on query-stream identity. The native
  counters record `pair_table_points` and `scan_window_len`. Improve-loop
  evidence: `research/ecbench_rootfind_improve_20261013/`.
- `rootfind-subtract:m=2|3` on a prime curve is `prime_fast`'s
  direct-subtraction root finding: the Semaev `S₃` root in the third
  variable is `x(R ∓ F_i)`, so every candidate of a target comes from one
  Montgomery batched inversion over the base's abscissae, each looked up
  in the base (`m = 2`) or probed against the negation-folded pair table
  (`m = 3`). It is `subtract`'s search, candidate for candidate in the same
  order, so the two arms issue the same queries and find the same
  relations; it is charged as `subtract` is (one addition per candidate
  point) plus one re-addition per hit, so their `S` differ by one addition
  per relation and nothing else. The difference the unit does not see, one
  inversion a target against `2|F_x|`, is in the `inversions` counter, the
  record's field-operation tally and wall time.

Read a root-finding arm against a Gröbner arm on the resource vector and on
wall time: the root finder's `S` is a complete count of group additions,
while a Gröbner arm's `S` is a lower bound (its word XORs are counted and
not charged), and the ratio of a complete count to a lower bound is not a
speedup. Against `mitm`, whose generic pair table does the same search in
the framework's own arithmetic, it is an *engineering* comparison in the
sense of AGENTS.md §3: the same work in a faster constant. `ecbench` has no
prime-curve Gröbner oracle (the note's Buchberger baseline terminates only
at `|F| = 6`, too small for a DLP), so the prime side is compared with the
existing `subtract` and `mitm` oracles.

**Calibration.** The unit has been checked against theory in
[`research/ecbench_calibration_20261002`](../../research/ecbench_calibration_20261002/README.md).
Over 2 384 verified runs, preregistered and replayed, every generic
method's mean `S` reaches its analysed constant within its 95 % interval
at `r ≈ 2^25.4` (prime) and `r = 2^39` (Koblitz): plain rho
`√(π/2)`, negation rho `√(π/4)`, signed-Frobenius rho `√(π/4m)` (at 1.006
of it), BSGS 1.5, 4/3 and 1, kangaroo 2. Round 1 missed two
preregistered bands for sample-size reasons, and that record is kept.

Two consequences to read every table with:

- **BSGS's `S` omits its memory.** BSGS stores about `√r` points and
  probes a hash table every step. Neither is a group operation, so a
  BSGS row can sit below rho in `S` and still be the slower and far
  larger method at any interesting size. The records keep `max_rss_kib`
  and the lookup counts. This is the known time–memory trade, never a
  finding.
- **The pair claw's `S` omits its memory too.** `claw.pair_table`, ported
  from `aburan28/cryptanalysis#175`, holds `M = table_scale·√(r/n)` classes.
  Every base logarithm is known, so it is a generic method (a randomised
  BSGS on signed Frobenius classes) whose mean is about
  `(c + 1/c)/√n` at `c = table_scale`, at best `2/√n` against rho's
  `√(π/4n)`. It is in the `claw` family, not `ic`, and it is not an
  index-calculus candidate for a `vs_rho` claim.
- **A lower bound is marked.** `cost.lower_bound` is set whenever
  anything is unpriced, and comparisons carry `bounded: true`. A rho with
  the negation map is a lower bound by this rule (its canonicalisations
  are counted, not charged), exactly as in the rest of the repository.
  The ratio of two such lower bounds is **not** a bound on the ratio of
  full costs. Read it as a counted-work diagnostic until both arms have
  complete, calibrated costs; do not infer an end-to-end speedup from it.

- **Full pricing is opt-in.** `ic.pipeline_full`,
  `rho.signed_frobenius_full` and `rho.signed_frobenius_strong_full` run
  their base method unchanged, then price every native counter (field
  multiplications, squarings and inversions, Artin–Schreier solves, hash
  probes, canonicalisations, row operations, word XORs) from
  `src/cryptanalysis/ecbench/full_prices.json`, one host's measured ratios to an affine
  addition, written by `ecbench price-units --curve a:n …`.
  `src/cryptanalysis/ecbench/full_price.rs` gives each counter an explicit
  rule. A counter without a rule stays unpriced, so the record is still
  marked a lower bound. The pinned `docs/ic/calibration.json` and every
  other method are unchanged. `ecbench trend --compare …` fits how a paired
  ratio moves with `ln r` and `ln n` across curves. First use:
  `research/ecbench_frobenius_orbit_vs_rho_20261014/`.

`ecbench resources --dir SESSION --out resources.json` reports, for every
measured attempt including failures, the native integer counters by phase,
charged GAE and clocks by phase, the target's exclusive online window,
method counters, fork-to-reap and solve time, CPU time, peak resident memory,
and solve-wide instructions and cycles when the host exposes a PMU. A PMU
sum is null unless **every** measured run in that arm has the counter.
These are separate units: a SAT conflict, table lookup and byte of memory
do not have universal GAE prices. The resource vector shows where the
counted `S` is incomplete. A full runtime comparison still requires
matched L2/L3 runs and the A/A noise gate in §8; the Mac's L0 timings
remain descriptive.

## 8. Statistics

- **The mean, not the median.** The floor `√(πr/2A)` is an expected
  value, so the comparable statistic is the mean `S`.
- **Paired, and a ratio of totals.** Runs are matched by (workload,
  round). Within a session the two arms of a round share an algorithm
  seed; across two sessions of one spec (a baseline binary and a
  candidate) they do too, and `same_seeds` says so. The ratio is
  `Σ S_B / Σ S_A` over the matched verified pairs, AGENTS.md §8's
  `baseline_total / candidate_total` in `S`. It is never a mean of
  ratios, which is biased whenever the arms spread differently (rho's
  cost varies with its seed, BSGS's does not).
- **Two-stage (cluster) bootstrap.** Intervals resample *workloads*, then
  pairs within each, and recompute the same ratio of totals. This
  matters: BSGS's cost is fixed by its target, so all its variance is
  between workloads, and an interval that held the workloads fixed would
  be a point. With one workload only the within-workload interval exists,
  and the comparison says `ci_method: within`. Runs without a partner
  make the comparison `partial`.
- **Per curve.** A ratio is reported per curve (one size each) as well as
  pooled. A claim about scaling reads the per-curve rows and fits
  exponents over at least four sizes (AGENTS.md §5).
- **An A/A in operations is exactly 1.** Two arms running one method on
  the same seeds do identical work, so the ops A/A is a check of
  determinism, not of noise. The noise in operations comes from targets
  and seeds, and the intervals above already carry it. The A/A arm is
  for wall time, where its interval is the session's noise floor.
- **Wall time is gated.** It is compared only within one session, only
  over pairs where both runs reached `isolation_required`, only with at
  least five such pairs. It is reported as the median per-pair ratio with
  a bootstrap interval, set beside the A/A interval (which needs five
  admitted pairs of its own), and marked `outside_noise` only when it
  excludes 1 and does not overlap the A/A interval. Anything less prints
  as `descriptive`, never as a result.
  This `compare` wall field uses the whole solve, including IC's reusable
  setup. For the primary one-target IC question, read the separate online
  windows and `ecbench claim`'s same-point interval; a cold wall ratio is
  never an online speedup.

A speedup in the sense of AGENTS.md §8 is still
`baseline_total_operations / candidate_total_operations`, over the whole
cold method with every phase priced. `ecbench compare` reports the ratio.
Classifying the change (advance, engineering, relabelling, accounting)
stays with the author.

## 9. Reproduction and the independent runner

**Exact replay.** Operation counts depend only on the code, the workload
and the seed, never on the host. `ecbench verify --replay N` re-executes
N measured runs (`--replay-all`, every deterministic one) and requires
the same answer, total, phase counts, counters, unpriced work and factor
base, bit for bit. Run on another machine, it is an independent check of
the figure. Records written before 2026-10-05 by an `ic.pipeline` arm
whose solver work was wall-priced (the descent-algebraic arms of
`research/ecbench_yield_sweep_20261004/sessions/koblitz`) carry gae
figures rounded to the binade of that wall term, because the term was
removed by subtraction; the audit reproduces that rounding from the
solver term and pre-removal total the record itself carries, reports the
replay as `identical (legacy gae rounding, …)`, and new records carry the
exact figure. The audit also recomputes every derived figure from the
record's own counts (`S = total/√r`, the floor, the ratio, the
lower-bound flag) and, for a session graded under the current rules,
regrades every run from its recorded observations. A figure edited by
hand and resealed is caught. The receipt names every file by SHA-256,
and its own SHA-256 is the replay certificate a claim cites. CI
re-audits every committed session under `research/ecbench_*/sessions`
on a Linux x86-64 runner with every run replayed.

**A method's counts are frozen by its id.** Changing what a registered
method counts breaks the replay of every committed session that used it,
and CI fails. That is by design: register a new method id
(`kangaroo.vow2`) instead. The records' hashes are unkeyed, so someone
who rewrites every observation and reseals every hash could forge a
quieter run. Git history and the replay on another host are the
witnesses against that, as with ICMS.

**Independent runner (isolab).** `ecbench isolab-job` writes an
`isolab.job/v1` that runs a spec on a lab worker:

```bash
./target/release/ecbench isolab-job --spec docs/ecbench/specs/smoke-prime.json --binary target/release/ecbench > job.json
```

isolab places the job on whole cores of one NUMA node (`smt: isolate`,
`numa_node: single`), moves threads and IRQs off, sets frequency policy,
gates on PSI, and runs `ecbench run --cpus inherit`. Inside, `ecbench`
pins its children to the placement it was given and grades every run
from its own observations, as anywhere else. After the clock stops, the
job's `verify` step runs `ecbench verify --replay 3 --exit-code` on the
worker. The session comes back as hashed artifacts.

- `--binary` (preferred) uploads a binary built for the worker's
  architecture. The job needs no network and the binary's hash is its
  identity.
- `--git-url URL --commit SHA` builds on the worker. The repository does
  not track `Cargo.lock`, so dependencies are resolved there and the job
  needs network. The job labels this, and every record still carries the
  binary hash a comparison checks.

### Claims

A **`vs_rho` claim** is the report `docs/ic/boundary_targets.json`
defines for the primary IC question: one index-calculus run against one
strong-rho run (`rho.signed_frobenius_strong`) on one public target,
compared on their one-target online windows (§3), keyed by
`(candidate_id, workload_id, run_id)` in the repository's IC1 convention.
The ledger calls its one-run ratio `online_speedup`; a checker `PASS`
certifies the report's pairing, timing fields and provenance conditions,
not a statistically established runtime gain. That requires at least five
admitted matched pairs, an A/A noise interval, and a paired 95% interval
excluding no improvement (§8). The cold end-to-end operation comparison
is a separate requirement.

```bash
./target/release/ecbench claim build --dir SESSION --ic ic --rho rho-strong --workload W... --out claim.json --exit-code
```

```bash
./target/release/ecbench claim check --report claim.json
```

`build` takes the IC and rho records of one workload and round (by
default the first measured round both arms ran) and assembles every
field the schema requires, from the session's own files:

- **The candidate identity** is the tournament's
  (`identity.py`, ported in `claim.rs` and pinned against its values):
  the curve record and EC1 id, the factor base rebuilt and inventoried
  under the cofactor convention (usable base `[h]B` without the identity,
  columns its sign-and-Frobenius orbits), every stage of the method
  resolved, and the IC sources hashed into the binary at compile time. A
  claim is therefore built only by the binary that measured the session.
  The claim also carries the complete canonical candidate and workload
  records beside their hashes, so a reader can audit the compact IDs.
  `ecbench`'s Koblitz bases fold the raw lifted points, so a raw orbit
  whose cofactor image is the identity (the 2-torsion point above
  `x = 0`) is a column with logarithm zero; the record discloses such
  columns under `factor_base.construction.solver_columns` rather than
  counting them.
- **The workload identity** needs a public target: a planted one is
  refused, as `identity.py` refuses it.
- **The windows** are the records' online windows, the IC one split into
  the five exclusive phases, both arms' start and stop events named.
- **Wall isolation.** A `vs_rho` wall comparison requires both runs to
  earn at least L2, even if the session spec lowers `isolation_required`
  for diagnostic runs. A stricter L3 spec still requires L3. The builder
  writes the effective floor and the native checker validates it, so an
  L0 session cannot become a passing wall claim by supplying a replay
  receipt. Operation-count comparisons remain available at every level.
- **The envelope** is the session's (binary, host class, CPU, NUMA node,
  one thread, timeout), identical for both arms by construction.
- **Independence.** The schema requires `independent_validation` and a
  replay certificate for each arm. The certificate is the SHA-256 of an
  `ecbench verify --replay-all` receipt that reproduced both runs;
  `--independent-receipt FILE --pointer WHERE` supplies it. The auditor
  must have a different measurement environment **and** a different
  hardware profile class. The latter omits mutable compiler, kernel and
  governor facts: changing `rustc` on `PATH` on one Mac was enough to
  change `ECBENV2` while staying on the same machine. A receipt from the
  session's own hardware class, one lacking the new auditor hardware
  class, one for other bytes, or one that did not reproduce both runs is
  refused. These self-reported fingerprints cannot by themselves attest
  physical independence; retain the external runner's job and custody
  receipt with the pointer. Without an admitted receipt the
  claim is built with null certificates, the checker fails it on exactly
  those fields, and the verdict says it is not yet a claim.
- **Late receipt attachment.** If the measuring binary is no longer
  available when the other-host receipt arrives, `ecbench claim attach`
  reconstructs the entire previously saved diagnostic report from the frozen
  session and the current claim sources. It requires exact equality before
  adding the receipt; a changed source, target, time, candidate, or other
  report field is refused. The new report retains the original binary hash.
  `claim build` still requires the measuring binary itself. Archive both the
  original report and the attached one.
- **Admissibility.** The verdict states the speedup `rho / IC` only when
  both runs earned the spec's isolation level; otherwise it is marked
  descriptive. The levels are in the report either way.

`check` is `boundary_autolab.py claim-check --stage vs_rho` natively: the
same required fields and aliases, the same validations in the same order,
the same output object. CI builds a claim from
`docs/ecbench/specs/claim-koblitz.json` and requires the two checkers to
agree on it.

## 10. The database

[`schema.sql`](schema.sql) is the schema, compiled into the binary.

```bash
./target/release/ecbench db sql research/ecbench_smoke_20261002/sessions/* | sqlite3 -bail ecbench.db
```

The database is an **index**, rebuilt from session directories,
comparison files, factor-base dumps and claim files; it is never the
record of a measurement. Loading is idempotent in one way only:
an insert tolerates a conflict on its table's primary key when it is the
same row again, and nothing else. A second row claiming an existing run
id is a UNIQUE error. A key arriving with a different identity (a curve
slug on another ICV1 string, a workload, method, factor-base or host id
on another hash) is stopped by the schema's identity triggers.

| table | holds |
|---|---|
| `curves`, `curve_representations`, `curve_constructions` | ICV1 slug, ICV1 string, family, sizes, floor; the generator measured and its EC1 alias; the constructor calls that produced it |
| `workloads` | one curve, generator and target, with the known answer |
| `algorithms` | every method configuration (`ECM1h…`): rho, BSGS, kangaroo, IC |
| `factor_bases`, `factor_base_points` | every factor base a run used, by `FB1h…`; with its points when an `ecbench fb` dump is loaded |
| `hosts`, `sessions`, `arms` | the host class and capsule, the session's build and reservation, its arms |
| `runs`, `phases`, `run_counters`, `isolation_blockers` | one row per execution, its phases, its counters, and why it did not reach a higher level |
| `run_solver` | the decomposition solver's statistics for an algebraic or SAT IC run (the record's `solver` block): the system's variables, equations and semi-regular degree, the solving degree reached, Macaulay rows, columns, degree and rank where the solver reports them, SAT variables, clauses, conflicts, decisions and propagations, and every solver-specific counter verbatim. Absent for table oracles. Informational: outside `S` beyond the relations phase's `solver_*` counters, and outside the replay comparison |
| `comparisons` | every saved comparison |
| `claims` | every loaded `vs_rho` claim, keyed by its IC1 `run_id`, with the checker's status at load time; a claim rebuilt with an independent receipt replaces its earlier form |

Views: `arm_workload_summary` (the AGENTS.md §2 table, per workload),
`method_by_curve` (each method's mean `S` and ratio to the floor per curve,
across sessions), `ic_phase_split` (each IC run's phase shares with its
factor base), `ic_yield` (the yield ledger: one row per IC run with its
curve, target, factor base, oracle, trials, relations, `yield` =
relations / trials, lookups, lift failures, matrix rank and rows, and the
`run_solver` statistics; `ic.pipeline`'s `relations` and
`ic.shared_rank`'s `hits` are read as the same count).

```sql
SELECT fb_family, oracle, count(*), round(avg(yield), 4), round(avg(lookups * 1.0 / relations), 1)
FROM ic_yield WHERE status = 'verified' GROUP BY fb_family, oracle ORDER BY avg(yield) DESC;
```

The same ledger, joined with the curve registry and the tournament's
rounds, is the Yield view of the lab browser (`docs/browser/`,
AGENTS.md §7c). A yield is a stage diagnostic: it ranks bases and oracles
on the same targets, and only the whole-pipeline `S` decides speed.

```sql
SELECT method, curve_slug, verified_runs, round(mean_s, 3), round(mean_ratio_to_floor, 2)
FROM method_by_curve ORDER BY curve_slug, mean_s;
```

## 11. Adding a method, a curve family or a factor base

- **A method:** an entry in `methods::registry()` and an arm in
  `methods::solve`. Once a session that uses it is committed, its counts
  are frozen (§9). A changed algorithm is a new id. Charge every group operation through `CountedGroup`,
  count everything else under a name ending `_uncharged`, never read the
  planted scalar, and return the recovered value unverified. Give it a
  unit test that recovers known logarithms on a prime and a Koblitz curve.
  The `ecbench-extend` skill walks through it.
- **A curve construction:** a variant of `workload::CurveSpec` that calls
  an existing constructor. Register any new curve before citing its slug
  (AGENTS.md §11). `conductor_gap` (`{"kind": "conductor_gap", "ell": 263,
  "t": 241, "level": "crater" | "floor", "seed": 1}`) builds CM pairs over
  `GF((t² + 19ℓ²)/4)` of prime order `q + 1 − t` whose endomorphism
  conductor is 1 at the crater and the prime `ℓ` at the floor, certified
  without computing an isogeny; it is the construction behind
  `research/conductor_gap_cm_pairs_20261010/`.
- **A factor base or oracle:** add it to `ic_framework` as a plug-in;
  `ic.pipeline` reaches it through its name. Dump it with `ecbench fb` to
  store its points.

## 12. Limits

Stated so that nothing here is read as more than it is:

- **A claim's identity covers Koblitz curves with `5 ≤ n ≤ 61`.** The IC1
  candidate identity is the tournament's (`identity.py`), whose adapter
  admits those curves and no others; a prime-curve IC run has an online
  window and no claim. A claim's independence rests on the receipt's
  recorded host class and file hashes: `ecbench` checks that the receipt
  is a passing audit of these exact files, run elsewhere, that
  reproduced both runs, and a reader checks the receipt itself at the
  pointer the claim cites.
- **Wide curves run two methods.** Prime curves need `p < 2^62`. Binary
  curves run every method for `m ≤ 62`. Koblitz curves with `64 ≤ m ≤ 127`
  (`koblitz_wide`, over the registry's modulus, so the m = 83 gate curve
  `icv1-f2m83-tm6151469093347-debefd74` of AGENTS.md §8a builds) run only
  `rho.signed_frobenius_strong` and `claw.pair_table`. Both are the one-word
  code made generic over the word, and the committed sessions' replays pin
  that the one-word walk did not change. Every other method refuses a wide
  curve. A subgroup order past `2^64` is written as a decimal string
  (`canonical::compat_u128`) and gets the `_v2` target laws
  (`uniform_scalar_sha256_v2`, `hash_to_subgroup_v2`). Of the registered
  wide Koblitz curves only m = 83 and m = 97 have a prime `#E/h`, so the
  cheapest wide solve is m = 83 itself: about `2^37` strong-rho steps,
  some 20 single-core hours at the 0.49 µs per step measured on Apple
  silicon. The claw's cost-minimising table there is about `2^37` entries,
  out of reach of one host's memory. `p256_factor_base` only constructs
  and inventories a factor base on the registered wide prime curve; it
  does not make that curve runnable by `ic.pipeline` and produces no
  operation-count or speed claim.
- **NUMA binding has met a two-node kernel, not two-socket hardware.** In
  a two-node QEMU guest the policy read back as `bind:<node>` and every
  anonymous page sat on the bound node, and that test found and fixed a
  wrong read-back (`research/ecbench_numa_vm_20261002`). Remote-memory
  latency, the confounder the binding exists for, is not emulated, so a
  two-socket host remains the test of its effect.
- **Hardware counters need a PMU.** Instructions and cycles are read
  through `perf_event_open` on Linux. macOS, the QEMU guest and most cloud
  VMs expose no hardware counters, and their records carry the reason
  instead of a count.
- **A separate Callgrind solve count is available on Linux x86-64.** Export
  a frozen child input with `ecbench profile-input --spec SPEC --seq N`, run
  `ECBENCH_CALLGRIND_SOLVE=1 valgrind --tool=callgrind --cache-sim=no
  --branch-sim=no --callgrind-out-file=PREFIX ecbench exec < INPUT.json`,
  then read it with `ecbench callgrind-ir --prefix PREFIX`. The parser counts
  every numbered part after the pre-solve reset through the post-solve reset,
  including internal phase dumps, and rejects missing markers. This is a
  complete user-space **implementation** instruction count for the one solve;
  it has its own `callgrind.Ir` unit. It is not a PMU count or an isolated
  native wall-time result. Keep its binary, Valgrind version, CPU feature
  dispatch, input, child output, and raw parts together. The
  [n37 protocol](../../research/ecbench_callgrind_solve_20261004/PROTOCOL.md)
  freezes a same-target IC/rho census using this path.
- **A target-only Callgrind interval is opt-in for the shared-rank IC and
  strong signed-Frobenius rho.** Set `ECBENCH_CALLGRIND_TARGET=1` alongside
  `ECBENCH_CALLGRIND_SOLVE=1`, then run `ecbench callgrind-online-ir
  --prefix PREFIX`. The parser requires one target interval inside the
  complete solve and reports pre-online, online and post-online Ir with
  their sum checked against the complete solve. It rejects profiles
  missing either boundary. This is simulated instruction attribution,
  not an isolated online wall claim; the
  [preregistered n37 gate](../../research/ecbench_n37_online_ir_20261004/PROTOCOL.md)
  defines its first comparison.
- **Rho's distinguished-point table is not counted.** Its stores happen
  once per distinguished point, a vanishing fraction of steps. BSGS and
  the kangaroo count their table operations.
- **One thread.** Parallel collision search is out of scope until a
  multi-thread run can be graded (a parallel change must be measured at
  one thread and at many, AGENTS.md §10).

## 13. GPU stage runner

`ecbench gpu-stage` runs the repository's existing CUDA ECC driver
(`gpu/ecc/bench`) with a device self-test before a bounded benchmark. Build
the CUDA binary on a GPU host with `make -C gpu/ecc bench`, then run, for
example:

```bash
ecbench gpu-stage --binary gpu/ecc/bench --out /tmp/ecbench-gpu-rho-1 \
  --mode rho --timeout-seconds 600 -- --iters 256 --fold 2
```

The command creates a new output directory and writes `selftest.stdout`,
`selftest.stderr`, `benchmark.stdout`, `benchmark.stderr`, and
`gpu-stage.json`, plus a private byte copy of the binary. The receipt binds
the binary SHA-256, `CUDA_VISIBLE_DEVICES`, exact argument
vector, exit code, timeout, timestamps, and output hashes. A failed device
self-test prevents the benchmark. The directory is never overwritten.
`--solve` is excluded from this stage command so a large-group solver does
not run indefinitely under a throughput request.

This entry point runs real CUDA code, but its current fixed curves and
targets do not share `ecbench run`'s workload identity or independently
replayed answer and operation ledger. The stage receipt therefore stays
separate from `records.jsonl`, `ecbench compare`, and the L0–L3 CPU score.
It has no GPU contention sensor or admitted CPU/GPU interval yet. To make
a GPU ECDLP method fully comparable, the CUDA driver must accept the
frozen EC Bench workload and seed, emit exact phase counts and recovered
answer, expose GPU device and contention telemetry, and pass independent
answer replay under the same record contract. The specific acceptance test
and remaining device-validation gap are tracked in
[the GPU matched-workload problem](../../problems/ecbench-gpu-matched-workload-contract.md).

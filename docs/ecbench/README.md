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

## 1. Quick start

```bash
cargo build --release --bin ecbench
```

```bash
./target/release/ecbench plan --spec docs/ecbench/specs/smoke-prime.json
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
    `recovery_check`.
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
| Build | a debug build or another binary is another experiment | the binary's SHA-256, the commit and dirty state of the worktree the binary sits in, `rustc -Vv`, compiler flags in the environment; on Linux each child execs `/proc/self/exe`, so a rebuild mid-session cannot swap the code (elsewhere a changed binary stops the session); a debug build loses L1 | recorded, gated (L1) |
| Host | the same binary differs by 17–20 % across hosts | the host capsule's stable facts hash to the class id; wall-clock ratios are formed only inside one session | recorded, controlled |
| What each method counts | one tool counts additions, another wall time, another instructions | one `CountedGroup` ledger for every method; native work the unit does not charge is counted as `*_uncharged` and makes the total a lower bound (§7) | controlled |
| Pricing drift | per-run measured ratios repriced identical runs by up to 8 % | IC native work is priced only at the repository's pinned ratios (`docs/ic/calibration.json`); a unit without one stays unpriced, never host-measured; an algebraic solver's wall-priced work is taken back out | controlled |
| Selective reporting | dropping a failure, rerunning until it looks good | every execution writes a sealed record before the next starts; the output directory is claimed with an atomic `mkdir` under the lock, so a session queued for a used directory fails instead of writing over it; an interrupted session is marked, never deleted | controlled |
| Interruption | a killed runner leaves narrowed masks behind and a child solving on an unprotected core | SIGINT, SIGTERM, SIGHUP and SIGQUIT are taken by one thread, which kills the child's process group; the session is marked `interrupted`, the eviction restored, and the runner exits `128 + signal`; the child dies with the runner (`PR_SET_PDEATHSIG`) | controlled |
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
| `bsgs.*` | baby steps, the giant stride, giant steps | table inserts and lookups |
| `kangaroo.vow` | jump-table set-up, starts and restarts, every jump | table inserts and lookups |
| `ic.pipeline` | factor base, oracle set-up (pair tables), relation trials, linear algebra, verification, each a phase; native work (lookups, row operations, square roots, Artin–Schreier solves, …) at the pinned ratio where the repository has one | native work with no pinned ratio for this curve; algebraic-solver operations |

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
- **A lower bound is marked.** `cost.lower_bound` is set whenever
  anything is unpriced, and comparisons carry `bounded: true`. A rho with
  the negation map is a lower bound by this rule (its canonicalisations
  are counted, not charged), exactly as in the rest of the repository.

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
the figure. The audit also recomputes every derived figure from the
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
- **The envelope** is the session's (binary, host class, CPU, NUMA node,
  one thread, timeout), identical for both arms by construction.
- **Independence.** The schema requires `independent_validation` and a
  replay certificate for each arm. The certificate is the SHA-256 of an
  `ecbench verify --replay-all` receipt made on a host of another class
  that reproduced both runs; `--independent-receipt FILE --pointer WHERE`
  supplies it, and a receipt from the session's own class, for other
  bytes, or that did not reproduce both runs is refused. Without one the
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
  (AGENTS.md §11).
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
- **Word-size curves only.** The counted group types hold `GF(p)` with
  `p < 2^62` and `GF(2^m)` with `m ≤ 62`. The m = 83 confidence gate
  (AGENTS.md §8a) needs a wide-word group type before `ecbench` can run
  it. `p256_factor_base` only constructs and inventories a factor base on the
  registered wide curve; it does not make that curve runnable by
  `ic.pipeline` and produces no operation-count or speed claim.
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

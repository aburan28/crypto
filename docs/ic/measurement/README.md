# ICMS v1: the index-calculus measurement standard

ICMS makes two index-calculus (IC) measurements comparable or says why they
are not. It does not replace the repository's accounting rules. It is the
machinery that enforces them:

- **Operation counts** stay the primary metric (AGENTS.md §2, §6).
- **Wall time** goes through `tools/isolated_bench.py` (§10).
- **The primary comparison** uses one target (AGENTS.md "Primary comparison uses one target").
- **Results** are keyed by candidate, workload and run (cryptanalysis AGENTS.md).

Everything here is implemented in [`tools/icms/`](../../../tools/icms) and
tested by `python3 -m unittest discover -s tools/icms/tests`. The published
page is [`docs/ic-measurement.html`](../../ic-measurement.html).

## Why this exists

Before ICMS, an IC figure in these repositories could change meaning between
two tables while keeping its name. The audit behind this standard found,
with file paths, in [`audit.json`](audit.json):

- **Several quantities called S.** `ic bench` counts group-addition
  equivalents with pinned ratios. `ic price` divides wall time by a batched
  addition measured on the host. The CI gate uses callgrind instructions.
  One leaderboard heads all of them "S".
- **Wall time inside "operation counts".** Every shipped polynomial solver
  in the IC framework is priced from host wall time (`priced_by: measured`),
  so an algebraic row's S moves with the machine, its load and its thread
  count.
- **No shared window.** Cold end to end, one-target online, whole process,
  amortised over k targets: one thread reported 0.286, 5.03 and 110.79 for
  the same n = 41 pipeline under three windows.
- **Different factor-base sizes under one name.** `abscissae` is
  "distinct x with a point" in one producer and "every x the construction
  admits" in another.
- **A thin host record.** Nothing recorded microcode, governor, turbo, SMT,
  sysctls, the kernel command line, steal time or the toolchain. Host
  effects of 17 to 20 percent and A/A spreads of ±20 percent were measured
  and left unexplained.

## The objects

```
spec (YAML)  ──icms run──►  session/                 ──icms compare──►  comparison.json
                            ├─ session.json          (plan, reservation, preflight, hashes)
                            ├─ capsule.json          (the host)
                            ├─ records.jsonl         (one record per execution)
                            ├─ specs/                (byte copies of the specs used)
                            └─ exec/NNNN/            (stdout, stderr, producer files)
```

| object | schema | identity |
|---|---|---|
| spec | `icms.spec/v1`, [`schema/spec.v1.json`](schema/spec.v1.json) | `ICS1h` + 12 hex of the canonical spec without label, description, role or notes |
| workload | the cryptanalysis AGENTS.md workload record | `W` + 12 hex: curve, targets, law, seed, cold target count |
| host class | `icms.capsule/v1` | `ENV1h` + 12 hex of the stable facts in `env_class()` |
| record | `icms.record/v1` | `ICR1h` + 12 hex of the record |
| session | `icms.session/v1` | `ICSESS1h` + 12 hex |

Canonical JSON is sorted, compact and ASCII, as both repositories already
hash. A float never enters an identity: the spec's identity view renders it
as its shortest round-trip decimal string.

## 1. The specification

One YAML file describes one configuration on one instance, and the
measurement it runs under:

```yaml
icms: icms.spec/v1
label: "K_1/GF(2^23): binary subspace dim 8, pair-table m=2"
role: candidate                     # candidate | baseline | reference
instance:
  curve: {regime: koblitz, degree: 23}
  workload: {targets: 1, law: known_answer, seeds: [123212651130]}
factor_base:
  family: binary_subspace           # a name from registry.json
  params: {dimension: 8}            # the family's own parameters
  quotient: [negation]              # folded symmetries; [] is one column per point
  large_primes: {mode: none}        # none | single | double, with bound and cycle finding
  # partition: {parts: 2, scheme: per_summand, assignment: [0, 1]}
decomposition:
  arity: 2
  method: mitm                      # subtract | enumerate | pair_table | mitm | mitm_frobenius | algebraic | symmetrised
  # encoding: expanded_semaev       # | symmetrized | chained_s3 | riemann_roch | orbit_coordinates
  # splitting: none                 # | chained_s3 | balanced_resultant | mitm_table
  # solver: {name: sat-cdcl, options: {sat_xor_encoding: native, sat_conflict_budget: 200000}}
relations: {collector: walk, stop: first_log}
linear_algebra: {method: incremental_gauss}
reference: {rho: rho.measured_matched, floor: floor.generic, rho_runs: 8}
accounting: {unit: crypto.S.gae_pinned}
measurement:
  window: cold_end_to_end
  repetitions: 5
  warmup: 1
  interleave: abab
  isolation_required: L2
  timeout_seconds: 120
execution: {adapter: crypto.ic_bench, threads: 1, cpus: 1}
```

**It is flexible where configurations differ, and strict where measurements
must agree.**

- **New variants need no schema change.** A family, unit, window or reference
  is a name in [`registry.json`](registry.json), so a new factor-base family is
  one registry entry. Family parameters live in `params`. Anything not yet named
  goes under `extensions` as `namespace.knob`, and it is part of the identity.
- **Solver options that change search must be pinned.** WDSat's XG mode
  changes conflict counts by two orders of magnitude, so the registry lists
  the options each solver must pin. A spec that leaves one out is refused.
- **Defaults are written out before hashing**, so an omitted default and an
  explicit one give the same `spec_id`.
- **One target, one workload.** The primary window takes exactly one target.
  Each seed in `seeds` is a separate workload; more seeds give independent
  targets, and repetitions re-run the same seed.

Validation runs in three layers, and any failure refuses the spec: the JSON
Schema, the registry vocabulary, and policy rules across sections. The
adapter then refuses anything its producer cannot do faithfully. For example,
`ic bench` refuses `relations.stop: full_rank`, because it stops when the
target column is pinned.

## 2. The window

A figure always names its window, from [`registry.json`](registry.json)
`windows`:

| window | primary | what it covers |
|---|---|---|
| `online_one_target` | yes | first target-dependent computation to verified recovery, one target, five exclusive phases (`target_query`, `target_pdp`, `target_relation_check`, `target_descent`, `recovery_check`) |
| `cold_end_to_end` | no | the whole method, cold; every phase charged; the window S is defined over |
| `whole_process` | no | exec to exit of the measured child, timed by the runner; a practicality figure |
| `stage_only` | no | one phase or solver stage: a stage diagnostic, never a speedup |

The phase names are those of `src/cryptanalysis/ic_measurement.rs`, plus
`isogeny`, which the cryptanalysis `T_cold` charges. A producer whose clocks
match the five online phases reports the online window as `exact`. One that
can only assemble it from coarser phases reports it as `derived`. `ic bench`
writes the target into every relation, so its online window is `derived`. A
derived window is usable for operation counts and refused as a primary
wall-clock speedup.

## 3. Units

Every unit is a separate registry entry, so a comparison across units is
refused rather than made silently.

| unit | what it is | deterministic |
|---|---|---|
| `crypto.S.gae_pinned` | `ic bench` group-addition equivalents with pinned ratios | only when no calibration ratio was measured and no solver was priced from wall time |
| `crypto.S.price_batched_add` | `ic price` wall ns / ns of a batched addition | no |
| `crypto.S.callgrind_ir` | callgrind instructions / (targets·√r) | yes, for a fixed binary and valgrind version |
| `cryptanalysis.rps` | exact counters weighted by frozen per-class weights | yes, per `calibration_id` |
| `count.group_additions`, `count.s3_solves`, `count.solver_conflicts` | plain counts | yes (conflicts: per solver, version and options) |
| `time.wall_ns`, `time.cpu_ns` | time | no |

Each record states whether its total was deterministic. It also states why
not when it was not, for example a solver priced from host wall time or a
native counter with no pinned price that was charged zero.

## 4. References and floors

The rho reference moved three times between 2026-09-21 and 2026-10-01. A
record therefore names its reference by id: `rho.negation`,
`rho.signed_frobenius`, `rho.automorphism`, `rho.measured_matched`, or
`rho.batch_kuhn_struik`, which is a diagnostic for k targets only and never
the one-target reference. A comparison refuses two different references.

## 5. Structural metrics

Every record carries the sizes that decide whether two configurations
solved the same problem. Unknown is `null`, never zero.

- **Factor base**: `usable_points` is the fb count of an IC1 candidate id.
  The record also carries `signed_points`, `abscissae_with_points`,
  `abscissae_allowed`, `columns`, `orbit_representatives`,
  `points_per_column` and `nominal_dimension`. These are never substituted
  for one another.
- **Decomposition system**: `n_vars`, `n_equations`, `degrees`,
  `semi_regular_degree`, `solving_degree`, `anf_monomials`, `cnf_variables`,
  `cnf_aux_variables`, `cnf_clauses`, `xor_rows` and `literals`.
- **Solver**: `conflicts`, `decisions`, `propagations`, `restarts`, `calls`,
  `op_unit`, `priced_by`, `budget_exceeded`, and the engine's own `extra`.
- **Relations and linear algebra**: `targets_tried`, `relations_found`,
  `hit_rate`, `rows`, `rank`, `dependent`, `work` and `work_unit`.

An adapter also checks declared structure against observed structure. A spec
that declares `quotient: [negation]` on a run whose report shows one point
per column produces a `consistency` entry, and that entry refuses every
comparison of the record.

## 6. The host capsule

`python3 tools/icms capsule` prints this; every session saves it whole.

- **Stable** facts are hashed into `env_class_id`, and two wall-clock arms
  must share it:
  - the CPU's vendor, model, family, stepping and microcode;
  - the hash of the full flag set, plus the features AGENTS.md §10 names
    (popcnt, avx2, avx512*, pclmulqdq, vpclmulqdq, gfni, aes, and the Arm64
    equivalents);
  - the cache topology and core and SMT siblings of every CPU;
  - frequency policy: the governor per CPU, `no_turbo`, `boost` and EPP;
  - SMT control and state;
  - the kernel release, its command-line flags (`isolcpus`, `nohz_full`,
    `rcu_nocbs`, `mitigations`, ...), the clocksource and transparent huge
    pages;
  - twenty sysctls (`kernel.randomize_va_space`, `kernel.nmi_watchdog`,
    `kernel.numa_balancing`, `kernel.sched_autogroup_enabled`,
    `vm.swappiness`, ...);
  - the virtualisation type, cgroup limits, memory, OS and libc;
  - the toolchain (`rustc -Vv`, cargo, gcc, clang, Python and its
    executable hash).
- **Volatile** facts are sampled before and after every run and never
  hashed: load, PSI for cpu, memory and io, per-CPU jiffies including
  hypervisor steal, current frequencies and temperatures.

A fact the kernel does not expose is `null`. The Firecracker containers
these repositories use expose no cpufreq directory, report microcode `0x1`,
and do not support SMT control. On such a host the frequency checks are
unknown, and unknown never passes a gate.

## 7. Isolation levels

Each run earns the highest level whose checks, and all lower levels'
checks, passed (`tools/icms/gates.py`):

| level | requires |
|---|---|
| **L0 recorded** | the run executed under a saved capsule. Enough for operation counts, which do not depend on contention. |
| **L1 pinned** | child CPU mask read back equal to the reservation; CPU 0 excluded; other threads evicted (`isolated_bench.evict`); the tree's CPU time / wall ≤ declared threads |
| **L2 quiet** | L1, and during the window: run-queue delay of the child ≤ 0.5 % (from `/proc/<pid>/task/*/schedstat`), zero hypervisor steal on the pinned CPUs, no foreign runnable task on them in any sample, no memory stall, ≤ 50 preemptions per second, and a quiet `isolated_bench` preflight |
| **L3 isolated** | L2, and the host configured for measurement: pinned CPUs in `isolcpus` or `nohz_full`, performance governor with turbo off, SMT off or siblings reserved, bare metal |

A spec can tighten a threshold but never loosen one. A wall-clock figure is
admitted only at or above the spec's `isolation_required` (default L2).

What the reference containers earn is recorded, not assumed. On the 4-vCPU
Firecracker VM, runs pin and evict cleanly (L1). They often fail L2 on steal,
which the guest can see, and they cannot reach L3: it is a VM, and frequency
is unobservable. Two failures this gate catches that `isolated_bench`'s
`contended` flag does not:

- **Self-oversubscription.** Four Rayon threads on one pinned CPU spent 75
  percent of their runnable time waiting for each other.
- **Steal.** Hypervisor steal on the pinned CPU is invisible to a
  process-tick count.

## 8. How a session runs

`python3 tools/icms run A.yaml B.yaml A.yaml --cpus 3 --out <new dir>`:

1. **Validate before anything runs.** Every spec is validated and checked by
   its adapter, and any refusal stops the session.
2. **Never overwrite.** The output directory must not exist. The specs are
   byte-copied in, and the capsule is captured.
3. **Lock, check and preflight.** The `isolated_bench` lock is taken, the CPU
   set checked (SMT siblings, never every CPU), and the machine must be
   quiet. `--allow-busy` records a refused preflight and continues, and every
   run then stays below L2.
4. **Evict.** Other movable threads leave the reserved CPUs for the whole
   session and are restored afterwards.
5. **Interleave.** Each arm runs its warm-ups, then the rounds: every arm
   once per round, in ABAB order or a seeded shuffle. Listing a spec twice
   adds an A/A arm, the session's noise floor.
6. **Run each execution measured:**
   - pinned between fork and exec, with a built environment;
   - only `PATH`, `HOME`, the locale and toolchain roots pass through;
   - `RAYON_NUM_THREADS`, `OMP_NUM_THREADS`, `PYTHONHASHSEED` and `LC_ALL`
     are pinned, plus whatever the spec declares;
   - engine knobs found in the operator's shell (`KIC_*`, `IC_*`, `GAUDRY_*`)
     are dropped and listed;
   - sampled every 50 ms for schedstat, threads and foreign runnable tasks,
     with steal and PSI totals read around it;
   - left unreaped at exit long enough to read its own schedstat;
   - killed as a process group at the timeout.
7. **Write a record per execution**, failures, timeouts and warm-ups
   included.

ICMS is the in-process equivalent of `isolated_bench.py reserve`: it imports
the tool's lock, preflight, CPU checks and eviction, and records the tool's
hash.

## 9. Comparing

`python3 tools/icms compare <session> --a 0 --b 1 --declare factor_base.params.dimension`

1. **Classify every differing spec field:**
   - **declared**, which is what the comparison is about;
   - **forbidden**, which defines the problem: instance, window, unit,
     reference, stop rule, threads or CPUs;
   - **confound**, which is anything else undeclared.

   A forbidden difference or a confound refuses the comparison. So do
   different workloads, and records whose report contradicts their spec.
2. **Operations.** The median of each arm's verified, completed runs in the
   spec's unit and window gives `ratio_b_over_a`. A host-dependent unit is
   labelled valid on this host class only. A run that did not complete makes
   the total unknown, not zero.
3. **Wall time.** It needs both arms at their declared level, paired by
   round, and at least five rounds. It is reported as the median per-round
   ratio with a percentile-bootstrap 95 % interval, beside the A/A interval.
   Below the level it is printed as descriptive only, never as a result.

A speedup in the sense of AGENTS.md §8 is still
`baseline_total_operations / candidate_total_operations` over the whole
cold method with every phase priced. ICMS reports the ratio. Classifying the
change as advance, engineering, relabelling or accounting stays with the
author (§3).

## 10. Version pinning

Every session records:

- the repository commit, dirty flag, dirty paths and the sha256 of the
  uncommitted diff;
- the sha256 of every ICMS source file, the registry and `isolated_bench.py`;
- per producer, the binary sha256, the compiled-in `docs/ic/calibration.json`
  hash, and the `Cargo.lock` present at build;
- the full toolchain versions in the capsule.

Gaps this standard records but does not close:

- **No pinned toolchain.** The repository has no `rust-toolchain.toml`, and
  CI installs a floating `stable`.
- **Unpinned dependency versions.** `Cargo.lock` is gitignored, so they are
  pinned only in the build tree a session saw.
- **Wrong commit in `ic` reports.** An `ic` report's `commit` is the working
  directory's HEAD at start-up, not the commit the binary was built from. Use
  the binary hash.

## 11. Adapters

| adapter | producer | status |
|---|---|---|
| `crypto.ic_bench` | `ic bench` (IC framework) | implemented; tested on frozen reports |
| `command` | any program that writes `$ICMS_METRICS_PATH` | implemented; how a new producer joins without an adapter of its own |
| `cryptanalysis.ic_bench` | cryptanalysis `experiments/ic-bench/bench.py` | see the adapter's docstring for the current state |
| `autoresearcher.index_calculus` | crypto-autoresearcher `index_calculus sweep` | see the adapter's docstring for the current state |

taskq can carry a session unchanged: submit `python3 tools/icms run ...` as
the task's argv on the `cpu-exclusive` queue. The session directory then
lands in `TASKQ_OUTPUT_DIR` with its own capsule.

## Files

- [`registry.json`](registry.json): phases, windows, units, references,
  factor-base families and size fields, PDP metrics, and the solver options
  that must be pinned.
- [`schema/spec.v1.json`](schema/spec.v1.json): the specification schema.
- [`specs/`](specs): specifications.
- [`sessions/`](sessions): frozen sessions the page cites.
- [`audit.json`](audit.json): the comparability findings that motivated
  this standard, each with its evidence path.
- [`tools/icms/`](../../../tools/icms): the implementation and its tests.

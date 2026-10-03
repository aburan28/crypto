---
name: ecbench-measure
description: Measure and compare ECDLP methods (Pollard rho variants, BSGS, kangaroo, index calculus) with the repository's native ecbench harness - write a spec, run an isolated interleaved session, audit it with exact replays, compare arms with paired bootstrap intervals, tabulate in S against the floor and matched rho, load the database, and land the evidence in a PR. Use for any new cross-method ECDLP measurement or speed claim.
---

# Measure with ecbench

Read `AGENTS.md` (the boundary/table/ratio rule, §§1–8, §10, §11, and the IC
measurement rules) and [`docs/ecbench/README.md`](../../../docs/ecbench/README.md)
before measuring. ecbench is the native harness for new cross-method work. Do
not write a Python harness or extend ICMS for it. Historical autolab, tournament
and ICMS evidence keeps its own protocols.

## 1. State the question before running anything

Write, in the note or PR that will carry the result:

- the **boundary**: the curve's generic floor `√(π/2A)` in `S` (ecbench records it per
  run) and the **reference** arm (`rho.negation` on a prime curve,
  `rho.signed_frobenius` on a Koblitz curve);
- the **hypothesis**, the success and stop conditions, and what is inadmissible
  (changing the workloads, seeds, unit or accounting after seeing results);
- the **sizes**: at least four curve sizes if any claim is about scaling (§5).

## 2. Write the spec

Copy `docs/ecbench/specs/smoke-prime.json` and edit it. Rules:

- Curves come from constructors (`prime_search`, `prime_roster`, `koblitz`,
  `binary_random`, `prime_explicit`). Run `ecbench plan --spec S` and check every
  slug is `registered=true`. Register new curves (AGENTS.md §11) before citing
  them in prose.
- `targets_per_curve` sets independent one-target workloads. Prefer more targets
  to more rounds: BSGS's cost is fixed by its target, and intervals resample
  workloads. Use 8 or more targets per curve for a claim, and 2 for a smoke run.
- Name every parameter that changes the search. `ic.pipeline` refuses to run
  without `factor_base` and `oracle`. Defaults are written out before hashing.
- Include the reference arm, the unmodified baseline, and the candidate. Add a
  `control` arm that repeats one method when wall time matters: it is the A/A
  noise floor.
- `rounds` ≥ 5 if wall time will be reported. Keep `warmup: 1`.
- `isolation_required: L2` (default) is what a wall-clock figure needs. Operation
  counts need no level.

## 3. Run

```bash
cargo build --release --bin ecbench
```

```bash
./target/release/ecbench run --wait --spec SPEC.json --out OUT_DIR
```

- Linux: `--cpus auto` reserves a whole core away from CPU 0. Run as root for
  eviction to reach other users' threads (otherwise runs stop at L0). Use
  `--cpus inherit` under isolab or a prepared cgroup.
- macOS: there is no affinity control, so runs are L0. That is fine for operation
  counts and never for wall time.
- Never run builds, tests or other agents' jobs during a timed session. Wrap your
  own heavy work in `python3 tools/isolated_bench.py busy -- CMD`, which takes the
  same lock.
- A busy host is refused. Wait, or pass `--allow-busy` and accept runs below L2.
- `OUT_DIR` must not exist; it is claimed atomically under the lock. Never delete
  or edit a session. Ctrl-C or a SIGTERM stops a session cleanly: the child is
  killed, evicted threads are restored and the session is marked `interrupted`. An
  interrupted session stays as evidence; start a new directory.

## 4. Audit

```bash
./target/release/ecbench verify --dir OUT_DIR --replay 12 --out OUT_DIR.audit.json
```

The audit must be `OK`, with every replay `identical`. Use `--replay-all` for a
claim: every deterministic run is re-executed, every derived figure is recomputed
from its counts, and a session graded under the current rules is regraded. Cite
the receipt's SHA-256 (printed) as the replay certificate. A failed audit means no
result: find out why.

## 5. Compare and tabulate

```bash
./target/release/ecbench compare --dir OUT_DIR --a REFERENCE_ARM --b CANDIDATE_ARM --save
```

```bash
./target/release/ecbench table --dir OUT_DIR
```

Read it honestly:

- The ratio is `Σ S_B / Σ S_A` over (workload, round) pairs, the totals ratio of
  AGENTS.md §8. `ops (ok)` with an interval excluding 1 is a measured difference in
  `S` on these curves. `incomplete` (a measured run did not verify) is not
  admissible. `partial` (runs without a partner) needs explaining. `bounded` means
  unpriced work: say so.
- For a baseline binary against a candidate binary, run the same spec into two
  sessions and pass `--b-dir`. Pairs then share seeds (`same_seeds: true`), and
  wall time is refused across sessions by design.
- Read the **per-curve** rows. A pooled ratio across sizes is not a scaling claim.
- `ci_method: within` (one workload) says nothing about other targets.
- Wall time is a result only when `admitted`, and it means something only when
  `outside_noise` is true. Otherwise quote it as descriptive or not at all.
- BSGS below rho in `S` is the time–memory trade, not a finding (README §7).
- Classify the change (advance, engineering, relabelling, accounting) by AGENTS.md
  §3, from the ratio to the floor and to the reference, not from the headline.

## 6. Store

```bash
./target/release/ecbench db sql OUT_DIR | sqlite3 -bail ecbench.db
```

The database is an index rebuilt from sessions. Keep sessions in the PR, not the
`.db` file. Useful views: `method_by_curve`, `arm_workload_summary`,
`ic_phase_split`.

## 7. Land it

Per AGENTS.md, open or update a PR in this task:

- the spec, the session directory under `research/<topic>_<date>/sessions/`, the
  audit receipt, the saved comparisons, and a README with the hypothesis, the
  host (`ecbench host`), the levels earned, the table and the decision;
- the scoreboard and the thread's note in the same PR when the round measures an
  IC variant (§7). Mark extrapolations as such;
- failures, timeouts and regressions included.

CI (`.github/workflows/ecbench.yml`) re-audits every committed session with
replays on Linux. A session whose counts do not reproduce there fails the PR.

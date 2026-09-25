---
name: ic-run-tournament
description: Prepare, run, resume and inspect bounded end-to-end index-calculus candidate tournaments in the crypto repository, preserving frozen inputs, full costs and independently verified ECDLP results.
---

# Run an IC tournament

For the repo-local development loop, platform checks and backend panels, use
[ic-autolab](../../../../.agents/skills/ic-autolab/SKILL.md) and
[AUTOLAB.md](../../AUTOLAB.md). New instruction rounds retain a frozen diverse
portfolio (`--selection-width 6 --exploration-slots 1`) before selection locks
one challenger. Native development screens cannot promote a winner.

Locate the repository, read `AGENTS.md`, then
`research/ic_candidate_tournament_20260915/OPERATIONS.md`. Use existing session
authorization and resource limits. A bounded local run does not imply starting
an indefinite service or provisioning remote hardware.

## Prepare a new round

Check the runner's `--help`, the source state, free disk space, amd64 architecture,
Rust dependencies and Valgrind version. The current instruction protocol is pinned
to Valgrind 3.22.0. A new profiler/ISA needs a new calibrated protocol.

```bash
python3 research/ic_candidate_tournament_20260915/tournament.py prepare \
  --out /absolute/path/to/new-round --profile pilot
```

Use `--candidates FILE` for a prepared registry and `--source-root PATH` for the
intended source tree. For a task about native time or rho parity, pass
`--require-native-progress`; this adds the paired native-time improvement gate.
Preserve the previous contract's target count with `--targets N` (default one).
Changing the target count creates a separate workload panel, with all setup
charged to each complete cold job. The runner copies and hashes source dependencies, builds in
isolation, freezes fixtures/configurations, and writes a frozen evaluator copy.
The pilot confirms on 60 public-target fixtures across five curve cells with three
repetitions per arm. It is a bounded configuration result, not a family-wide claim.

Run the evaluator path printed by `prepare`:

```bash
python3 /absolute/path/to/new-round/evaluator/tournament.py run --round /absolute/path/to/new-round
```

Stages are `aa`, `smoke`, `development`, `selection`, `confirmation`, `replay`.
`--stage NAME` runs one stage after its prerequisites. Keep logs. Each profile is
paired with a fresh native run, within the frozen CPU, memory and process budgets.

## Continue and inspect

Use the frozen evaluator's `status --round PATH`. Reissuing `run` reuses completed
receipts only after verification; it never overwrites evidence. An interrupted
trial directory without a receipt is retained and rejected. Investigate it and
start a new campaign rather than deleting it to obtain a lucky completion.

The A/A control must pass before comparisons. An incomplete baseline, timeout,
invalid relation, missing phase or unmatched target cannot produce a win.
Neither can a target recovered without its descent certificate: every IC arm
must report the relation `[a]G + [b]Q = Σ P_i` each logarithm came from, and
the checker verifies it. Only index-calculus algorithms are candidates; rho is
the reference. Keep
all such outcomes visible. Do not raise limits midway through a sealed campaign.

## Finish

Run `verify --round PATH` with the frozen evaluator, inspect `decision.json`, and
produce the research report and canonical scoreboard update using the reporting
instructions in `OPERATIONS.md`. A round is complete when its evidence, verdict
and scoreboard agree.

The current score is complete **user-space guest instructions**, including all
worker phases. It excludes kernel/device work and the external audit. Do not
rename this as curve additions, use profiled wall time as native performance, or
claim a new ECDLP exponent. Cold and amortized workloads need separate contracts.

### Factor-base policy panels

The default `--comparison-kind fixed-support` rejects different paired supports.
Use `--comparison-kind factor-base-policy` only for a predeclared policy comparison
on identical public ECDLP targets. This permits support differences across arms
while still rejecting unstable support within an arm/case. Record actual B and
rank per arm and retain all setup costs. See the policy section of `OPERATIONS.md`.

Keep rejected smoke candidates in the admission table, report and scoreboard.
`report.py` writes `admission.json` separately from frozen final measurements.
A successful evidence audit includes faithfully retained rejected receipts; it
does not mean every candidate solve was admitted.

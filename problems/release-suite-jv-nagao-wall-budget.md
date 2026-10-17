# Parallel release-suite Joux–Vitse wall-budget failure

- **Status:** investigating
- **Kind:** observed validation timeout; underlying cause unknown
- **Affected work:** [volcano execution validation](../research/full_volcano_graphwide_ic_20261017/execution/README.md)
- **Next action owner:** library test/F4 backend maintainer; unassigned
- **Last updated:** 2026-10-17

## Problem

The full release library suite cannot be called green: the unchanged
`cryptanalysis::jv_cover::tests::nagao_finds_planted_six_sums_and_agrees_with_the_oracle`
times out in a parallel full run while passing in isolation. Do not present
the isolated pass as a passing full suite or patch the scientific assertions
to accept missing decompositions.

## Evidence and uncertainty

The complete restricted run (`cargo test --release --locked --offline --lib
-- --test-threads=2`, RAYON_NUM_THREADS=1) has 4,225 passes, ten failures and
117 existing ignored tests. Nine failures are explicit local-socket permission
denials. Repeating the whole `pollard_collab` module with local networking
permitted passes all 37 tests. Logs are in
[`execution/evidence/`](../research/full_volcano_graphwide_ic_20261017/execution/evidence/):
`merged-release-tests-restricted.log`, `merged-network-tests-permitted.log`.

The remaining failure has a deterministic planted six-sum reference. Trial 0
completes with F4 work 715,704 and F4 wall time about 88.56 s; trial 1 reports
`timed_out: true`, F4 work zero, about 90.10 s and no decompositions, whereas
the MITM reference has one. The isolated exact test passes all six trials in
2.55 s with the same built release test binary; see
`merged-jv-isolated-permitted.log`. The permitted parallel rerun reproduces the
timeout while `the_whole_registry_builds_and_round_trips` is also running.
The final permitted full run has 4,234 passes, that one failure and 117
existing ignored tests; see `merged-release-tests-permitted.log`.

Concurrency, lock waiting, process restrictions and resource contention are
possible explanations, not established causes. The Joux–Vitse source/assertion
was not changed by the volcano patch, but a failing base-branch full run has
not been newly reproduced here. This evidence does not prove a wrong
decomposition, a solver speed result, or failure of the separate prime-field
volcano control/pilot/audits. The planned serial restart was not launched:
automatic permission review of stopping/restarting the active run timed out;
the current run was allowed to finish instead.

## What would unblock it

- [ ] Run the complete unchanged release suite serially on an otherwise idle
  host with local networking permitted, retaining all tests and existing
  solver budgets, and obtain a full passing receipt.
- [ ] Reproduce the parallel timeout on the base branch and instrument waiting
  versus solver execution before attributing a cause or changing budgets.

A serial pass establishes that validation configuration, not thread safety or
qualified performance. A parallel failure remains open until explained and
addressed without skipping tests or weakening exact-answer assertions.

## Next discriminating step

After the scientific worker is paused at a checkpoint, run the existing
cached binary through `cargo test --release --locked --offline --lib --
--test-threads=1`, with local networking permitted. Preserve the full result;
if it still fails, investigate that specific failure before another full run.
Do not run this heavy suite alongside measured campaign sessions.

## Resolution and history

2026-10-17: recorded both the full-run timeout and passing isolated control.
Unresolved; no all-suite pass or qualified runtime conclusion is asserted.

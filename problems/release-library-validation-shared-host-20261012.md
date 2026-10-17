# Release library validation under shared-host contention

- **Status:** resolved
- **Kind:** resource and sandbox validation limit, resolved by a
  controlled complete rerun
- **Affected work:** [isogeny and torsion fiber study](../research/isogeny_torsion_fibers_20261012/REPORT.md),
  [first validation receipts](../research/isogeny_torsion_fibers_20261012/validation/README.md),
  [paired-round validation](../research/isogeny_torsion_fibers_20261012/paired_round/validation/README.md)
- **Next action owner:** repository CI only if the final PR head needs
  separate confirmation
- **Last updated:** 2026-10-12

## Problem

The first complete `cargo test --release --lib` run failed under shared
host load and restricted loopback access: 3,996 passed, 13 failed and
113 ignored. Four existing `jv_cover` tests recorded timeouts or too few
accepted runs; nine `pollard_collab` tests received OS error 1 at socket
setup. A second complete run with loopback access and two test threads
passed: **4,009 passed, 0 failed, 113 ignored**. The local validation
gate is met on the source tree used for this PR; the failed first run
remains part of the record.

## Evidence and uncertainty

- **Observed:** [raw attempt logs](../research/isogeny_torsion_fibers_20261012/validation/README.md)
  for default, four-thread, serial and isolated runs. The default and
  serial attempts were interrupted; the four-thread attempt contains
  three `FAILED` lines but no final suite summary.
- **Observed in the paired round:** the [complete release log](../research/isogeny_torsion_fibers_20261012/paired_round/validation/release_lib.log)
  has the final 3,996/13/113 result. Two `jv_cover` assertions show
  `timed_out: true`; two require more accepted runs than survived.
  All nine `pollard_collab` failures show `PermissionDenied` and
  `Operation not permitted` at socket setup. The [outside-sandbox
  replay](../research/isogeny_torsion_fibers_20261012/paired_round/validation/loopback_replay.log)
  passed its selected socket test. The isolated [Nagao replay](../research/isogeny_torsion_fibers_20261012/paired_round/validation/jv_cover_isolated.log)
  passed 1/1 in 11.77 seconds after failing with a timeout in the parallel
  suite. The other three [isolated replay logs](../research/isogeny_torsion_fibers_20261012/paired_round/validation/README.md)
  also passed individually on this head. The entire affected
  [networking group](../research/isogeny_torsion_fibers_20261012/paired_round/validation/pollard_collab_outside_sandbox.log)
  passed 37/37 outside the restricted sandbox. A [host-load snapshot](../research/isogeny_torsion_fibers_20261012/paired_round/validation/host_load_during_lib.log)
  was retained during the run.
- **Controlled complete run:** the [second full log](../research/isogeny_torsion_fibers_20261012/paired_round/validation/release_lib_unsandboxed_2threads.log)
  records `RUST_TEST_THREADS=2 cargo test --release --lib` outside the
  restricted sandbox and its final 4,009/0/113 pass in 2,505.94 seconds
  of test execution. The four former `jv_cover` failures and affected
  networking tests passed within this complete run. No library or
  research-runner source changed after it began; later changes were to
  reports and receipts.
- **Derived or proved:** the study's standalone census, map check, and
  Lean proof checks passed; those do not replace this repository gate.
- **Suspected explanation:** host contention contributed to the
  `jv_cover` timeout behavior; isolated reruns distinguish that from
  an actual solver regression. Sandbox socket denial is directly
  observed for the `pollard_collab` failures.
- **Missing capability:** none for this local validation gate. Final-head
  repository CI is a separate review gate.
- **What this does not establish:** a library regression caused by this
  branch, a mathematical result, or how a future run behaves under
  different host load.

## What would unblock it

- [x] Run `cargo test --release --lib` to completion with logged host
  conditions and retain the full result. This completed attempt failed;
  it does not satisfy the passing gate.
- [x] Run the full suite to a pass on the unchanged source tree with
  loopback access and limited local concurrency; retain the exact output.
- [x] Capture the first run's assertions, rerun the four solver tests
  individually, rerun the affected networking group with socket access,
  then rerun the complete suite without weakening its assertions.

A complete pass satisfies the local release-library requirement for the
tested head. It does not validate the study's mathematical claims beyond
their separate proofs and exact controls.

## Recheck trigger

Reopen this record if the final PR head's required CI fails, a library
or runner source edit follows the controlled pass, or a reproducible
failure appears under the same two-thread, loopback-accessible conditions.
The authoring session leaves PR CI to the repository's separate workflow.

## Resolution and history

2026-10-12: opened after three incomplete full-suite attempts and three
passing isolated replays. The paired round then completed a full suite
with 13 failures, separated socket denial from recorded `jv_cover`
timeouts, and retained passing narrow replays. A second full release
suite with loopback access and two test threads passed 4,009/0/113;
this resolves the local gate while preserving the failed receipt.

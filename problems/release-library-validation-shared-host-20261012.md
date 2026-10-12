# Release library validation under shared-host contention

- **Status:** open
- **Kind:** resource and validation limit; cause of full-suite behavior
  not proved
- **Affected work:** [isogeny and torsion fiber study](../research/isogeny_torsion_fibers_20261012/REPORT.md),
  [local validation receipts](../research/isogeny_torsion_fibers_20261012/validation/README.md)
- **Next action owner:** repository CI or an isolated validation runner
- **Last updated:** 2026-10-12

## Problem

The required `cargo test --release --lib` gate has not produced a complete
passing result for this PR. The shared macOS host had load above 100 for
much of the run. A four-thread attempt recorded three failures in existing
`jv_cover` tests with solver wall caps, then was stopped. Each of those
tests passed in a separate one-thread rerun. The branch does not edit
library source, so it cannot by itself explain a changed library result;
that fact does not certify the complete suite.

## Evidence and uncertainty

- **Observed:** [raw attempt logs](../research/isogeny_torsion_fibers_20261012/validation/README.md)
  for default, four-thread, serial and isolated runs. The default and
  serial attempts were interrupted; the four-thread attempt contains
  three `FAILED` lines but no final suite summary.
- **Derived or proved:** the study's standalone census, map check, and
  Lean proof checks passed; those do not replace this repository gate.
- **Suspected explanation:** concurrent runs on a heavily loaded host
  hit the `jv_cover` solver's 120-second per-call wall caps. The exact
  assertions for the failed parallel tests were not captured before
  the suite was stopped, so this remains a hypothesis.
- **Missing capability:** a complete release library suite on a host
  with controlled load, or an otherwise conclusive diagnosis and fix
  for the three failures.
- **What this does not establish:** a library regression caused by this
  branch, a mathematical obstruction, or a passing full-suite result.

## What would unblock it

- [ ] Run `cargo test --release --lib` to completion on the final PR head
  with logged host conditions and retain the full result.
- [ ] If any test fails, capture its assertion and reproduce it under
  controlled load; fix an actual defect without weakening its oracle or
  correctness gate, then rerun the affected checks and the full gate.

A complete pass satisfies the local release-library requirement for the
tested head. It does not validate the study's mathematical claims beyond
their separate proofs and exact controls.

## Next discriminating step

- **Action:** execute the full release library suite on an isolated
  runner or in repository CI on this PR head.
- **Inputs and controls:** exact committed head, Rust/Cargo versions,
  test thread count, host-load receipt, and the three isolated test names.
- **Execution limits:** use the runner's ordinary test limits; retain a
  timeout as a timeout. Do not silently skip or expand a test wall cap.
- **Outcomes:** a complete pass resolves this gate; a reproducible
  isolated failure becomes a code defect to fix; a load-sensitive-only
  failure remains a host-isolation issue; incomplete execution leaves
  this problem open.
- **Recheck trigger:** first applicable PR validation run on the final head.

## Resolution and history

2026-10-12: opened after three incomplete full-suite attempts and three
passing isolated replays. Unresolved.

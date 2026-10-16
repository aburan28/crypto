# JV F4 controls complete too few cases in the full library suite

- **Status:** open
- **Kind:** validation gap; runtime cause unknown
- **Affected work:** [genus-two integration validation](../research/binary_genus2_cover_20261014/REPORT.md)
- **Next action owner:** unassigned
- **Last updated:** 2026-10-14

## Problem

On base 275beb93f3bf33ade6c19f19e8b336027a00b297, the full optimized library
suite reports too few completed comparison cases in two existing JV F4 tests.
This blocks an unqualified claim of a green whole-library gate. It does not
invalidate the separately checked binary genus-two cover.

The two unchanged tests are
`cryptanalysis::jv_cover::tests::replayed_f4_agrees_with_the_full_solver_and_the_oracle`
and `cryptanalysis::jv_cover::tests::stopped_f4_agrees_with_the_full_solver_and_the_oracle`,
using their existing small-context parameters (53,7). Their acceptance requires
at least 30 and 25 completed cases respectively. Do not reduce these thresholds
or mark the tests ignored to get a pass.

## Evidence and uncertainty

- **Observed:** `cargo test --release --lib --locked --offline` completed
  with 4087 passed, 2 failed, 117 ignored. The failed assertions report
  17 completed replay cases and 19 completed stopped cases. The full output
  is [root-library-tests.log](../research/binary_genus2_cover_20261014/runs/run-20261014T1442Z/root-library-tests.log).
- **Observed:** the replay test passes in an isolated execution of the same
  built test binary; [raw diagnostic](../research/binary_genus2_cover_20261014/runs/run-20261014T1442Z/jv-base-failure-replay.log).
- **Observed:** both tests pass together with `--test-threads=1` on that same
  binary; [serial diagnostic](../research/binary_genus2_cover_20261014/runs/run-20261014T1442Z/jv-serial-diagnostic.log).
- **Source binding:** this change modifies no root library source,
  Cargo.toml or Cargo.lock; the empty `root-library-diff.txt` records that
  comparison with the pinned base. The initial missing compile-time records
  fixture was materialized and that compile failure was resolved separately.
- **Suspected explanation:** resource contention or shared runtime state is
  possible, since the tests filter incomplete/timed-out cases. This is a
  hypothesis; rejection reasons for each omitted case are not in the failure
  output, and no cause is established.
- **Not established:** a failure of the solver/oracle equality on an accepted
  case, a defect caused by the new binary cover, or a green full-suite gate.

## What would unblock it

- [ ] Record the omitted-case reasons for these unchanged deterministic
  workloads, distinguish shared-state and resource-limit effects, and retain
  accepted-case oracle equality checks and original completion thresholds.
- [ ] Run the full original library command successfully with the corrected
  setup or a reviewed fix, preserving the failed baseline receipt.

Passing a smaller diagnostic alone does not meet the whole-suite exit criterion.

## Next discriminating step

Compare the two exact tests serially with their isolated receipts and the
full-suite receipt on the same built binary; retain all outcomes. If the
serial pair passes, instrument the omitted-case reasons before altering any
budget or shared state. If it fails, inspect the first incomplete case under
the existing parameters. No full-suite pass is inferred from the diagnostic.

## Resolution and history

2026-10-14: opened with the complete failing library run and passing isolated
replay control. Serial-pair diagnostic is retained alongside the run. Unresolved.

# Current ecbench integration does not compile

- **Status:** open
- **Kind:** implementation gap
- **Affected work:** [backend frontier validation](../docs/ecbench/backend-frontier-report/REPORT.md); root library and ecbench builds
- **Next action owner:** unassigned; ecbench integration maintainer
- **Last updated:** 2026-10-17

## Problem

The root library at `c6d64436a33359320ad9d5b42c4acebee2968b9a` cannot complete
release validation because existing ecbench/quartic integration declarations are
inconsistent. The new backend report therefore cannot be validated through a
root-library executable on this base. It is delivered as an independently
buildable native Rust tool that shares the harness's read-only pricing and
statistics source; no affected solver or curve implementation is changed.

## Evidence and uncertainty

- **Observed:** `cargo build --locked --release --bin ecbench` fails with 25
  compiler errors; `cargo test --locked --release --lib` fails with 26, before
  running tests. [Build log](../docs/ecbench/backend-frontier-report/validation/root-build.log)
  and [library test log](../docs/ecbench/backend-frontier-report/validation/root-lib-tests.log).
- Errors include two `Instance::Fp4` declarations, missing `PrimePublicFixture`,
  unexported `fbm`, missing `OnbWeightBase` imports and `sparsest`/
  `echelon_by_lead_with` definitions. The affected root files are unchanged by
  the backend report PR. Native tool validation is recorded separately.
- **Suspected explanation:** incomplete integration of independently added
  harness variants; no specific historical merge is identified as the cause.
- **Missing capability:** a compiling, tested root integration on this base.
- This does not show that the algorithms or mathematics are impossible, nor
  does it invalidate a frozen session's independently bound historical audit.

## What would unblock it

- [ ] Reconcile the variant declarations, exports, imports and missing helper
  interfaces without removing intended methods or changing frozen accounting.
- [ ] Pass the root release library suite and affected ecbench/variant tests,
  then independently replay affected frozen sessions.

These conditions restore root integration validation. They do not supply absent
large-field backend timings or establish a new performance claim.

## Next discriminating step

Build the recorded base in a clean checkout, repair the named integration
interfaces on a focused branch, and run the listed release checks plus affected
session audits. Stop and preserve diagnostics on a compile/test failure. A pass
would permit root-integrated report validation; a failure leaves this record
open. Recheck when the root integration changes.

## Resolution and history

2026-10-17: discovered during backend leaderboard work. The standalone report
path removes this build dependency from leaderboard generation, while root
validation remains unresolved. No cryptographic implementation was repaired or
reverted in this change.

# Native F5 controller implementation checks

Physical macOS ARM64 local checks, 2026-10-03 UTC, under the native busy lock:

- Four initial controller controls passed: strict scope/budgets/unknown fields,
  mathematical-only input matching the accepted geometry/log table,
  rejection of ancestor Cargo configuration, and simulated exact transport
  mutations including input/environment/output/claim/drain/exit status.
- Broad regression: 68 icprog tests, 7 worker/transport tests and 6 integration
  checks passed. Two legacy integration checks are deferred locally to CI,
  which installs/fetches their dependencies and runs `--include-ignored`.
- The transport simulation starts no worker. It labels its registration as
  a test fixture, uses a placeholder binary identity, and exercises only the
  data checker. Its historical report fails the new seed's mathematical gate;
  it cannot stand in for a new native execution.
- Local Clippy found no new controller warning; this host's older compiler
  reports existing unrelated library and unknown-lint warnings. CI's pinned
  compiler performs its own strict all-target lint check.
- All five controller controls passed after adding validation-only mode. The
  fifth verifies that validation-only, null or undeclared capsule mode cannot
  dispatch or admit a scientific job.

The implementation is not a scientific execution, ordinary yield result or
new native F5 runtime admission. Actual source freeze, durable preregistration,
sole invocation and frozen independent audit are separately pending. The local
template is not an executable registration. No old registration is resumed.

## Offline build validation

The first validation-only freeze completed from committed source
`92b24110d49de660c8e4451f907a58b5c49a09bd`. Its validation seal is
`421ac32d7833673ae10e414a398ccb887dcf02556d88a135804ad73869eb61bf`.
It inventories 5,953 immutable files. Both worker and checker were built from
the retained offline source/dependency tree, and the worker's compiled source
and build identity matched its descriptor. The only worker call supplied
`--build-identity`, with no job input. No scientific solver ran.

`build-validation-v1` retains the exact original config, mathematical-only job,
preparation, host context, registration/seal, Cargo.lock and all build/identity
logs and receipts. Its compact publication checks the original capsule and
identity-only receipt before copying. This build evidence does not publish the
complete capsule archive or establish archive custody of an actual scientific
registration. The full local capsule remains non-executable. Actual durable
preregistration and sole scientific dispatch remain separate pending gates.

The portable data replay checks the external registration seal, all retained
inputs, the complete original build-receipt set, dependency lockfile and compiled
source/build identity. It never executes an archived binary and does not verify
custody of the omitted full archive. Linux and macOS CI run this same replay and
retain their own receipt:

```sh
/tmp/ic-native-busy busy -- target/debug/icprog f5-control-replay-validation \
  --publication research/ic_candidate_tournament_20260915/goal_20260924/native-f5-control-v1/build-validation-v1 \
  --validation-registration-sha256 421ac32d7833673ae10e414a398ccb887dcf02556d88a135804ad73869eb61bf \
  --out /tmp/native-f5-build-validation-replay-NEW.json
```

Use a new receipt path. The replay rejects a wrong external seal, modified
inputs or original receipts, omitted receipts despite a regenerated outer
inventory, and runtime/performance claims added to this build-only record.

After adding portable replay and the explicit external-seal requirement for
original audit, all 73 local `icprog` tests passed, including nine controller
controls. The built CLI also replayed the committed bundle successfully;
`local-macos-arm64-build-validation-replay-v1.json` retains that data-only result.
Formatting and diff checks passed, and local Clippy found no new implementation
warning. Cross-platform final-head CI remains the publication acceptance gate.

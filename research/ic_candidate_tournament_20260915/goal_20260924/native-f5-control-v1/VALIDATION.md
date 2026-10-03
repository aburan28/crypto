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

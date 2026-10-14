# Local validation receipts

The study changes standalone Rust research/checker code, a Cargo example,
proof modules, data and reports. No library source under `src/` changed
since the preceding [controlled complete library pass](../../paired_round/validation/release_lib_unsandboxed_2threads.log).
The earlier pass recorded 4,009 passed, zero failed, 113 ignored on
the same library source; its [obstruction record](../../../../problems/release-library-validation-shared-host-20261012.md)
describes the host controls. This directory retains this round's runs
without rewriting the earlier receipt.

## Release library

The first `cargo test --release --lib` invocation was interrupted during
compilation before a terminal test result. Its
[stdout](cargo_release_lib.stdout) and [stderr](cargo_release_lib.stderr)
remain an incomplete infrastructure receipt.

The second exact `cargo test --release --lib` invocation completed under
restricted workspace loopback access. Its [stdout](cargo_release_lib_v2.stdout)
and [stderr](cargo_release_lib_v2.stderr) record **3,996 passed, 13 failed,
113 ignored**, in 1,348.08 seconds of test execution after a 19-minute
release build. Nine `pollard_collab` tests failed with `Operation not
permitted` during socket setup. Four `jv_cover` tests failed under this
shared-host run; the earlier controlled full pass on unchanged library
source passed those cases. Neither group belongs to the new standalone
study code. This is a failed full-suite receipt, not a passing gate.

The controlled complete rerun used loopback-accessible execution and
`env RUST_TEST_THREADS=2 cargo test --release --lib`. Its
[stdout](cargo_release_lib_controlled.stdout) and
[stderr](cargo_release_lib_controlled.stderr) record **4,009 passed,
zero failed, 113 ignored** in 1,243.74 seconds of test execution.
All four formerly failing `jv_cover` tests and all nine formerly failing
`pollard_collab` socket tests passed in this full run. This is the
passing local release-library gate on the unchanged library source.
The failed unrestricted result above remains part of the evidence.

## Study-specific validation

- `screen.rs` compiled natively, executed all 12 frozen rows with
  assertions, and replayed byte-identically; see [RUN.md](../RUN.md).
- `equivariant.rs` compiled natively, checked 51,492 rational points and
  replayed byte-identically after the retained missing-temporary-binary
  failure; see the [run report](../reports/run-20261014T0551Z/README.md).
- The native identity-certificate example passed two direct `rustc --test`
  tests and replayed both records with 64 sealed plus 64 extra checker-seed
  points. `cargo test --release --example isogeny_torsion_identity` also
  passed both tests; [stdout](cargo_example.stdout) and
  [stderr](cargo_example.stderr) are retained. Cargo regeneration matched
  both certificate JSON files byte-for-byte; see
  [CERTIFICATES.md](../CERTIFICATES.md).
- The pinned Lean project built both modules. The quotient-support theorem
  has no axioms; the order-three linear-action consequence uses standard
  Lean foundations. See [PROOF_MAP.md](../PROOF_MAP.md).
- The result PDF and both run PDFs were rendered and visually inspected.

The checked run logs have SHA-256 digests
`7d8c507e3a043f12c8f0a6c2e4ef30b51a9722b8346aa657ab272e2ad351c877`
for the unrestricted full stdout,
`31b054c8f62d32c3555921a2def688d3345add4e108e5d3ff764ba9a2a75b0f5`
for the controlled full stdout, and
`09d698477ddcfcbf00d4773c00cb4006b922326baa78a189f90b40bbf44805c7`
for the Cargo example stdout. The output files are retained as separate
receipts. No whole-method ECDLP timing is inferred from validation runtime.
The study's [artifact manifest](../ARTIFACTS.sha256) covers these
receipts; verify it from the study directory with
`shasum -a 256 -c ARTIFACTS.sha256`.

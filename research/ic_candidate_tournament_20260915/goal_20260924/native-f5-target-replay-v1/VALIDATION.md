# Native F5 target replay validation

Local physical macOS ARM64 validation, 2026-10-03 UTC, through the shared native
busy lock:

- `cargo test --locked --bin icprog --bin prepared_sat_worker --test icprog`:
  64 icprog tests, 7 worker/transport tests and 6 integration checks passed.
  Two legacy integration checks require Valgrind or an archived Git commit;
  the harness CI installs/fetches them and runs with `--include-ignored`.
- Six target-verifier controls cover retained independent mathematics,
  changed query/witness/log/scalar/engine/exhaustion/phase data, a false
  negative, incomplete-target unknown phases, a missing phase key, and
  create-only portable replay.
- Two added transport controls cover exact 256 KiB job delivery with declared
  environment and a nonreading child given a 1 MiB job under a 30 ms deadline.
  The latter retains partial-delivery and confirmed process-group drain.
- The final local CLI receipt in `local-macos-arm64-replay-v3.json` reports
  `PASS_NATIVE_PREPARED_F5_TARGET_MATHEMATICS`, all three historical target
  attempts and independently verified scalar 24886. It executes no child.
- `git diff --check` passed. Existing unrelated library warnings remain.
- Local Clippy found one unnecessary temporary vector in the new verifier;
  it was replaced with a borrowed single-element slice before final validation.
  `local-macos-arm64-replay.json` preserves the original pre-lint checker
  receipt. Both receipts replay the same unchanged producer data and scalar.
- The first Linux CI wrapper failure is retained in `first-linux-ci-failure.txt`.
  Its direct mathematical controls passed. Checker self-identification now
  permits a bounded 512 MiB debug/test executable with file-specific failure
  diagnostics; its new local six-test run passed. Earlier checker receipts
  remain preserved. Final Linux CI is a separate required gate.

The local receipt identifies the actual checker binary and source. Linux/macOS
CI receipts identify their own builds; binary hashes need not match across
platforms. All must reconstruct the same target queries, statuses, scalar,
preparation and historical phase costs. No new solver result, natural yield,
source-bound runtime admission, calibrated timing or speedup is claimed.
CI status is determined from the final PR head, not this local validation note.

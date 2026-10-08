# Implementation validation

This validates a target-free producer API on the disclosed synthetic n17
instance. It is not a registered execution, natural-yield measurement, complete
one-target solve or performance comparison. The full goal remains active.

At the final local source, all seven `prepared_ordinary::tests` controls passed
with zero failures. They replay retained inputs and mock transport failures;
they execute neither MatrixF5 search nor a native CryptoMiniSat process. In
particular, the 216 historical attempts derive all 29 logs without importing
the known-log table into the preparation implementation. The SAT model test
checks the historical source assignment and a deliberately changed assignment.
Those controls do not estimate the proposed new panel's success rate.
All eight independent `ordinary_preparation::tests` controls also passed,
including the new producer-to-auditor schema, full-base and query-law check
using mock timeouts. A deliberately changed producer scalar is rejected.

The raw final test and Clippy output is retained losslessly as
`local-tests-and-clippy-v4-20261004.txt.gz`; decompress with `gzip -cd`.
The earlier three outputs remain alongside it, including the redundant closure
warning that was corrected before the final validation. `local-build-receipt.txt`
records the earlier source; `local-build-receipt-v4.txt` records the final source,
lockfile and log hashes and the toolchain. None of these
receipts attests source-bound scientific execution or an isolated CPU host.

Final command, run from this owned worktree through the shared busy launcher:

```sh
/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy -- sh -c '
  rustfmt --edition 2021 src/cryptanalysis/prepared_ordinary.rs src/bin/icprog/ordinary_preparation.rs &&
  CARGO_TARGET_DIR=/Volumes/SSD990/llm/tmp/ic-native-preparation-build-cache-20261004 cargo test --locked --offline --lib prepared_ordinary::tests -- --test-threads=1 &&
  CARGO_TARGET_DIR=/Volumes/SSD990/llm/tmp/ic-native-preparation-build-cache-20261004 cargo test --locked --offline --bin icprog ordinary_preparation::tests -- --test-threads=1 &&
  CARGO_TARGET_DIR=/Volumes/SSD990/llm/tmp/ic-native-preparation-build-cache-20261004 cargo clippy --locked --offline --lib --bin icprog --tests
'
```

The command exited zero. Local Clippy reported no warning in the new module,
but this is **not a clean strict Clippy pass**: local Rust 1.93.1 does not know
one configured newer Clippy lint, and existing modules emit unrelated warnings.
The repository's pinned CI toolchain and strict all-target checks are the
acceptance gate. Rust formatting and `git diff --check` passed locally.
Linux/macOS correctness and exact-head CI remain pending when this document
is committed; any completed receipts must identify the tested source head.

The owned Cargo cache is separate from the immutable scientific capsules. No
old registration was invoked, restored, resumed or refilled. The ignored
development Cargo.lock is hashed here but is not a scientific build freeze.
After terminal validation, disposable incremental objects from this task's
owned cache were removed to recover disk space; logs and receipts are retained.
The future executor still needs accepted, pinned source/assets/binaries and
an independently audited one-use registration.

This source-only change adds no IC1 result, new curve, ecbench session, timing
ratio or tournament round. It updates the autolab skill and status links; the
scoreboard, leaderboard and lab-browser measurements remain unchanged. Source
controls must not appear there as a new measured winner.

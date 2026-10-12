# Local validation receipts

Date: 2026-10-12. Base revision:
`e3090cbc98c657106fef10020da8ccc19be49a3f`.
This branch changes the fiber report and adds standalone Rust study files;
it does not change `src/` or `Cargo.toml`.

## Passed

- `rustc --edition 2021 -O -D warnings` compiled `census.rs`,
  `map_check.rs`, and `plot.rs`. The exact run commands and checks are in
  the [run receipt](../reports/run-20261012T1947Z/README.md).
- `rustfmt --edition 2021 --check` passed for all three Rust sources.
- `lake build` passed with the pinned Lean toolchain; the exact output
  and `#print axioms` results are in [build.log](../lean/build.log).
- `typst compile` built both PDFs; each opened as a one-page letter PDF
  with the diagram and graph visible on visual inspection.
- The three `jv_cover` tests flagged in the concurrent library run each
  passed when run alone with `RUST_TEST_THREADS=1`, `CARGO_NET_OFFLINE=true`
  and `--exact --nocapture`: [six-sum test](jv_cover_exact.stdout)
  (12.71 s), [stopped F4](jv_cover_stopped.stdout) (50.08 s), and
  [replayed F4](jv_cover_replayed.stdout) (88.05 s). Their stderr
  transcripts are retained beside the stdout files.

## Incomplete or failed gates

- `CARGO_NET_OFFLINE=true cargo test --release --lib` compiled the test
  binary in 15 min 13 s, then was stopped after 1,366 reported passes as
  many long tests ran concurrently on a host with load above 100.
  [stdout](cargo_release_lib_interrupted.stdout) and
  [stderr](cargo_release_lib_interrupted.stderr) are preserved. Exit 130
  means interruption, not a passing suite.
- `RUST_TEST_THREADS=4 CARGO_NET_OFFLINE=true cargo test --release --lib`
  was stopped after 1,693 reported passes and three existing `jv_cover`
  tests marked `FAILED`. The named tests all passed in the isolated
  reruns above. The full attempt's [stdout](cargo_release_lib_threads4_failed.stdout)
  and [stderr](cargo_release_lib_threads4_failed.stderr) are preserved.
  The suite did not reach a final summary; the exact failure assertions
  are unavailable in this captured output. Contention is supported by
  the isolated reruns, but the complete suite remains unverified.
- `RUST_TEST_THREADS=1 CARGO_NET_OFFLINE=true cargo test --release --lib`
  was stopped after four passes when an existing ElGamal test remained
  active under host load around 95. [stdout](cargo_release_lib_serial_interrupted.stdout)
  and [stderr](cargo_release_lib_serial_interrupted.stderr) are preserved.
- `python3 scripts/check_curve_names.py --diff HEAD` exited 1 with
  174 repository-wide issues. Its [complete output](curve_names_diff_head.log)
  has no paths changed by this branch. It reports existing unregistered
  slugs and an alias-table mismatch; this result is not a clean checker
  pass.

The full release library gate is unresolved. The observed failures are
not evidence that the added standalone study or corrected proof is wrong;
they also do not establish that the whole repository suite passes. The
[validation problem](../../../problems/release-library-validation-shared-host-20261012.md)
states the exact exit condition. No CI result was read or inferred here.

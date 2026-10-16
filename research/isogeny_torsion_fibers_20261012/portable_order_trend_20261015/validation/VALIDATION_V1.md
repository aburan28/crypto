# Local validation receipt, 2026-10-15

The required `cargo test --release --lib` was run from this branch with
`CARGO_TARGET_DIR=target/order-batch-check`. It compiled in 21m 12s, ran
4,200 tests and exited 101: **4,070 passed, 13 failed, 117 ignored** in
937.29s. This is a failed suite, not a green check. Nine
`pollard_collab::{cairn,net}` tests failed because this sandbox denied
local socket operations (`Os { code: 1, kind: PermissionDenied }`). Four
`jv_cover` tests failed after their time budgets were exhausted while
other release suites were contending for the shared host. The four are
`nagao_finds_planted_six_sums_and_agrees_with_the_oracle`,
`replayed_f4_agrees_with_the_full_solver_and_the_oracle`,
`stopped_f4_agrees_with_the_full_solver_and_the_oracle`, and
`the_solver_loses_no_rational_solution_at_small_p`. The first run's full
terminal stream was not redirected to a durable file; its final summary
and panic excerpts were retained in the task transcript. This branch
changes no file in `src/cryptanalysis/jv_cover.rs` or
`src/cryptanalysis/pollard_collab/`.

The exact follow-up checks were:

| command | result |
| --- | --- |
| `cargo test --release --lib cryptanalysis::pollard_collab:: -- --test-threads=1` with approved local socket access | 37 passed, zero failed, 1.98s |
| `RUST_TEST_THREADS=1 cargo test --release --lib cryptanalysis::jv_cover::tests::nagao_finds_planted_six_sums_and_agrees_with_the_oracle -- --exact --nocapture` | 1 passed, 23.05s; six planted trials all returned checked decompositions |
| same exact command for `stopped_f4_agrees_with_the_full_solver_and_the_oracle` | 1 passed, 56.24s |
| same exact command for `replayed_f4_agrees_with_the_full_solver_and_the_oracle` | 1 passed, 55.59s |
| same exact command for `the_solver_loses_no_rational_solution_at_small_p` | 1 passed, 11.93s |
| `cargo test --release` on the nine touched example targets | all nine compiled and executed; each has zero unit tests |
| `cargo test --release --example different_j_host_performance`, followed by a release build and local artifact-shaped self-check | compiled; analyzer matched 252/252 arithmetic rows, 5,355 scalar outputs and 37,061 relation witnesses against the committed Mac records |
| synthetic tamper of one copied `total_fp_muls` value with the replay hash recomputed | analyzer exited 101 at `arithmetic drift 409:433 (1, 0, "source")`; the original committed run was not changed |
| `rustfmt --check`, `actionlint`, `bash -n`, and `git diff --check` | passed |
| native three-cell fiber and batch checkers | 252 new child exits zero, 5,355 new scalar outputs and 37,061 relation witnesses replayed; see `../RESULTS_V2.md` |

A second full-suite run with approved socket access and four test threads
was started. It again reported the `nagao` and `replayed_f4` Jacobian-cover
tests as failed under concurrent load. It was deliberately interrupted
after those failures, so it has **no final complete-suite count**. The
[raw partial log](full-release-threads4-partial.log) is 144,667 bytes,
SHA-256 `7e64a450178c7f8bfab7f81108a027ce41a9af9fb6fe5b1db7d6f002a4c6c110`.
No test was skipped or disabled to obtain a green claim. The independent
Linux PR checks must be reviewed on their actual head and reported
separately. This validation result does not alter the arithmetic replay
or promote a performance claim.

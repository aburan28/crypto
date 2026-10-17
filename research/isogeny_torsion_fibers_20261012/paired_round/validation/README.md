# Local validation for the paired isogeny round

Host: shared Darwin arm64, kernel 25.6.0 (L0), 2026-10-12 local;
Cargo 1.93.1 and rustc 1.93.1 (Homebrew). A
[load snapshot](host_load_during_lib.log) was retained while the release
library suite compiled. The frozen v4 run and
same-host replay both completed; see the [run receipt](../reports/run-20261013T0200Z/README.md).

| check | command | result |
|:--|:--|:--|
| release example build and full v4 execution | `cargo build --release --example isogeny_fiber_pair_20261012`; `target/release/examples/isogeny_fiber_pair_20261012` | passed; 48 verified arm records |
| independent arithmetic, automorphism and isomorphism verification | `rustc --edition 2021 -O -D warnings research/isogeny_torsion_fibers_20261012/paired_round/verify.rs -o /private/tmp/isogeny-paired-verifier`; run it on `raw-v4.tsv` | passed; 48 arm records, 89,764 cubic-automorphism group-law pairs, 44,882 scaling-isomorphism pairs, and pointwise `φ=ψ_u∘[m]`; [receipt](independent_arithmetic_v4.log) |
| Rust formatting | `rustfmt --check research/isogeny_torsion_fibers_20261012/paired_round/verify.rs research/isogeny_torsion_fibers_20261012/paired_round/plot.rs examples/isogeny_fiber_pair_20261012.rs` | passed; [log](rustfmt.log) |
| derived table and graph replay | `rustc --edition 2021 -O -D warnings research/isogeny_torsion_fibers_20261012/paired_round/plot.rs -o /private/tmp/isogeny-paired-plot`; run on `raw-v4.tsv` with temporary outputs, then `cmp` summary, phase summary and SVG against the committed v4 files | all three `cmp` commands passed byte for byte; original files were not rewritten |
| Lean conditional identities and concrete roots/scalings | `lake build` in `research/isogeny_torsion_fibers_20261012/lean/` | passed; [receipt](../../lean/cubic_full_build.log) includes axiom printouts. `lean4checker` is unavailable. Full curve group-action formalization remains open. |
| curve ID selftest | `python3 scripts/curve_id.py selftest` | passed; [log](curve_id_selftest.log). This is registry tooling, not cryptographic research execution. |
| curve names in changed paths | `python3 scripts/check_curve_names.py --diff HEAD` | repository-wide check exits nonzero with 174 existing registry issues; no named model from this round appears in the [staged-path log](curve_names_diff_staged.log). The earlier [HEAD comparison](curve_names_diff_head.log) is retained. |
| touched example test target | `cargo test --release --example isogeny_fiber_pair_20261012` | passed after 28m56s of compilation; the example has zero unit tests; [log](example_test.log) |
| required release library suite | `cargo test --release --lib` | **failed, exit 101:** 3,996 passed, 13 failed, 113 ignored; [full log](release_lib.log) |
| local socket replay outside the restricted sandbox | `cargo test --release --lib cryptanalysis::pollard_collab::cairn::tests::the_client_speaks_enough_http -- --exact --test-threads=1` | passed, 1/1; [log](loopback_replay.log) |
| isolated `jv_cover` timeout replay | `cargo test --release --lib cryptanalysis::jv_cover::tests::nagao_finds_planted_six_sums_and_agrees_with_the_oracle -- --exact --test-threads=1` | passed, 1/1 in 11.77s; [log](jv_cover_isolated.log) |
| isolated `jv_cover` long solver replay | `cargo test --release --lib cryptanalysis::jv_cover::tests::the_solver_loses_no_rational_solution_at_small_p -- --exact --test-threads=1` | passed, 1/1 in 20.52s; [log](jv_cover_isolated_solver.log) |
| isolated stopped F4 replay | `cargo test --release --lib cryptanalysis::jv_cover::tests::stopped_f4_agrees_with_the_full_solver_and_the_oracle -- --exact --test-threads=1` | passed, 1/1 in 29.68s; [log](jv_cover_isolated_stopped.log) |
| isolated replayed F4 replay | `cargo test --release --lib cryptanalysis::jv_cover::tests::replayed_f4_agrees_with_the_full_solver_and_the_oracle -- --exact --test-threads=1` | passed, 1/1 in 31.81s; [log](jv_cover_isolated_replayed.log) |
| affected local networking group outside restricted sandbox | `cargo test --release --lib cryptanalysis::pollard_collab:: -- --test-threads=1` | passed, 37/37; [log](pollard_collab_outside_sandbox.log) |
| full suite with loopback access and limited local concurrency | `RUST_TEST_THREADS=2 cargo test --release --lib` outside restricted sandbox | **passed:** 4,009 passed, 0 failed, 113 ignored; test execution 2,505.94s; [full log](release_lib_unsandboxed_2threads.log) |

The release example test target contains no `#[test]` cases; its value is
compilation of the final source under Cargo's test configuration. The
48-arm execution and independent verifier provide the substantive
correctness checks for this round. The full library gate failed: four
`jv_cover` tests reported timeouts or too few accepted runs under shared
host load, and nine `pollard_collab` tests reported OS error 1,
`Operation not permitted`, at socket setup in the restricted sandbox.
The 37-test affected networking group passed outside that sandbox, and
all four failing `jv_cover` tests passed individually. The controlled
full suite then passed on the same source tree. The first full-suite
failure remains a separate, valid record of the restricted environment;
its result is not rewritten by the later pass. No library, example,
or verifier source changed after the controlled suite began; only
reports and validation receipts were updated. The [problem record](../../../../problems/release-library-validation-shared-host-20261012.md)
tracks this resolution and its environment limit.

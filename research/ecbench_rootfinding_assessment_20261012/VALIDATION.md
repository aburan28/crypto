# Validation log

The run-producing source and specs were committed as `4ffb2c91600a8e82dc7e8790a253e1f0e5b2d523` before execution. The current branch later merged `origin/main`, which changed only two unrelated examples. The study itself changes no Rust source.

## Ecbench and artifact checks

| Check | Outcome |
| --- | --- |
| `CARGO_BUILD_JOBS=4 cargo build --release --locked --bin ecbench` | passed; binary SHA-256 `be956553a5a2cd1227e8f1add8dd24e11dffba5068c799234b55e44ffc3a1b61` |
| `ecbench run` for both frozen specs | both sessions completed; 128/128 verified executions each, 24 measured per arm |
| `ecbench audit --replay-all` for both sessions | both passed; 128 sealed records and 96 identical measured-run replays per session; zero problems |
| `jq -s -f analysis.jq` on each sealed `records.jsonl` | both summaries generated; all four arms in each summary have 24 verified measured rows |
| `typst compile` for `report.typ` and `run_report.typ` | passed; results report 2 pages, run report 1 page; pages visually inspected |
| `git diff --cached --check` before the evidence commit | passed with PDF and PNG paths marked binary by directory-local `.gitattributes` |

## Repository library test

`CARGO_BUILD_JOBS=4 cargo test --release --locked --lib` compiled and ran all 4,123 library tests. Its first, highly parallel sandbox run ended `3997 passed; 13 failed; 113 ignored`. Nine `pollard_collab` network tests failed with `PermissionDenied (Operation not permitted)` while binding local sockets in the sandbox. Four `jv_cover` solver tests failed under the heavily contended parallel run; one receipt explicitly recorded `timed_out: true`, and the others did not satisfy their expected counts. These failures are retained as validation observations, not omitted from the report.

Focused reruns of the failing groups:

| Command | Outcome |
| --- | --- |
| `cargo test --release --locked --lib cryptanalysis::jv_cover::tests::the_solver_loses_no_rational_solution_at_small_p -- --exact --test-threads=1` | 1 passed |
| `cargo test --release --locked --lib cryptanalysis::jv_cover::tests:: -- --test-threads=1` | 14 passed, 1 ignored; includes all four previously failing tests |
| `cargo test --release --locked --lib cryptanalysis::pollard_collab:: -- --test-threads=1` with local loopback access | 37 passed; includes all nine previously failing tests |

A serial whole-suite rerun with loopback access was started and stopped after the four `jv_cover` failures had passed again. It was an optional repetition of the already passing 3,997 tests; it did not produce a whole-suite pass receipt. The publication gate is therefore supported by one complete full-suite run plus passing reruns of every failing test, with the environmental conditions disclosed. CI remains an independent check after the PR opens.

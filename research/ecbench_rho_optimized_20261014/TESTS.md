# Validation

Measured source commit: `2fd065a6c32533a81c1936121e17daa03c4448e8`.

| Check | Result |
|---|---|
| optimized/legacy frozen n=61 profile | pass; scalar, 1,591,865 walk steps, 1,698,247 canonicalizations, and all charged counters identical |
| six-curve ECBench session | pass; 864/864 records verified and 720/720 measured replays identical |
| historical factor-base session | pass; replay-all 20/20 identical |
| historical F4/F5/F6 session | pass; replay-all 48/48 identical |
| historical root-finder session | pass; eight sampled rho rows identical |
| raw Callgrind manifest | pass; every retained part matches `callgrind/SHA256SUMS` |
| every study JSON document | pass with `jq empty` |
| report plot source | pass; Python AST parse succeeds |
| report PDFs | pass; both are two-page letter PDFs and were visually inspected |
| benchmark copies | pass; PDFs and SVGs are byte-equal to the study copies |
| measured-source `rustfmt --edition 2021 --check` on changed Rust files under Rust 1.94.1 | pass |
| measured-source `cargo clippy --locked --all-targets -- -D warnings` under Rust 1.98.1 | pass; warning-free |
| measured-source `cargo test --locked --release --lib` under Rust 1.98.1 | pass; 4,060 passed, 0 failed, 114 ignored; 302.43 s |
| measured-source `cargo test --locked --release --test ecbench` under Rust 1.98.1 | pass; 15 passed, 0 failed; 75.75 s |
| measured-source `cargo test --locked --release --test ecbench_large_prime` under Rust 1.98.1 | pass; 5 passed, 0 failed; 0.05 s |
| post-main-sync `rustfmt --edition 2021 --check` on touched Rust files under Rust 1.94.1 | pass |
| post-main-sync `cargo clippy --locked --all-targets -- -D warnings` under Rust 1.98.1 | pass; warning-free |
| post-main-sync `cargo test --locked --release --lib` under Rust 1.98.1 | pass; 4,083 passed, 0 failed, 114 ignored; 218.89 s |
| final-tree focused ECBench comparison unit tests under Rust 1.98.1 | pass; 2 passed, 0 failed (4,195 filtered out) |
| post-main-sync `cargo test --locked --release --test ecbench` under Rust 1.98.1 | pass; 15 passed, 0 failed; 75.84 s |
| post-main-sync `cargo test --locked --release --test ecbench_large_prime` under Rust 1.98.1 | pass; 5 passed, 0 failed; 0.05 s |
| leaderboard path-check tests | pass; 9/9 |
| leaderboard workflow path coverage | pass; push and pull-request filters cover all 215 builder inputs |
| leaderboard and lab-browser freshness checks | pass; ranked rows unchanged and the rho audit is current |
| synchronized curve-cover catalog | pass; 2,100/2,100 registry curves verified and generated cover artifacts refreshed |
| authored/generated publication whitespace check | pass; raw Callgrind files excluded because Valgrind emits significant trailing spaces and their byte hashes are retained |

The result PDF SHA-256 is
`3e9595054b63fca608f541c49c1133563946bbe0ce246c26faf660622a8a5d7c`;
the run-report PDF SHA-256 is
`3c3a7fe614b51aa043d7b365144c538e8ff3f8f178fbb63fdf83b9f6af2b1171`.

The native-wall acceptance gate remains untested, not waived: Linux PSI was
unavailable, so the preregistered L2 session refused. The retained L0 session
is used only for exact work accounting and descriptive diagnostics.

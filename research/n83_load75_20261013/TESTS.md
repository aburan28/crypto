# Verification commands and status

| Command / check | Status | Evidence |
| --- | --- | --- |
| `rustfmt --check examples/koblitz_orbit_dlp_fast.rs` | PASS | Frozen source format |
| `cargo test --release --locked --offline --example koblitz_orbit_dlp_fast` with pinned `CARGO_TARGET_DIR` | PASS, 12 tests | Focused release test output from pre-run gate |
| `cargo build --release --locked --offline --example koblitz_orbit_dlp_fast` with pinned `CARGO_TARGET_DIR` | PASS | `fc4b62da8e5dc7980c3e2a27e774166d594b5661641833252e8ec8ae6894fa40` binary |
| `sh -n run.sh` and `sh -n run_panel_v4.sh` | PASS | Shell parse gates |
| `/usr/bin/time -p true` | PASS | `timer_preflight.log` |
| `python3 -m py_compile run_with_usage.py` and memory-child preflight | PASS | `usage_preflight.json` |
| `cargo check --manifest-path analyzer/Cargo.toml --locked --offline` | PASS | `analyzer_primary_build.log`, `analyzer_panel.log` |
| `cargo run --manifest-path analyzer/Cargo.toml --locked --offline -- . primary_v2` | PASS: two replayed arms, K2202 interrupted | `analysis_primary_v2.json` |
| `cargo run --manifest-path analyzer/Cargo.toml --locked --offline -- . panel` | PASS: two TIMEOUT rows, no selection | `analysis_panel.json` |

The producer's complete one-target rows and the panel's censored rows are
kept distinct. The live K2202 process was stopped under the recorded local
memory guard. No Tier 0/1 receipt was produced; the AWS request was rejected
before allocation, and RunPod credentials were unavailable in this session.

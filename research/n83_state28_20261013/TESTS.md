# Verification commands and status

| Command or check | Status | Evidence |
| --- | --- | --- |
| `rustfmt --check examples/koblitz_orbit_dlp_fast.rs` | PASS | Frozen source formatting |
| `cargo test --release --locked --offline --example koblitz_orbit_dlp_fast` with pinned target dir | PASS, 13 tests | `example_tests.log` |
| `cargo build --release --locked --offline --example koblitz_orbit_dlp_fast` with pinned target dir | PASS | `producer_build.log`, binary SHA-256 in `PROTOCOL.md` |
| `sh -n research/n83_state28_20261013/run.sh` | PASS | Parse gate before runs |
| `python3 -m py_compile research/n83_state28_20261013/analyze.py` | PASS | Checker parse gate |
| Generic scalar replay on the frozen public point and subgroup generator | PASS | Pre-run checker control |
| `python3 research/n83_state28_20261013/analyze.py research/n83_state28_20261013 control` | PASS, one full replay and one TIMEOUT | `analysis_control.json` |
| `python3 research/n83_state28_20261013/analyze.py research/n83_state28_20261013 primary` | PASS, two resource-guard statuses | `analysis_primary.json` |
| Pinned Lean 4 v4.34.1 theorem build and `#print axioms` | PASS, conditional numeric results | `proof_check.json`, `proof_check.log` |

All execution attempts retain their own launch, status, output and resource
receipts under `runs/`. The two primary resource stops are censored for the
factor-base ranking.

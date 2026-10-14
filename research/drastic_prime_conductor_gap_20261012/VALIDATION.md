# Local validation record

Host: the working clone at `crypto-merge` (macOS, arm64). Recorded
2026-10-12 UTC after the six sessions were sealed and audited.

- `cargo test --release --lib`: 4009 passed, 0 failed, 113 ignored
  (183.64 s). Log retained at the session start time in the task record.
- `cargo build --release --example drastic_prime_gap_pairs`,
  `drastic_gap_specs`, `drastic_gap_analyze`: build clean; the examples
  are main-only (no test functions) and are exercised end-to-end by this
  study (search → specs → analysis reproduced byte-identically on rerun
  of the analyzer).
- `ecbench verify --replay-all` on all six sessions: audit OK, every
  measured record replayed with identical answers and counts. See
  `sessions/*.audit.json`; exact commands in `run_sessions.sh`.
- `python3 scripts/site/test_build.py`: 48 tests, OK.
- `python3 scripts/build_ic_leaderboard.py --check`: current.
- The pinned full-library run had 4009 passed and 0 failed; no fixture
  or N83 sampling failure appeared in this environment.

Environment limitations, retained honestly:

- This host has no CPU isolation topology service, so ecbench recorded
  every run at L0 with `--allow-busy` (the runner prints this and
  degrades gracefully; counted operation totals are unaffected, wall
  figures are not interpreted anywhere).
- Sage was unavailable; certification uses the class-number-one
  conductor argument with order verification by scalar multiplication,
  exactly as the prior round recorded (`sage_used: false`).
- The Conductor coordination service was unreachable at session time
  (connection reset); the check was attempted before editing, as the prior
  round also recorded.

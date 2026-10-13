# Local validation record

Host: the working clone at `crypto-merge` (macOS, arm64). Recorded
2026-10-12 UTC after the seven round-3 sessions were sealed and audited.

- `cargo test --release --lib`: 4011 passed, 0 failed, 113 ignored
  (182.56 s), including the new `quartic_classes_name_the_class_and_solve`
  unit test and the prior rounds' `cubic_classes_name_the_class_and_solve`.
- `cargo check --all-targets`: clean.
- `ecbench verify --replay-all` on all seven sessions: audit OK, every
  measured record replayed with identical answers and counts; A/A exact
  pairs with zero mismatches per session (128/128/64/64 rho pairs).
- The paired target law round-trips through the runner-child protocol
  (`ChildInput.target_kind` now serialises `planted_shared`; the workload
  rebuild check that caught the gap is the session's own first failure
  mode and is retained as the design's guard).
- `python3 scripts/site/test_build.py`: 48 tests, OK.
- `python3 scripts/build_ic_leaderboard.py --check`: current.
- Registry rebuilt (2095 curves — the four new D = -4 pair curves and
  the D = -19 gap-sweep pair registered), covers, leaderboard, browser
  data, panel index and the dashboard regenerated; the round-3 panel is
  registered as retained.

Environment limitations, retained honestly:

- This host has no CPU isolation topology service, so ecbench recorded
  every run at L0 with `--allow-busy`; counted totals are unaffected
  and wall figures are descriptive only. Tier-1 Linux isolation
  (`tools/isolated_bench.py` needs PSI and `taskset`) is the recorded
  step-3 blocker.
- Sage unavailable: class-number-one certification with
  scalar-multiplication order verification, as all prior rounds recorded.
- The Conductor coordination service was unreachable at session time;
  the check was attempted before editing.

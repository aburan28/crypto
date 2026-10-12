# Local validation record

Host: the working clone at `crypto-merge` (macOS, arm64). Recorded
2026-10-12 UTC after the five round-2 sessions were sealed and audited.

- `cargo test --release --lib`: 4008 passed, 2 failed, 113 ignored in
  the full-suite run; both failures are `jv_cover` solver-agreement
  tests that pass their isolated retry (`cargo test --release --lib
  jv_cover`: 14 passed, 0 failed, 1 ignored), the same retry-passing
  class the prior round's validation record documents.
- New unit test `cubic_classes_name_the_class_and_solve` passes: the
  six-class canonicalisation maps every point to the representative its
  multiplier names, the zeta orbit collapses to one class, and
  `rho_walk_with` on `CubicClasses` solves planted targets and reports
  `A = 6`.
- `cargo check --all-targets`: clean.
- `ecbench verify --replay-all` on all five sessions: audit OK, every
  measured record replayed with identical answers and counts; A/A exact
  pairs with zero mismatches per session. See `sessions/*.audit.json`.
- `python3 scripts/site/test_build.py`: 48 tests, OK.
- `python3 scripts/build_ic_leaderboard.py --check`: current.
- Registry rebuilt (2091 curves — the four new D = -3 pair curves
  registered), covers, leaderboard, browser data, panel index and the
  dashboard regenerated; the renderer's `retainedPanels` now includes
  the three conductor-gap panels so future renders preserve them.

Environment limitations, retained honestly:

- This host has no CPU isolation topology service, so ecbench recorded
  every run at L0 with `--allow-busy`; counted totals are unaffected
  and wall figures are descriptive only.
- Sage unavailable: certification is the class-number-one conductor
  argument with scalar-multiplication order verification, as both prior
  rounds recorded.
- The Conductor coordination service was unreachable at session time;
  the check was attempted before editing.
- Round-1 regression fixed in passing: the round-1 dashboard render had
  dropped the `conductor-gap-ic37-20261012` panel (and the round-1
  panel itself was rendered out before its commit); both are restored
  and registered as retained panels here.

# Agent instructions

## Standing workflow rules (from the owner)

- **Always commit finished work and open or update a PR.** Landed rungs,
  ledger promotions, evidence directories, and code changes land as commits
  on a feature branch. If the branch already has an open PR, extend that
  PR and update its title and description to match the new scope instead of
  opening a second one. Do not leave finished work uncommitted.
- **Verify before pushing:** run `cargo test --release --lib`, the touched
  examples' tests, and the relevant Python suites (for example
  `research/sat_factor_base_review_20260908/autolab/test_boundary_autolab.py`)
  and fix failures first.
- **Claim hygiene** (see `docs/ic/boundary_targets.json`): public synthetic
  / known-answer fixtures only; fail-closed evidence; no key-recovery,
  asymptotic sub-rho, or deployed-curve-security claims; multi-target and
  amortized results stay secondary; ledger rows promote only with complete
  measurement fields and independent replay.
- **Conductor:** run `conductor check --summary "…" --scope path:…` before
  editing; report scope expansion before editing outside reserved paths. Do
  not publish chat transcripts or secrets as task metadata.

## Pointers

- Boundary ledger (machine): `docs/ic/boundary_targets.json`
- Boundary scoreboard (human): `docs/ic/BOUNDARY_TARGETS.md`
- Autolab runner: `research/sat_factor_base_review_20260908/autolab/`
- Current ladder frontiers: `RESEARCH_KOBLITZ_INDEX_CALCULUS.md`,
  `RESEARCH_ECC2K130_IC_FEASIBILITY.md`

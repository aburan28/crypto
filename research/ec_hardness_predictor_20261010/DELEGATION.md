# Delegation brief: execute PLAN.md phases 0–3 (then 4–6)

You are the executing agent. Everything you need is in
[`PLAN.md`](PLAN.md); this file is the order of operations and the rules.
Do not rely on any conversation that produced the plan.

## Repositories

- `crypto` worktree: `/Volumes/SSD990/crypto/elliptic-curve-complexity-ml-4b4228`,
  branch `feat/elliptic-curve-complexity-ml-4b4228`. Rust producers live here.
  Read `AGENTS.md` and `CLAUDE.md` first. Run all commands from this
  directory; never `cd` to the main checkout; never bare `git stash`.
- `ml-cryptanalysis`: `/Volumes/SSD990-2/ml-cryptanalysis` (Python learner,
  branch `hardness-predictor`, create it from `main`). Read `docs/PLAN.md`.
- Data goes only under `/Volumes/SSD990-2/ec-hardness/` (`SSD990` is full).

## Before editing

1. `conductor check --summary "ec hardness predictor: <phase>" --scope path:<paths>`.
   At plan time the control plane timed out; if it still does, say so in
   your report and keep edits inside the reserved paths listed in
   `PLAN.md` §12 (last row).
2. Baselines green: `cargo test --release --lib` (crypto, long; run in
   background), `python -m pytest -q` (ml-cryptanalysis).

## Order of work

Phase 0 → 1 → 2 → 3 of `PLAN.md` §11, each with its gate, each committed
with a PR (ready for review, not draft; extend the existing PR instead of
opening a second). Phases 4–6 only after Phase 3's gate passes and after
reporting. Do not start a fleet or GPU run without reporting the
`taskq`/`isolab` availability check first.

## Rules that override speed

- Preserve the plan's parameters (trial counts, `ℓ_fb`, `m`, seeds, tiers).
  If a budget is exceeded, subsample as §6.5/§11 Phase 3 say and record it;
  do not silently shrink.
- Every label row carries its ceiling columns (§4.4) and `status`; failed
  and timed-out rows stay.
- Every exact invariant from Rust is cross-checked against `csd` at toy size
  (Phase 1 gate) before any label is produced.
- Planted controls (§8.4) must pass before any blind result is read.
- No hardness statement about P-256, secp256k1 or ECC2K-130; no "speedup",
  "breakthrough" or exponent language; `S` stays unset.
- Report per `AGENTS.md`: requirement → verified complete / implemented but
  unverified / partial / blocked / not attempted, with commands, commits,
  wall time and the exact remaining gap. Distinguish tool limits from
  mathematical obstructions.

## Report format at the end of each phase

```
Phase N — <name>
gate: pass | fail (evidence path)
requirements table (see above)
commits / PR links
measured numbers that replace a §6.5 estimate
blockers and the next phase's first command
```

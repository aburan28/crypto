# IC boundary autolab

Agent-facing runner for pushing
[`docs/ic/boundary_targets.json`](../../../../docs/ic/boundary_targets.json)
(`schema_version` 2). Public synthetic / known-answer fixtures only.

This replaces the earlier local KIC autolab scripts whose sources were lost
(only `__pycache__` remains under sibling `autolab_*` directories). The Rust
producers `koblitz_rank_fixture` and `koblitz_rho_fixture` are the measurement
engines; this Python control plane owns locks, run directories, ledger pinning,
and fail-closed measurement-schema checks.

## Agent skill

Personal skill (not shipped in this tree): `~/.agents/skills/ic-boundary-autolab/SKILL.md`.
Install or copy into a Cursor skills path if agents should auto-discover it.

## Quick start

```bash
# From repo root
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py plan
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py preflight

# Setup smoke (n=13, 1 fixture) — verifies deps + producers + schema wiring
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py \
  launch --beat smoke.koblitz.vs_rho.n13

# Priority beat: Koblitz whole-process wall-clock at n=37 (start with 1 fixture)
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py \
  launch --beat koblitz.vs_rho.n37_wall

# Same beat, amortized 1024-fixture panel (long)
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py \
  launch --beat koblitz.vs_rho.n37_wall --fixtures 1024

# Alternate priority: charged crossover attempt at n=41
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py \
  launch --beat koblitz.vs_rho.n41_charged

# Scaffold only (writes commands without executing)
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py \
  launch --beat koblitz.vs_rho.n37_wall --prepare-only
```

## Commands

| Command | Purpose |
|---------|---------|
| `plan` | Ledger agent priorities + beat ids + copy-paste launch lines |
| `preflight` | Ledger schema v2, cargo/rustc, producer sources, host pin |
| `launch --beat …` | Lock → run dir → build producers → measure → claim draft + schema check |
| `status [--run-id]` | Print `runs/<id>/state.json` (default: `runs/current.json`) |
| `verify [--run-id]` | Rehash review manifest |
| `claim-check --report PATH --stage STAGE` | Fail-closed measurement-schema validation |

## Layout

```
autolab/
  protocol.json          # beat registry + producer contract
  boundary_autolab.py    # CLI
  README.md
  test_boundary_autolab.py
  runs/
    current.json         # pointer to active run
    autolab.lock         # exclusive lock (pid)
    <run-id>/
      state.json
      inputs/            # pinned protocol + ledger snapshot
      artifacts/         # preflight, claim_draft, claim_check, candidate
      logs/              # producer stdout/stderr
      receipts/          # whole-process wall timings
```

## Measurement schema (fail closed)

Every beat claim must include **all** required fields for its stage from
`docs/ic/boundary_targets.json` → `measurement_schema`, plus global provenance.
`claim-check` and `launch` both enforce this. Incomplete reports get
`SCHEMA_INCOMPLETE` and must not promote a ledger row.

For `vs_rho` the required set includes `timing_class`, `ic_cost`, `rho_cost`,
`automorphism_discount`, `all_stages_charged_same_series`, `verdict`,
`claim_boundary`, and `independent_replay_pointer`.

## Claiming a ledger beat

1. Freeze the public fixture (curve `n`/`a`, seeds, fixture counts, caps).
2. `launch` the beat; keep `artifacts/claim_draft.json` + receipts.
3. Ensure `claim_check.status == PASS` (fill any missing fields first).
4. Independent recomputation on a second process / author.
5. PR that updates `current` → `history`, sets a new `next_target`, and links
   evidence under `research/` or `docs/ic/runs/`.

Do **not** combine best-of-breed component costs from different runs into a
synthetic `vs_rho` win.

## Dependencies

- Rust toolchain (`cargo`, `rustc`)
- Python 3.10+ (3.13 fine)
- Producer sources:
  - `examples/koblitz_rank_fixture.rs`
  - `examples/koblitz_rho_fixture.rs`
- Optional: CryptoMiniSat at `KIC_AUTOLAB_CMS` or `/opt/homebrew/bin/cryptominisat5`
  (needed for SAT-hybrid experiments, not for the direct/rho wall-clock path)

## Related historical dirs

Sibling directories (`autolab_n37_*`, `autolab_n41_*`, `autolab_beats_rho`, …)
retain run artifacts / bytecode from earlier campaigns. Prefer this control
plane for new boundary pushes.

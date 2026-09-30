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

# Optional one-target legacy diagnostic (n=13); the online-speedup gate stays blocked
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py \
  launch --beat smoke.koblitz.vs_rho.n13

# One-target legacy diagnostic at n=37; this cannot claim an online speedup
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py \
  launch --beat koblitz.vs_rho.n37_wall --fixtures 1

# Alternate one-target diagnostic at n=41
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py \
  launch --beat koblitz.vs_rho.n41_charged --fixtures 1

# Scaffold only (writes commands without executing)
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py \
  launch --beat koblitz.vs_rho.n37_wall --prepare-only
```

## Commands

| Command | Purpose |
|---------|---------|
| `plan` | Ledger priorities, target-count gate, and primary-speedup eligibility |
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
      receipts/          # whole-process/cold-start timing diagnostics
  evidence/
    <date>-<topic>/      # bundles promoted out of runs/ for a ledger change
```

`runs/` is gitignored in full. Treat it as scratch: anything a ledger row cites
has to be promoted into `evidence/` first. See "Claiming a ledger beat" below.

## Measurement schema (fail closed)

Every beat claim must include **all** required fields for its stage from
`docs/ic/boundary_targets.json` → `measurement_schema`, plus global provenance.
`claim-check` and `launch` both enforce this. Incomplete reports get
`SCHEMA_INCOMPLETE` and must not promote a ledger row.

The primary `vs_rho` claim is one previously unseen target, paired on the exact
same point and resource envelope. Time IC after reusable setup at its first
target-dependent operation, and rho at its first target-dependent walk. Stop
after independently verified scalar recovery. Exclude launch, input loading,
and known-answer fixture construction. Report candidate, workload, and run IDs
with manifest hashes; use the canonical run-id convention. Report
`rho_online_wall_ms / ic_online_wall_ms`, the exact online stages, and rho's
worker, walk/collision, and distinguished-point memory policy. The five
exclusive IC online phase costs must sum to its online total. Every panel point
is a separate one-target workload row.

The current producers do not expose these online intervals. Their whole-process
wall and operation-counted readings are diagnostic only, and `plan` reports the
primary speedup as blocked until online instrumentation and verification fields
are available. The launcher accepts exactly one fixture and rejects batches;
old multi-target results stay historical.

### Legacy producer timing

The current producers report whole-process wall time or operation-counted
totals. They do not identify the point where target-dependent work begins after
reusable setup, nor do they emit the exclusive online phase record. Keep their
`full_algorithm_charged_total_ms`, setup, and per-target values as labeled
diagnostics only. Do not divide a shared total by the target count, sum rows
that repeat shared setup, or use either value to claim single-target speedup.

### Historical target mode

The old beats pin `target_mode=independent`. Preserve that configuration when
reproducing those diagnostics; the target-mode sweep is not a replacement for a
paired single-target online result. See
[`evidence/20260912-koblitz-vs-rho-no-crossover/sweep_target_modes.py`](evidence/20260912-koblitz-vs-rho-no-crossover/sweep_target_modes.py).

## Claiming a ledger beat

1. Freeze exactly one public target (curve `n`/`a`, point, seeds, and resource envelope).
2. `launch` the beat; keep `artifacts/claim_draft.json` + receipts.
3. Keep current-producer outputs labeled diagnostic-only; do not fill absent
   online fields with estimates. Claim-check can pass only after the producers
   emit the required intervals, matched target, resource policy, and verification.
   A promotable comparison must also include independent replay certificate
   SHA-256 digests and identical, nonempty IC/rho resource-envelope records; a
   boolean flag alone is insufficient.
4. Independently recompute a complete result on a second process / author.
5. Promote the bundle: copy `artifacts/`, `inputs/`, `receipts/` and
   `state.json` into `evidence/<date>-<topic>/`. Leave `logs/` behind — it is
   tens of megabytes and its per-relation records carry factor-base point
   coordinates, target point keys and walk coefficients. Reduce any cost
   component you need out of it into a committed aggregate first.
6. Promote only a complete one-target online comparison. Historical diagnostic
   evidence must stay labeled and cannot set the primary `next_target`.

`runs/*` is gitignored, so a ledger row citing a path under `runs/` is dangling
the moment the directory is cleaned. This is not hypothetical: the `vs_rho`
record retracted on 2026-09-12 cited
`autolab_n37_direct_retry/…/verdict.json`, which appears in no commit, so it was
never independently replayable.

Do **not** combine best-of-breed component costs from different runs into a
synthetic `vs_rho` win, and do not subtract a cost component from one arm
without applying the same policy to the other.

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

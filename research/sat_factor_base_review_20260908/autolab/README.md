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
| `launch-panel --beat … [--targets L] [--blocks B] [--k K \| --k-candidates …] [--cpu N] [--resume RUN_ID]` | Multi-target batch panel beats only (see below); `launch` rejects them |
| `launch-single --beat … [--workloads W] [--tune-workloads T] [--k K \| --k-candidates …] [--cpu N] [--resume RUN_ID]` | One-target online panel beats only (see below) |
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

### Batch panel beats

Beats with `"launch_mode": "batch_panel"` (currently
`koblitz.compact_orbit.n61_panel`: compact-orbit shared-log DLP against the
frozen KS v2 batched rho at a=0 n=61) run through `launch-panel`. One run:

1. generates the tune and eval corpora (`KIC_RHO_GENERATE_ONLY=1`, same
   derivation as the walking rho) and checks they are disjoint;
2. tunes K on the tune corpus of the same L and picks the lowest whole-process
   wall, skipping (and recording) candidates whose estimated IC RSS exceeds
   available memory; `--k` fixes K instead;
3. runs `B` paired blocks in alternating order on the eval corpus, recording
   wall, CPU, RSS, load, available memory and swap per run, and the retired
   instruction count on macOS;
4. checks every target in both arms, then replays every IC and rho record with
   the pure-Python checkers;
5. writes an `end_to_end_dlp` claim draft (`claim_draft.json` /
   `claim_check.json`) and a supplementary `vs_rho` draft
   (`claim_draft_vs_rho.json` / `claim_check_vs_rho.json`).

`L < targets_minimum` (1,024) is refused. The `vs_rho` draft fails closed by
design: the batch producers emit no single-target online fields, and none are
filled in. A panel can reach `PENDING_INDEPENDENT_VALIDATION` under
`end_to_end_dlp`, never a `vs_rho` promotion. With `--cpu N` every producer is
pinned with `taskset`; wrap the whole command in
`tools/isolated_bench.py reserve` (Linux) for AGENTS.md section-10 timings, as
`.github/workflows/compact-orbit-n61-isolated.yml` does.

`--resume RUN_ID` continues an interrupted panel in place. It requires the
same git head, rebuilds the producers and requires identical binary hashes,
and regenerates both corpora byte-for-byte. Finished K rows and blocks are
kept. A partially written tune or block log is renamed to `*.abortedN` and
rerun.

The memory skip uses `ic_rss_model` (`compact_orbit_rss`): 24 B per touched
state (K²·n) plus a 16 B-slot root table sized to `next_pow2(4·K²·n)`, plus
fixed bytes and a headroom margin. It matches the measured peaks within
+0.3% to +4% for K=400..1,200. The old flat 94 B per state under-predicted
above K=800 (7.7 GiB estimated vs 10.0 GiB measured at K=1,200), and it
missed the table doubling at K≈1,482.

### One-target online panel beats

Beats with `"launch_mode": "single_target_panel"` (currently
`koblitz.compact_orbit.n61_single_target`) run through `launch-single`
(`single_target_panel.py`). Every workload is one previously unseen public
target with no known scalar, solved by a fresh process in each arm:

- IC: `examples/koblitz_orbit_dlp_fast_online.rs`, a versioned copy of the
  frozen compact-orbit producer. It builds the base, ranks, indexes and solves
  the linear algebra outside the online interval. It then times the target
  query, PDP probing, relation checks, descent and recovery check as exclusive
  phases that sum to `online_ms`.
- Rho: `examples/koblitz_rho_batch_ks_strong_online.rs`, a versioned copy of
  strong rho, run at rung 3. Its online timer starts after Q is built
  (`KIC_RHO_TARGET_POINT`). It reports walk + collision and the `[d]G = Q`
  check separately.

The frozen originals are untouched, so older evidence keeps verifying.

One run:

1. derives disjoint tune and eval targets by hash-to-curve
   (`target_law.domain`);
2. tunes K on the tune workloads (lowest median whole-process wall);
3. checks that each online producer's untimed output equals its frozen
   original;
4. runs the eval workloads in alternating order, recording resources;
5. writes a replay certificate per arm and workload: an out-of-process Python
   recomputation of the digest, the target, the subgroup membership and
   `[d]G = Q`, plus for IC the relation points and their sum;
6. writes one `vs_rho` claim and claim-check per workload under
   `artifacts/claims/`, `claim_draft.json` / `claim_check.json` for the first
   eval workload, and `single_target_summary.json` with medians and bootstrap
   intervals.

The online speedup is `rho_online_wall_ms / ic_online_wall_ms`. The cold
ratio, process start to verified recovery and including IC setup, is reported
next to it. Both arms share one resource envelope (one worker, one thread, no
memory cap). `independent_validation` comes from the same-host Python replay,
as in `research/ic_single_target_20260930`. An independent-host replay is
still needed before any ledger change.

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

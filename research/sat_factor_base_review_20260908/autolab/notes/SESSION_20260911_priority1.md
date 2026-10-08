# IC boundary session 2026-09-11 (priority #1)

Public synthetic / known-answer only. No ledger promotion.

## Runs

| beat | run_id | status | claim_check | ic_cost | rho_cost | 20% gate |
|---|---|---|---|---:|---:|---|
| smoke n13 | 20260911T032158Z-388c18a6de | PENDING_IV | PASS | (smoke) | (smoke) | n/a |
| n53_probe | 20260911T032209Z-6625ece8b6 | PENDING_IV | PASS | 1.816e6 wall_ms | 8.02e3 | FAIL (IC >> rho; construction probe) |
| n37_wall | 20260911T035251Z-d07cfbf020 | PENDING_IV | PASS | 387.2 | 287.8 | FAIL (ratio ~1.35) |
| n41_charged | 20260911T035344Z-46db668412 | PENDING_IV | PASS | 19413 | 376 | FAIL (ratio ~52) |

## Interpretation

- **n=53**: producers exit 0 past prior ceiling; schema-complete draft. Not a wall crossover.
- **n=37 wall**: IC slower than rho (~1.35×); no ≥20% win.
- **n=41 charged**: IC much slower than rho (~52×); no ≥20% win.

## Next (ledger priorities)

1. Independent validation of the three drafts above (required before any promotion discussion).
2. Priority #2: balanced factor_base at n=53 under resource gate (may need new beat if not in protocol).
3. Or dig into n37 wall gap (why IC wall > rho) before 1024-fixture panels.

Do **not** launch `--fixtures 1024` until a 1-fixture path shows a plausible win.

## Priority #2 — n=53 factor_base (manual)

- Path: `runs_manual/n53_factor_base_20260911/`
- `|F|=19928`, `K=188` orbits, frobenius+negation closed
- `construction_wall_ms` (producer setup) ≈ 10305 ms
- `retained_bytes` ≈ 89.9 MiB; peak RSS ≈ 0.62 GiB (under 16 GiB cap)
- claim-check **PASS** (`claim_draft.json`); verdict `N53_FACTOR_BASE_MEASURED_UNDER_CAP`
- Not a ledger promotion (needs independent validation / PR)

Note: this run also collected RANK_PLUS_32 (whole process ~1030s wall); construction metric uses producer setup total, not full collection wall.

### Construction-only replay (`KIC_SUMMARY_ONLY=1`; `KIC_EXACT_COVERAGE_ONLY` was a no-op)

- Path: `runs_manual/n53_factor_base_construction_only_20260911/`
- Same `|F|=19928`, `K=188`, base hash matches full run
- Producer `total_setup_ms` ≈ 9051; wall-clock ≈ 1157.5s (support/index dominated in wall, not fully reflected in setup fields)
- Peak RSS ≈ 0.62 GiB; claim-check PASS → `N53_FACTOR_BASE_CONSTRUCTION_ONLY_UNDER_CAP`
- Protocol beat: `koblitz.factor_base.n53`

### Control-plane wiring (2026-09-11)

- `boundary_autolab.py` now launches `stage=factor_base` as **direct-only** (no rho).
- `draft_factor_base_claim` + schema validation against `measurement_schema.factor_base`.
- Unit tests cover omit-rho commands and factor_base claim PASS.
- prepare-only run: `20260911T044648Z-d796b6530f`
- Full control-plane launch: `20260911T044704Z-ac2c948295`
  - status `PENDING_INDEPENDENT_VALIDATION`, claim-check **PASS**
  - `|F|=19928`, `K=188`, base_hash matches manual (`e2c0276a…`)
  - `construction_wall_ms` ≈ 19486; retained ≈ 89.9 MiB; peak RSS ≈ 1.89 GiB (under 16 GiB)
  - whole-process wall ≈ 1047 s
  - verdict `DRAFT_PENDING_INDEPENDENT_VALIDATION` (not a ledger promotion)

## Priority #3 — n=31 dim-16 m=2 decomposition

- Tooling: `examples/symmetrised_oracle_bench --only 0 31 2`
- Protocol beat: `koblitz.decomposition.n31_m2_dim16`
- Smoke (1 target, ~148 s wall): **dim V=16**, x-chained F4 32 vars / 31 eqs / deg 2, found@70359 ms, **FFD=3**, gate ok; sym F4 FFD=4 gate ok; disagreements=0
- Smoke claim-check **PASS** → `runs_manual/n31_decomp_m2_dim16_20260911/claim_draft_smoke.json`
- 8-target panel (~612 s wall, peak RSS ≈ 197 MiB): x-chained **5 found / 3 refuted / 0 inconc**, medians 20.9 s found / 75.9 s refuted, **FFD=3**, gate ok; sym gate ok FFD=4; disagreements=0
- Panel claim-check **PASS** → `runs_manual/n31_decomp_m2_dim16_20260911/claim_draft_panel8.json`
- Ledger gates (median≤1h, FFD, disagreements=0, dim16 m=2 quadratic) **met in draft**; not a promotion (needs independent validation)
- Independent replay of panel8 started (`panel_targets8_replay.*`); compare on wake

### Independent replay (panel8)
- Replay EXIT 0; verdict mix / FFD / gate **identical**; found/refuted medians within ~2%
- Receipt: `replay_receipt_panel8.json` → `MATCH_WITHIN_NOISE`

## Priority #4 — binary decomposition ℓ=8 (probe)

- `semaev_decomp_bench 8`: pairs@ℓ=8 = **0.068 s** (8.9× vs triples/gen); still Θ(2^(2ℓ)) baseline
- Claim draft: `runs_manual/binary_decomp_l8_20260911/claim_draft_pairs_baseline.json`
- Ledger gates **not met**: no sub-2^(2ℓ) oracle, **FFD unmeasured**
- Next: SAT (or other) vs pairs at ℓ=8 + subspace FFD ladder
- Full-field FFD sweep n=3..7: **FFD=3** everywhere → `ffd_fullfield_summary.json` (does not close l=8 subspace FFD gate)
- Priority #4 ledger beat remains **blocked** on sub-2^(2ℓ) DoF (pairs is the 2^(2ℓ) baseline; SAT historically slower)

## Priority #5 — prime j=0 end_to_end_dlp 16-bit

- Tooling: `examples/j0_known_answer_bits`
- 14-bit control: IC≈363 ms, ρ≈5.1 ms, agree ✓
- **16-bit**: IC≈1851 ms, ρ≈20.5 ms, ic_agrees_rho ✓, matches truth ✓
- claim-check → see claim_check_16bit.json; not a promotion (needs independent validation + fuller stage timers)
- Independent replay: second 16-bit draw also ic_agrees_rho ✓ → `replay_receipt_16bit.json`

## Priority #6 — relation_yield probe (binary Koblitz n=19)

- eta-sweep via rank_fixture SUMMARY: **no coverage fields** (instrumentation miss)
- `koblitz_factor_base_yield 19 1 256`: m=2 mean coverage=0.0947, max=0.3086, identities=12
- claim-check PASS → `runs_manual/relation_yield_koblitz_n19_20260911/claim_draft_yield.json` (first publication probe; CI still weak)

## Priority #7 — FFD backfill

- pairs ℓ=8: FFD **inapplicable** (allowed) + full-field S3 FFD=3 pointer; claim-check PASS
- n31 m=2: already complete
- **Prime Fp FFD (gbrl)**: `research/gbrl/src/ffd.rs`; ladder fb=3..6 → FFD 6/8/10/12 (mean 9.0)
  - claim-check PASS → `runs_manual/prime_m3_ffd_20260911/`
- **Prime ≥14-bit m=3 + witnesses + FFD**: planted 3-sum on y²=x³+x+1 / F_1000003
  - `|E|=1000727` (20 bits ≥14), FFD=8, GB 98 steps / ~14 ms
  - planted vanishes + group-sum OK; recovered planted up to perm; 6 group-verified witnesses
  - claim-check **PASS** → `runs_manual/prime_m3_witness_ffd_20260911/`
  - Not a promotion (toy non-j=0 curve; single planted draw)
- Open: binary subspace SAT FFD; multi-draw FFD stats; j=0 curve variant
- Inventory: `runs_manual/ffd_backfill_20260911/backfill_inventory.json`

## Loop status (2026-09-11T16:00Z)

| Pri | Beat | Result |
|---|---|---|
| 1 | vs_rho n37/41/53 | measured, no 20% win |
| 2 | factor_base n53 | PASS draft + control-plane |
| 3 | n31 m2 dim16 decomp | PASS + replay MATCH |
| 4 | binary sub-2^(2ℓ) | **blocked** (pairs=baseline) |
| 5 | j0 16-bit e2e | PASS + replay agree ρ |
| 6 | relation_yield n19 | PASS probe (weak CI) |
| 7 | FFD backfill | Fp FFD + ≥14-bit planted m=3+witnesses **PASS drafts**; SAT FFD remains |

Next DoF: binary subspace SAT FFD, or stop loop (p4 hard-blocked).

# FFD/DoR backfill inventory (2026-10-12)

Priority **#7** from the IC boundary autolab skill / `BOUNDARY_TARGETS.md`
("backfill FFD/DoR on algebraic decomposition claims that only cite wall-clock
/ conflicts"). Status **`PENDING_INDEPENDENT_VALIDATION`**. **No ledger
promotion.** Public synthetic / known-answer evidence only.

`docs/ic/boundary_targets.json` `agent_priorities` currently lists ranks 1–6
only; this inventory follows the skill's explicit rank-7 backfill item and the
scoreboard wording that asks to fill FFD/DoR where only wall or conflicts are
cited.

## Ledger gaps found

| Ledger path | Ledger status | What already fills it | Still open |
|---|---|---|---|
| `regimes.binary…decomposition…ffd.subspace_sat_pairs_ladder` | `unmeasured` (conflicts/wall only) | **pairs:** `binary.decomposition.l8_pairs_baseline` → FFD `inapplicable` (claim-check PASS). **Subspace algebraic:** `runs_manual/binary_s4_ffd_l8_20260920/` → FFD measured min=max=mean=7 on x1-fixed (16 draws); claim-check PASS | No registered beat translates historical **SAT** conflict/wall rows into Macaulay FFD |
| `regimes.prime…decomposition.history[0]` (S3 2-decomp pre-S4) | `unmeasured` (wall/α fits only) | Current S4 3-decomp record: FFD `inapplicable` (CZ root finding). Adjacent measured FFD: `prime_m3_ffd_20260911`, `prime_m3_witness_ffd_20260911` (m=3, not S3 2-decomp) | No `protocol.json` beat for S3 2-decomp FFD/DoR; do **not** invent one |

Machine-readable map: [`inventory.json`](inventory.json).

## Claim-check re-runs (this host, 2026-10-12)

Replayed with:

```bash
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py \
  claim-check --report PATH --stage decomposition
```

| Draft | Result file | Status |
|---|---|---|
| `runs_manual/binary_decomp_l8_20260911/claim_draft_pairs_baseline.json` | `binary_decomp_l8_20260911.json` | PASS |
| `runs_manual/binary_s4_ffd_l8_20260920/claim_draft_s4_subspace_ffd.json` | `binary_s4_ffd_l8_20260920.json` | PASS |
| `runs_manual/prime_m3_ffd_20260911/claim_draft.json` | `prime_m3_ffd_20260911.json` | PASS |
| `runs_manual/prime_m3_witness_ffd_20260911/claim_draft.json` | `prime_m3_witness_ffd_20260911.json` | PASS |

No new producer measurement was launched: every fill above is an existing
registered draft / manual evidence directory already on `main`. The control
plane has no launchable `relation_yield` invent-beat and no S3-2decomp-FFD
beat in `protocol.json`.

## Exact next registered step (if more FFD is required)

1. **Binary SAT conflict rows:** only after a producer is registered that
   builds a Boolean Macaulay system for the SAT-frontier Semaev instances and
   a beat is added to `protocol.json` — then `claim-check --stage decomposition`.
2. **Prime S3 2-decomp history cell:** register an FFD harness for that arm
   (or document `ffd.status=inapplicable` with reason if it is pure
   root-finding without a GB system) before re-measuring; do not invent a
   beat outside `protocol.json`.
3. Ledger row edits stay out of scope here (`PENDING_INDEPENDENT_VALIDATION`).

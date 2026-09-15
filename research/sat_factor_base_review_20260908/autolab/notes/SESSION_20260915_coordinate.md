# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T14:05:01Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, instances **1**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| orbit_factorized planted n13→53 | edge-free SAT through n53 |
| unrestricted n13 | finds models; n19+ cliff remains |
| pair-support nogoods | sound; do not close n19 cliff |
| **factor_base n53** `20260915T122253Z-229f0ee245` | claim_check PASS; **CONSTRUCTION_INDEPENDENTLY_REPLAYED**; no ledger promotion |
| pairing/seed/branch sweeps | n19 cliff persists |
| **parallel nogood install** | rayon scan+expand; **n53 planted+nogoods SAT** (~52s scan, 39856 clauses, pair_table=0) |
| lazy_roots+nogoods n53 | CDCL-heavy; aborted >13min — prefer full S3 |
| formula accounting | n13/19/53 incl. n53+parallel_nogoods |

## Priority signals

1. Edge-free planted extraction through n53 stands **with** compressed pair-support nogoods.
2. factor_base n53 construction independently replayed — still not a ledger promotion.
3. Unrestricted n19 remains open; half-trace lazy roots hurt large-n CDCL.

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Phase-save / restart heuristics for unrestricted n19 (pair_table=0, full S3 + nogoods).
2. Optional binary l=8 pairs baseline; do not promote ledger rows.
3. Re-checkout IC branch if worker flipped to rho; never mutate AWS.

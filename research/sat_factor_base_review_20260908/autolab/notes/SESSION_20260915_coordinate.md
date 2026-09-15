# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T13:32:43Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, instances **2**
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
| pair_then_pair / reversed | n13 ok; **n19 still 0 models** |
| alternate pairing∈{0,1,2} + seeds{1,2,3,5} | **n19 cliff persists** |
| formula/memory accounting | n13/19/53 orbit sizes logged; **pair_table=0**; planted n53 pair_then_pair SAT |
| n53+nogoods install | CPU-bound >6min (engineering note) |
| n31 decomp | prior IV MATCH_WITHIN_NOISE — skip re-run |

## Priority signals

1. Edge-free planted extraction through n53 stands (incl. pair_then_pair).
2. factor_base n53 construction independently replayed — still not a ledger promotion.
3. Unrestricted n19 remains open after pairing/seed/branch-order sweeps.

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Optimize compressed-nogood install at n53, or phase-save/restart heuristics for n19 (pair_table=0).
2. Optional: binary l=8 pairs baseline; do not promote ledger rows.
3. Re-checkout IC branch if worker flipped to rho; never mutate AWS.

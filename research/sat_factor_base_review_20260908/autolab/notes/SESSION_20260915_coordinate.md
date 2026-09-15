# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T14:12:46Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, instances **2**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| orbit_factorized planted n13→53 | edge-free SAT through n53 |
| parallel nogood install | n53 planted+nogoods SAT (~52s, 39856 clauses, pair_table=0) |
| **factor_base n53** | claim_check PASS; CONSTRUCTION_INDEPENDENTLY_REPLAYED; no ledger promotion |
| **phase-restart shots** (`KIC_PHASE_RESTART_SHOTS`) | n13 ok; **n19 still 0 models** (seeds 1/7/11, 80k–120k) |
| binary l=8 pairs draft | claim-check PASS (schema only); not remeasured; no ledger promotion |
| unrestricted n19 cliff | remains open |

## Priority signals

1. Edge-free planted n53 with compressed nogoods stands.
2. factor_base n53 independently replayed — not a ledger promotion.
3. Unrestricted n19 survives phase-scramble multi-shot restarts.

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Try smaller-eta orbit bases or alternate decision priorities for unrestricted n19 (pair_table=0).
2. Optional fresh binary l=8 pairs remount only if asked; do not promote ledger.
3. Re-checkout IC branch if worker flipped to rho; never mutate AWS.

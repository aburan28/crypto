# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T13:13:22Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, instances **2**
- Running tagged workers observed: multiple `g7.2xlarge` + `g7e.2xlarge` (+ vivado/f2/ingest)
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| orbit_factorized planted n13→53 | edge-free SAT through n53 |
| unrestricted n13 | finds models; n19+ cliff remains |
| pair-support nogoods | sound; do not close n19 cliff |
| **factor_base n53** `20260915T122253Z-229f0ee245` | claim_check PASS; **CONSTRUCTION_INDEPENDENTLY_REPLAYED** (F=19928,K=188,hash match); no ledger promotion |
| lazy_pair_roots static fix | empty body → `constrain_s3_root`; theory opt-in |
| **pair_then_pair branch order** | planted n13 SAT; unrestricted n13 finds ≥1 model @100k; **n19 still 0 models @150k** |

## Priority signals

1. Edge-free planted extraction through n53 stands.
2. Ledger beat **factor_base n53** construction independently replayed — still not a ledger promotion.
3. Unrestricted n19 search still open; `KIC_ORBIT_BRANCH_ORDER=pair_then_pair` does not close the cliff under tested budgets.

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Alternate pairing / phase-block heuristics for unrestricted n19 (keep pair_table=0).
2. Optional ledger beats: n31 decomp distribution; do not promote factor_base without explicit ledger edit after IV review.
3. Re-checkout IC branch if worker flipped to rho; never mutate AWS.

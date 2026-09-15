# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15 ~12:27Z)

- Slots **12**, alive **10**, retired 2, errors 0
- Aggregate **~106.4 B it/s**, ~0.346% of 2^60.9
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| orbit_factorized planted n13→53 | edge-free SAT through n53 |
| unrestricted n13 | finds models; n19+ cliff remains |
| pair-support nogoods | sound; do not close n19 cliff |
| **factor_base n53** `20260915T122253Z-229f0ee245` | **claim_check PASS**, PENDING_IV; F=19928, K=188, ~2.2s, ~74.9MB |
| lazy_pair_roots static fix | empty body → `constrain_s3_root`; theory opt-in |
| lazy-static + nogoods unrestricted n13 | UNKNOWN 0 models @150k (half-trace selectors hurt CDCL) |

## Priority signals

1. Edge-free planted extraction through n53 stands.
2. Ledger beat **factor_base n53** has a schema-PASS draft awaiting independent validation.
3. Unrestricted n19 search still open; full S3 + compressed nogoods preferred over half-trace selectors for CDCL.

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Independent validation replay of factor_base n53 draft.
2. Search heuristics / pair-then-pair decomposition for unrestricted n19 (keep pair_table=0).
3. Re-checkout IC branch if worker flipped to rho; never mutate AWS.

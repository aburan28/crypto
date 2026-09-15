# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T15:44:28Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, instances **2**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| planted orbit S5 n13→53 (full S3) | edge-free SAT; +parallel nogoods |
| reusable pair-support certs | export/reload; n53 scan ~48s→0 |
| **multi-target cert harvest** | n19 planted 5/5 group-valid; n53 planted 3/3; cert-shared digests; pair_table=0 |
| n13 natural + cert | edge-free loaded; hits=0/5 (no promotion) |
| unrestricted n19 cliff | still open |

## Priority signals

1. Edge-free planted extraction through n53 stands under full S3.
2. Reusable pair-support certificates now amortize across multi-target planted harvests without edge tables.
3. Unrestricted n19 cliff remains open.

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Positive regular-support without joins, or structural n19 rewrite.
2. Keep full S3 for planted n53 edge-free relations.
3. Never mutate AWS from this agent.

# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T15:30:14Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, instances **2**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| planted orbit S5 n13→53 (full S3) | edge-free SAT; +parallel nogoods at n53 |
| factor_base n53 | independently replayed; **no ledger promotion** |
| unrestricted n19 cliff | still open under heuristic/combo sweeps |
| **reusable pair-support cert** | export/reload at n13/19/53; planted SAT preserved; **n53 scan 48s→0**; pair_table=0 |
| binary l8 / j0-16 IVs | MATCH / PASS known-answer; **no promotion** |

## Priority signals

1. Edge-free **full-S3** planted extraction through n53 stands.
2. Reusable compressed relative-Frobenius pair-support certificates now exist **without** edge/pair tables (ledger #1/#2 direction).
3. Unrestricted n19 cliff remains open under tested encodings.

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Wire certs into multi-target harvest / positive regular-support without joins.
2. Keep full S3 for planted n53 edge-free relations.
3. Never mutate AWS from this agent.

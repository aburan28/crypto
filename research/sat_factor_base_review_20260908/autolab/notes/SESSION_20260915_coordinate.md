# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T14:44:57Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, instances **2**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| planted orbit S5 n13→53 (full S3) | edge-free SAT; +parallel nogoods at n53 |
| factor_base n53 | independently replayed; **no ledger promotion** |
| unrestricted n19 cliff | persists under eta/branch/phase levers |
| **lazy pair roots planted** | n13 SAT; **n23 UNKNOWN @100k** (full S3 SAT @61); n53 search aborted >8min |
| lazy formula size | ~2× vars/clauses vs full S3; pair_table=0 |

## Priority signals

1. Edge-free **full-S3** planted extraction through n53 stands (with/without nogoods).
2. Lazy pair roots are **not** currently a viable path to group-valid n≥23 planted relations under tested budgets.
3. Unrestricted n19 cliff remains open.

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Keep full S3 for n53 edge-free relations; treat lazy roots as accounting/experimental only until search recovers.
2. New unrestricted-n19 lever or divert to other ledger beats without promotion.
3. Never mutate AWS from this agent.

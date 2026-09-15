# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T18:00:18Z)

- ASG `ecc2k130-g7`: desired/InService **6/6**
- ASG `ecc2k130-workers`: desired **4**, InService **0**, instances **0**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| planted n53 full-S3 + cert | SAT, pair_table=0 |
| **intermediates_then_pair + pair_sum + relative-positive** | **n19 unrestricted/natural SAT group-valid (draft)** |
| ablations (drop either lever) | UNKNOWN @80k |
| n23 unrestricted same stack | UNKNOWN @100k |
| CMS5 on prior combo (pair_then_pair) | INDETERMINATE @150k |
| j0-16 3-shot IV | PASS |

## Priority signals

1. Edge-free planted n53 stands under full S3.
2. n19 unrestricted cliff **draft-broken** under intermediates_then_pair+pair_sum+relative-positive (needs external replay); n23 still open.
3. Keep full S3 for planted n≥19 / n53 (lazy unsuitable).
4. Skip n37/n53 symmetrised mid-dim (blocked).

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Independent CMS/external replay of the n19 intermediates_then_pair win; push n23/n53 budgets.
2. Keep claim hygiene fail-closed — draft ≠ ledger promotion ≠ vs_rho.
3. Never mutate AWS from this agent.

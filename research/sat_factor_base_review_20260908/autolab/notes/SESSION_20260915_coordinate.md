# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T16:58:00Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, InService **0**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| planted n53 full-S3 + cert (tick sanity) | SAT, pair_table=0, invalid_lifts=0 |
| lazy + relative-positive | fails planted n19 — keep full S3 |
| **eta=1/16 + relative-positive** | planted SAT (smaller formula); natural UNKNOWN @100k |
| eta=1/16 + relative + symmetry_break | natural UNKNOWN @80k |
| n41 dim21 @budget5k | still inconclusive |

## Priority signals

1. Edge-free planted n53 stands under full S3 + reusable nogood cert.
2. Smaller η shrinks n19 relative-positive formulas but does not close the unrestricted cliff.
3. Lazy roots remain unsuitable for planted n≥19 / n53.
4. Skip n37/n53 symmetrised mid-dim (blocked).

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Non-η structural n19 rewrite, or much larger n41 F4 budget.
2. Keep full S3 for planted n53 edge-free relations.
3. Never mutate AWS from this agent.

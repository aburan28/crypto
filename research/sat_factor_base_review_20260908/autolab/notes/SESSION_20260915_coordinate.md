# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T18:23:17Z)

- ASG `ecc2k130-g7`: desired/InService **6/6**
- ASG `ecc2k130-workers`: desired **4**, InService **0**, instances **0**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| intermediates_then_pair stack | n19 draft-SAT (multi-seed) |
| **n19 IV @300k** | planted 1–5,7,11 + natural 1–3 all SAT group-valid |
| **n23 seed1 @300k** | **SAT @168284** (replay identical); seeds 2/3 + natural UNKNOWN |
| CMS without branch heuristic | still INDETERMINATE @200k |

## Priority signals

1. Edge-free planted n53 stands under full S3.
2. n19 unrestricted cliff **draft-broken with budget-robust multi-seed IV**; n23 seed1 draft-SAT, hit-rate incomplete.
3. Keep full S3 for planted n≥19 / n53 (lazy unsuitable).
4. Skip n37/n53 symmetrised mid-dim (blocked).

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Raise n23 budgets / more seeds; chart path to n53; keep fail-closed hygiene.
2. Draft ≠ ledger promotion ≠ vs_rho ≠ key recovery.
3. Never mutate AWS from this agent.

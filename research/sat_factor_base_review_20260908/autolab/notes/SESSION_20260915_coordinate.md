# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T17:43:30Z)

- ASG `ecc2k130-g7`: desired/InService **6/6**
- ASG `ecc2k130-workers`: desired **4**, InService **0**, instances **0**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| planted n53 full-S3 + cert | SAT, pair_table=0 |
| eta=1/16 relative-positive | planted SAT; natural UNKNOWN |
| CMS5 on eta=1/16 relative-positive | natural INDETERMINATE @300k; planted SATISFIABLE |
| pair_sum_trie intermediate domain | planted-unit SAT; unrestricted n19 UNKNOWN @100k |
| **pair_sum_trie + relative-positive combo** | planted-unit SAT (14 conf); unrestricted/natural UNKNOWN |
| **n31 third-seed FFD (424242)** | FFD x=3/sym=4 stable; gate ok; x found 4/4; 0 inconclusive |
| n41 dim21 @budget20k | still inconclusive (effort 10005) |

## Priority signals

1. Edge-free planted n53 stands under full S3.
2. n19 unrestricted cliff persists under internal CDCL and external CMS at tested budgets.
3. Keep full S3 for planted n≥19 / n53 (lazy unsuitable).
4. Skip n37/n53 symmetrised mid-dim (blocked).

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Deeper non-domain n19 rewrite, or factor_base/relation_yield / j0 IV beats.
2. Keep full S3 for planted n53 edge-free relations.
3. Never mutate AWS from this agent.

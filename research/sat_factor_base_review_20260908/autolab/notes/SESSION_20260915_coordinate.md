# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T16:45:00Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, InService **1**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| planted orbit S5 n13→53 (full S3) | edge-free SAT |
| relative-frame positive | planted SAT; n19 cliff open |
| n41 dim21 @budget2k | constructible; inconclusive |
| **lazy + relative-positive** | n13 planted SAT; **n19 planted/natural UNKNOWN** — do not use lazy for n53 |
| **n41 dim21 @budget5k** | still inconclusive (effort 2504, ~349s) |

## Priority signals

1. Edge-free planted n53 requires **full S3** (lazy+relative does not carry planted n19).
2. n41 mid-dim cell needs substantially larger F4 budgets before claim-grade panels.
3. Unrestricted n19 cliff remains open.
4. Skip n37/n53 symmetrised mid-dim (geometrically blocked).

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. New structural n19 rewrite (not lazy), or n41 with much larger node budget / fewer arms.
2. Keep full S3 for planted n53 edge-free relations.
3. Never mutate AWS from this agent.

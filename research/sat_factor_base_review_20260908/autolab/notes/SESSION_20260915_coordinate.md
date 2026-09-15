# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T16:23:00Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, instances **2**, InService **1**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| planted orbit S5 n13→53 (full S3) | edge-free SAT |
| relative-frame positive | planted SAT; n19 cliff open |
| n31 distributional redraw | FFD/gate stable (not promotion) |
| **n41 dim21 resource probe** | constructible (~2.10M FB); 1-target F4 @budget2000 inconclusive ~314s |
| n19 relative+binary+fieldbits @60k | UNKNOWN |

## Priority signals

1. Edge-free planted extraction through n53 stands under full S3.
2. n41 is the correct geometric successor to n31 mid-dim; needs larger budgets before claim-grade panels.
3. Unrestricted n19 cliff remains open.
4. Do not burn budget on n37/n53 symmetrised mid-dim (blocked).

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Higher-budget n41 multi-target F4 panel, or new structural n19 rewrite.
2. Keep full S3 for planted n53 edge-free relations.
3. Never mutate AWS from this agent.

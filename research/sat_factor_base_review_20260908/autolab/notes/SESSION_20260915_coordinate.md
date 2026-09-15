# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T16:04:00Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, instances **2**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| planted orbit S5 n13→53 (full S3) | edge-free SAT; +parallel nogoods |
| relative-frame positive | planted SAT; n19 clauses shrink vs absolute; cliff open |
| n37 quadratic cell | geometrically blocked |
| **n31 distributional redraw** | seed 66142 ×4 F4: FFD 3/4 stable, gate ok |
| symmetrised dim geometry map | n41/n47 mid-dim candidates; n37/n53 blocked |
| n19 shift_then_reps + relative-positive @40k | UNKNOWN |

## Priority signals

1. Edge-free planted extraction through n53 stands under full S3.
2. n31 mid-dim cell FFD/gate stable under a second-seed distributional redraw (not promotion).
3. Unrestricted n19 cliff remains open.
4. Advance quadratic cell past n31 should target n41/n47 geometrically, not n37/n53.

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Structural n19 rewrite, or cautious n41 dim21 resource probe, or larger n31 distribution.
2. Keep full S3 for planted n53 edge-free relations.
3. Never mutate AWS from this agent.

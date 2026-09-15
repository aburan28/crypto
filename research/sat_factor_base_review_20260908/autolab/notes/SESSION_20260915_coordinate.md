# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T15:52:25Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, instances **2**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| planted orbit S5 n13→53 (full S3) | edge-free SAT; +parallel nogoods |
| reusable pair-support certs | export/reload; n53 scan ~48s→0 |
| multi-target cert harvest | n19 planted 5/5; n53 planted 3/3; pair_table=0 |
| absolute positive regular support | planted SAT; unrestricted n19 still UNKNOWN; ~10× clauses |
| **relative-frame positive** | mux+canonical; planted SAT; n19 clauses 90031→53857; cliff still open |
| n37 quadratic cell | **geometrically blocked** (Phi_37 deg 36 over F2 → only dim=1) |

## Priority signals

1. Edge-free planted extraction through n53 stands under full S3.
2. Relative-frame positive is a real structural shrink vs absolute expansion but does not close unrestricted n19.
3. n37 symmetrised dim~19 cell is unavailable; do not burn budget there.
4. Unrestricted n19 cliff remains open.

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Other structural n19 rewrite (or distributional n31 / usable-dim decomp n).
2. Keep full S3 for planted n53 edge-free relations.
3. Never mutate AWS from this agent.

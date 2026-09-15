# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T15:05:33Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, instances **2**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| planted orbit S5 n13→53 (full S3) | edge-free SAT; +parallel nogoods at n53 |
| factor_base n53 | independently replayed; **no ledger promotion** |
| unrestricted n19 cliff | persists under eta/branch/phase/karatsuba/natural/cms/field-bits |
| lazy pair roots planted | n13 SAT; **n23 UNKNOWN @100k**; n53 abort >8min |
| karatsuba mul | n13 OK; **n19 still 0 @80k** |
| natural targets | n13 0@40k / n19 0@60k |
| **CMS emit_xor** | planted n13 SAT; unr n13 SAT@250k; **n19 INDET @300k** |
| **field_bits + PHASE_INIT** | n13 OK; **n19 still 0 @80k** (phase 0 and 1) |
| n53 planted reconfirm | valid=1, invalid_lifts=0, pair_table=0 |

## Priority signals

1. Edge-free **full-S3** planted extraction through n53 stands — group-valid public-synthetic relation under planted units.
2. Lazy pair roots are **not** currently a viable path to group-valid n≥23 planted relations under tested budgets.
3. Unrestricted n19 cliff remains open; search-heuristic levers (incl. CMS / field-bits / phase-init) do not close it.

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Keep full S3 for n53 edge-free relations; treat lazy roots as accounting/experimental only.
2. Prefer a **structural** encoding change for unrestricted n19, or divert to other ledger beats without promotion.
3. Never mutate AWS from this agent.

# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T15:19:03Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, instances **2**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| planted orbit S5 n13→53 (full S3) | edge-free SAT; +parallel nogoods at n53 |
| factor_base n53 | independently replayed; **no ledger promotion** |
| unrestricted n19 cliff | persists under broad heuristic + binary-rep/field-bits combo |
| lazy pair roots planted | n13 SAT; **n23 UNKNOWN @100k**; n53 abort |
| binary l8 pairs IV | MATCH_WITHIN_NOISE; **no promotion** |
| **binary-rep + field_bits** | n13 OK; **n19 still 0 @80k** (pair_table=0) |
| **j0 16-bit IV (3-shot)** | **PASS_KNOWN_ANSWER_AGREES_RHO**; not vs_rho crossover; **no promotion** |

## Priority signals

1. Edge-free **full-S3** planted extraction through n53 stands.
2. Unrestricted n19 needs a deeper structural encoding rewrite; search/combo levers are exhausted under tested budgets.
3. Prior PASS drafts continue to independently validate without promotion (binary l8, j0-16).

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Prototype reusable pair-support without materializing an edge/pair table (ledger priority #1/#2), or divert to another non-promotion IV.
2. Keep full S3 for planted n53 edge-free relations.
3. Never mutate AWS from this agent.

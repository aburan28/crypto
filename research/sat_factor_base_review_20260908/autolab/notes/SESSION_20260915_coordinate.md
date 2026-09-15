# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T15:14:37Z)

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
| karatsuba / natural / CMS / field_bits | do not close n19 |
| **binary l8 pairs IV** | **MATCH_WITHIN_NOISE** (0.068→0.066s); claim-check PASS; **no promotion** |
| n53 planted witness compact | valid=1, invalid_lifts=0, pair_table=0 |

## Priority signals

1. Edge-free **full-S3** planted extraction through n53 stands — group-valid public-synthetic relation under planted units.
2. Lazy pair roots / search-heuristic levers do **not** currently unlock unrestricted n19 or group-valid n≥23 planted lazy paths under tested budgets.
3. Binary l8 pairs baseline independently remeasured (MATCH_WITHIN_NOISE); still not a sub-2^(2l) claim and **not** ledger promotion.

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Prefer a **structural** encoding change for unrestricted n19 (heuristics exhausted).
2. Optional: other ledger beats without promotion (n37 wall smoke, prime j0 IV) if S5 blocked.
3. Never mutate AWS from this agent.

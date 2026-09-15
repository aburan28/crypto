# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T14:20:28Z)

- ASG `ecc2k130-g7`: desired/alive **6/6**; `ecc2k130-workers`: desired **4**, instances **2**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| planted orbit S5 n13→53 | edge-free SAT; +parallel nogoods at n53 |
| factor_base n53 | claim_check PASS; independently replayed; **no ledger promotion** |
| phase-restart / pairing / seed sweeps | **n19 unrestricted cliff persists** |
| smaller eta (1/16) | n19 base shrinks (1 col); planted SAT; unrestricted still 0 models |
| `shift_then_reps` branch order | n13 finds models; n19 still 0 |
| n23 | planted SAT; unrestricted UNKNOWN @60k |
| n17 | not constructible for a=0 |

## Priority signals

1. Edge-free planted extraction through n53 with compressed nogoods stands.
2. factor_base n53 independently replayed — not a ledger promotion.
3. Unrestricted search cliff is real by n19 under all tested search levers (pair_table=0).

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. New encoding/decision lever beyond polarity/eta/branch-order, or divert to other ledger beats without promotion.
2. Re-checkout IC branch if worker flipped to rho; never mutate AWS.

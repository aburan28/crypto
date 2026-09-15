# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15 ~12:20Z)

- Bucket `s3://ecc2k130-590183823895`
- Slots 8, alive **6**, retired 2, errors 0
- Aggregate **~86.9 B it/s**, ~0.34% of 2^60.9
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111` (keep IC work off the rho PR branch;
sibling may flip the shared worktree — re-checkout this branch each tick).

## Runs this session

| Beat / probe | Path / run_id | Result |
|---|---|---|
| smoke `vs_rho` n13 | `runs/20260915T110615Z-cbfe798888` | PASS / PENDING_IV |
| `vs_rho` n37_wall 1fx | `runs/20260915T110708Z-6bac5fa80e` | PASS draft; IC ≫ ρ wall |
| relative Frobenius pair support n13→53 | `runs_manual/relative_pair_stats_20260915/` | compression = n; 0 equivariance misses |
| j0 16-bit e2e | `runs_manual/prime_j0_e2e_16bit_20260915/result.json` | ic_agrees_rho ✓ |
| orbit_factorized planted ladder n13→53 | `runs_manual/orbit_factorized_s5_20260915/` | **edge-free SAT** through n53 |
| orbit_factorized n13 unrestricted 1M | same | SAT, **8** valid models, then UNKNOWN |
| orbit_factorized n19/n23 unrestricted | same | **UNKNOWN**, 0 models |
| **pair-support nogoods** (one_hot) | `KIC_ORBIT_PAIR_SUPPORT_NOGOODS=1` | planted n13 **SAT**; unrestricted n13 finds models; **n19 still UNKNOWN** (76 nogoods, 0 models @200k) |
| `koblitz.factor_base.n53` | `runs/20260915T122253Z-229f0ee245` | ACTIVE measurement (`koblitz_rank_fixture` live) |

## Priority #1 signal

1. Edge-free **planted** extraction works through **n=53** (0 pair table).
2. Unrestricted search cliff remains at **n19+**.
3. Compressed relative-Frobenius **pair-support nogoods** land and are sound
   (planted SAT), but **do not** close the unrestricted n19 cliff alone.

Claim boundary: planted units ≠ unrestricted n53 search ≠ vs_rho ≠ ledger promotion.

## Next ticks

1. Attach relative-support witnesses / selector-native lazy roots (not only
   exceptional nogoods) before four x-factors are fixed; re-test n19.
2. Harvest factor_base n53 run `20260915T122253Z-229f0ee245` when measurement ends.
3. Independent validation of prior PASS drafts before any ledger promotion.
4. Re-checkout IC branch if worker flipped to rho; never mutate AWS.

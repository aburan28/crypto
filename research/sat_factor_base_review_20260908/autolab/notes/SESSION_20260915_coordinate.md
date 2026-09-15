# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15 ~11:15Z tick)

- Bucket `s3://ecc2k130-590183823895`
- Slots 8, alive **6**, retired 2, errors 0
- Aggregate **87.073 B it/s**, checkpointed 7.159e15 iters ≈ 0.333% of 2^60.9
- ETA ~285d at current rate
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111` (keep IC work off the rho PR branch;
sibling may flip the shared worktree — re-checkout this branch each tick).

## Runs this session

| Beat / probe | Path / run_id | Result |
|---|---|---|
| preflight | — | ok |
| smoke `vs_rho` n13 | `runs/20260915T110615Z-cbfe798888` | PASS / PENDING_IV |
| `vs_rho` n37_wall 1fx | `runs/20260915T110708Z-6bac5fa80e` | PASS draft; IC 638 vs ρ 194 wall (no 20% win) |
| relative Frobenius pair support n13/19/23/37/41/53 | `runs_manual/relative_pair_stats_20260915/` | compression = n every rung; 0 equivariance misses |
| j0 16-bit e2e | `runs_manual/prime_j0_e2e_16bit_20260915/result.json` | ic_agrees_rho ✓, ic_matches_truth ✓ |
| orbit_factorized planted ladder n13→53 | `runs_manual/orbit_factorized_s5_20260915/extraction_panel.json` | **edge-free SAT** with planted x+chain units; n53: 1 valid tuple, 0 bad lifts, 0 pair table, 1053 conflicts / 11.1s solve |
| orbit_factorized n13 unrestricted 1M | same dir | SAT with **8** valid tuples, 0 bad lifts, then UNKNOWN (censored) |

## Priority #1 signal

1. `relative_pair_stats`: n-fold support compression through n=53 (accounting only).
2. **Edge-free planted extraction**: `factorized_s5_over_frobenius_orbit_factor_base`
   recovers group-valid planted relations at n=13,19,23,37,41,**53** with
   `pair_table_entries=0` and `pair_selector_variables=0`.
3. Unrestricted n13 (no planted units) already yields multiple group-valid
   models under a 1M-conflict budget (censored UNKNOWN after 8 models).

Claim boundary: planted units certify formula soundness / planted recovery,
not unrestricted n53 search, SAT advantage, or vs_rho crossover. No ledger
promotion.

## Next ticks

1. Unrestricted (no planted units) budgets at n19–n37; then n53 if stable.
2. Feed relative-Frobenius pair-support into branching / lazy roots so
   unrestricted search prunes before four factors are fixed.
3. Independent validation of prior PASS drafts before any ledger promotion.
4. Re-checkout `cursor/ic-boundary-experiments-d111` if the shared worker
   is back on the rho branch; never mutate AWS.

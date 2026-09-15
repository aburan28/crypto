# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| [Rho challenge spot instances](https://cursor.com/agents/bc-01a0a48f-2eba-7a64-8666-ddc1a66e3a89) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| [Experiment coordination](https://cursor.com/agents/bc-01a0a4b9-9938-7c24-8a14-465cbab7d111) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only)

- Bucket `s3://ecc2k130-590183823895`
- ASG `ecc2k130-workers` desired 8 / in-service 1 (Oregon spot); plus 3 us-east-1 on-demand
- Aggregate ~57.8 B it/s, ~0.33% of expected 2^60.9 work, 4 alive slots
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111` (IC dirty tree moved off the rho PR branch).

## Runs this session

| Beat / probe | Path / run_id | Result |
|---|---|---|
| preflight | — | ok |
| smoke `vs_rho` n13 | `runs/20260915T110615Z-cbfe798888` | PASS / PENDING_IV |
| `vs_rho` n37_wall 1fx | `runs/20260915T110708Z-6bac5fa80e` | PASS draft; IC 638 vs ρ 194 wall (no 20% win) |
| relative Frobenius pair support n13/19/23/37/41/53 | `runs_manual/relative_pair_stats_20260915/` | compression = n at every rung; n53: 150 reps, 53×, 47.7 MiB vs 1.88 GiB expanded LB, 0 discrepancies |
| j0 16-bit e2e | `runs_manual/prime_j0_e2e_16bit_20260915/result.json` | ic_agrees_rho ✓, ic_matches_truth ✓ (~2.1 s IC / 6.7 ms ρ) |
| orbit_factorized S5 n13 planted | `runs_manual/orbit_factorized_s5_20260915/` | 0 pair-table / 0 edge selectors; 100k conflicts → UNKNOWN |

## Priority #1 signal

`relative_pair_stats` measures reusable pair-level support over representative +
relative-Frobenius shift states without materializing the expanded endpoint
table. Observed state / root-entry compression ratio equals `n` on every
successful rung; equivariance sample checks report 0 discrepancies.

This is **support accounting only** — not relation extraction, not a SAT win,
not a vs_rho crossover.

## Next ticks

1. Wire orbit-factorized S5 + lazy pair roots toward a group-valid n=53
   extraction without the pair/edge table; keep memory + formula accounting.
2. Independent validation of prior PASS drafts before any ledger promotion.
3. Keep polling rho campaign read-only; escalate only if alive slots drop to 0.

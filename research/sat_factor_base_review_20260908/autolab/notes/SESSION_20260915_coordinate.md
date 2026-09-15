# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T19:05:00Z)

- ASG `ecc2k130-g7`: desired/InService **6/6**
- ASG `ecc2k130-workers`: desired **4**, InService **0**, instances **0**
- Do **not** run `fleet.sh` / `infra.sh` from this agent

## Branch

`cursor/ic-boundary-experiments-d111`

## Runs this session (highlights)

| Beat / probe | Result |
|---|---|
| intermediates_then_pair stack | n19 draft-SAT (multi-seed) |
| **n19 IV @300k** | planted 1–5,7,11 + natural 1–3 all SAT group-valid |
| **n23 seed1 @300k** | **SAT @168284** (replay identical) |
| **n23 @1M hit-rate** | planted seed2 **SAT @590623**; seeds 3/5 UNKNOWN @1M; natural seed1 UNKNOWN @1M; natural seed2 running |
| **Bound lift 23→41** | `pair_sum_trie` + relative-positive relative mode; absolute positive still ≤23 |
| **n37 cert export** | succeeded (`pair_table=0`, edge_selectors=0, 12 orbit reps) |
| **n37 planted-units reload** | restarted clean-env via nohup (prior run polluted by `KIC_DIMACS_PATH`); in flight |
| CMS without branch heuristic | still INDETERMINATE @200k |

## Priority signals

1. Edge-free planted n53 stands under full S3.
2. n19 unrestricted cliff **draft-broken with budget-robust multi-seed IV**.
3. n23 hit-rate incomplete: planted ~2/4 SAT (seed1@300k, seed2@1M); seeds 3/5 + natural1 UNKNOWN @1M.
4. Admitted ladder to n37 unblocked; cert export proves encoding constructs edge-free.
5. Keep full S3 for planted n≥19 / n53 (lazy unsuitable).
6. Skip n37/n53 symmetrised mid-dim (blocked).

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Collect natural seed2 @1M + n37 planted-units clean reload.
2. If n37 planted units SAT, push unrestricted n37 / n41 / n53 under same stack.
3. Draft ≠ ledger promotion ≠ vs_rho ≠ key recovery.
4. Never mutate AWS from this agent.

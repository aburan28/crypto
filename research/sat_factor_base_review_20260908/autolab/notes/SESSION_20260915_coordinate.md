# Experiment coordination — 2026-09-15

Public synthetic / known-answer only. No ledger promotion. No AWS mutations
(rho fleet owned by sibling agent).

## Parallel agents

| Agent | Role | Status |
|---|---|---|
| Rho challenge spot instances (`bc-01a0a48f…`) | G7e ASG / ECC2K-130 DP collection | RUNNING — leave AWS alone |
| Experiment coordination (`bc-01a0a4b9…`) | IC boundary autolab + priority probes | this session |

### Rho campaign snapshot (read-only, 2026-09-15T20:45:00Z)

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
| **n23 @1M hit-rate** | planted seed2 **SAT @590623**; seeds 3/5 UNKNOWN; natural seed1 UNKNOWN; **natural seed2 SAT @846240** (first natural n23 @1M) |
| **Bound lift 23→41** | `pair_sum_trie` + relative-positive relative mode; absolute positive still ≤23 |
| **n37 cert export** | succeeded (`pair_table=0`, edge_selectors=0, 12 orbit reps) |
| **n41 pairthen chain planted** | **SAT @389** group-valid pair_table=0 |
| **unrestricted n37 eta1/16** | UNKNOWN @50k; 200k in flight; eta1/2@100k aborted |
| **bound lift →53** | pair_sum/relative-positive relative mode 41→53 |
| **n37 chain units breakthrough** | Missing `KIC_PLANTED_CHAIN_ROOT_UNITS=1`. With it: n37 pairthen SAT@93; fullstack SAT@92 group-valid pair_table=0; n53 control SAT@3631. Unrestricted@100k in flight |
| **n37 planted-units diagnostics** | n23 fullstack control SAT@115; n37 lean UNKNOWN@100k; eta1/2 UNKNOWN@1k; eta1/16 UNKNOWN@5k (50k in flight); pair_table=0 |
| **n37 formula accounting** | edge-free: vars≈428k, +clauses≈6.89M, pair_table=0; planted units UNKNOWN@1k (~94s/1k); 100k aborted (~2.6h); no-positive@20k diag in flight |
| CMS without branch heuristic | still INDETERMINATE @200k |

## Priority signals

1. Edge-free planted n53 stands under full S3.
2. n19 unrestricted cliff **draft-broken with budget-robust multi-seed IV**.
3. n23@1M draft hit-rate: planted 2/4 incl. prior seed1; natural 1/2 (seed2 SAT@846240). Seeds 3/5 + natural1 UNKNOWN.
4. Admitted ladder to n37 unblocked; cert export proves encoding constructs edge-free.
5. Keep full S3 for planted n≥19 / n53 (lazy unsuitable).
6. Skip n37/n53 symmetrised mid-dim (blocked).

Claim boundary: drafts ≠ ledger promotion ≠ vs_rho ≠ key recovery.

## Next ticks

1. Collect unrestricted n37 @100k; if SAT push natural n37 + n41/n53 under winning stack with chain units for planted sanity.
2. If n37 planted units SAT, push unrestricted n37 / n41 / n53 under same stack.
3. Draft ≠ ledger promotion ≠ vs_rho ≠ key recovery.
4. Never mutate AWS from this agent.

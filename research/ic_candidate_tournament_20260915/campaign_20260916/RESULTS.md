# Continued IC tournament results

Batch (16-target) parity verdict: **True**; selected implementation **combined_descent**.
Single-target parity verdict: **True**; selected implementation **tiny_batch1** (round-0006).

Completed four new tournaments with **6,288 profiled trials**, each paired with a fresh native run. Every listed round passed its independent artifact/correctness audit. The 36-trial batch screen is separate development evidence.

## Profiled instruction cost relative to matched rho

| Round | Targets per cold job | Selected candidate | Instructions / rho | Paired 95% interval |
|---|---:|---|---:|---|
| [round-0002](../runs/round-0002/REPORT.md) | 1 | batch16 | 13.0349 | not calculated in original protocol |
| [round-0003b](../runs/round-0003b/REPORT.md) | 1 | combined_batch8 | 7.1494 | [6.563380843095217, 7.844932268218591] |
| [round-0004](../runs/round-0004/REPORT.md) | 1 | folded_lift_batch4 | 3.0429 | [2.688274518239304, 3.4833783647230367] |
| [round-0005-batch16](../runs/round-0005-batch16/REPORT.md) | 16 | combined_descent | 0.7195 | [0.6374471292425982, 0.8219836239443673] |
| [round-0006](../runs/round-0006/REPORT.md) | 1 | tiny_batch1 | 0.8264 | [0.8032467836918035, 0.8488166193357406] |

## Native process time relative to matched rho

| Round | Targets per cold job | Selected candidate | Time / rho | Paired 95% interval |
|---|---:|---|---:|---|
| round-0002 | 1 | batch16 | diagnostic only | old polling-based timing excluded |
| round-0003b | 1 | combined_batch8 | 2.7579 | [2.408587066800431, 3.0866033786867497] |
| round-0004 | 1 | folded_lift_batch4 | 1.5213 | [1.4825577742485183, 1.5518235985269213] |
| round-0005-batch16 | 16 | combined_descent | 0.7796 | [0.6991620583548833, 0.8736653449888768] |
| round-0006 | 1 | tiny_batch1 | 0.9766 | [0.957057247571661, 0.9946203513598987] |

## Single-target parity (round-0006)

[Pre-registered](ROUND6.md) before measurement; seed 2026091606; 1,584 profiled trials with paired native runs, all verified; incumbent round-0004's `folded_lift_batch4`. The selected `tiny_batch1` runs the same bounded algorithm — the identical factor-base point set, three summands, a folded pair-sum table, a full-rank solve with every column logarithm certified in the group, a `[a]G + [b]Q` descent and the unchanged general-arithmetic final check — in single-word arithmetic: a Euclidean field inverse, a normal-basis orbit key, density-sized table rows, walked probes, no thread pool, and a directly written report ([patch](round6-tiny.patch), [winner record](WINNER-single-target.json)).

Confirmation and replay, one target per cold job, five cells:

| Comparison | Instructions | Paired 95% interval | Native time | Paired 95% interval |
|---|---:|---|---:|---|
| winner / incumbent (confirmation) | 0.2702 | [0.2320, 0.3073] | 0.6511 | [0.6236, 0.6801] |
| winner / rho (confirmation) | 0.8264 | [0.8032, 0.8488] | 0.9766 | [0.9571, 0.9946] |
| winner / rho (replay) | 0.8264 | [0.8032, 0.8488] | 0.9854 | [0.9641, 1.0062] |

Per-cell winner/rho instruction ratios 0.788–0.850 and native ratios 0.947–1.012; all within the 1.10 parity margin on both stages. The winner uses fewer instructions than rho with the whole interval below one. In native process time the point estimates are below one but the replay interval reaches 1.006, so the native result is parity, not a demonstrated speed advantage. The controls in the same round: `combined_descent` (the batch-panel winner's source, 0.949 of the incumbent), `fast_report` (serialisation only, 0.922) and `tiny_fulltable` (the table rule disabled, 0.2785 of the incumbent against 0.2711 for `tiny`; the rule's gain is on n13a0, where it builds two of seven rows).

Where the ratio comes from: on n13a0 rho's job is about 6.2 million instructions, of which curve construction, target hashing, general-arithmetic final verification and startup — shared by both arms — are about 4.3 million and the rho solve about 1.9 million, most of it the jump-table scalar multiplications and the general-arithmetic collision finish rather than the walk. The IC arm's own work has to fit in that solve plus ten percent of the job; on that cell it now takes about 0.7 million. The shipped rho's per-job cost on these cells is therefore mostly fixed setup; a rho specialised the same way has not been measured and would be cheaper. This is an implementation result in fixed-compiler Valgrind amd64 instructions and matched native time on five small Koblitz cells, class engineering; no arithmetic-complexity, family-wide or cryptographic-size claim follows.

## Interpretation

Every ratio uses a fresh matched rho run in the same round. The 16-target panel charges all setup once to the complete job and solves every target; it is separate from the single-target result, and no ratio combines the two panels. Rho uses the existing per-target solver API on the same constructed curve. Additional cross-target rho optimizations, and a rho specialised like the round-0006 winner, have not been measured here.

The parity rule was declared before measurement: both candidate/rho upper paired 95% limits and every curve-cell ratio must be at most 1.10, in instructions and native time, on confirmation and replay. The full decisions contain the replay evidence.

These are implementation improvements in fixed-compiler Valgrind amd64 guest instructions and matched native wall time. Kernel/device and external-audit work are outside the instruction count. No arithmetic-complexity, broader-family, or cryptographic-size claim follows. The K-instruction rank floor is deliberately weak.

The successful mechanisms are exact arithmetic substitutions, folded pair-table construction, smaller relation batches, cheaper descent initialization/checking and, in round-0006, a single-word implementation of the whole pipeline. Every returned scalar still passes general worker verification and the independent Python checker.

## Reproduce and review

- [16-target winner source/configuration](WINNER.json) and [cumulative source patch](WINNER.patch).
- [Single-target winner record](WINNER-single-target.json), its [patch against the round-0005 winner source](round6-tiny.patch), the [serialisation-only](round6-fast_report.patch) and [full-table](round6-tiny_fulltable.patch) controls, and [round6_candidates.py](round6_candidates.py) to rebuild the trees.
- [Operation plan](PLAN.md), [single-target successor plan](ROUND4.md), [batch plan](ROUND5.md), and [single-target parity pre-registration](ROUND6.md).
- [Operating guide](../OPERATIONS.md) and [controller tests](controller-tests.json).
- [Arithmetic/table equivalence tests](preflight-next-tests.log) and [installed skill validation](skill-validation.json).
- [Retained failed build attempt](../runs/round-0003/prepare_failure.json); no measurements came from it.

Production library defaults were not changed. The frozen source, worker, configurations, raw profiles, receipts and replay evidence are retained under each linked round.

# Disclosed-point F4/F5 dimension-6 pilot — same-budget replay

Status: **replay of the canonical series, not a new registration**.
`promotion_eligible=false`. The merged record is
[../RESULT.md](../RESULT.md) (PR #961). This directory repeats
`max_trials=1`, `node_budget=4096`, and the five inventory points. Its build
SHA-256 matches that series. Its worker binary hash does not, because the
worker was rebuilt from the same source commit. Do not rerun `run_pilot.py`.
Comparison kind: `factor-base-policy`. Competitive `S`, rho ratios, and
scoreboard rows stay unset. Valgrind 3.22.0 is absent on this host, so
`instruction_count` is null on every receipt. `wall_ns` is process
practicality only.

## What was run

| Item | Binding |
| --- | --- |
| Worker commit | `765c3c5f19032bd852163805f257c56babef2040` |
| rustc | `rustc 1.94.1 (e408947bf 2026-03-25)` |
| Host | Linux x86_64, Intel Xeon, 4 logical CPUs, `RAYON_NUM_THREADS=1` |
| Worker SHA-256 | `46325885cc72e9c4928f180f6a687e85d198887ef93da4488db471920d7ef62b` |
| Build SHA-256 | `de3cb8b896f31f03668f1d0eb14302fef2b1e0bc0be303e7f4021a70dd085335` |
| Inventory | `standard-subspace-d6-inventory-control.json` SHA-256 `be3053b4b637311e8255f294807fde412510fbda2464e678a6a16fdea397d403` |
| Encoder panel | `panel.json` SHA-256 `26b7bcc0c82f16e8d7a1b1ad7ebb3cb036b0eb926211d8b16d60b47c844f32c6` |
| Static preflight | `PASS_STATIC_LAYOUT_ONLY` (0 cells above `MAX_VARS=64`) |
| Budget | `f4` and `f5`, `max_trials=1`, `batch_trials=1`, `groebner_degree=3`, `node_budget=4096`, dense relation LA, 180s wall cap, 8 GiB address-space cap |
| Points | The five public targets already stored in the inventory control |
| Summary SHA-256 | `04230b42ce3a6392b9a2681ac7ab385e43db41aba5e2617466597246a65a2e4c` |

`generic_solver_feasibility.py --require-pass` ran against the detached
worker checkout before any job. No new target was sampled.

## Table

One unit: verified complete one-target IC recovery under the frozen budget.
Every row is incomplete, so the speedup column is unknown.

| Cell | Solver | Engine | PDP outcome | `unsupported` | Node budget exhausted | Accepted relations | Folded columns | Geometric points | Wall ms (practicality) |
| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| n17a1 | f4 | MatrixF4 degree 3 | incomplete | false | true | 0 | 29 | 63 | 6152.20 |
| n17a1 | f5 | MatrixF5 degree 3 | incomplete | false | true | 0 | 29 | 63 | 4349.90 |
| n19a0 | f4 | MatrixF4 degree 3 | incomplete | false | true | 0 | 27 | 65 | 11422.81 |
| n19a0 | f5 | MatrixF5 degree 3 | incomplete | false | true | 0 | 27 | 65 | 7303.17 |
| n23a0 | f4 | MatrixF4 degree 3 | incomplete | false | true | 0 | 33 | 75 | 19897.29 |
| n23a0 | f5 | MatrixF5 degree 3 | incomplete | false | true | 0 | 33 | 75 | 11777.67 |
| n23a1 | f4 | MatrixF4 degree 3 | incomplete | false | true | 0 | 23 | 53 | 20089.48 |
| n23a1 | f5 | MatrixF5 degree 3 | incomplete | false | true | 0 | 23 | 53 | 11669.86 |
| n31a0 | f4 | MatrixF4 degree 3 | incomplete | false | true | 0 | 27 | 69 | 126161.83 |
| n31a0 | f5 | MatrixF5 degree 3 | incomplete | false | true | 0 | 27 | 69 | 35753.46 |

Each job exited 2 with report `status=incomplete`. Each ordinary query recorded
exactly 4096 reductions, which is the frozen node budget. Collector dispatch
was `Groebner` / `trial-keyed-sample` with no pair table. Final relation LA
attempts were 0. `relation_la`, `matrix_build`, target descent, and recovery
phase times are null, so the eleven-phase scientific total and any `S` ratio
stay **unknown**. Folded column counts match the inventory control. Usable
subgroup `B` remains that control's census; this solve report does not replace it.

Class: **engineering / diagnostic**. The encoder floor `m·ℓ+(m−2)·n ≤ 64` at
`ℓ=6` was already satisfied before these jobs. Nothing here is an advance, a
family qualification, or an ECC2K-130 transfer.

## Decision

1. Dimension-6 `standard_subspace` clears the Semaev encoder cap. Both `f4` and `f5` enter the named matrix engines on all five disclosed cells.
2. Under `max_trials=1` and `node_budget=4096`, no cell produced an accepted relation or a verified one-target recovery.
3. Do not rerun this budget. Do not redispatch seeds `2026092901` or `2026092902`. Do not open a fresh competitive F4/F5 panel from these reused points.

Raw jobs, stdout, stderr, and `report.json` files are under `runs/`. The
controlled build receipt is under `build-receipt/` because this tournament
ignores `**/build/` and `**/worker`. The local executable and source tar
remain on the measurement host; `build-receipt/local-artifacts.json` records
their SHA-256 and byte counts.

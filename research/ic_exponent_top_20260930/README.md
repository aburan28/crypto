# Ledger §23: the top-end exponent, in range

> **Withdrawn before anything ran.** AGENTS.md's IC measurement rules
> (#1076) made one previously unseen target the primary comparison, and
> bar a ratio to batch rho at `k = 32` as a headline. This declaration's
> primary figure was that ratio. Nothing below ran. §23 was declared
> again as a single-target round in
> [`../ic_single_target_20260930/`](../ic_single_target_20260930/PROTOCOL.md).
> This file is kept unchanged below as the record.

A measurement round (AGENTS.md §3, accounting) on the Koblitz collection
thread of [ledger §20](../ic_exponent_20260926/)–[§22](../ic_descent_20260930/).
It was declared in [`PROTOCOL.md`](PROTOCOL.md) before anything ran, and is
written up in `research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md`
§23.

**Sizes.** Every Koblitz curve the library can build (`n ≤ MAX_N = 63`)
with `r ≥ 2^36`. That is the ladder's four largest, plus
`icv1-f2m57-tm747311035-c1f545af` (`2^38.0`) and `icv1-f2m59-tm943548413-98844ecc` (`2^44.5`).

**What is measured.** Each size's recipe is re-swept on this round's
binary, then priced on eight sets, with every timed process isolated.
The exponent is fitted over the six sizes.

## Files

| file | what it is |
|:--|:--|
| `PROTOCOL.md` | the declaration: sizes, recipe, sets, accounting, isolation, targets, stop rules |
| `predict.py`, `prediction.json` | §20's frozen model and law, carried to the six sizes; each size's sweep grid |
| `run.py` | the runs, in the declared order; resumable, never overwrites (`IC=… IC_RUNS=…`) |
| `analyse.py` | every figure the note and the page quote |
| `curve_records.json`, `curve_ids.py`, `curve_ids.json` | the six curves' exact field and curve records, and their EC1 aliases and full UIDs |

The runs, the analysis and the host manifest are added with the results.

## Reproducing

    cargo build --release --bin ic                    # at the declaration commit; kept outside the tree
    cd research/ic_exponent_top_20260930
    IC=<ic> IC_COMMIT=<commit> IC_RUNS=runs python3 run.py all
    IC_RUNS=runs python3 analyse.py > analysis.json

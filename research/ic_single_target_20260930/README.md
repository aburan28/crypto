# Ledger §23: one unseen point, online and cold

The Koblitz collection thread's first single-target round, under
AGENTS.md's IC measurement rules (#1076). It was declared in
[`PROTOCOL.md`](PROTOCOL.md) before anything ran. It is written up in
`research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md` §23.

§23 was first declared at `k = 32`, in
[`../ic_exponent_top_20260930/`](../ic_exponent_top_20260930/). That
declaration was withdrawn before anything ran, when the rule made one
target the primary comparison.

**What is measured.** On each of the six Koblitz curves the library can
build with `r ≥ 2^36`, 64 public hash-to-curve points are measured, one
process each:
- the index calculus's reusable setup;
- its online interval, split into five exclusive phases;
- one-target rho on the same point;
- both answers, replayed outside the process.

Every row is keyed by the canonical identities and checked by the
repository's claim checker.

## Files

| file | what it is |
|:--|:--|
| `PROTOCOL.md` | the declaration: question, arms, boundaries, sizes, recipe, rows, accounting, isolation, targets, stop rules |
| `predict.py`, `prediction.json` | the declared prediction, from §22's frozen records and formulas |
| `make_params.py` | one parameter file per target: §20's recipe, one public target |
| `run.py` | the runs, in the declared order, each process isolated; resumable, never overwrites (`IC=… IC_RUNS=…`) |
| `claims.py` | identities (`identity.py`), independent replays (`oracle.py`), a `vs_rho` claim per row and its check |
| `analyse.py` | every figure the ledger and the page quote |
| `curve_records.json`, `curve_ids.json` | the six curves' exact records, EC1 aliases and UIDs (as in the withdrawn declaration) |

The runs, claims, manifests, analysis and host manifest are added with
the results.

## Reproducing

    cargo build --release --bin ic      # at the declaration commit; kept outside the tree
    cd research/ic_single_target_20260930
    IC=<ic> IC_COMMIT=<commit> IC_RUNS=runs python3 run.py all
    IC_RUNS=runs python3 claims.py
    IC_RUNS=runs python3 analyse.py > analysis.json

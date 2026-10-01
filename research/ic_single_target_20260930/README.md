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
| `render_rows.py` | the ledger's Markdown rows and the page's HTML rows, printed from `analysis.json` |
| `curve_records.json`, `curve_ids.json` | the six curves' exact records, EC1 aliases and UIDs (as in the withdrawn declaration) |
| `analysis.json` | `analyse.py`'s output: per-size figures with bootstrap intervals, fits, A/A, diagnostics, prediction against measurement |
| `runs/host.json` | the host manifest: commit, `rustc`, CPU and features, cores, memory, binary SHA-256, isolation tool hash |
| `runs/pin/` | the pin: the batch pricer on §20's M1 files, counts against §22's record |
| `runs/k<a>n<n>/` | per target: the parameter file, each process's report (`*.price.json`), isolation record and stderr; `*-retry1.*` is the rerun of a contended attempt, which is kept |
| `runs/run.log`, `runs/refusals.log`, `runs/psi-waits.log` | the run's log, every refused start, and the PSI wait before each process |
| `claims/k<a>n<n>/` | per row: the `vs_rho` claim and its independent replay; `claims/summary.json` counts the checker's verdicts |
| `manifests/{candidates,workloads,rho}/` | the canonical candidate, workload and rho manifests the claims cite |

## Results

All 410 processes completed, two were contended and were run again (the
contended attempts are kept), and none failed. All 408 claims pass
`validate_claim`; both arms' scalars replay in `oracle.py` and agree.
The pin reproduced §22's counts.

| curve | log₂ r | online speedup [95%] | canonical step | cold ratio | S, set-up | break-even | over BL model |
|:--|--:|--:|--:|--:|--:|--:|--:|
| `K_1/GF(2^47)` | 36.6 | 9.77× [7.83, 12.5] | 4.34× | 4.95× | 4.20 | 7 | 2.4× |
| `K_0/GF(2^57)` | 38.0 | 13.5× [10.7, 17.3] | 5.97× | 1.98× | 8.08 | 15 | 3.3× |
| `K_0/GF(2^41)` | 39.0 | 9.38× [7.43, 12.1] | 4.57× | 8.04× | 5.51 | 10 | 2.5× |
| `K_0/GF(2^53)` | 44.3 | 17.5× [13.6, 22.6] | 9.37× | 14.7× | 7.37 | 16 | 2.0× |
| `K_1/GF(2^59)` | 44.5 | 16.9× [13.0, 22.1] | 8.90× | 15.7× | 9.00 | 18 | 2.8× |
| `K_0/GF(2^61)` | 47.2 | 8.84× [6.76, 11.9] | 4.68× | 31.6× | 15.10 | 37 | 8.0× |

Online, the index calculus is faster than rho on the same point at all six
sizes. Cold, with its reusable set-up, rho is faster at all six, so by
AGENTS.md §2 the method is not faster than rho. The online gain sits
2.0–8.0× above a generic walk given the same set-up as precomputation
(Bernstein–Lange; a model, no such walk ran). At `2^38.0` both arms pay a
48–58 ms curve construction, which is why its cold ratio is low; without
it the six sizes read 6.3–32.5×. The ledger §23 has the readings, the
graded targets and the class (accounting: no algorithm changed).

## Reproducing

    cargo build --release --bin ic      # at the declaration commit; kept outside the tree
    cd research/ic_single_target_20260930
    IC=<ic> IC_COMMIT=<commit> IC_RUNS=runs python3 run.py all
    IC_RUNS=runs python3 claims.py
    IC_RUNS=runs python3 analyse.py > analysis.json

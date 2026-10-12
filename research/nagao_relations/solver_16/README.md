# solver_16 — the Riemann–Roch and S′4 encodings under Weil descent, in front of XOR-native SAT, on ECC2K-130

Successor to [`RESEARCH_ECC2K130_WDSAT.md`](../../notes/ecc2k130/RESEARCH_ECC2K130_WDSAT.md)
(#498, toy fields, circuit-IR export) and to `solver_10`–`solver_15`
([`RESEARCH_ECC2K130_RR_SOLVER_PANEL.md`](../../notes/ecc2k130/RESEARCH_ECC2K130_RR_SOLVER_PANEL.md),
the algebraic RR solvers on the real curve). The write-up is
[`RESULTS.md`](RESULTS.md); the frozen scope is [`contract.json`](contract.json),
committed with its prediction before execution.

| file | role |
|:--|:--|
| `weil.py` | symbolic Weil descent over F₂ in a polynomial basis, `V = {deg x < d}`; Trimoska's S′4 layout and the Nagao norm form with `b` eliminated; WDSat ANF writer and static-config sizing; CryptoMiniSat CNF-XOR adapter with deterministic conflict bracketing; independent witness verifier |
| `test_weil.py` | 7 tests: exhaustive tiny-field oracles for both encodings against the group law, RR ⊆ S4 zero set, reproduction of Trimoska's published `Xn15l5-1-S.anf` equation for equation, WDSat FIND_ALL enumeration against the oracle |
| `contract.json` | frozen hypothesis, prediction, boundaries, panel, gates, falsifier |
| `run.py` | campaign driver; appends to `raw.jsonl`, resumable |
| `targets.json` | every minted target, seed 20260920161 |
| `raw.jsonl` | one record per solver cell, append-only, with every witness and its verification |
| `summarize.py` → `summary.json` | aggregates, fits, matched ratios; cites, never re-solves |
| `pilot.py`, `pilot.json`, `pilot_calibration.log` | sizing pilots — **development data**, other seeds, not evidence |

Run from anywhere with Python 3.10+, gcc, `pip install pycryptosat==5.15.0`, and the
pinned WDSat source fetched by
`research/index_calculus_baseline_20260914/pilot/build_pilot.py --cache-dir /tmp/ec-baseline-cache`:

```sh
python3 research/nagao_relations/solver_16/test_weil.py
python3 research/nagao_relations/solver_16/run.py        # ~3 h on 4 cores; resumable
python3 research/nagao_relations/solver_16/summarize.py
```

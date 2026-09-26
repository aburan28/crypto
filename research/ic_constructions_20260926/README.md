# Ledger §21: the Koblitz workflow's constructions, built once

An engineering round (AGENTS.md §3) on the Koblitz collection thread of
[ledger §20](../ic_exponent_20260926/). Declared in
[`PROTOCOL.md`](PROTOCOL.md) before any candidate code existed; written
up in `research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md`
§21, and drawn on `docs/index-calculus-scoreboard.html`.

**The change** (candidate `33936731`): a `FrobeniusFactorBase` builds its
projected signed-orbit map, point index and decomposition-check classes
once and shares them with every consumer, and the single-word map
projects one point per Frobenius orbit instead of one per point.

## Files

| file | what it is |
|:--|:--|
| `PROTOCOL.md` | the declaration: change, inputs, procedure, measures, targets, stop rules |
| `inputs.sha256` | the 36 frozen §20 parameter files this round prices, by hash |
| `run.py` | the runs, in the declared order; resumable, never overwrites |
| `analyse.py` | every figure the note and the page quote → `analysis.json` |
| `render_rows.py` | the table rows, rendered from `analysis.json` (`html` or `md`) |
| `host.json` | host manifest and both binaries' sha256 |
| `runs/main/k{a}n{n}/M{j}/r{i}-{baseline,candidate}.price.json` | the 360 main-comparison reports, with stderr |
| `runs/control1/` | `ic workflow` against `ic price` on the candidate, `M1` at every size |
| `runs/rho/` | batch rho re-priced on the candidate at `n = 41`, `M1`, against §20's counts |
| `runs/threads/` | four Rayon threads, `n = 53` and `n = 61`, `M1`, three ABAB rounds |
| `runs/constructions/` | per-constructor prices on both binaries (`examples/koblitz_construction_prices.rs`) |
| `runs/run.log` | the run's own log |
| `unit_shift/` | the unit's code in both binaries, normalised (`disasm.sh` → `*.s`, `disasm.txt`): ledger §21.4 |

## Reproducing

The two binaries are built from their commits with the rustc in
`host.json` and kept outside the tree:

    git checkout e7022b75 && cargo build --release --bin ic                      # baseline ic
    git checkout 305d5078 && cargo build --release --example koblitz_construction_prices
    git checkout 33936731 && cargo build --release --bin ic --example koblitz_construction_prices

(`305d5078` only adds the pricer example and the declaration to
`e7022b75`; its library is `e7022b75`'s.)  Then, with each binary's path
in the environment:

    cd research/ic_constructions_20260926
    IC_BASELINE=… IC_CANDIDATE=… PRICES_BASELINE=… PRICES_CANDIDATE=… python3 run.py all
    python3 analyse.py > analysis.json
    python3 render_rows.py md

Every single-thread process runs with `RAYON_NUM_THREADS=1` under
`taskset -c 2`, one at a time; the thread check uses four threads under
`taskset -c 0-3`.

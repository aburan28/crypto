# Ledger §22: the Koblitz descent's fixed cost per target

An engineering round (AGENTS.md §3) on the Koblitz collection thread of
[ledger §20](../ic_exponent_20260926/) and [§21](../ic_constructions_20260926/).
Declared in [`PROTOCOL.md`](PROTOCOL.md) before any candidate code existed;
written up in `research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md`
§22.

## Files

| file | what it is |
|:--|:--|
| `PROTOCOL.md` | the declaration: change, inputs, procedure, measures, targets, stop rules |
| `inputs.sha256` | the 36 frozen §20 parameter files this round prices, by hash (the same as §21's) |
| `probe/k{a}n{n}-M1.baseline.json` | the probe on the unmodified library: each target's descent, part by part |
| `probe/k{a}n{n}-M1.logs.json` | the log table the probe used, from `ic workflow` on the same file (baseline binary) |
| `run.py` | the runs, in the declared order; resumable, never overwrites |
| `analyse.py` | every figure the note and the page quote → `analysis.json` |
| `render_rows.py` | the table rows, rendered from `analysis.json` (`html` or `md`) |
| `host.json` | host manifest and all four binaries' sha256 |
| `runs/main/k{a}n{n}/M{j}/r{i}-{baseline,candidate}.price.json` | the 360 main-comparison reports, with stderr |
| `runs/control1/` | `ic workflow` against `ic price` on the candidate, `M1` at every size |
| `runs/rho/` | batch rho re-priced on the baseline at `n = 41`, `M1`, against §20's counts |
| `runs/threads/k0n41-M1/` | four Rayon threads, three ABAB rounds |
| `runs/probe/` | the descent probe on both libraries, six sizes, same log tables |
| `runs/run.log` | the run's own log |

## The probe

`examples/koblitz_descent_prices.rs` rebuilds a parameter file's curve,
base and pair table, loads the log table, regenerates the 32 targets and
times each descent. It asserts every target is recovered, and records
its trials per target, which must equal the pricer's for the same file.

    ic workflow --params P --dir D                     # the baseline ic
    RAYON_NUM_THREADS=1 taskset -c 2 \
      target/release/examples/koblitz_descent_prices P D/logs.json folded 5

It was built at the declaration commit (the library is `main` at
`57e7ce3a`) and has sha256 `75e59a86…`; the baseline `ic` has sha256
`07e527a4…`. Both are kept outside the tree.

## Reproducing the comparison

    git checkout 57e7ce3a && cargo build --release --bin ic                          # baseline ic
    git checkout f231d7d1 && cargo build --release --example koblitz_descent_prices  # baseline probe
    git checkout 20ff4765 && cargo build --release --bin ic --example koblitz_descent_prices
    cd research/ic_descent_20260930
    IC_BASELINE=… IC_CANDIDATE=… PROBE_BASELINE=… PROBE_CANDIDATE=… python3 run.py all
    python3 analyse.py > analysis.json
    python3 render_rows.py md

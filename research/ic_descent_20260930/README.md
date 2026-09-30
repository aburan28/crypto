# Ledger §22: the Koblitz descent's fixed cost per target

An engineering round (AGENTS.md §3) on the Koblitz collection thread of
[ledger §20](../ic_exponent_20260926/) and [§21](../ic_constructions_20260926/).
Declared in [`PROTOCOL.md`](PROTOCOL.md) before any candidate code existed;
written up in `research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md`
§22.

**v1 and v2.**
- **v1** (`runs/`, `analysis.json`) was pinned with `taskset` alone. AGENTS.md
  §10 does not accept that for a timed number.
- **v2** (`runs-isolated/`, `analysis-isolated.json`) ran every timed step
  again through `tools/isolated_bench.py`, under PROTOCOL.md v2, which was
  declared before the rerun. v2 is the evidence.
- v1 is kept as it ran and is never pooled with v2. Ledger §22.7 sets the
  two side by side.

## Files

| file | what it is |
|:--|:--|
| `PROTOCOL.md` | the declaration: change, inputs, procedure, measures, targets, stop rules |
| `inputs.sha256` | the 36 frozen §20 parameter files this round prices, by hash (the same as §21's) |
| `probe/k{a}n{n}-M1.baseline.json` | the probe on the unmodified library: each target's descent, part by part |
| `probe/k{a}n{n}-M1.logs.json` | the log table the probe used, from `ic workflow` on the same file (baseline binary) |
| `run.py` | the runs, in the declared order; resumable, never overwrites (`IC_ISOLATE=1 IC_RUNS=…` for v2) |
| `curve_records.json` | the nine curves' exact field and curve records, from `examples/koblitz_curve_records.rs` |
| `curve_ids.py`, `curve_ids.json` | their EC1 aliases and full curve UIDs (docs/curve-identities.md) |
| `analyse.py` | every figure the note and the page quote → `analysis.json` (v1) or, with `IC_RUNS=runs-isolated`, `analysis-isolated.json` (v2) |
| `analysis-isolated.json` | v2's figures: clean pairs only, the A/A noise floor, and the isolation's accounting |
| `render_rows.py` | the table rows, rendered from `analysis.json` or `analysis-isolated.json` (`html` or `md`) |
| `host.json` | v1's host manifest and all four binaries' sha256 |
| `runs-isolated/` | v2, laid out like `runs/`, each report beside its `*.isolation.jsonl` record |
| `runs-isolated/host.json` | v2's host manifest: the binaries, the byte-identical copy, the isolation tool's sha256, the SMT siblings |
| `runs-isolated/aa/` | the A/A: the baseline against its copy, `M1`, nine sizes, five rounds |
| `runs-isolated/refusals.log` | every start the tool refused, with its reason, and the retry after 15 s |
| `runs-isolated/run.log` | v2's own log, including its one stop and resume between processes |
| `runs/main/k{a}n{n}/M{j}/r{i}-{baseline,candidate}.price.json` | v1: the 360 main-comparison reports, with stderr |
| `runs/control1/` | Control 1, `ic workflow` against `ic price` on the candidate, `M1` at every size; counts only, so it stands for v2 too |
| `runs/rho/` | v1: batch rho re-priced on the baseline at `n = 41`, `M1`, against §20's counts |
| `runs/threads/k0n41-M1/` | v1: four Rayon threads, three ABAB rounds (v2: `runs-isolated/threads/k0n41-M1-t3/`, three threads) |
| `runs/probe/` | v1: the descent probe on both libraries, six sizes, same log tables |
| `runs/run.log` | v1's own log |

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

The binaries, kept outside the tree:

    git checkout 57e7ce3a && cargo build --release --bin ic                          # baseline ic
    git checkout f231d7d1 && cargo build --release --example koblitz_descent_prices  # baseline probe
    git checkout 20ff4765 && cargo build --release --bin ic --example koblitz_descent_prices
    cp <baseline ic> <baseline ic copy>                                               # v2's A/A arm

v2, the evidence:

    cd research/ic_descent_20260930
    IC_ISOLATE=1 IC_RUNS=runs-isolated IC_BASELINE=… IC_BASELINE_COPY=… IC_CANDIDATE=… \
      PROBE_BASELINE=… PROBE_CANDIDATE=… python3 run.py all
    IC_RUNS=runs-isolated python3 analyse.py > analysis-isolated.json
    python3 render_rows.py md analysis-isolated.json

v1, as it ran:

    cd research/ic_descent_20260930
    IC_BASELINE=… IC_CANDIDATE=… PROBE_BASELINE=… PROBE_CANDIDATE=… python3 run.py all
    python3 analyse.py > analysis.json
    python3 render_rows.py md

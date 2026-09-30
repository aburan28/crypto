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

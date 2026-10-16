# Counted Cantor against the cost model's 10 g² multiplications

2026-10-14. Result: [`REPORT.pdf`](REPORT.pdf). Plan: [`PROTOCOL.md`](PROTOCOL.md), frozen before the
full run (hash in `manifest.json`).

| file | what |
|:--|:--|
| `PROTOCOL.md` | question, claims, instrument, controls, decision rule, what the dry run showed |
| `rows.csv` | 10400 counted results: one row per genus, prime, operation and trial |
| `summary.csv` | per-genus means, the fit with bootstrap interval, local exponents, verdict counts (computed natively, `summary` mode) |
| `manifest.json` | commits, hashes, arguments, timings |
| `plot_figs.py`, `build_report.py` | layout only: figures and `REPORT.typ` from the recorded CSVs |
| `ratio_by_genus.svg`, `verdict_fragility.svg`, `REPORT.typ`, `REPORT.pdf` | graphs and report |

Reproduce, from the repository root, after `cargo build --release --example jacobian_op_count`:

```sh
E=./target/release/examples/jacobian_op_count; D=research/jacobian_op_count_20261014
$E run --out $D/rows.csv --seed 20261014 --trials 200 --genera 2,3,4,5,6,8,10,12,14,16,20,24,32 --primes 1009,65521
$E summary --rows $D/rows.csv --jacobian docs/curves/jacobian-transfers.json --out $D/summary.csv
```

# Leaving a blocked twist tower: 3-isogeny walks and the class ceiling

2026-10-14. Result: [`REPORT.pdf`](REPORT.pdf). Plan: [`PROTOCOL.md`](PROTOCOL.md),
frozen before the catalog and null runs (hash in `manifest.json`).

| file | what |
|:--|:--|
| `PROTOCOL.md` | question, claims, controls, decision rule, runtime estimate |
| `control_out.txt` | positive controls A (root finder, N = 20) and B (walk, N = 18): both PASS |
| `catalog_towers.csv` | one row per frozen blocked twist tower (13): reach and escape by depth |
| `null_towers.csv` | 2000 random blocked starts per (N, l), 12 configurations, seed 20261014 |
| `summary.csv` | counts and Clopper–Pearson intervals, computed natively (`summary` mode) |
| `class_ceiling.csv` | exact class counts for N = 8 … 22 (`classes` mode) |
| `manifest.json` | commits, hashes, arguments, timings |
| `plot_figs.py`, `build_report.py` | layout only: figures and `REPORT.typ` from the recorded CSVs |
| `escape_by_config.svg`, `class_ceiling.svg`, `REPORT.typ`, `REPORT.pdf` | graphs and report |

Reproduce, from the repository root, after `cargo build --release --example twist_walk_blocked_towers`:

```sh
E=./target/release/examples/twist_walk_blocked_towers; D=research/twist_walk_blocked_towers_20261014
$E control > $D/control_out.txt
$E catalog --registry docs/curves/registry.json --out $D/catalog_towers.csv
$E null --registry docs/curves/registry.json --samples 2000 --seed 20261014 --out $D/null_towers.csv
$E summary --catalog $D/catalog_towers.csv --null $D/null_towers.csv --out $D/summary.csv
$E classes --min-n 8 --max-n 22 --out $D/class_ceiling.csv
```

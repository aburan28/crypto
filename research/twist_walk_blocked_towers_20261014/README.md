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
| `pilot_towers.csv`, `pilot_towers_summary.csv` | additive, after the freeze: the twist tally that selected the 13 frozen towers, recomputed natively (`tally` mode); its 13 twist-blocked rows equal `catalog_towers.csv` exactly |
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
$E tally --jacobian docs/curves/jacobian-transfers.json --out $D/pilot_towers.csv   # also writes pilot_towers_summary.csv
```

`tally` was added after `PROTOCOL.md` was frozen. The protocol took the 13
towers from a Python pilot in the experiment-proposals work; `tally` reproduces
that selection natively (179 binary towers on 91 curves: none blocked as
recorded, 113 with Tr(a) = 1, 13 twist-blocked on 5 curves) so the claim no
longer depends on the pilot script. It changes no frozen choice, and the
protocol file is untouched.

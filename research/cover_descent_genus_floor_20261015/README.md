# `cover_descent` rows repriced at the GHS genus floor

2026-10-15. Result: [`REPORT.pdf`](REPORT.pdf). This is a correction to
`docs/curves/jacobian-transfers.json`, made after the odd-characteristic
descent study ([`../odd_char_descent_genus3_20261014/`](../odd_char_descent_genus3_20261014/REPORT.md))
showed that no GHS cover of `E/GF(p^k)` with `F_p(j) = GF(p^k)` has genus below 9,
49 and 41 at `k` = 6, 10 and 14. Until now every `cover_descent` row was priced at
genus `k`. Rows are now priced at the floor, labelled with that genus, and keep their
old class as `class_at_genus_k`.

| file | what |
|:--|:--|
| `REPORT.pdf`, `REPORT.typ` | result table, verdict changes, graph, scope |
| `route_changes.csv` | `cover_descent` rows repriced, by `k`, class before and after |
| `verdict_changes.csv` | curve verdicts, old file to new file |
| `before_after.svg`, `before_after.png` | the graph |
| `tabulate.py`, `plot_fig.py`, `build_report.py` | extraction and layout only; the classes are computed natively by `jacobian_transfers build` |

Headline, on the committed file: of the 846 `cover_descent` rows that were
`formula_advantage`, 599 still are (576 at `k` = 6, 23 at `k` = 14), 41 are now
`formula_marginal` and 206 are `boundary` (all of `k` = 10, most of `k` = 14). Curve
verdicts go from 848 / 50 / 1201 to 601 / 67 / 1431 (`formula_advantage` /
`formula_marginal` / `boundary`). The floor is a lower bound for GHS covers only,
and the rows stay `proposed`: no cover is exhibited.

Reproduce, from the repository root:

```sh
cargo build --release --bin jacobian_transfers
./target/release/jacobian_transfers check                       # the committed file is current
git show 3c5d435ee:docs/curves/jacobian-transfers.json > /tmp/old.json   # the file before this change
D=research/cover_descent_genus_floor_20261015
python3 $D/tabulate.py /tmp/old.json docs/curves/jacobian-transfers.json $D
python3 $D/plot_fig.py $D/route_changes.csv $D/before_after.svg
python3 $D/build_report.py $D/route_changes.csv $D/verdict_changes.csv $D/REPORT.typ
```

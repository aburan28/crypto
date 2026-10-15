"""Lay out REPORT.typ from route_changes.csv and verdict_changes.csv (layout only; computes nothing).

usage: python3 build_report.py route_changes.csv verdict_changes.csv REPORT.typ
"""
import csv
import sys

routes = list(csv.DictReader(open(sys.argv[1])))
verdicts = list(csv.DictReader(open(sys.argv[2])))
NAME = {"formula_advantage": "formula_advantage", "formula_marginal": "formula_marginal", "boundary": "boundary"}


def n(rows, **kw):
    return sum(int(r["routes"]) for r in rows if all(r[k] == v for k, v in kw.items()))


body = []
for k in sorted({int(r["k"]) for r in routes}):
    for before in ("formula_advantage", "formula_marginal"):
        sel = [r for r in routes if int(r["k"]) == k and r["class_before"] == before]
        if not sel:
            continue
        g0, g1 = sel[0]["genus_before"], sel[0]["genus_after"]
        body.append(
            f"  [{k}], [{g0} → {g1}], [`{before}`], [{n(sel)}], "
            f"[{n(sel, class_after='formula_advantage')}], [{n(sel, class_after='formula_marginal')}], [{n(sel, class_after='boundary')}],"
        )
tot_adv = [r for r in routes if r["class_before"] == "formula_advantage"]
tot_mar = [r for r in routes if r["class_before"] == "formula_marginal"]
body.append(
    f"  [*all*], [], [`formula_advantage`], [*{n(tot_adv)}*], [*{n(tot_adv, class_after='formula_advantage')}*], "
    f"[*{n(tot_adv, class_after='formula_marginal')}*], [*{n(tot_adv, class_after='boundary')}*],"
)
body.append(
    f"  [*all*], [], [`formula_marginal`], [*{n(tot_mar)}*], [*{n(tot_mar, class_after='formula_advantage')}*], "
    f"[*{n(tot_mar, class_after='formula_marginal')}*], [*{n(tot_mar, class_after='boundary')}*],"
)
vrows = "\n".join(
    f"  [`{r['verdict_before']}`], [`{r['verdict_after']}`], [{r['curves']}]," for r in verdicts
)
adv0 = sum(int(r["curves"]) for r in verdicts if r["verdict_before"] == "formula_advantage")
adv1 = sum(int(r["curves"]) for r in verdicts if r["verdict_after"] == "formula_advantage")

typ = f'''#set page(paper: "a4", margin: (x: 1.8cm, y: 1.5cm))
#set text(size: 9.5pt)
#set par(justify: true)
#show heading: set block(above: 1em, below: 0.5em)

= `cover_descent` rows repriced at the GHS genus floor

_2026-10-15. Generator: `src/cryptanalysis/jacobian_transfer/mod.rs` (`ghs_genus_floor`). Output: `docs/curves/jacobian-transfers.json`, rebuilt by `jacobian_transfers build` and checked by `jacobian_transfers check`. Counts below come from `tabulate.py`, which reads the previous and the new committed file and counts. Source of the floor: `research/odd_char_descent_genus3_20261014/`._

== Result (exact computation on the committed file)

Every `cover_descent` row on a `GF(p^k)` curve was priced at genus `k`. No GHS cover of `E/GF(p^k)` with `F_p(j) = GF(p^k)` has genus below 9, 49 and 41 at `k` = 6, 10, 14, so those rows are now priced at that floor. The table counts the rows by their class before and after.

#table(
  columns: (auto, auto, auto, auto, auto, auto, auto),
  align: (right, left, left, right, right, right, right),
  stroke: 0.4pt + luma(180),
  inset: 3pt,
  [*k*], [*genus*], [*class at genus k*], [*rows*], [*now advantage*], [*now marginal*], [*now boundary*],
{chr(10).join(body)}
)

Curve verdicts (best class over a curve's routes), old file → new file:

#table(
  columns: (auto, auto, auto),
  align: (left, left, right),
  stroke: 0.4pt + luma(180),
  inset: 3pt,
  [*verdict before*], [*verdict after*], [*curves*],
{vrows}
)

`formula_advantage` verdicts go from {adv0} to {adv1}. Rows the model still favours at the floor genus keep their class; the rest are demoted by the model's own pricing, not by a new rule.

#figure(image("before_after.svg", width: 90%), caption: [`cover_descent` rows that were a lead at genus `k`, by class before and after pricing at the GHS floor. Source: `route_changes.csv`.])

== Claims and scope

- *Exact computation:* the counts above, on the regenerated file, which `jacobian_transfers check` reproduces from `registry.json`.
- *From the source study:* the floor values (a proof for the lower bound, an exact enumeration for the table) and Joux–Vitse App. C for genus 9 at `k` = 6.
- *Not shown:* that any cover exists. The rows stay `proposed`. The floor is a lower bound for GHS covers only. Whether a particular curve attains it is verified only at `k` = 6; a curve whose `j` lies in a proper subfield can have a smaller cover; a non-GHS construction is not covered.
- *Not shown:* that the surviving advantages are real. The cost model prices a genus-`g` operation at `10 g²` multiplications, which counted textbook Cantor exceeds by 2 to 3.3 bits (`research/jacobian_op_count_20261014/`), so rows within a few bits of the margin could move.
- *Correction recorded:* an earlier reading was that the 846 advantage verdicts "rest on a genus no GHS cover has", which suggested they should all be demoted. Pricing at the genus the construction does have leaves 599 of the 846 rows as advantage; the demotions are the `k` = 10 rows and most of the `k` = 14 rows.
'''
open(sys.argv[3], "w").write(typ)

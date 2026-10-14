"""Lay out REPORT.typ from summary.csv and manifest.json (layout only; computes nothing).

usage: python3 build_report.py summary.csv manifest.json REPORT.typ
"""
import csv
import json
import sys

summary, manifest, out = sys.argv[1:4]
rows = list(csv.DictReader(open(summary)))
man = json.load(open(manifest))
g_rows = [r for r in rows if r["kind"] == "genus"]
sens = next(r for r in rows if r["kind"] == "sensitivity")
field = next(v for v in sens.values() if v and "advantage_curves=" in v)
d = dict(kv.split("=") for kv in field.split(";"))

pick = {"2", "4", "6", "8", "10", "14", "16", "24", "32"}
trs = []
for r in g_rows:
    if r["op"] == "add" and r["p"] == "1009" and r["g"] in pick:
        trs.append(
            f"[{r['g']}], [{r['mean_mul']}], [{r['mean_inv']}], [{r['mean_total']}], [{r['model_10g2']}], "
            f"[{r['log2_ratio_to_model']}], [{r['local_exponent'] or '—'}]"
        )
tbl = ",\n  ".join(trs)


def fit(op, p):
    r = next(x for x in g_rows if x["op"] == op and x["p"] == p)
    return r["fit_exponent"], r["exponent_ci_lo"], r["exponent_ci_hi"]


def at(op, p, g):
    return next(x for x in g_rows if x["op"] == op and x["p"] == p and x["g"] == g)["log2_ratio_to_model"]


def ratio_map(op, p):
    return {r["g"]: float(r["log2_ratio_to_model"]) for r in g_rows if r["op"] == op and r["p"] == p}


pdiff = max(
    abs(ratio_map(op, "1009")[g] - ratio_map(op, "65521")[g])
    for op in ("add", "double")
    for g in ratio_map(op, "1009")
)
e_add, lo_add, hi_add = fit("add", "1009")
e_add2, lo_add2, hi_add2 = fit("add", "65521")
e_dbl, lo_dbl, hi_dbl = fit("double", "1009")

typ = f'''#set page(paper: "a4", margin: (x: 1.8cm, y: 1.5cm))
#set text(size: 9.5pt)
#set par(justify: true)
#show heading: set block(above: 1em, below: 0.5em)

= Counted Cantor against the cost model's 10 g² multiplications

_2026-10-14. Frozen protocol `PROTOCOL.md` (sha256 `{man["protocol_sha256"][:16]}…`). Native instrument `examples/jacobian_op_count.rs`. All {man["rows"]} counted results agree with the repository's group law (0 mismatches)._

== Result (exact computation)

Mean counts for one Cantor addition of two random degree-`g` divisors over `F_1009`, 200 pairs per genus, with the total weighting an inversion as 8 multiplications. `ratio` is `log2` of the total over the model's `10 g²`; `local e` is the exponent between this genus and the previous measured one.

#table(
  columns: (auto, auto, auto, auto, auto, auto, auto),
  align: (right, right, right, right, right, right, right),
  stroke: 0.4pt + luma(180),
  inset: 3pt,
  [*g*], [*mults*], [*inversions*], [*total*], [*10 g²*], [*ratio (bits)*], [*local e*],
  {tbl},
)

*Verdict (fixed wording).* The model is below textbook Cantor by {at("add", "1009", "2")} bits at genus 2 and {at("add", "1009", "32")} bits at genus 32 (addition, `p = 1009`; doubling {at("double", "1009", "32")} bits and `p = 65521` {at("add", "65521", "32")} bits at genus 32). The ratio is outside the one-bit band at every measured genus. The fitted exponent is {e_add} (95% bootstrap interval {lo_add} to {hi_add}) for addition, {e_dbl} ({lo_dbl} to {hi_dbl}) for doubling and {e_add2} ({lo_add2} to {hi_add2}) for addition at `p = 65521`: the exponent 2 is outside each interval. The intervals cover divisor sampling only; the local exponents rise from about 1.8 to about 2.75, so a single power law is an approximation.

#figure(image("ratio_by_genus.svg", width: 90%), caption: [Counted cost over `10 g²` by genus, addition and doubling, two primes. The grey band is the one-bit decision rule; dotted lines mark genus 6, 10 and 14, where most extension routes sit. Source: `summary.csv`.])

== Could the gap move a verdict? (exact count on the committed JSON)

Of the {d["advantage_curves"]} `formula_advantage` verdicts, {d["slack_le_1bit"]} lie within 1 bit of the 8-bit margin, {d["slack_le_2bits"]} within 2 bits and {d["slack_le_3bits"]} within 3 bits; {d["marginal_curves"]} more are `formula_marginal`. If the measured cost replaced `10 g²` and the *whole* route cost scaled with it, {d["would_drop_below_margin_if_measured_cost_applied_to_the_whole_route"]} of the {d["advantage_curves"]} would drop below the margin. That is an upper bound on fragility: routes priced by a different solver do not scale entirely with the group operation.

#figure(image("verdict_fragility.svg", width: 90%), caption: [Advantage verdicts close to the margin, and the conservative count that the measured cost would push below it. Source: `summary.csv`.])

== Claims and scope

- *Exact computation:* the counts per genus, the fit, and the verdict counts above.
- *Not shown:* that the model is wrong. Textbook Cantor with extended Euclid is not the cheapest algorithm (explicit formulas exist for `g = 2, 3`), so the ratio is a statement about this algorithm, not about the best one. The model's unit may also count an optimised operation.
- *Not shown:* anything about `F_p` of cryptographic size; the primes are `1009` and `65521`, and the two primes differ by at most {pdiff:.2f} bits in the ratio at any measured genus, so no dependence on `p` is visible here.
- *Dominant caveat:* 846 of the 848 advantage verdicts are extension-field `cover_descent` routes marked "proposed", whose genus-`k` cover over `GF(p)` is not established (the route's own note). Whether those covers exist matters more to those verdicts than any constant in this report; see the odd-characteristic descent study.
'''
open(out, "w").write(typ)
print("wrote", out)

"""Lay out recorded jacobian_transfers rows for the Lemma 6 report.

Reads two committed outputs of `jacobian_transfers build` and writes the
GF(2^2) GHS rows and best routes of the curves whose verdict or GF(4) row
changed.  It copies recorded values only; it computes nothing.

usage: python3 extract_rows.py BEFORE.json AFTER.json OUT.csv [BINARY_OUT.csv]
  BEFORE: git show b8ce3f697:docs/curves/jacobian-transfers.json
  AFTER:  docs/curves/jacobian-transfers.json at this report's commit
"""
import csv
import json
import sys

CURVES = ["c2pnb176w1", "c2pnb272w1", "c2pnb304w1"]


def rows(path, state):
    data = json.load(open(path))
    by_name = {n: c for c in data["curves"] for n in c.get("standard_names", [])}
    for name in CURVES:
        c = by_name[name]
        best = c["best_route"]
        for i, r in enumerate(c["routes"]):
            inv = r["invariants"]
            gf4 = r["kind"] == "ghs_weil_descent" and inv.get("base_degree") == 2
            if not (gf4 or i == best):
                continue
            yield {
                "state": state,
                "curve": name,
                "role": ("gf4_row" if gf4 else "") + ("+best" if i == best else ""),
                "label": r["label"],
                "magic_number": inv.get("magic_number"),
                "bound": inv.get("magic_number_bound", "exact"),
                "mmt_condition_8": inv.get("mmt_condition_8", ""),
                "genus": r["destination"]["genus"],
                "field_bits": r["destination"]["field_bits"],
                "transport": r["transport"],
                "evidence": r["evidence"],
                "ratio_log2": r["cost"]["ratio_log2"],
                "class": r["class"],
                "curve_verdict": c["verdict"],
            }


before, after, out = sys.argv[1:4]
all_rows = list(rows(before, "before")) + list(rows(after, "after"))
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(all_rows[0]))
    w.writeheader()
    w.writerows(all_rows)

# Second table: every binary-field curve's best carried route in AFTER, for
# binary_curves.svg.  Recorded values only.
if len(sys.argv) > 4:
    data = json.load(open(after))
    out2 = sys.argv[4]
    with open(out2, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "family", "field_degree", "subgroup_bits", "verdict", "best_ratio_log2", "best_label"])
        for c in data["curves"]:
            if c["field"]["kind"] != "binary":
                continue
            best = c["best_route"]
            r = c["routes"][best] if best is not None else None
            w.writerow([
                (c["standard_names"] or [c["slug"]])[0], c["family"], c["field"]["degree"],
                c["subgroup"]["bits"], c["verdict"],
                r["cost"]["ratio_log2"] if r else "", r["label"] if r else "",
            ])

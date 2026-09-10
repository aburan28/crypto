#!/usr/bin/env python3
"""Turn experiments/20_residual_walk_panel.json into the markdown tables
used by RESEARCH_RESIDUAL_WALKS.md.

    python3 scripts/summarize_residual_walk_panel.py experiments/20_residual_walk_panel.json
"""
import json
import math
import sys
from collections import defaultdict

TAGS = ["A", "B", "C1", "C2", "R"]


def mean(xs):
    xs = list(xs)
    return sum(xs) / len(xs) if xs else float("nan")


def fmt(v, digits=1):
    if isinstance(v, float):
        if math.isnan(v):
            return "-"
        if abs(v) >= 1000:
            return f"{v:,.0f}"
        return f"{v:.{digits}f}"
    return str(v)


def recovery(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[(r["bits"], r["factor_base"], r["tag"])].append(r)
    keys = sorted({(b, f) for (b, f, _) in groups})
    print("### P1 — operations to recover d\n")
    print("Mean over seeds. `samples/pred` divides the residual count by "
          "`√(2n(B+1))` (A, B, C1, C2) or by `√(πn/2)` (R). `ops/rho` divides "
          "total group operations by `√(πn/2)`; `×R` divides them by plain "
          "rho's *measured* total on the same instances.\n")
    print("| bits | B | tag | seeds | samples | samples/pred | ops/sample | total ops | ops/rho | ×R | full | trivial | correct |")
    print("|---:|---:|:--|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|")
    for (bits, fb) in keys:
        r_ops = mean(x["total_ops"] for x in groups.get((bits, fb, "R"), []))
        for tag in TAGS:
            g = groups.get((bits, fb, tag))
            if not g:
                continue
            pred = "predicted_rho_steps" if tag == "R" else "predicted_birthday_samples"
            print("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
                bits, fb, tag, len(g),
                fmt(mean(x["samples"] for x in g)),
                fmt(mean(x["samples"] / x[pred] for x in g), 2),
                fmt(mean(x["total_ops"] / x["samples"] for x in g)),
                fmt(mean(x["total_ops"] for x in g)),
                fmt(mean(x["total_ops"] / x["predicted_rho_steps"] for x in g)),
                fmt(mean(x["total_ops"] for x in g) / r_ops if r_ops else float("nan")),
                fmt(mean(x["full_decompositions"] for x in g)),
                fmt(mean(x["collisions_trivial"] for x in g)),
                "yes" if all(x["correct"] for x in g) else "NO",
            ))
    print()
    print("Breakdown of the operation count (means, same runs):\n")
    print("| bits | B | tag | setup | walk | replay | verify | LA ms | wall ms | table |")
    print("|---:|---:|:--|---:|---:|---:|---:|---:|---:|---:|")
    for (bits, fb) in keys:
        for tag in TAGS:
            g = groups.get((bits, fb, tag))
            if not g:
                continue
            print("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
                bits, fb, tag,
                fmt(mean(x["setup_ops"] for x in g)),
                fmt(mean(x["walk_ops"] for x in g)),
                fmt(mean(x["replay_ops"] for x in g)),
                fmt(mean(x["verify_ops"] for x in g)),
                fmt(mean(x["linear_algebra_ms"] for x in g)),
                fmt(mean(x["wall_ms"] for x in g)),
                fmt(mean(x["table_entries"] for x in g)),
            ))
    print()


def fixed_budget(rows):
    if not rows:
        return
    r0 = rows[0]
    print("### P2 — relation yield at a fixed budget\n")
    print(f"n = {r0['n']} (2^{math.log2(r0['n']):.1f}), B = {r0['factor_base']}, "
          f"budget = {rows[0]['total_ops']:,} group operations, no early stop.\n")
    print("| tag | samples | accepted | independent | rank | solved | ops at solve | full | trivial |")
    print("|:--|---:|---:|---:|---:|:--|---:|---:|---:|")
    for r in rows:
        print("| {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
            r["tag"], f"{r['samples']:,}", f"{r['accepted']:,}", r["relations_independent"], r["rank"],
            "yes" if r["correct"] else "no",
            f"{r['ops_at_solve']:,}" if r["ops_at_solve"] else "-",
            r["full_decompositions"], r["collisions_trivial"]))
    print()


def rank_cap(rows):
    if not rows:
        return
    r0 = rows[0]
    print("### P3 — rank of the r-adding residual walk vs r\n")
    print(f"n = {r0['n']}, B = {r0['factor_base']}, budget {r0['total_ops']:,} ops, no early stop.\n")
    print("| r | samples | verified | independent | dependent | rank | full hits | rank ≤ r+1+full | solved | ops at solve |")
    print("|---:|---:|---:|---:|---:|---:|---:|:--|:--|---:|")
    for r in rows:
        cap = r["multipliers"] + 1 + r["full_decompositions"]
        print("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
            r["multipliers"], f"{r['samples']:,}", r["relations_verified"], r["relations_independent"],
            r["relations_dependent"], r["rank"], r["full_decompositions"],
            "yes" if r["rank"] <= cap else "NO",
            "yes" if r["correct"] else "no",
            f"{r['ops_at_solve']:,}" if r["ops_at_solve"] else "-"))
    print()


def filter_trap(rows):
    if not rows:
        return
    base = rows[0]
    print("### P4 — rejecting residuals outside x < M (strategy A)\n")
    print(f"n = {base['n']}, p = {base['p']}, B = {base['factor_base']}. "
          "Predicted cost ratio against the unfiltered run is `√(p/M)`.\n")
    print("| M | samples | accepted | accepted/pred | total ops | measured ratio | predicted ratio | correct |")
    print("|:--|---:|---:|---:|---:|---:|---:|:--|")
    for r in rows:
        m = r["filter_bound"]
        if m is None:
            label, pred_ratio, acc_pred = "none (M = p)", 1.0, r["predicted_birthday_samples"]
        else:
            label = f"p / 2^{round(math.log2(r['p'] / m))}"
            pred_ratio = math.sqrt(r["p"] / m)
            acc_pred = math.sqrt(2 * m * (r["factor_base"] + 1))
        print("| {} | {} | {} | {} | {} | {} | {} | {} |".format(
            label, f"{r['samples']:,}", f"{r['accepted']:,}", fmt(r["accepted"] / acc_pred, 2),
            f"{r['total_ops']:,}", fmt(r["total_ops"] / base["total_ops"], 2), fmt(pred_ratio, 2),
            "yes" if r["correct"] else "no"))
    print()


def dps(rows):
    if not rows:
        return
    r0 = rows[0]
    print("### P5 — distinguished points on the collision-preserving walks\n")
    print(f"n = {r0['n']}, B = {r0['factor_base']}.\n")
    print("| tag | dp bits | samples | stored | stored/samples | total ops | replay ops | walks | abandoned | correct |")
    print("|:--|---:|---:|---:|---:|---:|---:|---:|---:|:--|")
    for r in rows:
        print("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
            r["tag"], r["dp_bits"], f"{r['samples']:,}", f"{r['table_entries']:,}",
            fmt(r["table_entries"] / r["samples"], 4), f"{r['total_ops']:,}", f"{r['replay_ops']:,}",
            r["walks"], r["abandoned_walks"], "yes" if r["correct"] else "no"))
    print()


def mitm(rows):
    if not rows:
        return
    print("### P6 — meet in the middle on 4-decompositions\n")
    print("| bits | n | B | pair table | coincidences | targets | relations | independent | total ops | ops/relation | ops/√(πn/2) | pred. matches/target | correct |")
    print("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|")
    for r in rows:
        rho = math.sqrt(math.pi * r["n"] / 2)
        print("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
            r["bits"], r["n"], r["factor_base"], f"{r['pair_table']:,}", r["pair_coincidences"], r["targets"],
            r["relations_verified"], r["relations_independent"], f"{r['total_ops']:,}",
            fmt(r["ops_per_independent_relation"] or float("nan")), fmt(r["total_ops"] / rho),
            fmt(r["predicted_matches_per_target"], 2), "yes" if r["correct"] else "no"))
    print()


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "experiments/20_residual_walk_panel.json"
    with open(path) as f:
        panel = json.load(f)
    recovery(panel["recovery"])
    fixed_budget(panel["fixed_budget"])
    rank_cap(panel["rank_cap"])
    filter_trap(panel["filter_trap"])
    dps(panel["distinguished_points"])
    mitm(panel["mitm"])


if __name__ == "__main__":
    main()

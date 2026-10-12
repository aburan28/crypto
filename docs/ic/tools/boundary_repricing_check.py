#!/usr/bin/env python3
"""Check that a repricing repriced and nothing else.

usage: boundary_repricing_check.py <before.json> <after.json> [--out <json>]

Round 4 changed the unit's conversion ratios and nothing else, so every
native counter of every row must be identical between the two runs and
only the converted numbers may move.  That is the claim, and a claim of
the form "this is accounting, not a gain" is worth checking rather than
asserting: a repricing that quietly changed a trial count would be a
change to the method wearing an accountant's coat.

For each row paired by `(regime, instance, variant, repeat)` this
compares every field of `group_ops` and every key of `native`, in every
phase, and reports:

  * how many rows are identical in all of them, and the first field of
    the first row that is not;
  * the distribution of `S_after / S_before`, which is the size of the
    correction the repricing applied;
  * the rows where it is largest, since those are the ones whose earlier
    figures a reader should not compare across runs.

Wall-clock fields are excluded: they are a practicality note, they never
priced anything, and they are expected to differ.
"""
import json, math, sys, statistics
from collections import OrderedDict

PHASES = ("factor_base", "relations", "linear_algebra", "verify")


def rows_of(run):
    out = {}
    for inst in run["ledger"]["instances"]:
        seen = OrderedDict()
        for v in inst["variants"]:
            seen.setdefault(v["name"], []).append(v)
        out[(inst["regime"], inst["curve"]["name"])] = (inst, seen)
    return out


def counters(v):
    """Every count in a row, with nothing converted and no wall clock."""
    flat = {"trials": v["trials"], "relations_found": v["relations_found"],
            "independent": v["independent"], "dependent": v["dependent"],
            "rank": v["rank"], "recovered": v.get("recovered"),
            "verified": v["verified"], "signed_points": v["signed_points"],
            "columns": v["columns"], "summands": v["summands"]}
    for p in PHASES:
        if p not in v:
            continue
        for k, x in v[p]["group_ops"].items():
            flat[f"{p}.group_ops.{k}"] = x
        for k, x in v[p]["native"].items():
            flat[f"{p}.native.{k}"] = x
    return flat


before = json.load(open(sys.argv[1]))
after = json.load(open(sys.argv[2]))
out_path = next((sys.argv[i + 1] for i, a in enumerate(sys.argv) if a == "--out"), None)

A, B = rows_of(before), rows_of(after)
identical, mismatches, ratios, per_row = 0, [], [], []
for key in B:
    if key not in A:
        continue
    (inst_a, rows_a), (inst_b, rows_b) = A[key], B[key]
    for name, runs_b in rows_b.items():
        runs_a = rows_a.get(name)
        if not runs_a:
            continue
        for i in range(min(len(runs_a), len(runs_b))):
            ca, cb = counters(runs_a[i]), counters(runs_b[i])
            bad = [k for k in set(ca) | set(cb) if ca.get(k) != cb.get(k)]
            if bad:
                mismatches.append((key, name, i, sorted(bad)[:4]))
            else:
                identical += 1
            r = runs_b[i]["s"] / runs_a[i]["s"] if runs_a[i]["s"] else float("nan")
            if not math.isnan(r):
                ratios.append(r)
                per_row.append((r, key[0], key[1], name))

print(f"rows compared: {identical + len(mismatches)}")
print(f"identical in every counter: {identical}")
print(f"rows whose counters moved: {len(mismatches)}")
for key, name, i, bad in mismatches[:8]:
    print(f"   {key[0]} {key[1]} {name} repeat {i}: {bad}")

if ratios:
    lo, hi = min(ratios), max(ratios)
    print()
    print(f"S after / S before: median {statistics.median(ratios):.4f}, "
          f"range {lo:.4f} to {hi:.4f}")
    print(f"rows repriced by more than 2%: {sum(1 for r in ratios if abs(r - 1) > 0.02)} of {len(ratios)}")
    print("\nlargest repricings:")
    for r, regime, inst, name in sorted(per_row, key=lambda x: -abs(math.log(x[0])))[:10]:
        print(f"   {r:7.4f}  {regime:8s} {inst[:26]:26s} {name}")

if out_path:
    doc = OrderedDict([
        ("what", "Check that a repricing of the boundary ledger's unit changed only converted "
                 "numbers: every native counter and group-operation count of every row compared "
                 "between the two runs, and the distribution of the correction it applied."),
        ("before", sys.argv[1]),
        ("after", sys.argv[2]),
        ("rows_compared", identical + len(mismatches)),
        ("rows_identical_in_every_counter", identical),
        ("rows_whose_counters_moved", len(mismatches)),
        ("counter_mismatches", [{"regime": k[0], "instance": k[1], "variant": n,
                                 "repeat": i, "fields": b} for k, n, i, b in mismatches]),
        ("s_after_over_before_median", round(statistics.median(ratios), 6) if ratios else None),
        ("s_after_over_before_min", round(min(ratios), 6) if ratios else None),
        ("s_after_over_before_max", round(max(ratios), 6) if ratios else None),
        ("rows_repriced_above_2_percent", sum(1 for r in ratios if abs(r - 1) > 0.02)),
        ("largest_repricings", [{"ratio": round(r, 6), "regime": g, "instance": i, "variant": n}
                                for r, g, i, n in sorted(per_row, key=lambda x: -abs(math.log(x[0])))[:40]]),
    ])
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    print(f"\nwritten {out_path}")

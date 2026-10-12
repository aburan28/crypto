#!/usr/bin/env python3
"""How much of a round's move is the round, and how much is the draw?

usage: boundary_repeat_spread.py <run.json> [<run.json> …] [--min-repeats N]

Every row of a boundary-ledger run is a mean over repeats that differ
only in the target they were handed.  A two-summand relation search
closes when the factor-base graph first carries a cycle, and how long
that takes is a property of the draw, not of the code — so the same row,
same instance and same binary can cost a factor of several more on one
repeat than on another.  That spread is the resolution of every
comparison built on these runs: a round whose speedup is inside it has
not been measured, however cleanly the means separate.

For each row this prints `max/min` of the per-repeat total operation
count, and summarises the distribution over rows, split by summand
count, because the `m = 2` rows are the ones that close by cycles.
"""
import json, math, sys
from collections import OrderedDict

paths = [a for a in sys.argv[1:] if not a.startswith("--")]
min_repeats = int(next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--min-repeats=")), "2"))


def rows_of(inst):
    seen = OrderedDict()
    for v in inst["variants"]:
        seen.setdefault(v["name"], []).append(v)
    return seen


def quantile(xs, q):
    if not xs:
        return float("nan")
    xs = sorted(xs)
    i = q * (len(xs) - 1)
    lo, hi = math.floor(i), math.ceil(i)
    return xs[lo] if lo == hi else xs[lo] + (xs[hi] - xs[lo]) * (i - lo)


rows, per_m = [], {}
for path in paths:
    run = json.load(open(path))
    for inst in run["ledger"]["instances"]:
        for name, runs in rows_of(inst).items():
            if len(runs) < min_repeats:
                continue
            ops = [r["total_gae"] for r in runs]
            spread = max(ops) / min(ops)
            rows.append(OrderedDict([
                ("run", path.split("/")[-1]),
                ("regime", inst["regime"]),
                ("instance", inst["curve"]["name"]),
                ("variant", name),
                ("m", runs[0]["summands"]),
                ("repeats", len(runs)),
                ("spread_max_over_min", round(spread, 3)),
            ]))
            per_m.setdefault(runs[0]["summands"], []).append(spread)

print(f"| repeats | rows | m | median max/min | 75th | 90th | max |")
print(f"|:--|--:|--:|--:|--:|--:|--:|")
for m in sorted(per_m):
    xs = per_m[m]
    print(f"| ≥{min_repeats} | {len(xs)} | {m} | {quantile(xs, 0.5):.2f} | {quantile(xs, 0.75):.2f} "
          f"| {quantile(xs, 0.9):.2f} | {max(xs):.2f} |")
allx = [x for xs in per_m.values() for x in xs]
print(f"| ≥{min_repeats} | {len(allx)} | all | {quantile(allx, 0.5):.2f} | {quantile(allx, 0.75):.2f} "
      f"| {quantile(allx, 0.9):.2f} | {max(allx):.2f} |")

print("\nWidest rows:\n")
print("| run | regime | instance | variant | m | repeats | max/min |")
print("|:--|:--|:--|:--|--:|--:|--:|")
for r in sorted(rows, key=lambda r: -r["spread_max_over_min"])[:12]:
    print(f"| {r['run']} | {r['regime']} | {r['instance']} | {r['variant']} | {r['m']} | {r['repeats']} "
          f"| {r['spread_max_over_min']:.2f} |")

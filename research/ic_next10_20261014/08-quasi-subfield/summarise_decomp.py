#!/usr/bin/env python3
"""Tabulate arm D (decomp_*.jsonl) as Markdown and decomp_summary.json.
Usage: python3 -I summarise_decomp.py <dir>"""
import glob
import json
import math
import sys

D = sys.argv[1]


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 2
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


rows = []
for p in sorted(glob.glob(f"{D}/decomp_*.jsonl")):
    for line in open(p):
        if line.strip():
            rows.append(json.loads(line))
rows.sort(key=lambda r: (r["n"], r["k"], r["base"].startswith("random"), r["base"]))
out = []
print("| n | k | base | V on curve | yield (95% CI) | mean ordered decomps | refutation degrees (6 targets) | median word XORs |")
print("|--:|--:|---|--:|---|--:|---|--:|")
for r in rows:
    lo, hi = wilson(r["decomposed"], r["yield_targets"])
    degs = [t["refute_degree"] if t["refute_degree"] is not None else f">{r['d_max']}" for t in r["refute"]]
    steps = [t["steps"] for t in r["refute"]]
    xs = sorted(t["word_xors"] for t in r["refute"])
    med = xs[len(xs) // 2] if xs else None
    cap = any("CAP" in s for s in steps)
    print(f"| {r['n']} | {r['k']} | {r['base']} | {r['v_on_curve']}/{r['v_size']} | "
          f"{r['yield']:.3f} ({lo:.3f}–{hi:.3f}) | {r['mean_ordered_decomps']:.2f} | "
          f"{' '.join(map(str, degs))}{' (size cap hit)' if cap else ''} | {med} |")
    out.append(dict(n=r["n"], k=r["k"], base=r["base"], yield_=r["yield"], ci=[lo, hi],
                    mean_decomps=r["mean_ordered_decomps"], refute_degrees=degs, steps=steps,
                    median_word_xors=med))
json.dump(out, open(f"{D}/decomp_summary.json", "w"), indent=1)

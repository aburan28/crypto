#!/usr/bin/env python3
"""Plot recorded A/A (or A/B) timing rows: ns per pair by rung and round.

Lays out values the harness recorded (``pdp_speedup time --out``); computes
nothing new beyond the per-arm median drawn as a line.  Usage:

  plot_aa.py times.jsonl out.svg --a A --b B --title "..."
"""
import argparse, json, statistics
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("times")
ap.add_argument("out")
ap.add_argument("--a", default="A")
ap.add_argument("--b", default="B")
ap.add_argument("--title", default="")
ap.add_argument("--ylabel", default="ns per pair (single core)")
args = ap.parse_args()

rows = [json.loads(l) for l in open(args.times) if l.strip()]
by = {}
for r in rows:
    tag, rest = r["label"].rsplit("-l", 1)
    rnd = int(rest.split("-r")[1])
    by.setdefault((r["l"], tag), []).append((rnd, r["ns_per_pair"]))

ls = sorted({k[0] for k in by})
fig, ax = plt.subplots(figsize=(7.2, 3.8))
colors = {args.a: "#1f77b4", args.b: "#d95f02"}
for ti, tag in enumerate([args.a, args.b]):
    xs, meds = [], []
    for i, l in enumerate(ls):
        pts = sorted(by.get((l, tag), []))
        off = (ti - 0.5) * 0.18
        ax.scatter([i + off] * len(pts), [p[1] for p in pts], s=22, color=colors[tag],
                   label=f"{tag} (each round)" if i == 0 else None, zorder=3)
        if pts:
            m = statistics.median(p[1] for p in pts)
            ax.hlines(m, i + off - 0.08, i + off + 0.08, color="black", lw=1.4, zorder=4,
                      label="median" if (i == 0 and ti == 0) else None)
ax.set_xticks(range(len(ls)))
ax.set_xticklabels([f"l={l}\n(n={3*l})" for l in ls])
ax.set_ylabel(args.ylabel)
ax.set_title(args.title, fontsize=9)
ax.grid(axis="y", alpha=0.3)
ax.legend(fontsize=8, loc="upper left")
n_rounds = max(len(v) for v in by.values())
fig.text(0.01, 0.01, f"recorded rows: {len(rows)}; rounds per arm per rung: {n_rounds}", fontsize=7)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(args.out)

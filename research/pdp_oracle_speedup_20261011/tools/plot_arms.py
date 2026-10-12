#!/usr/bin/env python3
"""Plot recorded timing rows by arm: per-round points and per-arm medians against the rung.

  plot_arms.py times.jsonl out.svg [--metric ns_per_target|ns_per_pair] [--title T] [--arms a,b,c]
Contended runs are drawn hollow (read from the isolated_bench records next to the times file).
"""
import argparse, json, os, statistics
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("times"); ap.add_argument("out")
ap.add_argument("--metric", default="ns_per_target")
ap.add_argument("--title", default="")
ap.add_argument("--arms", default="")
ap.add_argument("--iso", default="")
a = ap.parse_args()
iso_path = a.iso or os.path.join(os.path.dirname(a.times), "iso.jsonl")
cont = {}
if os.path.exists(iso_path):
    for l in open(iso_path):
        if l.strip():
            r = json.loads(l); cont[r["label"]] = bool(r["run"]["contended"])
rows = [json.loads(l) for l in open(a.times) if l.strip()]
tags = {}
for r in rows:
    tag = r["label"].rsplit("-l", 1)[0]
    tags.setdefault(tag, []).append(r)
order = [t for t in (a.arms.split(",") if a.arms else tags) if t in tags]
ls = sorted({r["l"] for r in rows})
fig, ax = plt.subplots(figsize=(7.4, 4.0))
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
for ti, tag in enumerate(order):
    for i, l in enumerate(ls):
        pts = [r for r in tags[tag] if r["l"] == l]
        off = (ti - (len(order) - 1) / 2) * 0.13
        for r in pts:
            hollow = cont.get(r["label"], False)
            ax.scatter([i + off], [r[a.metric]], s=18, facecolors="none" if hollow else colors[ti % 10],
                       edgecolors=colors[ti % 10], zorder=3, label=tag if (i == 0 and r is pts[0]) else None)
        if pts:
            m = statistics.median(r[a.metric] for r in pts)
            ax.hlines(m, i + off - 0.05, i + off + 0.05, color="black", lw=1.2, zorder=4)
ax.set_yscale("log")
ax.set_xticks(range(len(ls))); ax.set_xticklabels([f"l={l}\n(n={3*l})" for l in ls])
ax.set_ylabel({"ns_per_target": "ns per refuted target", "ns_per_pair": "ns per own pair"}.get(a.metric, a.metric))
ax.set_title(a.title, fontsize=9); ax.grid(axis="y", alpha=0.3); ax.legend(fontsize=7, ncol=2)
fig.text(0.01, 0.01, f"{len(rows)} recorded runs; hollow = marked contended by isolated_bench; black bars = per-arm medians", fontsize=7)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(a.out)

#!/usr/bin/env python3
"""Plot recorded gate counts per rung (pdp_speedup gates): bit-sliced group-law
task vs the baseline formulation derived from counted operations.
  plot_gates.py gates.jsonl out.svg
"""
import json, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

rows = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
ls = [r["l"] for r in rows]
f3 = [r["task_fused"] for r in rows]
f1 = [r["f1_gcd_derived"]["gates"] for r in rows]
fig, ax = plt.subplots(figsize=(7.2, 3.6))
w = 0.38
xs = range(len(ls))
ax.bar([x - w / 2 for x in xs], f1, w, label="baseline formulation (gcd), derived from counted operations", color="#d95f02")
ax.bar([x + w / 2 for x in xs], f3, w, label="group-law formulation, one task (counted by the gate-counting lane)", color="#1f77b4")
ax.axhline(1.8e5, color="k", ls="--", lw=1)
ax.text(len(ls) - 0.5, 1.8e5, " 1.8e5 (plan: bit-slicing cannot reach 5x)", ha="right", va="bottom", fontsize=7)
ax.axhline(9e4, color="gray", ls=":", lw=1)
ax.text(len(ls) - 0.5, 9e4, " 9e4 (plan: 10x budget at n=33)", ha="right", va="bottom", fontsize=7)
ax.set_xticks(list(xs))
ax.set_xticklabels([f"l={r['l']}\n(n={r['n']})" for r in rows])
ax.set_ylabel("gates per evaluation (fused ternary-logic count)")
ax.set_yscale("log")
ax.grid(axis="y", alpha=0.3)
ax.legend(fontsize=7, loc="upper left")
ax.set_title("Gates per pair-evaluation of the root test, by formulation", fontsize=9)
fig.tight_layout()
fig.savefig(sys.argv[2])

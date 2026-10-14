#!/usr/bin/env python3
"""Plot recorded per-pair operation counts (pdp_speedup count) against the rung.

Lays out values the instrumented build recorded; computes only the division by pairs
that the harness itself printed under "per_pair".
  plot_counts.py counts.jsonl out.svg
"""
import json, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

rows = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
ls = [r["l"] for r in rows]
series = [("pclmul", "carry-less multiplies (PCLMULQDQ)"), ("mul", "field multiplications"),
          ("sqr", "field squarings"), ("reduce", "stand-alone reductions"),
          ("clmul", "stand-alone unreduced products")]
fig, ax = plt.subplots(figsize=(7.2, 3.6))
for key, label in series:
    ax.plot(ls, [r["per_pair"][key] for r in rows], marker="o", label=label)
ax.set_xlabel("rung l (n = 3l)")
ax.set_ylabel("counted operations per pair")
ax.set_xticks(ls)
ax.set_yscale("log")
ax.grid(alpha=0.3)
ax.legend(fontsize=7)
ax.set_title("Baseline pairs-and-solve: counted field operations per pair (instrumented build)", fontsize=9)
fig.tight_layout()
fig.savefig(sys.argv[2])

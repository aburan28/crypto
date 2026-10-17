"""Lay out recorded values from experiments/koblitz_m3/*.jsonl as SVG.
Plots recorded numbers only; computes no result."""
import json, sys
from pathlib import Path
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[4]
here = Path(__file__).parent
cmp_rows = [json.loads(l) for l in open(root / "experiments/koblitz_m3/solver_compare_n17.jsonl")]
full = [json.loads(l) for l in open(root / "experiments/koblitz_m3/full_solves_e3.jsonl")]

# Figure 1: per-call seconds vs l, n = 17, by solver (timeouts marked)
fig, ax = plt.subplots(figsize=(6.4, 4.0))
style = {"SymmetricE3": ("e3 enumeration", "o-"), "SymmetricGroebner": ("symmetric Groebner", "s-"), "Chain": ("S3 chain (baseline)", "^-")}
names = {"e3": "SymmetricE3", "sym": "SymmetricGroebner", "chain": "Chain"}
for key, (label, fmt) in style.items():
    pts = [(r["l"], r["us_per_call"] / 1e6) for r in cmp_rows if r.get("m3_solver") == key and "us_per_call" in r]
    if pts:
        ls, ys = zip(*sorted(pts))
        ax.plot(ls, ys, fmt, label=label)
    for r in cmp_rows:
        if names.get(r.get("m3_solver")) == key and "timed_out_after_s" in r:
            ax.plot([r["l"]], [r["timed_out_after_s"] / r["probes_requested"]], "x", color="red", markersize=10)
            ax.annotate("timeout: >75 s/call", (r["l"], r["timed_out_after_s"] / r["probes_requested"]), textcoords="offset points", xytext=(-90, 6), fontsize=8)
ax.set_yscale("log")
ax.set_xlabel("l (dim V)")
ax.set_ylabel("seconds per decomposition call (measured)")
ax.set_xticks([3, 4, 5])
ax.set_title("m = 3, n = 17, a6 = 1, monomial V, 4 probes per point")
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(here / "per_call_n17.svg")

# Figure 2: full-solve wall seconds, n = 17 and 19, by V
fig, ax = plt.subplots(figsize=(6.4, 3.6))
for i, (n, geo) in enumerate([(17, False), (17, True), (19, False), (19, True)]):
    ys = [r["seconds"] for r in full if r["n"] == n and r["geometric"] == geo]
    ax.scatter([i] * len(ys), ys, color="C0" if not geo else "C1")
ax.set_xticks(range(4))
ax.set_xticklabels(["n=17 mono", "n=17 geom", "n=19 mono", "n=19 geom"])
ax.set_yscale("log")
ax.set_ylabel("wall seconds per full solve (measured)")
ax.set_title("m = 3 full IC solves with e3 enumeration; each point one curve; all BSGS-verified")
ax.title.set_fontsize(9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(here / "full_solves.svg")
print("ok")

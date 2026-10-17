"""Draws reports/<run>/table.tsv (from summarize.py); no new computation."""
import sys, csv
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
run = sys.argv[1]
rows = list(csv.DictReader(open(f"{run}/table.tsv"), delimiter="\t"))
cells = [r["cell"] for r in rows if r["engine"] == "f4"]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for ax, e in zip(axes, ("f4", "f6p")):
    rr = [r for r in rows if r["engine"] == e]
    x = range(len(rr))
    ax.axhline(1, color="#888", lw=0.8)
    lo = [float(r["C_scatter/runs"]) - float(r["C_scatter/runs_lo"]) for r in rr]
    hi = [float(r["C_scatter/runs_hi"]) - float(r["C_scatter/runs"]) for r in rr]
    ax.errorbar([i - 0.12 for i in x], [float(r["C_scatter/runs"]) for r in rr], yerr=[lo, hi], fmt="o", capsize=3, color="#1f6fb4",
                label="scatter / runs (the non-BMI2 path before / after)")
    ax.plot([i + 0.12 for i in x], [float(r["C_bmi2/runs"]) for r in rr], "s", color="#d1495b", label="bmi2 / runs")
    ax.set_xticks(list(x)); ax.set_xticklabels(cells, rotation=45, ha="right", fontsize=8)
    ax.set_title({"f4": "inherited F4", "f6p": "F6-IC + pair index"}[e], fontsize=9)
axes[0].set_ylabel("median CPU ratio (> 1 = runs cheaper)")
axes[0].legend(fontsize=7.5, loc="upper left")
fig.suptitle("Eager rewrite through shift runs, x86 proxy, 9 cells × 5 rounds (run-20261016T2335Z)", fontsize=9)
fig.tight_layout(); fig.savefig(f"{run}/fig_runs.svg"); fig.savefig(f"{run}/fig_runs.png", dpi=150)

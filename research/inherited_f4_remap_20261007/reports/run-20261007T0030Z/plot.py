"""Draw the recorded values of the inherited-F4 remap round as SVG.

Plots only numbers already written to ../../results/: the per-input
Callgrind ratios in run-1/ir.tsv (column `ratio`, written by run.sh), and the
raw instruction counts in main-7666b708-check.tsv.  Nothing is computed
beyond converting instructions to units of 10^9 for the axis.
Usage: python plot.py   (from this directory)
"""
import csv
import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt

RES = "../../results"


def arm(input_id):
    # Arm from the frozen spec layout: ladder seq = 1 (mod 3) is inherited F4,
    # 0 (mod 3) is F6-IC; holdout seq even is inherited F4, odd is F6-IC.
    s, q = input_id.split("-")
    q = int(q)
    if s.startswith("L"):
        return "inherited F4" if q % 3 == 1 else "F6-IC"
    return "inherited F4" if q % 2 == 0 else "F6-IC"


rows = list(csv.DictReader(open(f"{RES}/run-1/ir.tsv"), delimiter="\t"))
groups = ["L13", "H13", "H13b", "L19", "H19", "H19b"]
fig, ax = plt.subplots(figsize=(8.5, 4.2))
colors = {"inherited F4": "#2a6fb0", "F6-IC": "#d0702a"}
x = 0
ticks, labels = [], []
for g in groups:
    sel = [r for r in rows if r["input"].split("-")[0] == g]
    start = x
    for r in sorted(sel, key=lambda r: int(r["input"].split("-")[1])):
        ax.bar(x, float(r["ratio"]), color=colors[arm(r["input"])], width=0.8)
        x += 1
    ticks.append((start + x - 1) / 2)
    labels.append(f"{g}\n(n={len(sel)})")
    x += 1.5
ax.axhline(0.97, color="#444", linestyle="--", linewidth=1)
ax.text(x - 1, 0.972, "preregistered threshold 0.97 (pooled, confirmatory m=19)", ha="right", fontsize=8)
ax.axhline(1.0, color="#999", linewidth=0.8)
ax.set_ylim(0.85, 1.01)
ax.set_xticks(ticks, labels, fontsize=8)
ax.set_ylabel("candidate / baseline Callgrind Ir (whole process)")
ax.set_title("Inherited-F4 remap round, run-1: per-input instruction ratio\nbaseline a3652edd, one Callgrind run per binary and input (deterministic)", fontsize=10)
for a, c in colors.items():
    ax.bar(0, 0, color=c, label=a)
ax.legend(fontsize=8, loc="lower left")
fig.tight_layout()
fig.savefig("ir_ratio_per_input.svg")

check = [r for r in csv.DictReader((l for l in open(f"{RES}/main-7666b708-check.tsv") if not l.startswith("#")), delimiter="\t")]
fig, ax = plt.subplots(figsize=(7.5, 3.8))
series = [("round2_baseline_ir", "round-2 baseline a3652edd", "#999999"),
          ("round2_candidate_ir", "round-2 candidate (not merged)", "#2a6fb0"),
          ("main_7666b708_ir", "main 7666b708e (#1646, #1671)", "#3a9a5a")]
for i, r in enumerate(check):
    for j, (k, lab, c) in enumerate(series):
        ax.bar(i * 4 + j, int(r[k]) / 1e9, color=c, label=lab if i == 0 else None)
ax.set_xticks([i * 4 + 1 for i in range(len(check))], [r["input"] for r in check])
ax.set_yscale("log")
ax.set_ylabel("Callgrind Ir, 10^9 (log scale)")
ax.set_title("Superseding check 2026-10-11: same frozen inputs, three binaries\n(main's counted outputs differ from the round-2 baseline; not a gated A/B)", fontsize=10)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("main_check_ir.svg")

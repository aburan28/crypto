"""Draws reports/<run>/table.tsv (from summarize.py); no new computation."""
import sys, csv
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
run = sys.argv[1]
rows = list(csv.DictReader(open(f"{run}/table.tsv"), delimiter="\t"))
f6 = [r for r in rows if r["engine"] == "f6p"]
f4 = [r for r in rows if r["engine"] == "f4"]
x = range(len(f6))
fig, ax = plt.subplots(figsize=(9, 4.4))
ax.axhline(1, color="#888", lw=0.8)
lo = [float(r["C_off/gate"]) - float(r["C_off/gate_lo"]) for r in f6]
hi = [float(r["C_off/gate_hi"]) - float(r["C_off/gate"]) for r in f6]
ax.errorbar([i - 0.15 for i in x], [float(r["C_off/gate"]) for r in f6], yerr=[lo, hi], fmt="o", color="#1f6fb4", capsize=3,
            label="F6-IC CPU off / gate (median; bar = worst–best round pairing)")
ax.plot([i + 0.15 for i in x], [float(r["U_off/gate"]) for r in f6], "s", color="#d1495b", label="F6-IC stage unit off / gate (exact)")
ax.plot([i for i in x], [float(r["C_off/gate"]) for r in f4], "^", color="#999", mfc="none", label="inherited F4 CPU off / gate (median)")
ax.set_xticks(list(x)); ax.set_xticklabels([r["cell"] for r in f6], fontsize=8)
ax.set_ylabel("ratio (> 1 = gate-first cheaper)")
ax.set_title("Gate-first, x86 (BMI2, eager rewrite): 9 cells, 5 rounds × 3 reps (run-20261016T2232Z)", fontsize=9)
ax.legend(fontsize=7.5, loc="upper left")
fig.tight_layout(); fig.savefig(f"{run}/fig_gate.svg"); fig.savefig(f"{run}/fig_gate.png", dpi=150)
fig, ax = plt.subplots(figsize=(9, 4.0))
ax.axhline(1, color="#888", lw=0.8)
for k, (col, lab, c) in enumerate([("U_off/gate", "gate", "#1f6fb4"), ("U_off/rref", "rref root", "#e09f3e"), ("U_off/gaterref", "gate + rref root", "#2a9d8f")]):
    ax.plot([i + (k - 1) * 0.2 for i in x], [float(r[col]) for r in f6], "o", color=c, label=f"F6-IC stage unit off / {lab}")
ax.set_xticks(list(x)); ax.set_xticklabels([r["cell"] for r in f6], fontsize=8)
ax.set_ylabel("stage-unit ratio (exact; > 1 = cheaper)")
ax.set_title("F6-IC stage units per arm, run-20261016T2232Z", fontsize=9)
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{run}/fig_arms.svg"); fig.savefig(f"{run}/fig_arms.png", dpi=150)

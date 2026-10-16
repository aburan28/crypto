"""Plot the recorded largest Macaulay-matrix width (max_cols) against b per arm,
from runs/*.jsonl only. Rendering aid only (AGENTS.md "no Python")."""
import glob, json, os, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
runs = os.path.join(here, "..", "..", "runs")
ARMS = ["A0", "N", "A1z", "A1g", "A4", "A11"]
COL = {"A0": "#2a78d6", "N": "#eb6834", "A1z": "#1baf7a", "A1g": "#eda100", "A4": "#e87ba4", "A11": "#008300"}
MK = {"A0": "o", "N": "s", "A1z": "^", "A1g": "v", "A4": "D", "A11": "P"}
OFF = {a: (i - 2.5) * 0.06 for i, a in enumerate(ARMS)}

rows = []
for f in sorted(glob.glob(os.path.join(runs, "*-m*.jsonl"))):
    rows += [json.loads(l) for l in open(f) if l.strip()]

plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#888", "axes.labelcolor": "#222",
                     "xtick.color": "#444", "ytick.color": "#444"})
fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0), sharey=True)
for ax, m in zip(axes, [2, 3]):
    for a in ARMS:
        ok = sorted([(r["b"], r["max_cols"]) for r in rows
                     if r["arm"] == a and r["m"] == m and r["status"] == "ok" and r["target"] == 0])
        to = [(r["b"], r["max_cols"]) for r in rows
              if r["arm"] == a and r["m"] == m and r["status"] == "timeout" and r["target"] == 0 and "max_cols" in r]
        if ok:
            ax.plot([b + OFF[a] for b, _ in ok], [c for _, c in ok], color=COL[a], lw=1.5,
                    marker=MK[a], ms=6, mec="white", mew=1, label=a)
        for b, c in to:
            ax.plot([b + OFF[a]], [c], marker=MK[a], ms=8, mfc="none", mec=COL[a], mew=1.5, ls="none")
    ax.set_yscale("log")
    ax.set_title(f"m = {m}")
    ax.set_xlabel("b = log2|V| per summand (bits)")
    ax.grid(True, color="#e6e6e6", lw=0.6)
axes[0].set_ylabel("largest F4 matrix width (max_cols, columns)")
axes[0].legend(ncol=3, fontsize=8, frameon=False, loc="lower right")
fig.suptitle("Largest F4 matrix width vs factor-base bits, p1 = 786433, target 0 of each cell (run-20261015T1657Z)", fontsize=9.5)
fig.text(0.5, 0.005, "Filled: completed systems. Open: stopped at the 300 s F4 budget (width reached so far, a lower bound). "
        , ha="center", fontsize=7.5, color="#555")
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
fig.savefig(os.path.join(here, "cols_vs_b.svg"))
fig.savefig(os.path.join(here, "cols_vs_b.png"), dpi=150)

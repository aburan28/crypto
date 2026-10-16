"""Plot recorded D against b per arm, from runs/*.jsonl only (no computation
beyond reading the rows). Rendering aid only (AGENTS.md "no Python")."""
import glob, json, os, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
runs = os.path.join(here, "..", "..", "runs")
ARMS = ["A0", "N", "A1z", "A1g", "A4", "A11"]
COL = {"A0": "#2a78d6", "N": "#eb6834", "A1z": "#1baf7a", "A1g": "#eda100", "A4": "#e87ba4", "A11": "#008300"}
MK = {"A0": "o", "N": "s", "A1z": "^", "A1g": "v", "A4": "D", "A11": "P"}
OFF = {a: (i - 2.5) * 0.06 for i, a in enumerate(ARMS)}  # horizontal jitter so overlaps stay visible

rows = []
for a in ARMS:  # the six arms' own files, not the F1-F4 files merged later
    for m in [2, 3]:
        rows += [json.loads(l) for l in open(os.path.join(runs, f"{a}-m{m}.jsonl")) if l.strip()]

plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#888", "axes.labelcolor": "#222",
                     "xtick.color": "#444", "ytick.color": "#444"})
fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0), sharey=True)
for ax, m in zip(axes, [2, 3]):
    for a in ARMS:
        ok = sorted([(r["b"], r["D"]) for r in rows if r["arm"] == a and r["m"] == m and r["status"] == "ok"])
        bad = [r for r in rows if r["arm"] == a and r["m"] == m and r["status"] != "ok"]
        if not ok and not bad:
            continue
        xs = [b + OFF[a] for b, _ in ok]
        ax.plot(xs, [d for _, d in ok], color=COL[a], lw=1.5, marker=MK[a], ms=6,
                mec="white", mew=1, label=a)
        for r in bad:
            if r.get("D") is not None:  # timeout: open marker at the partial D (a lower bound)
                ax.plot([r["b"] + OFF[a]], [r["D"]], marker=MK[a], ms=8,
                        mfc="none", mec=COL[a], mew=1.5, ls="none")
            else:  # memory abort: no D was recorded; a cross on the top band
                ax.plot([r["b"] + OFF[a]], [9.0], marker="x", ms=7, color=COL[a], mew=1.8, ls="none")
                ax.annotate(f"{a} abort\n(memory)", (r["b"] + OFF[a], 9.0), xytext=(-8, 0),
                            textcoords="offset points", ha="right", va="center", fontsize=6.5, color="#555")
    ax.set_title(f"m = {m}", color="#222")
    ax.set_xlabel("b = log2|V| per summand (bits)")
    ax.grid(color="#e5e5e5", lw=0.6)
    ax.set_axisbelow(True)
axes[0].set_ylabel("D = solving_degree_max (f4_fp, grevlex)")
axes[0].legend(frameon=False, ncol=3, fontsize=8, loc="lower right")
for ax in axes:
    ax.set_ylim(1.5, 9.6)
fig.suptitle("Solving degree vs factor-base bits, p1 = 786433, 2 planted targets per cell (run-20261015T1657Z)",
             fontsize=9.5, color="#222")
fig.text(0.5, 0.005, "Filled: completed systems (both targets plotted; overlapping points coincide). "
         "Open: F4 budget timeout, at its partial D (a lower bound). Cross at top: memory abort, no D recorded.",
         ha="center", fontsize=7.5, color="#555")
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
fig.savefig(os.path.join(here, "d_vs_b.svg"))
fig.savefig(os.path.join(here, "d_vs_b.png"), dpi=150)
print("wrote d_vs_b.svg/png from", len(rows), "rows")

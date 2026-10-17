"""Draws reports/<run>/table.tsv files (from summarize.py); no new computation."""
import sys, csv
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
runs = sys.argv[1:-1]; out = sys.argv[-1]
CELLS = [("form0", "1"), ("form0", "4"), ("form2", "1"), ("form2", "4")]
ARMS = [("off", "all three off (old defaults)", "#444"), ("nopack", "direct pack off", "#1f6fb4"),
        ("nounpack", "direct unpack off", "#e09f3e"), ("noreuse", "table reuse off", "#2a9d8f")]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
for ax, case in zip(axes, ("f5_n20_m20_d4", "f5_n24_m24_d4")):
    ax.axhline(1, color="#888", lw=0.8)
    for ri, run in enumerate(runs):
        rows = [r for r in csv.DictReader(open(f"{run}/table.tsv"), delimiter="\t") if r["case"] == case]
        for ai, (arm, lab, col) in enumerate(ARMS):
            for ci, (form, th) in enumerate(CELLS):
                rr = [r for r in rows if r["arm"] == arm and r["form"] == form and r["threads"] == th]
                if not rr: continue
                x = ci + (ai - 1.5) * 0.17 + ri * 0.05
                ax.plot(x, float(rr[0]["arm/on"]), "o" if ri else "^", color=col, mfc=col if ri else "none",
                        label=f"{lab} (run {ri + 1})" if ci == 0 else None, ms=5)
    ax.set_xticks(range(4)); ax.set_xticklabels(["Reduced, 1 t", "Reduced, 4 t", "Selective, 1 t", "Selective, 4 t"], fontsize=8)
    ax.set_title(case, fontsize=9)
axes[0].set_ylabel("median wall ratio arm / defaults-on (> 1: on is faster)")
axes[1].legend(fontsize=6.5, loc="upper right")
fig.suptitle("matrix F5 switches, 4-vCPU Xeon; run 1 = 3 rounds (open), run 2 = 15 rounds (filled)", fontsize=9)
fig.tight_layout(); fig.savefig(f"{out}/fig_f5.svg"); fig.savefig(f"{out}/fig_f5.png", dpi=150)

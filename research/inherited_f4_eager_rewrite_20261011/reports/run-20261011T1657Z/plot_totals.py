# Lays out the recorded arm totals from raw/*.jsonl; computes nothing new
# beyond the median of the three recorded rounds shown beside the points.
import glob, json, os, statistics
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
vals = {}
for path in sorted(glob.glob(os.path.join(here, "raw", "*.jsonl"))):
    name = os.path.basename(path)
    cell = "n=23 dim 8 (275 pts, 4 targets)" if name.startswith("n23") else "n=31 dim 9 (553 pts, 1 target)"
    mode = "lazy (KIC_F4_EAGER=0)" if "eager0" in name else "eager (default)"
    for line in open(path):
        if '"total_ms"' in line:
            d = json.loads(line)
            arm = "inherited F4" if d["arm"] == "f4" else "F6-IC + pair index"
            vals.setdefault((cell, arm, mode), []).append(d["total_ms"])

cells = sorted({k[0] for k in vals})
arms = ["inherited F4", "F6-IC + pair index"]
modes = ["lazy (KIC_F4_EAGER=0)", "eager (default)"]
colors = {modes[0]: "#8a8f98", modes[1]: "#2a6fdb"}
fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), sharey=False)
for ax, cell in zip(axes, cells):
    for i, arm in enumerate(arms):
        for j, mode in enumerate(modes):
            ys = vals[(cell, arm, mode)]
            x = i * 3 + j
            ax.bar(x, statistics.median(ys), color=colors[mode], width=0.8,
                   label=mode if i == 0 else None)
            ax.scatter([x] * len(ys), ys, color="black", s=10, zorder=3)
            ax.text(x, max(ys) + 25, f"{statistics.median(ys):.0f}", ha="center", fontsize=8)
    ax.set_xticks([0.5, 3.5])
    ax.set_xticklabels(arms)
    ax.set_title(cell, fontsize=10)
    ax.set_ylabel("probe total_ms (wall, ms)")
    ax.set_ylim(0, 1800)
axes[0].legend(fontsize=8, loc="upper right")
fig.suptitle("f6_ic_probe, same binary (6bdf93f4b), run-20261011T1657Z: bars = median of 3 rounds, dots = each round", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(here, "totals.svg"))

# Lays out summary.json (examples/pdp_shape_replay_summary.rs): each cell's
# default ÷ replay word operations, by ladder. No result is computed here.
import json, os
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
S = json.load(open(os.path.join(here, "..", "summary.json")))
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "plum": "#8e5bd1", "ink": "#0b0b0b",
     "ink2": "#52514e", "grid": "#e6e5e1", "surface": "#fcfcfb", "wash": "#9aa3ad"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C["grid"], "axes.labelcolor": C["ink2"],
                     "xtick.color": C["ink2"], "ytick.color": C["ink2"], "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": C["surface"], "axes.facecolor": C["surface"]})
cols = {"frozen (pilot)": C["wash"], "chain": C["blue"], "chain-holdout": C["orange"], "holdout": C["aqua"]}
fig, ax = plt.subplots(figsize=(8.4, 3.9))
x, labels, pos = 0, [], []
for L in S["ladders"]:
    for c in L["cells"]:
        v = c["default_over_replay"]
        ax.bar(x, v - 1, bottom=1, color=cols[L["ladder"]], width=0.7,
               label=L["ladder"] if c is L["cells"][0] else None)
        ax.text(x, max(v, 1) + 0.006, f"{v:.3f}", ha="center", fontsize=6.5, color=C["ink"])
        labels.append(f"{c['curve']} m={c['m']}" + (f" base {c['factor_index']}" if c["factor_index"] else "")
                      + (f" t{c['first_target']}" if c["first_target"] else ""))
        pos.append(x)
        x += 1
    x += 0.6
ax.set_xticks(pos, labels, rotation=55, ha="right", fontsize=6.5)
ax.axhline(1, color=C["ink2"], lw=0.8)
ax.set_ylim(0.7, 1.22)
ax.set_ylabel("default ÷ replay (word ops)\n> 1: replay cheaper")
ax.grid(axis="y", color=C["grid"], lw=0.6)
ax.legend(frameon=False, fontsize=7.5, ncol=4, loc="upper right")
ax.set_title("root shape replay in inherited F4: every registered cell costs more (verdicts identical)", fontsize=9)
fig.tight_layout()
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(os.path.join(here, f"ratio_by_cell.{ext}"), **kw)

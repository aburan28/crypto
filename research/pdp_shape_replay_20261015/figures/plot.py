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

# Lead 3: decompose_by_e3 against the pair table, per cell (raw/e3_vs_table.jsonl).
rows = [json.loads(l) for l in open(os.path.join(here, "..", "raw", "e3_vs_table.jsonl"))]
fig, ax = plt.subplots(figsize=(7.2, 3.6))
xs = range(len(rows))
ax.bar([x - 0.2 for x in xs], [r["e3_field_muls_per_target"] for r in rows], width=0.4, color=C["orange"], label="decompose_by_e3 (counted)")
ax.bar([x + 0.2 for x in xs], [r["table_muls_per_target"] for r in rows], width=0.4, color=C["blue"], label="pair table, 5|F| − 3")
for x, r in zip(xs, rows):
    ax.text(x - 0.2, r["e3_field_muls_per_target"] * 1.3, f"{r['e3_over_table']:.0f}×", ha="center", fontsize=7, color=C["ink"])
ax.set_yscale("log")
ax.set_ylim(10, 1e8)
ax.set_xticks(list(xs), [f"n={r['n']}, ℓ={r['l']}" for r in rows], fontsize=7.5, rotation=30, ha="right")
ax.set_ylabel("F_(2^n) multiplications per target")
ax.grid(axis="y", color=C["grid"], lw=0.6)
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.set_title("m = 3 by enumerating e₃ against the root-finding table (16 uniform targets a cell)", fontsize=9)
fig.tight_layout()
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(os.path.join(here, f"e3_vs_table.{ext}"), **kw)

"""Layout only: plots values recorded in data/summary_p*.json.

python research/smith_class_invariant_20261015/figures/plot.py
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
plt.rcParams.update({"font.size": 8, "svg.fonttype": "none"})

ps = [7, 11]
S = {p: json.loads((DATA / f"summary_p{p}.json").read_text()) for p in ps}
NAMES = [
    ("X1_cycle_type", "X1 cycle type"),
    ("X2_sylow2", "X2 2-Sylow"),
    ("X3_j4_profile", "X3 J[4] profile\n(frozen; invalid)"),
    ("X3plus_j4_sets", "X3+ J[4] sets\n(frozen; invalid)"),
    ("X3c_j4_profile_corrected", "X3c J[4] profile\n(corrected)"),
    ("X3cplus_j4_descent_sets_corrected", "X3c+ descent sets\n(corrected)"),
]

fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.4), sharey=False)
for ax, p in zip(axes, ps):
    s = S[p]
    c = s["C"]
    rows = [("unconditioned", s["null_unconditioned"], None)]
    for key, lab in NAMES:
        cd = s["candidates"][key]
        rows.append((lab, cd["null"], cd["placebo"]))
    for i, (lab, nl, pl) in enumerate(rows):
        y = len(rows) - 1 - i
        ax.plot([nl["lo"], nl["hi"]], [y, y], color="#3a6ea5", lw=3, solid_capstyle="butt")
        ax.plot(nl["median"], y, "o", color="#3a6ea5", ms=3.5)
        if pl is not None:
            ax.plot([pl["lo"], pl["hi"]], [y - 0.28, y - 0.28], color="#c58b2b", lw=1.6)
    ax.axvline(c, color="#b23b3b", lw=1, ls="--")
    ax.text(c, len(rows) - 0.45, f" observed C = {c:.4f}", color="#b23b3b", fontsize=7, va="bottom")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=6.5)
    ax.set_xlabel("class share C (999 permutations; bar = central 95%)")
    ax.set_title(f"p = {p}")
    ax.set_ylim(-0.7, len(rows) - 0.1)
axes[1].set_yticklabels([])
fig.text(
    0.01,
    0.005,
    "blue: null conditioned on (stratum, invariant), dot = median; orange: L-class placebo of equal cardinality. "
    "Source data/summary_p*.json.",
    fontsize=6,
)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(HERE / "class_invariant_nulls.svg")
fig.savefig(HERE / "class_invariant_nulls.png", dpi=150)

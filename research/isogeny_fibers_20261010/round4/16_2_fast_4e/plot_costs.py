"""Plot recorded values from sample_n131.json / sample_n233.json (no computation
beyond reading and arranging recorded numbers).  Writes costs.svg and
plot_data.json (the exact values plotted)."""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
METHODS = [
    ("binary_ld_rP", "[r]P, binary LD"),
    ("tnaf_ld_rP", "[r]P, tau-NAF LD"),
    ("closed_form_chain", "closed form, H by squarings"),
    ("closed_form_table", "closed form, H by byte table"),
]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, MUTED, GRID = "#1f1f1e", "#6b6a63", "#e4e3dd"

data = {}
for n in (131, 233):
    s = json.loads((HERE / f"sample_n{n}.json").read_text())
    ops = s["counted_ops_members_of_G"]
    bench = {b["set"]: {r["method"]: r for r in b["rows"]} for b in s["bench"]}
    data[n] = {
        "points_members_of_G": ops["points"],
        "mul": {m: ops[m]["mul_per_check"] for m, _ in METHODS},
        "sqr": {m: ops[m]["sqr_per_check"] for m, _ in METHODS},
        "wall_ns_median": {
            m: bench["members_of_G"][m]["contended_wall_ns_per_check_median"] for m, _ in METHODS
        },
        "bench_points": bench["members_of_G"]["closed_form_chain"]["points"],
        "bench_reps": bench["members_of_G"]["closed_form_chain"]["reps"],
    }
(HERE / "plot_data.json").write_text(json.dumps(data, indent=2))

plt.rcParams.update({"font.size": 9, "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED})
fig, axes = plt.subplots(1, 3, figsize=(11, 4.4))
panels = [
    ("mul", "field multiplications per check\n(counted, members of G)"),
    ("sqr", "field squarings per check\n(counted, members of G)"),
    ("wall_ns_median", "contended wall time per check, ns\n(median of reps, members of G)"),
]
width = 0.2
for ax, (key, title) in zip(axes, panels):
    for i, ((m, label), col) in enumerate(zip(METHODS, COLORS)):
        xs = [j + (i - 1.5) * (width + 0.02) for j in range(2)]
        vals = [max(data[n][key][m], 0.5) for n in (131, 233)]
        ax.bar(xs, vals, width, color=col, label=label, edgecolor="white", linewidth=1)
        for x, n in zip(xs, (131, 233)):
            v = data[n][key][m]
            txt = f"{v:.0f}" if v >= 1 else "0"
            ax.text(x, max(v, 0.5) * 1.15, txt, ha="center", va="bottom", fontsize=6.5, color=INK,
                    rotation=90)
    ax.set_yscale("log")
    ax.set_xticks([0, 1], ["n = 131\n(ECC2K-130 curve)", "n = 233\n(sect233k1)"])
    ax.set_title(title, fontsize=9)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_ylim(0.4, None)
    top = ax.get_ylim()[1]
    ax.set_ylim(0.4, top * 8)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=8, frameon=False)
fig.suptitle(
    "Lead 16.2: closed-form 4E test vs [r]P = O on y^2 + xy = x^3 + 1 "
    f"(counts: {data[131]['points_members_of_G']} / {data[233]['points_members_of_G']} members of G; "
    "a 0 bar is drawn at 0.5 on the log axis)",
    fontsize=9,
)
fig.text(0.5, 0.075, f"wall time: {data[131]['bench_points']} members of G x {data[131]['bench_reps']} reps per n, single thread, 4-core host shared with other agents (contended)", ha="center", fontsize=7, color=MUTED)
fig.tight_layout(rect=(0, 0.1, 1, 1))
fig.savefig(HERE / "costs.svg")
print("wrote costs.svg")

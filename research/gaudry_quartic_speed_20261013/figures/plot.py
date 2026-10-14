# Lays out the before/after records (raw/*.jsonl, written by
# gaudry_quartic_replay) and the S5 run (../gaudry_quartic_replay_20261013/
# raw/s5_p269.jsonl) as one figure: median wall time per phase of the
# full solve, before and after, at d = 4, 5, 6, and the S5 phases after.
import json, glob, os, statistics as st, collections
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
raw = os.path.join(here, "..", "raw")
s5 = os.path.join(here, "..", "..", "gaudry_quartic_replay_20261013", "raw", "s5_p269.jsonl")
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100", "plum": "#8e5bd1",
     "ink": "#0b0b0b", "ink2": "#52514e", "grid": "#e6e5e1", "surface": "#fcfcfb"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C["grid"], "axes.labelcolor": C["ink2"],
                     "xtick.color": C["ink2"], "ytick.color": C["ink2"],
                     "figure.facecolor": C["surface"], "axes.facecolor": C["surface"]})
PH = [("echelon_ms", "elimination", C["blue"]), ("nf_ms", "normal forms", C["orange"]),
      ("charpoly_ms", "char. polynomial", C["aqua"]), ("eigen_ms", "eigenvectors", C["plum"])]
rows = collections.defaultdict(list)
for f in glob.glob(os.path.join(raw, "*.jsonl")):
    arm, _, d = os.path.basename(f)[:-6].split("_")
    for l in open(f):
        rows[(arm, int(d[1:]))].append(json.loads(l)["baseline"])
s5rows = [json.loads(l)["baseline"] for l in open(s5)]
panels = []
for d in (4, 5, 6):
    panels.append((f"generic d = {d}", [(arm, {k: st.median(x[k] for x in rows[(arm, d)]) / 1e3 for k, _, _ in PH}) for arm in ("before", "after")]))
panels.append(("S5, d = 8 (3 systems)", [("before", None), ("after", {k: st.median(x[k] for x in s5rows) / 1e3 for k, _, _ in PH})]))
fig, axes = plt.subplots(1, 4, figsize=(8.0, 3.4))
for ax, (title, bars) in zip(axes, panels):
    bottom = [0.0, 0.0]
    for k, label, col in PH:
        v = [b[1][k] if b[1] else 0.0 for b in bars]
        ax.bar([0, 1], v, bottom=bottom, color=col, label=label, width=0.6)
        bottom = [a + b for a, b in zip(bottom, v)]
    for i, (arm, b) in enumerate(bars):
        t = bottom[i]
        ax.text(i, t * 1.02 if b else 0, (f"{t:.2f} s" if t < 10 else f"{t:.0f} s") if b else "not timed\n(~30 min)",
                ha="center", va="bottom", fontsize=7.5, color=C["ink"])
    ax.set_xticks([0, 1], [b[0] for b in bars])
    ax.set_title(title, fontsize=8.5)
    ax.set_ylim(0, max(bottom) * 1.18)
    ax.grid(axis="y", color=C["grid"], lw=0.6)
axes[0].set_ylabel("median wall of one full solve (s)")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, fontsize=8, ncol=4, loc="lower center")
fig.tight_layout(rect=(0, 0.07, 1, 1))
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(os.path.join(here, f"phase_wall.{ext}"), **kw)

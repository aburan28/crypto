# Lays out the recorded values in data.json as three SVG figures (and PNG
# renders).  Computes nothing new: axes are the recorded aggregates.
import json, os
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(here, "data.json")))
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100",
     "ink": "#0b0b0b", "ink2": "#52514e", "grid": "#e6e5e1", "surface": "#fcfcfb"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C["grid"], "axes.labelcolor": C["ink2"],
                     "xtick.color": C["ink2"], "ytick.color": C["ink2"], "axes.titlecolor": C["ink"],
                     "figure.facecolor": C["surface"], "axes.facecolor": C["surface"]})

def style(ax):
    ax.grid(True, color=C["grid"], linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(here, name + ".svg"))
    fig.savefig(os.path.join(here, name + ".png"), dpi=150)
    plt.close(fig)

STYLES = {"standard": (C["blue"], "o"), "random": (C["ink2"], "s"),
          "nb_window": (C["aqua"], "^"), "nb_break": (C["orange"], "v"),
          "subfield": (C["yellow"], "D")}

# ── Fig 1: D* vs factor-base dimension, n=11 and n=13 ──────────────────
fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0), sharey=True)
for ax, n in zip(axes, (11, 13)):
    for st in ("standard", "random", "nb_window", "nb_break"):
        pts = sorted([(c["dim"], c["d_star_mean"]) for c in D["m2_grid"]
                     if c["n"] == n and c["structure"] == st])
        if not pts: continue
        xs, ys = zip(*pts)
        col, mk = STYLES[st]
        ax.plot(xs, ys, marker=mk, color=col, label=st, ms=4)
    style(ax)
    ax.set_xlabel("factor-base subspace dimension n'")
    ax.set_title("n = %d" % n)
    ax.set_yticks([2, 3, 4])
axes[0].set_ylabel("D* (mean refutation degree)")
axes[0].legend(frameon=False, fontsize=8)
fig.suptitle("Degree of regularity rises ~1 degree per +2 dims; structure: closure beats windows",
             fontsize=9, color=C["ink"])
save(fig, "dstar_vs_dim")

# ── Fig 2: work-per-relation proxy vs dimension (log y) ────────────────
fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0))
for ax, n in zip(axes, (11, 13)):
    for st in ("standard", "random"):
        pts = sorted([(c["dim"], c["work_per_relation_proxy"]) for c in D["m2_grid"]
                     if c["n"] == n and c["structure"] == st
                     and c["work_per_relation_proxy"]])
        if not pts: continue
        xs, ys = zip(*pts)
        col, mk = STYLES[st]
        ax.plot(xs, ys, marker=mk, color=col, label=st, ms=4)
    style(ax)
    ax.set_xlabel("factor-base subspace dimension n'")
    ax.set_title("n = %d" % n)
    ax.set_yscale("log")
axes[0].set_ylabel("work / relation proxy (Macaulay rows x cols / frac)")
axes[0].legend(frameon=False, fontsize=8)
fig.suptitle("Work per relation has an interior optimum: enlarging the base past it is counterproductive",
             fontsize=9, color=C["ink"])
save(fig, "wpr_vs_dim")

# ── Fig 3: Shantz-Teske delta sweep ────────────────────────────────────
fig, ax = plt.subplots(figsize=(4.6, 3.0))
for key, (col, mk) in (("n11", (C["blue"], "o")), ("n13", (C["orange"], "s"))):
    pts = [(r["delta"], r["work_per_relation_proxy"]) for r in D["delta_sweeps_m2"][key]
           if r.get("work_per_relation_proxy")]
    xs, ys = zip(*pts)
    ax.plot(xs, ys, marker=mk, color=col, label=key, ms=4)
style(ax)
ax.set_xlabel("Shantz-Teske delta (dim = (n - delta)/m)")
ax.set_ylabel("work / relation proxy")
ax.set_yscale("log")
ax.legend(frameon=False, fontsize=8, title="standard base, m=2")
ax.set_title("Smaller bases via delta: 6.5x (n=11) and 11x (n=13) per-relation gains",
             fontsize=9, color=C["ink"])
save(fig, "delta_sweep")
print("figures written")

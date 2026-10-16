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

# ── Fig 1: D* vs factor-base dimension, n=11, 13, 17, 19 ──────────────
fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.6), sharey=True)
for ax, n in zip(axes.flat, (11, 13, 17, 19)):
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
for ax in axes[:, 0]:
    ax.set_ylabel("D* (mean refutation degree)")
axes[0, 0].legend(frameon=False, fontsize=8)
fig.suptitle("Degree of regularity rises ~1 degree per +2 dims; structure: closure beats windows",
             fontsize=9, color=C["ink"])
save(fig, "dstar_vs_dim")

# ── Fig 2: work-per-relation proxy vs dimension (log y) ────────────────
fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.6), sharey=False)
for ax, n in zip(axes.flat, (11, 13, 17, 19)):
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
for ax in axes[:, 0]:
    ax.set_ylabel("work / relation proxy")
axes[0, 0].legend(frameon=False, fontsize=8)
fig.suptitle("Work per relation has an interior optimum: enlarging the base past it is counterproductive",
             fontsize=9, color=C["ink"])
save(fig, "wpr_vs_dim")

# ── Fig 3: Shantz-Teske delta sweep, normalized to each n's delta=-1 ──
fig, ax = plt.subplots(figsize=(4.6, 3.0))
for key, (col, mk) in (("n11", (C["blue"], "o")), ("n13", (C["orange"], "s")),
                       ("n15", (C["aqua"], "^")), ("n17", (C["yellow"], "D")),
                       ("n19", (C["ink2"], "v"))):
    rows = [r for r in D["delta_sweeps_m2"][key] if r.get("work_per_relation_proxy")]
    if not rows: continue
    base = next(r["work_per_relation_proxy"] for r in rows if r["delta"] == min(x["delta"] for x in rows))
    xs = [r["delta"] for r in rows]
    ys = [r["work_per_relation_proxy"] / base for r in rows]
    ax.plot(xs, ys, marker=mk, color=col, label=key, ms=4)
style(ax)
ax.set_xlabel("Shantz-Teske delta (dim = (n - delta)/m)")
ax.set_ylabel("work / relation proxy")
ax.set_yscale("log")
ax.legend(frameon=False, fontsize=8, title="standard base, m=2")
ax.set_title("Per-relation gain from delta grows with n (each curve = 1 at its delta=-1)",
             fontsize=9, color=C["ink"])
save(fig, "delta_sweep")

# ── Fig 4: symmetric S4, the arity+size lever: D* vs (cell) ────────────
fig, ax = plt.subplots(figsize=(6.4, 3.0))
cells = [("n9 l2", "sym4_n9_l2"), ("n15 l2", "sym4_n15_l2"), ("n15 l3", "sym4_n15_l3"),
        ("n19 l3", "sym4_n19_l3"), ("n21 l3", "sym4_n21_l3")]
xpos = range(len(cells))
for i, (label, key) in enumerate(cells):
    for r in D["symmetric_s4"].get(key, []):
        col, mk = STYLES.get(r["structure"], (C["ink"], "o"))
        if r["d_star_mean"] is None: continue
        ax.plot([i], [r["d_star_mean"]], marker=mk, color=col, ms=6,
                label=r["structure"] if i == 0 else None)
style(ax)
ax.set_xticks(list(xpos))
ax.set_xticklabels([c[0] for c in cells])
ax.set_ylabel("D* (symmetric S4, m=3)")
ax.set_xlabel("field and base dim (m*l vs n gives delta; e.g. n=15 l=2 is delta=+9)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("Heavy delta on the symmetric S4 encoding drops D* by 2 degrees (5 -> 3)",
             fontsize=9, color=C["ink"])
save(fig, "sym4_delta")
print("figures written")
print("figures written")

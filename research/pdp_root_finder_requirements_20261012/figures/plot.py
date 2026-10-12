# Lays out the recorded values in data.json as two SVG figures (and PNG
# renders).  It computes nothing new: the only arithmetic is log2 of a
# recorded count divided by the recorded count of the smallest cell of the
# same line, so that word XORs and multiply-adds share one dimensionless axis.
import json, math, os
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(here, "data.json")))
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100",
     "ink": "#0b0b0b", "ink2": "#52514e", "grid": "#e6e5e1", "surface": "#fcfcfb", "wash": "#9aa3ad"}
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

# ---- Figure 2: measured solving degrees against the admissible degrees -----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
B = D["binary_f4_split_dreg"]
series_b = [("unsplit S4, m = 3", "unsplit_m3", C["blue"], "o"),
            ("S3 chain, m = 3", "chain_m3", C["orange"], "s"),
            ("two-way 2+2, m = 4", "two_way_m4", C["aqua"], "^"),
            ("three-way 2+1+1, m = 4", "three_way_m4", C["yellow"], "D")]
for label, key, col, mk in series_b:
    xs = [r[1] + (0.06 if key == "two_way_m4" else -0.06 if key == "three_way_m4" else 0) for r in B[key]]
    ys = [r[2] for r in B[key]]
    ax1.plot(xs, ys, linestyle="none", marker=mk, markersize=7, color=col, markeredgecolor=C["surface"],
             markeredgewidth=1.5, label=label)
adm = D["admissible_macaulay_degree"]
ax1.axhline(adm["binary_m4_n131"], color=C["ink2"], linestyle="--", linewidth=1.2)
ax1.text(2.72, adm["binary_m4_n131"] + 0.15, "admissible on ECC2K-130 at m = 4 (derived, README 6.4)",
         color=C["ink2"], fontsize=8)
ax1.set_xlim(2.6, 5.4); ax1.set_ylim(0, 10.2)
ax1.set_xticks([3, 4, 5]); ax1.set_xlabel("factor-base dimension ℓ")
ax1.set_ylabel("F4 solving degree (highest step degree at which the basis grew)")
ax1.set_title("Binary: Weil-descended Semaev systems, Koblitz curves n = 11…19 (split-dreg grid, g = 0)", fontsize=9)
ax1.legend(fontsize=7.5, loc="upper center", ncol=2, frameon=False)
style(ax1)

P = D["prime_f4_kummer_p1"]
for label, key, col, mk in [("m = 2", "m2", C["blue"], "o"), ("m = 3", "m3", C["orange"], "s"), ("m = 4", "m4", C["aqua"], "^")]:
    xs = [r[0] for r in P[key]]; ys = [r[1] for r in P[key]]
    ax2.plot(xs, ys, marker=mk, markersize=7, linewidth=2, color=col, markeredgecolor=C["surface"],
             markeredgewidth=1.5, label=label)
ax2.axhline(adm["prime_m4_256bit"], color=C["ink2"], linestyle="--", linewidth=1.2)
ax2.text(7.3, adm["prime_m4_256bit"] + 0.15, "admissible on a 256-bit prime curve at m = 4 (derived, README 6.4)",
         color=C["ink2"], fontsize=8)
ax2.set_xlim(7, 23); ax2.set_ylim(0, 10.2)
ax2.set_xticks([8, 10, 12, 14, 16, 18, 20, 22]); ax2.set_xlabel("tower variables N = m·t")
ax2.set_ylabel("F4 solving degree D (f4_fp_tower)")
ax2.set_title("Prime: Kummer-tower PKM systems at p₁, two targets per cell (PKM rounds 2–7)", fontsize=9)
ax2.legend(fontsize=8, loc="lower right", frameon=False)
style(ax2)
fig.suptitle("Measured F4 solving degrees on both tracks, against the degree the per-target budget admits at m = 4", fontsize=10)
save(fig, "solving_degree")

# ---- Figure 3: growth of F4's per-target work, indexed to the smallest cell -
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
# binary: unsplit m = 3 at n = 13 and n = 17, ℓ = 4 -> 5
rows = {(r[0], r[1]): r[3] for r in B["unsplit_m3"]}
for n, col in [(13, C["blue"]), (17, C["orange"])]:
    xs = [4, 5]; ys = [math.log2(rows[(n, l)] / rows[(n, 4)]) for l in xs]
    ax1.plot(xs, ys, marker="o", markersize=7, linewidth=2, color=col, markeredgecolor=C["surface"], markeredgewidth=1.5,
             label=f"F4 unsplit m = 3, n = {n}: {ys[1]:.1f} bits per rung")
rs = D["reference_slopes_bits_per_rung_of_ell_binary"]
ax1.plot([4, 5], [0, rs["enumeration_m3"]], color=C["wash"], linewidth=1.5, label="enumeration of pairs: 2 bits per rung (derived)")
ax1.plot([4, 5], [0, rs["pair_table_m3"]], color=C["wash"], linewidth=1.5, linestyle=":", label="pair table: 1 bit per rung (derived)")
ax1.set_xticks([4, 5]); ax1.set_xlim(3.8, 5.2); ax1.set_ylim(-0.5, 10.5)
ax1.set_xlabel("factor-base dimension ℓ (n fixed)")
ax1.set_ylabel("log₂(work per target ÷ work at ℓ = 4)")
ax1.set_title("Binary: word XORs per target, unsplit S4, g = 0 (split-dreg; two rungs: a ratio, not a fit)", fontsize=9)
ax1.legend(fontsize=7.5, loc="upper left", frameon=False)
style(ax1)

for label, key, col, mk in [("m = 2 (from N = 12)", "m2", C["blue"], "o"), ("m = 3 (from N = 9)", "m3", C["orange"], "s"), ("m = 4 (from N = 8)", "m4", C["aqua"], "^")]:
    pts = [(r[0], r[2]) for r in P[key] if r[2] is not None]
    n0, w0 = pts[0]
    xs = [p[0] - n0 for p in pts]; ys = [math.log2(p[1] / w0) for p in pts]
    slope = ys[-1] / xs[-1]
    ax2.plot(xs, ys, marker=mk, markersize=7, linewidth=2, color=col, markeredgecolor=C["surface"], markeredgewidth=1.5,
             label=f"F4 {label}: {slope:.1f} bits per tower variable")
sl = D["reference_slopes_bits_per_tower_variable"]
xx = [0, 10]
ax2.fill_between(xx, [0, 10 * sl["pair_table"]["m3"]], [0, 10 * sl["enumeration"]["m4"]], color=C["wash"], alpha=0.18, linewidth=0)
ax2.plot(xx, [0, 10 * sl["enumeration"]["m4"]], color=C["wash"], linewidth=1.5, label="enumeration, m = 4: 0.75 bits per variable (derived)")
ax2.plot(xx, [0, 10 * sl["pair_table"]["m3"]], color=C["wash"], linewidth=1.5, linestyle=":", label="pair table, m = 3: 0.33 bits per variable (derived)")
ax2.text(6.2, 2.4, "generic oracles:\n0.33–0.75 bits per variable", color=C["ink2"], fontsize=8)
ax2.set_xlim(-0.3, 10.3); ax2.set_ylim(-0.5, 21)
ax2.set_xlabel("tower variables added since the smallest cell, ΔN")
ax2.set_ylabel("log₂(multiply-adds ÷ multiply-adds at the smallest cell)")
ax2.set_title("Prime: F4 multiply-adds per target, Kummer towers at p₁ (PKM rounds)", fontsize=9)
ax2.legend(fontsize=7.5, loc="upper left", frameon=False)
style(ax2)
fig.suptitle("Growth of F4's per-target work against the growth of the generic oracles (both panels indexed to their smallest cell)", fontsize=10)
save(fig, "cost_growth")
print("ok")

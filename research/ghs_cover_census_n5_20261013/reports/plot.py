"""Render recorded values only. Every plotted number is a `group` or `mr` record written by
cover_census_n.rs (curve counts per minimal genus, structure counts per (m, r)); this script derives
nothing beyond dividing a recorded count by the recorded group size for bar heights.
Usage: python3 plot.py <run-results-dir> <outdir>"""
import json, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
def load(name):
    p = os.path.join(src, name + ".jsonl")
    return [json.loads(l) for l in open(p)] if os.path.exists(p) and os.path.getsize(p) > 0 else []

INK, MUTED, GRID = "#1a1a19", "#6b6a63", "#e4e3dc"
plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})
PAL = ["#1f6fd1", "#36a3a3", "#8bbf3a", "#e3a92b", "#e06a2e", "#b8456b", "#7a4fb5", "#555555"]

def stacked(ax, rows, groups, title, genera):
    col = {g: PAL[i % len(PAL)] for i, g in enumerate(genera)}
    for k, (key, label) in enumerate(groups):
        rec = next((r for r in rows if r["key"] == key), None)
        if rec is None:
            ax.text(k, 0.5, "none", ha="center", color=MUTED, fontsize=7); continue
        h = rec["curves_by_min_genus"]; n = rec["curves"]; y = 0.0
        for g in genera:
            c = h.get(str(g), 0)
            if c == 0: continue
            ax.bar(k, c / n, bottom=y, width=0.62, color=col[g], edgecolor="white", linewidth=0.8)
            if c / n > 0.07: ax.text(k, y + c / (2 * n), str(g), ha="center", va="center", color="white", fontsize=7)
            y += c / n
        ax.text(k, -0.24, f"n={n}", ha="center", va="top", color=MUTED, fontsize=7, transform=ax.get_xaxis_transform())
    ax.set_xticks(range(len(groups))); ax.set_xticklabels([g[1] for g in groups], fontsize=7.5)
    ax.set_ylim(0, 1.0); ax.set_title(title, color=INK, fontsize=8.5); ax.grid(axis="y", color=GRID, linewidth=0.7); ax.set_axisbelow(True)

# Figure 1: n = 3, minimum over the isogeny class (exhaustive) or within one 2-isogeny step (sample), by #E(K) mod 4
G4 = [("odd", "odd"), ("2mod4", "2 mod\n4"), ("0mod4", "0 mod\n4")]
n3 = [("R1-n3-q5", 5, "class"), ("R1-n3-q7", 7, "class"), ("A-n3-q11", 11, "class"), ("B-n3-q13", 13, "class"), ("C-n3-q17-sample200", 17, "2isog")]
gen3 = [3, 5, 7, 9, 11, 13, 17]
arms = [(a, q, s, [r for r in load(a) if r.get("rec") == "group" and r["by"] == "order_mod4" and r["scope"] == s]) for a, q, s in n3]
arms = [x for x in arms if x[3]]
if arms:
    fig, axes = plt.subplots(1, len(arms), figsize=(1.75 * len(arms) + 0.6, 3.0), sharey=True)
    axes = axes if len(arms) > 1 else [axes]
    for ax, (a, q, s, rows) in zip(axes, arms):
        stacked(ax, rows, G4, f"q={q}\n" + ("class min, all curves" if s == "class" else "≤1 step, sample"), gen3)
    axes[0].set_ylabel("fraction of curves")
    fig.legend(handles=[Patch(color=PAL[i], label=f"g={g}") for i, g in enumerate(gen3)], frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.07), ncol=len(gen3), fontsize=7.5)
    fig.text(0.5, -0.1, "n = 3: minimal genus by #E(K) mod 4 (q = 5, 7 reproduce PR #1759)", ha="center", color=INK)
    fig.savefig(f"{out}/n3_min_genus_by_mod4.svg", bbox_inches="tight"); fig.savefig(f"{out}/n3_min_genus_by_mod4.png", dpi=150, bbox_inches="tight")

# Figure 2: n = 5, minimal genus by (#E(K) mod 4, #E(K)[2])
G5 = [("odd/1", "odd\nE[2]=1"), ("2mod4/2", "2 mod 4\nE[2]=2"), ("0mod4/2", "0 mod 4\nE[2]=2"), ("0mod4/4", "0 mod 4\nE[2]=4")]
gen5 = [5, 17, 25, 33, 41, 49, 65, 97]
n5 = [("D-n5-q3-cosets", 3, "self", "q=3\nself, all curves"), ("D-n5-q3-cosets", 3, "class", "q=3\nclass min, all curves"),
      ("E-n5-q5-sample60", 5, "self", "q=5\nself, sample"), ("E-n5-q5-sample60", 5, "2isog", "q=5\n≤1 step, sample")]
arms = [(t, [r for r in load(a) if r.get("rec") == "group" and r["by"] == "mod4_E2" and r["scope"] == s]) for a, q, s, t in n5]
arms = [x for x in arms if x[1]]
if arms:
    fig, axes = plt.subplots(1, len(arms), figsize=(3.0 * len(arms) + 0.4, 3.1), sharey=True)
    axes = axes if len(arms) > 1 else [axes]
    for ax, (t, rows) in zip(axes, arms):
        stacked(ax, rows, G5, t, gen5)
    axes[0].set_ylabel("fraction of curves")
    fig.legend(handles=[Patch(color=PAL[i], label=f"g={g}") for i, g in enumerate(gen5)], frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.07), ncol=len(gen5), fontsize=7.5)
    fig.text(0.5, -0.2, "n = 5: minimal genus by (#E(K) mod 4, #E(K)[2]); derived floors 33 / 33 / 25 / 5", ha="center", color=INK)
    fig.savefig(f"{out}/n5_min_genus_by_2torsion.svg", bbox_inches="tight"); fig.savefig(f"{out}/n5_min_genus_by_2torsion.png", dpi=150, bbox_inches="tight")

# Figure 3: n = 5, structures by (r, m) summed over all censused models, against the derived allowed set
for arm, q in (("D-n5-q3-cosets", 3), ("E-n5-q5-sample60", 5)):
    mr = [r for r in load(arm) if r.get("rec") == "mr"]
    if not mr: continue
    fig, ax = plt.subplots(figsize=(6.4, 2.4))
    for e2, mk, c in ((1, "o", PAL[0]), (2, "s", PAL[4]), (4, "^", PAL[2])):
        pts = [r for r in mr if r["E2"] == e2]
        ax.scatter([r["r"] + {1: -0.18, 2: 0, 4: 0.18}[e2] for r in pts], [r["m"] for r in pts], s=[12 + 10 * len(str(r["structures"])) for r in pts], marker=mk, color=c, label=f"E(K)[2] has {e2} point(s)", alpha=0.85)
    ax.set_xlabel("r = |U| (branch points of the compositum)"); ax.set_ylabel("m"); ax.set_yticks([1, 2, 3, 4, 5]); ax.set_xticks(range(4, 21))
    ax.axhspan(1.5, 3.5, color="#f3d6d6", alpha=0.6, lw=0); ax.text(12, 2.5, "m = 2, 3 excluded (Lemma M)", ha="center", va="center", color="#9b3b3b", fontsize=7.5)
    ax.set_title(f"n = 5, q = {q}: every structure of every model by (r, m); marker area ~ log10(count)", fontsize=8.5, color=INK)
    ax.grid(color=GRID, linewidth=0.6); ax.set_axisbelow(True); ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.savefig(f"{out}/n5_m_r_q{q}.svg", bbox_inches="tight"); fig.savefig(f"{out}/n5_m_r_q{q}.png", dpi=150, bbox_inches="tight")

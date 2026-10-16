"""Lay out the recorded scores (summary.md, written by
`factor_base_shape_ladder score`) as one figure. Computes nothing: every
plotted value is a mean already in summary.md."""
import sys
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = sys.argv[1] if len(sys.argv) > 1 else "."
rows = []
for line in open(f"{here}/summary.md"):
    p = [c.strip() for c in line.strip().strip("|").split("|")]
    if len(p) == 10 and p[0] in ("G1", "G2", "G3", "G5") and p[4] != "-":
        rows.append(dict(goal=p[0], cell=p[1], arm=p[2], n=int(p[3]), mean=float(p[4]),
                         cens=int(p[5]), degs=[int(x) for x in p[7].split()]))

fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
titles = {"G1": "G1 shape at m = 3 (class-0)", "G3": "G3 modulus panel, coord window",
          "G2": "G2 twisted chains", "G5": "G5 Möbius pole, m = 2"}
colors = ["#2f6db3", "#d1495b", "#3a9d5d", "#e0a02a", "#7a5cb8", "#555555", "#2bb5b5", "#b35c2f"]
for ax, goal in zip(axes.flat, ["G1", "G3", "G2", "G5"]):
    rs = [r for r in rows if r["goal"] == goal]
    cells = sorted({r["cell"] for r in rs}, key=lambda c: (len(c), c))
    arms = []
    for r in rs:
        a = r["arm"] if goal != "G3" else f"modulus {len([x for x in arms if x.startswith('modulus')])}"
        if goal != "G3" and a not in arms:
            arms.append(a)
    if goal == "G3":
        for ci, c in enumerate(cells):
            cr = [r for r in rs if r["cell"] == c]
            for k, r in enumerate(cr):
                x = ci + (k - (len(cr) - 1) / 2) * 0.08
                ax.plot(x, r["mean"], "o", color=colors[ci % len(colors)], ms=5)
        ax.set_xticks(range(len(cells)), cells)
        ax.text(0.02, 0.02, "one dot per modulus", transform=ax.transAxes, fontsize=8, color="#444")
    else:
        if goal == "G5":
            arms = [a for a in arms if not a.startswith("beta_rand") or a == "beta_rand"]
        if goal == "G2":
            arms = sorted(arms, key=lambda a: (a.split("/")[0], a.split("/")[1] != "plain", a))
        w = 0.8 / max(1, len(arms))
        for ai, a in enumerate(arms):
            xs, ys, cs = [], [], []
            for ci, c in enumerate(cells):
                r = next((r for r in rs if r["cell"] == c and r["arm"] == a), None)
                if r:
                    xs.append(ci - 0.4 + w * (ai + 0.5))
                    ys.append(r["mean"])
                    cs.append(r["cens"])
            ax.plot(xs, ys, "o", color=colors[ai % len(colors)], label=a, ms=5,
                    mfc=colors[ai % len(colors)])
            for x, y, c in zip(xs, ys, cs):
                if c:
                    ax.annotate(f"{c} cens.", (x, y), fontsize=6, xytext=(2, 3), textcoords="offset points")
        ax.set_xticks(range(len(cells)), cells)
        ax.legend(fontsize=6, ncol=2, loc="center" if goal in ("G1", "G2") else "lower right")
    ax.set_title(titles[goal], fontsize=10)
    ax.set_ylabel("mean refutation degree (censored at bound)")
    ax.grid(axis="y", alpha=0.3)
fig.suptitle("Factor-base shape search G1/G2/G3/G5 — stage diagnostic, toy fields, no S", fontsize=11)
fig.tight_layout()
fig.savefig(f"{here}/shape_search_degrees.svg")
fig.savefig(f"{here}/shape_search_degrees.png", dpi=110)

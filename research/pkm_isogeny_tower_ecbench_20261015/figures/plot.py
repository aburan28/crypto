"""Lays out recorded values of research/pkm_isogeny_tower_ecbench_20261015
(tables/summary-*.tsv and tables/panel3-summary.tsv, written by
summarise.sh from the sealed sessions and panel 3).  No computation beyond
reading the tables.  Writes SVG and a 150 dpi PNG of each figure."""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
TAB = HERE.parent / "tables"
BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
INK, MUTED, GRID = "#1f2328", "#5a6270", "#d7dbe0"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "svg.fonttype": "none"})
SESS = ["p20", "p22", "p24", "p26", "p28", "q20", "q22", "q24"]


def rows(name):
    with open(TAB / name) as f:
        return list(csv.DictReader(f, delimiter="\t"))


ratios = rows("summary-ratios.tsv")
arms = rows("summary-arms.tsv")


def ratio(s, a, b):
    for r in ratios:
        if r["session"] == s and r["A"] == a and r["B"] == b:
            return float(r["ratio_B_over_A"]), float(r["lo95"]), float(r["hi95"])
    return None


def arm(s, a):
    for r in arms:
        if r["session"] == s and r["arm"] == a:
            return r
    return None


def save(fig, name):
    fig.savefig(HERE / f"{name}.svg")
    fig.savefig(HERE / f"{name}.png", dpi=150)
    plt.close(fig)


# 1. isogeny / plain and isogeny / kummer, with the aux share.
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.2, 5.6), sharex=True,
                              gridspec_kw={"height_ratios": [3, 2]})
x = range(len(SESS))
for off, (a, col, lab) in zip((-0.12, 0.12), [("ic-plain", BLUE, "isogeny ÷ plain base (equal size)"),
                                             ("ic-kummer-rootfind", ORANGE, "isogeny ÷ Kummer tower (same t)")]):
    xs, ys, lo, hi = [], [], [], []
    for i, s in enumerate(SESS):
        r = ratio(s, a, "ic-iso-rootfind")
        if r:
            xs.append(i + off); ys.append(r[0]); lo.append(r[0] - r[1]); hi.append(r[2] - r[0])
    ax.errorbar(xs, ys, yerr=[lo, hi], fmt="o", ms=6, color=col, capsize=3, lw=1.5, label=lab)
ax.axhline(1.0, color=MUTED, lw=1, ls="--")
ax.set_ylabel("paired ratio of S (95% cluster bootstrap)")
ax.grid(axis="y", color=GRID, lw=0.6)
ax.legend(frameon=False, loc="upper left")
ax.set_title("IC + prime-isogeny (rootfind-mitm:m=3), charged search, against equal-size bases", color=INK, fontsize=10)
share = [100 * float(arm(s, "ic-iso-rootfind")["mean_fb_share"]) for s in SESS]
tried = [arm(s, "ic-iso-rootfind")["aux_curves_tried"] for s in SESS]
ax2.bar(list(x), share, width=0.5, color=AQUA)
for i, (v, t) in enumerate(zip(share, tried)):
    ax2.text(i, v + 0.5, f"{v:.1f}%\n{t} curves", ha="center", va="bottom", fontsize=7.5, color=INK)
ax2.set_ylim(0, max(share) * 1.45)
ax2.set_ylabel("aux search share of S (%)")
ax2.grid(axis="y", color=GRID, lw=0.6)
ax2.set_xticks(list(x))
ax2.set_xticklabels([f"{s}\n{'v₂=12' if s[0]=='p' else 'v₂=1'}" for s in SESS])
ax2.set_xlabel("curve (bits of p; v₂(p−1)); 8 targets × 2 rounds per arm, L0, counts only")
fig.tight_layout()
save(fig, "ratios")

# 2. mean S per arm per curve.
fig, ax = plt.subplots(figsize=(7.2, 3.8))
ARMS = [("ic-iso-rootfind", BLUE, "IC + isogeny, rootfind-mitm"), ("ic-iso-mitm", AQUA, "IC + isogeny, mitm"),
        ("ic-kummer-rootfind", ORANGE, "IC + Kummer, rootfind-mitm"), ("ic-plain", YELLOW, "IC + plain base, rootfind-mitm"),
        ("rho-neg", MAGENTA, "rho, negation map")]
w = 0.16
for k, (a, col, lab) in enumerate(ARMS):
    xs, ys = [], []
    for i, s in enumerate(SESS):
        r = arm(s, a)
        if r:
            xs.append(i + (k - 2) * w); ys.append(float(r["mean_S"]))
    ax.bar(xs, ys, width=w * 0.9, color=col, label=lab)
ax.set_yscale("log")
ax.set_ylabel("mean S (GAE / √r), log scale")
ax.set_xticks(list(x)); ax.set_xticklabels(SESS)
ax.grid(axis="y", color=GRID, lw=0.6)
ax.legend(frameon=False, fontsize=7.5, ncol=3, loc="upper left")
ax.set_ylim(0.5, 2000)
ax.set_title("Mean S per arm (16 measured runs each; every run verified)", color=INK, fontsize=10)
fig.tight_layout()
save(fig, "s_by_arm")

# 3. panel 3: auxiliary curves tried against t.
p3 = rows("panel3-summary.tsv")
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.4))
for c, col, lab in [("p24", BLUE, "icv1-fp24-t2705-27cbb01b (v₂ = 12)"), ("q24", ORANGE, "icv1-fp24-t4451-d5368b2d (v₂ = 1)")]:
    rs = [r for r in p3 if r["curve"] == c]
    t = [int(r["t"]) + (-0.1 if c == "p24" else 0.1) for r in rs]
    m = [float(r["mean_tried"]) for r in rs]
    lo = [m_ - float(r["min_tried"]) for m_, r in zip(m, rs)]
    hi = [float(r["max_tried"]) - m_ for m_, r in zip(m, rs)]
    ax.errorbar(t, m, yerr=[lo, hi], fmt="o", ms=5, color=col, capsize=2, lw=1.2, label=lab)
    ax2.plot([int(r["t"]) for r in rs], [float(r["gae_per_curve_tried"]) for r in rs], "o-", ms=5, lw=2, color=col)
tt = list(range(4, 15))
ax.plot(tt, [2 ** v for v in tt], color=MUTED, ls="--", lw=1, label="2^t")
ax.set_yscale("log"); ax.set_xlabel("t"); ax.set_ylabel("aux curves tried (mean, min–max of 16 seeds)")
ax.grid(color=GRID, lw=0.6); ax.legend(frameon=False, fontsize=6.5, loc="upper left")
ax2.set_xlabel("t"); ax2.set_ylabel("charged GAE per curve tried")
ax2.grid(color=GRID, lw=0.6)
fig.suptitle("Panel 3: the auxiliary-curve search (ecbench fb, 16 seeds per t)", color=INK, fontsize=10)
fig.tight_layout()
save(fig, "panel3")

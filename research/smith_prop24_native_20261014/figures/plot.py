"""Lay out recorded values only (no computation of results).

Reads research/smith_prop24_native_20261014/data/*.json written by
`examples/smith_prop24.rs` and draws:

  figures/smith_prop24_results.svg (+ .png at 150 dpi)
      (a) Question 1: observed minus predicted fraction, with the recorded
          98.75% Clopper-Pearson interval, against Smith's Example 2;
      (b) Question 2: direct D, class-level C and ceiling U per prime.
  figures/class_share.svg (+ .png)
      post hoc, descriptive: the recorded distribution of the share of
      full_success curves among the tractable curves of each class.

Usage: python figures/plot.py  (from the study directory)
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#d9d8d4", "#ffffff"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
    }
)


def load(name):
    p = DATA / name
    return json.loads(p.read_text()) if p.exists() else None


def results_figure():
    q1 = load("random_p2147483647_summary.json")
    census = [s for s in (load(f"census_p{p}_summary.json") for p in (7, 11)) if s]
    fig, (a, b) = plt.subplots(
        1, 2, figsize=(7.4, 3.5), gridspec_kw={"width_ratios": [1.15, 1]}
    )
    fig.patch.set_facecolor(SURFACE)

    # (a) Question 1, in percentage points from the prediction.
    short = [
        "rational tractable\nsubgroup (curves)",
        "Prop. 2 passes\n(subgroups)",
        "Prop. 4 passes, given\nProp. 2 (subgroups)",
        "full_success\n(curves)",
    ]
    fr = q1["fractions"]
    ys = list(range(len(fr)))[::-1]
    margin = fr[0]["margin"] * 100
    a.axvspan(-margin, margin, color=GRID, alpha=0.45, lw=0, zorder=0)
    a.axvline(0, color=INK2, lw=0.8, zorder=1)
    for y, f in zip(ys, fr):
        pred = f["prediction"]
        v = (f["value"] - pred) * 100
        lo, hi = (f["cp_lo"] - pred) * 100, (f["cp_hi"] - pred) * 100
        a.errorbar(
            v, y, xerr=[[v - lo], [hi - v]], fmt="o", ms=5, color=BLUE,
            ecolor=BLUE, elinewidth=2, capsize=0, zorder=3,
            label="this run, N = 10⁶ (98.75% CP)" if y == ys[0] else None,
        )
        a.plot(
            (f["smith_example2"] - pred) * 100, y - 0.22, marker="D", ms=5,
            mfc=SURFACE, mec=ORANGE, mew=1.6, ls="none", zorder=3,
            label="Smith Example 2, N = 10⁶" if y == ys[0] else None,
        )
    a.set_yticks(ys, short)
    a.set_xlim(-1.25, 1.25)
    a.set_ylim(-0.7, len(fr) - 0.4)
    a.set_xlabel("observed − predicted (percentage points)")
    a.text(margin * 0.97, -0.55, "±1 pp margin", ha="right", va="center",
           color=INK2, fontsize=7)
    a.set_title("(a) Q1: p = 2³¹ − 1, random curves", loc="left", fontsize=9)
    h, lab = a.get_legend_handles_labels()
    a.legend(h[::-1], lab[::-1], loc="upper center", frameon=False, fontsize=7,
             bbox_to_anchor=(0.45, -0.24), ncol=1)
    a.grid(axis="x", color=GRID, lw=0.5)

    # (b) Question 2: D, C, U per prime.
    asym_d = census[0]["asymptotic_direct_prediction"] if census else 0.1857
    asym_u = census[0]["asymptotic_ceiling"] if census else 0.5016
    b.axvline(asym_d, color=INK2, lw=0.9, ls=":", zorder=1)
    b.axvline(asym_u, color=INK2, lw=0.9, ls="--", zorder=1)
    b.text(asym_d + 0.005, len(census) - 0.35, "Thm 2: 0.1857", ha="left",
           va="bottom", color=INK2, fontsize=7)
    b.text(asym_u - 0.005, -0.62, "S₈ ceiling 0.5016", ha="right",
           va="center", color=INK2, fontsize=7)
    for i, s in enumerate(census):
        y = len(census) - 1 - i
        b.plot([s["direct_full"], s["class_tractable"]], [y, y], color=GRID, lw=2,
               zorder=1, solid_capstyle="round")
        b.plot(s["direct_full"], y, "o", ms=7, color=BLUE, mec=SURFACE, mew=1.5,
               zorder=3, label="D: direct full_success" if i == 0 else None)
        b.plot(s["class_full"], y, "o", ms=7, color=ORANGE, mec=SURFACE, mew=1.5,
               zorder=3, label="C: class has a full_success member" if i == 0 else None)
        b.plot(s["class_tractable"], y, "o", ms=7, mfc=SURFACE, mec=INK, mew=1.4,
               zorder=3, label="U: class has a tractable member" if i == 0 else None)
        # Label placed mid-segment (layout only).
        b.text((s["direct_full"] + s["class_full"]) / 2, y - 0.12, f"r = {s['r']:.3f}",
               ha="center", va="top", color=INK, fontsize=7.5)
    b.set_yticks(range(len(census)), [f"p = {s['p']}" for s in census][::-1])
    b.set_ylim(-0.75, len(census) - 0.1)
    b.set_xlim(0.1, 0.56)
    b.set_xlabel("fraction of curves (curve + twist)")
    b.set_title("(b) Q2: per isogeny class", loc="left", fontsize=9)
    b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), frameon=False,
             fontsize=7, ncol=1)
    b.grid(axis="x", color=GRID, lw=0.5)
    fig.subplots_adjust(left=0.19, right=0.985, top=0.91, bottom=0.36, wspace=0.42)
    for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
        fig.savefig(HERE / f"smith_prop24_results.{ext}", facecolor=SURFACE, **kw)
    plt.close(fig)


def share_figure():
    desc = load("describe_classes.json")
    if not desc:
        return
    fig, ax = plt.subplots(figsize=(7.4, 2.5))
    fig.patch.set_facecolor(SURFACE)
    bins = desc[0]["bins"]
    width = 0.38
    for k, (d, color) in enumerate(zip(desc, (BLUE, ORANGE))):
        tot = sum(d["curves_by_bin"])
        xs = [i + (k - 0.5) * (width + 0.02) for i in range(len(bins))]
        ax.bar(xs, [c / tot for c in d["curves_by_bin"]], width=width, color=color,
               label=f"p = {d['p']} ({d['tractable_classes']} classes)", zorder=2)
    ax.set_xticks(range(len(bins)), bins, fontsize=7)
    ax.set_xlabel("share of a class's tractable curves that have full_success")
    ax.set_ylabel("share of curves in\ntractable classes")
    ax.set_title("Post hoc, descriptive: spread of full_success within classes",
                 loc="left", fontsize=9)
    ax.legend(frameon=False, fontsize=7.5)
    ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
    fig.tight_layout()
    for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
        fig.savefig(HERE / f"class_share.{ext}", facecolor=SURFACE, **kw)
    plt.close(fig)


if __name__ == "__main__":
    results_figure()
    share_figure()

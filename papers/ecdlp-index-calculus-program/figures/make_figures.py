"""Render the figures of the index-calculus programme paper from data.json.

Plots recorded values only (AGENTS.md: Python lays out what was recorded; it
computes no result).  Each figure is written as PDF (for LaTeX), SVG and a
150 dpi PNG from the same source.

    python3 make_figures.py            # render every figure
"""
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "data.json")))

plt.rcParams.update({
    "font.size": 9,
    "font.family": "serif",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "svg.fonttype": "none",
})

# One categorical palette for the whole paper (colour-blind-safe Okabe-Ito).
C = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#999999"]
INK = "#222222"


def save(fig, name):
    for ext in ("pdf", "svg"):
        fig.savefig(os.path.join(HERE, f"{name}.{ext}"), bbox_inches="tight")
    fig.savefig(os.path.join(HERE, f"{name}.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_timeline():
    w = D["weekly_commits"]
    weeks = w["weeks"]
    x = range(len(weeks))
    fig, ax = plt.subplots(figsize=(6.4, 2.8))
    ax.bar(x, w["total"], color="#dddddd", label="all commits", zorder=1)
    for i, (name, vals) in enumerate(w["themes"].items()):
        ax.plot(x, vals, marker="o", ms=3, lw=1.4, color=C[i], label=name, zorder=3)
    ax.set_xticks(list(x))
    ax.set_xticklabels(weeks, rotation=0, fontsize=7)
    ax.set_ylabel("commits per ISO week")
    ax.set_xlabel("2026 ISO week (W27 = 29 Jun; W37 = 7 Sep; W42 = 12 Oct, partial)")
    ax.legend(fontsize=6.5, ncol=2, frameon=False, loc="upper left")
    ax.set_title("Commit activity on main, 2026-07-01 to 2026-10-14 (keyword-matched themes)",
                 fontsize=8.5)
    save(fig, "fig-timeline")


def fig_rootfinding_oracle():
    rows = D["rootfinding_oracle"]["rows"]
    labels = [r["label"] for r in rows][::-1]
    vals = [r["ratio"] for r in rows][::-1]
    cols = [C[0] if r["kind"] == "per decomposition" else C[1] for r in rows][::-1]
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    ax.barh(range(len(vals)), vals, color=cols)
    ax.set_xscale("log")
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels(labels, fontsize=7)
    for i, v in enumerate(vals):
        ax.text(v * 1.08, i, f"{v:,}×", va="center", fontsize=6.5, color=INK)
    ax.set_xlabel("Gröbner baseline time / root-finding time (log scale)")
    ax.set_xlim(10, 2e5)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=C[0], label="oracle stage, per decomposition"),
                       Patch(color=C[1], label="complete DLP, CPU time")],
              fontsize=7, frameon=False, loc="lower right")
    ax.set_title("Root-finding oracle against Gröbner solvers (stage vs end-to-end)",
                 fontsize=8.5)
    save(fig, "fig-rootfind-oracle")


def fig_k13():
    d = D["k13_rootfind_progression"]
    rows = d["rows"]
    labels = [r["step"] for r in rows]
    vals = [r["S"] for r in rows]
    fig, ax = plt.subplots(figsize=(6.4, 2.8))
    ax.bar(range(len(vals)), vals, color=C[2])
    for i, v in enumerate(vals):
        ax.text(i, v + 0.6, f"{v:.2f}", ha="center", fontsize=7)
    refs = d["references"]
    styles = [(C[0], "-"), (C[5], "--"), (C[6], ":")]
    for (col, ls), r in zip(styles, refs):
        if r["S"] < 60:
            ax.axhline(r["S"], color=col, ls=ls, lw=1.2, label=f"{r['label']}: {r['S']:.2f}")
    f4 = [r for r in refs if r["S"] > 60][0]
    ax.plot([], [], color=C[5], ls="--", lw=1.2,
            label=f"{f4['label']}: {f4['S']:.1f} (off scale)")
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, rotation=18, ha="right", fontsize=6.5)
    ax.set_ylim(0, 52)
    ax.set_ylabel("mean cold S (complete count)")
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    ax.set_title("Koblitz E₀, n = 13: root-finding IC arms across four ecbench loops "
                 "(8 targets × 3 rounds each, L0)", fontsize=8)
    save(fig, "fig-k13-progression")


def fig_cold_gap():
    rows = sorted(D["cold_gap"]["rows"], key=lambda r: r["ratio"])
    fig, ax = plt.subplots(figsize=(6.4, 2.9))
    for i, r in enumerate(rows):
        ax.barh(i, r["ratio"], color=C[0] if not r["lower_bound"] else "white",
                edgecolor=C[0], hatch="////" if r["lower_bound"] else None)
        ax.text(r["ratio"] * 1.07, i, f"{r['ratio']:g}×", va="center", fontsize=6.5)
    ax.axvline(1, color=C[5], lw=1.2)
    ax.text(1.05, -0.9, "parity with matched rho", color=C[5], fontsize=6.5)
    ax.set_xscale("log")
    ax.set_xlim(0.8, 400)
    ax.set_ylim(-1.3, len(rows) - 0.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r["label"] for r in rows], fontsize=7)
    ax.set_xlabel("IC cold S / matched rho cold S (log scale; hatched = lower bound on IC cost)")
    ax.set_title("Cold one-target frontier: every measured IC arm is above rho", fontsize=8.5)
    save(fig, "fig-cold-gap")


def fig_online():
    rows = D["koblitz_online"]["rows"]
    fig, ax = plt.subplots(figsize=(6.4, 2.7))
    for r in rows:
        adm = r["status"].startswith("admitted")
        x = r["n"] + (0.8 if "dense" in r["status"] else 0)
        if adm:
            ax.errorbar(x, r["ratio"], yerr=[[r["ratio"] - r["lo"]], [r["hi"] - r["ratio"]]],
                        fmt="s", color=C[2], ms=6, capsize=3, label="admitted (isolated, 64 targets, 95% CI)")
        else:
            ax.plot(x, r["ratio"], "o", mfc="white", mec=C[1], ms=6,
                    label="exploratory (shared host / no isolation receipt)")
        ax.text(x + 0.6, r["ratio"], f"{r['ratio']:g}×", fontsize=6.5, va="center")
    cold = 1 / D["koblitz_online"]["cold_ratio_n61_ic_over_rho"]
    ax.plot(61, cold, "v", color=C[5], ms=6, label="n=61 cold (setup charged): 1/83.4")
    ax.axhline(1, color="#888", lw=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("Koblitz field degree n (GF(2^n))")
    ax.set_ylabel("rho online / IC online (same point)")
    h, l = ax.get_legend_handles_labels()
    uniq = dict(zip(l, h))
    ax.legend(uniq.values(), uniq.keys(), fontsize=6.5, frameon=False, loc="lower left")
    ax.set_xlim(58, 90)
    ax.set_title("Retained-setup online boundary on Koblitz curves (setup excluded unless marked)",
                 fontsize=8.5)
    save(fig, "fig-koblitz-online")


def fig_prime():
    g = D["gaudry_quartic"]
    p = D["pkm"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.4, 2.6), gridspec_kw={"width_ratios": [1, 1.4]})
    a1.bar([0, 1], [g["C4_S4"], g["C4_D4"]], color=[C[6], C[3]])
    a1.set_yscale("log")
    a1.set_xticks([0, 1])
    a1.set_xticklabels(["Gaudry S₄\n(x-line)", "FGHR D₄\n(Y-line)"], fontsize=7)
    a1.set_ylabel("C₄: F_p mults per k=4 solve")
    for i, v in enumerate([g["C4_S4"], g["C4_D4"]]):
        a1.text(i, v * 1.3, f"{v:.2e}", ha="center", fontsize=6.5)
    a1.set_ylim(1e9, 1e13)
    a1.set_title(r"$\mathbb{F}_{p^4}$, p = 269: " + "47.5\u201348.4\u00d7", fontsize=8)
    lo = [r - l for r, l in zip(p["ratio"], p["lo"])]
    hi = [h - r for r, h in zip(p["ratio"], p["hi"])]
    a2.errorbar(p["bits"], p["ratio"], yerr=[lo, hi], fmt="o-", color=C[0], capsize=3, ms=4)
    a2.set_yscale("log")
    a2.axhline(1, color=C[5], lw=0.8)
    a2.set_ylim(0.8, 300)
    a2.set_xlabel("prime curve size (bits)")
    a2.set_ylabel("IC S / rho-negation S")
    a2.set_title("PKM Kummer-tower base + root finding (95% CI)", fontsize=8)
    save(fig, "fig-prime")


def fig_engine():
    rows = D["engine"]["rows"]
    fig, ax = plt.subplots(figsize=(6.4, 2.5))
    vals = [r["ratio"] for r in rows][::-1]
    ax.barh(range(len(vals)), vals, color=C[4])
    ax.set_xscale("log")
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels([r["label"] for r in rows][::-1], fontsize=7)
    for i, v in enumerate(vals):
        ax.text(v * 1.08, i, f"{v:g}×", va="center", fontsize=6.5)
    ax.set_xlim(1, 600)
    ax.set_xlabel("before / after (time or counted work; log scale)")
    ax.set_title("Gröbner-engine engineering, selected measured ratios", fontsize=8.5)
    save(fig, "fig-engine")


if __name__ == "__main__":
    fig_timeline()
    fig_rootfinding_oracle()
    fig_k13()
    fig_cold_gap()
    fig_online()
    fig_prime()
    fig_engine()
    print("ok")

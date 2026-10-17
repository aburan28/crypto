"""Lay out the recorded values of the greedy factor-base yield run.

Plotting only: every number drawn is a value in data/rows.json or
data/summary.json, written by examples/koblitz_greedy_factor_base_yield.rs.
The one reduction is the expected-trial formula (U + 1) * T / covered,
which the census already records per base; it is read, not recomputed,
except for the derived optimum noted in the README.

    python3 graphs/plot_yield.py          # from this study directory
"""

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "legend.fontsize": 8.5})

HERE = pathlib.Path(__file__).resolve().parent
STUDY = HERE.parent
RUN = "run 2026-10-14, seeds 1-20, exhaustive subgroup targets, m = 2"
PANELS = [(0, 13), (0, 15), (1, 17)]
# Column range shown per panel: up to where coverage saturates.
XMAX = {(0, 13): 12, (0, 15): 8, (1, 17): 34}
COLORS = {"greedy": "#1f5fa8", "pool": "#6fa8dc", "random": "#444444", "prod": "#b45f06"}

rows = json.loads((STUDY / "data/rows.json").read_text())
summ = json.loads((STUDY / "data/summary.json").read_text())


def title(a, n):
    return f"$K_{a}/2^{{{n}}}$"


def series(a, n, arm):
    return sorted(
        [r for r in rows if (r["a"], r["n"]) == (a, n) and r["arm"] == arm],
        key=lambda r: (r["columns"], r["seed"]),
    )


def mean_by_u(a, n, arm, field):
    out = {}
    for r in series(a, n, arm):
        out.setdefault(r["columns"], []).append(r[field])
    return {u: sum(v) / len(v) for u, v in out.items()}


def sd_by_u(a, n, arm, field):
    import statistics

    out = {}
    for r in series(a, n, arm):
        out.setdefault(r["columns"], []).append(r[field])
    return {u: (statistics.stdev(v) if len(v) > 1 else 0.0) for u, v in out.items()}


def save(fig, stem):
    fig.savefig(HERE / f"{stem}.svg", format="svg")
    fig.savefig(HERE / f"{stem}.png", format="png", dpi=150)
    plt.close(fig)


def fig_coverage():
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 3.7), sharey=True)
    for ax, (a, n) in zip(axes, PANELS):
        s = [x for x in summ if (x["a"], x["n"]) == (a, n)]
        us = [x["columns"] for x in s]
        g = [x["greedy_coverage"] for x in s]
        rm = [x["random_mean_coverage"] for x in s]
        rs = [x["random_sd_coverage"] for x in s]
        pm = [x["pool_mean_coverage"] for x in s]
        ax.plot(us, g, "-o", ms=3.5, color=COLORS["greedy"], label="greedy (all orbits scored)")
        ax.plot(us, pm, "--", color=COLORS["pool"], label="greedy, 64 sampled per step")
        ax.plot(us, rm, "-", color=COLORS["random"], label="random orbits (mean of 20)")
        lo = [m - d for m, d in zip(rm, rs)]
        hi = [m + d for m, d in zip(rm, rs)]
        ax.fill_between(us, lo, hi, color=COLORS["random"], alpha=0.15, linewidth=0)
        prod = series(a, n, "production")
        if prod:
            ax.scatter(
                [r["columns"] for r in prod],
                [r["coverage"] for r in prod],
                marker="x",
                s=28,
                color=COLORS["prod"],
                label="shipped sampler (its own column counts)",
                zorder=5,
            )
        ax.set_title(f"{title(a, n)}, $r={s[0]['subgroup_order']}$, $T={s[0]['targets']}$")
        ax.set_xlabel("signed Frobenius orbits (columns U)")
        ax.set_xlim(0, XMAX[(a, n)])
        ax.set_ylim(0, 1.02)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("coverage: fraction of targets as a 2-sum")
    axes[1].legend(loc="lower right")
    fig.suptitle(
        "Coverage at equal unknowns. Bands: ±1 sd over the 20 random seeds. " + RUN,
        fontsize=9,
    )
    fig.tight_layout()
    save(fig, "fig1_coverage_vs_columns")


def fig_ratio():
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    palette = ["#1f5fa8", "#b45f06", "#38761d"]
    for color, (a, n) in zip(palette, PANELS):
        s = [x for x in summ if (x["a"], x["n"]) == (a, n)]
        us = [x["columns"] for x in s]
        ratio = [x["coverage_ratio"] for x in s]
        lo = [x["coverage_ratio_ci95"][0] for x in s]
        hi = [x["coverage_ratio_ci95"][1] for x in s]
        ax.plot(us, ratio, "-o", ms=3.5, color=color, label=title(a, n))
        ax.fill_between(us, lo, hi, color=color, alpha=0.18, linewidth=0)
    ax.axhline(1.0, color="black", lw=0.8, ls=":")
    ax.set_xlabel("signed Frobenius orbits (columns U), equal unknowns and equal points")
    ax.set_ylabel("coverage ratio: greedy / mean random")
    ax.set_xlim(0, 26)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    ax.set_title(
        "Coverage ratio with 95% bootstrap interval (resampling random seeds)\n" + RUN,
        fontsize=9,
    )
    fig.tight_layout()
    save(fig, "fig2_coverage_ratio")


def fig_trials():
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 3.5), sharey=False)
    for ax, (a, n) in zip(axes, PANELS):
        g = {r["columns"]: r["expected_trials"] for r in series(a, n, "greedy")}
        rm = mean_by_u(a, n, "random", "expected_trials")
        us = sorted(g)
        ax.plot(us, [g[u] for u in us], "-o", ms=3.5, color=COLORS["greedy"], label="greedy")
        ax.plot(
            us,
            [rm[u] for u in us],
            "-",
            color=COLORS["random"],
            label="random orbits (mean of 20)",
        )
        ug = min(us, key=lambda u: g[u])
        ur = min(us, key=lambda u: rm[u])
        ax.scatter([ug], [g[ug]], s=70, facecolors="none", edgecolors=COLORS["greedy"], zorder=5)
        ax.scatter([ur], [rm[ur]], s=70, facecolors="none", edgecolors=COLORS["random"], zorder=5)
        ax.text(
            0.97,
            0.95,
            f"greedy min {g[ug]:.2f} at U={ug}\nrandom min {rm[ur]:.2f} at U={ur}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8"),
        )
        ax.set_title(title(a, n))
        ax.set_xlabel("signed Frobenius orbits (columns U)")
        ax.set_xlim(0, XMAX[(a, n)])
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("expected trials $(U+1)\\,T/\\mathrm{covered}$")
    axes[0].legend(fontsize=8.5)
    fig.suptitle("Expected trials to collect U+1 relations, from the census. " + RUN, fontsize=9)
    fig.tight_layout()
    save(fig, "fig3_expected_trials")


def fig_cost():
    a, n = 1, 17
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    for arm, label, color in [
        ("greedy", "greedy (all orbits scored)", COLORS["greedy"]),
        ("greedy_pool64", "greedy, 64 sampled per step", COLORS["pool"]),
        ("random", "random orbits (uniform draw)", COLORS["random"]),
    ]:
        ms = mean_by_u(a, n, arm, "build_ms")
        us = sorted(ms)
        ax.plot(us, [ms[u] for u in us], "-o", ms=3.5, color=color, label=label)
    ax.set_yscale("log")
    ax.set_xlabel("signed Frobenius orbits (columns U)")
    ax.set_ylabel("base construction time, ms (log scale)")
    ax.set_title(
        "Construction wall time, includes orbit enumeration (shared by all arms)\n"
        + f"{title(a, n)}, mean over seeds for the sampled arms; "
        + RUN,
        fontsize=8.5,
    )
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    save(fig, "fig4_construction_time")


if __name__ == "__main__":
    fig_coverage()
    fig_ratio()
    fig_trials()
    fig_cost()
    print("wrote fig1..fig4 (svg + png)")

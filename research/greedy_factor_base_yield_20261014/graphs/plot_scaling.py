"""Lay out the recorded values of the sampled-selector scaling runs.

Plotting only: every point is a value in scaling/*/summary.json, written by
examples/koblitz_greedy_scaling.rs.  The one transform is the horizontal axis
of the pool figure, sqrt(2 ln pool), which is the expected best of `pool`
standard normal draws to leading order; it is an axis, not a fit.

    python3 graphs/plot_scaling.py          # from this study directory
"""

import json
import math
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "legend.fontsize": 8.5})

HERE = pathlib.Path(__file__).resolve().parent
SCALING = HERE.parent / "scaling"
RUN = "run 2026-10-14, 64 candidate orbits per step, exact pair-sum coverage"


def load(path):
    return json.loads((SCALING / path / "summary.json").read_text())


main = load("main")
large = load("k0n37-large")
# The large-base run replaces the main run's K0/2^37 points.
rows = [r for r in main if (r["a"], r["n"]) != (0, 37)] + large
curves = sorted({(r["a"], r["n"]) for r in rows}, key=lambda k: [r for r in rows if (r["a"], r["n"]) == k][0]["subgroup_bits"])
palette = plt.cm.viridis([i / max(1, len(curves) - 1) for i in range(len(curves))])


def label(a, n):
    bits = [r for r in rows if (r["a"], r["n"]) == (a, n)][0]["subgroup_bits"]
    return f"$K_{a}/2^{{{n}}}$, $r=2^{{{bits:.1f}}}$"


def save(fig, stem):
    fig.savefig(HERE / f"{stem}.svg", format="svg")
    fig.savefig(HERE / f"{stem}.png", format="png", dpi=150)
    plt.close(fig)


def fig_scaling():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))
    for color, (a, n) in zip(palette, curves):
        s = sorted([r for r in rows if (r["a"], r["n"]) == (a, n)], key=lambda r: r["columns"])
        ax1.errorbar(
            [r["columns"] for r in s],
            [r["ratio"] for r in s],
            yerr=[
                [r["ratio"] - r["ratio_ci95"][0] for r in s],
                [r["ratio_ci95"][1] - r["ratio"] for r in s],
            ],
            fmt="-o",
            ms=3,
            lw=1,
            capsize=2,
            color=color,
            label=label(a, n),
        )
        ax2.plot(
            [r["random_mean_coverage"] for r in s],
            [r["ratio"] for r in s],
            "-o",
            ms=3,
            lw=1,
            color=color,
        )
    for ax in (ax1, ax2):
        ax.axhline(1.0, color="black", lw=0.8, ls=":")
        ax.grid(alpha=0.25)
        ax.set_ylabel("coverage ratio, greedy / random")
    ax1.set_xscale("log")
    ax1.set_xlabel("signed Frobenius orbits (columns U)")
    ax1.legend(fontsize=7.5, loc="upper left")
    ax2.set_xscale("log")
    ax2.set_xlabel("coverage of the random base (fraction of targets)")
    fig.suptitle("Greedy over random at equal columns, past the dense selector's cap. " + RUN, fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig5_scaling_ratio")


def fig_pool():
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    pools = [8, 64, 512]
    data = {p: load(f"pool{p}") for p in pools}
    for color, (a, n) in zip(["#1f5fa8", "#b45f06", "#38761d"], [(1, 17), (0, 23), (1, 23)]):
        xs, ys, lo, hi = [], [], [], []
        for p in pools:
            r = [x for x in data[p] if (x["a"], x["n"], x["columns"]) == (a, n, 16)][0]
            xs.append(math.sqrt(2 * math.log(p)))
            ys.append(r["ratio"] - 1)
            lo.append(r["ratio"] - r["ratio_ci95"][0])
            hi.append(r["ratio_ci95"][1] - r["ratio"])
        ax.errorbar(xs, ys, yerr=[lo, hi], fmt="-o", capsize=3, color=color, label=label(a, n))
    ax.set_xticks([math.sqrt(2 * math.log(p)) for p in pools])
    ax.set_xticklabels([f"{p}\n({math.sqrt(2 * math.log(p)):.2f})" for p in pools])
    ax.set_xlabel("candidate orbits scored per step (and $\\sqrt{2\\ln \\mathrm{pool}}$)")
    ax.set_ylabel("coverage ratio $-\\,1$, greedy / random")
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    ax.set_title(
        "Gain against pool size, 16 columns, 6 seeds, 95% bootstrap bars\nrun 2026-10-14, exact pair-sum coverage",
        fontsize=9,
    )
    fig.tight_layout()
    save(fig, "fig6_pool_sensitivity")


if __name__ == "__main__":
    fig_scaling()
    fig_pool()
    print("wrote fig5, fig6 (svg + png)")

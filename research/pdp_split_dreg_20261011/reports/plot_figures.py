"""Lay out the recorded values of ANALYSIS.json as SVG figures.

Layout only (AGENTS.md "Implementation language: no Python"): every number
plotted is read from ANALYSIS.json, which the native `pdp_split_analyze`
wrote from the run records and the ecbench session. Nothing is computed
here beyond choosing what to draw.

    python plot_figures.py ANALYSIS.json OUT_DIR
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

# Reference palette slots 1-3 (dataviz skill, validated all-pairs, light).
PLAN_STYLE = {
    "unsplit": ("#2a78d6", "o", "unsplit S4 (m=3)"),
    "split-2+1": ("#eb6834", "s", "S3 chain (m=3: 2+1)"),
    "split-2+1+1": ("#eb6834", "D", "S3 chain (m=4: 2+1+1, three-way)"),
    "split-2+2": ("#1baf7a", "^", "two-way split (m=4: S3 + S4 link)"),
}
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "svg.fonttype": "none",
    }
)


def cell_label(c):
    return f"m={c['m']}\nK{c['a']}/2^{c['n']}\nℓ={c['ell']}"


def ordered_cells(cells):
    return sorted(cells, key=lambda c: (c["m"], c["ell"], c["n"]))


def per_g_panels(analysis, field, ylabel, title, out, log=False):
    cells = ordered_cells(analysis["cells"])
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True)
    for ax, g in zip(axes, (0, 2, 4)):
        for xi, c in enumerate(cells):
            vs = [v for v in c["variants"] if v["g"] == g]
            offs = {p: k for k, p in enumerate(sorted({v["plan"] for v in vs}))}
            for v in vs:
                color, marker, _ = PLAN_STYLE[v["plan"]]
                ys = v[field] if isinstance(v[field], list) else [v[field]]
                dx = (offs[v["plan"]] - 0.5) * 0.25
                ax.scatter(
                    [xi + dx] * len(ys),
                    ys,
                    s=36,
                    color=color,
                    marker=marker,
                    edgecolor="white",
                    linewidth=1.0,
                    zorder=3,
                )
        ax.set_xticks(range(len(cells)))
        ax.set_xticklabels([cell_label(c) for c in cells], fontsize=7)
        ax.set_title(f"g = {g} guessed unknowns", fontsize=9, color=INK)
        if log:
            ax.set_yscale("log")
    axes[0].set_ylabel(ylabel)
    handles = [
        plt.Line2D([], [], color=s[0], marker=s[1], linestyle="", label=s[2])
        for s in PLAN_STYLE.values()
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8)
    fig.suptitle(title, fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    fig.savefig(out)
    plt.close(fig)


def projection_figure(analysis, out):
    cells = [c for c in ordered_cells(analysis["cells"]) if c["rho_strong"]]
    fig, ax = plt.subplots(figsize=(10, 4.2))
    xs = []
    for xi, c in enumerate(cells):
        xs.append(cell_label(c))
        # Rho references on this curve: neutral horizontal ticks, labelled.
        ax.hlines(c["rho_strong"]["mean_s"], xi - 0.4, xi + 0.4, color=INK2, linewidth=2)
        ax.hlines(
            c["rho_signed_frobenius"]["mean_s"],
            xi - 0.4,
            xi + 0.4,
            color=INK2,
            linewidth=2,
            linestyle=(0, (2, 2)),
        )
        plans = sorted({v["plan"] for v in c["variants"]})
        for k, plan in enumerate(plans):
            color, marker, _ = PLAN_STYLE[plan]
            for v in c["variants"]:
                if v["plan"] != plan or not v["projection"]:
                    continue
                ax.scatter(
                    xi + (k - 0.5) * 0.3,
                    v["projection"]["s_relations"],
                    color=color,
                    marker=marker,
                    s=30 + 10 * v["g"],
                    edgecolor="white",
                    linewidth=1.0,
                    zorder=3,
                )
        if not any(v["projection"] for v in c["variants"]):
            ax.text(xi, c["rho_strong"]["mean_s"] * 30, "no projection:\nλ = 0 in 256 targets"
                    if c["lambda"] == 0 else "no projection:\ngate failed", ha="center",
                    fontsize=7, color=INK2)
    ax.set_yscale("log")
    ax.set_xticks(range(len(cells)))
    ax.set_xticklabels(xs, fontsize=7)
    ax.set_ylabel("S (group-addition equivalents / √r), log scale")
    handles = [
        plt.Line2D([], [], color=s[0], marker=s[1], linestyle="", label=s[2])
        for s in PLAN_STYLE.values()
    ] + [
        plt.Line2D([], [], color=INK2, linewidth=2, label="strong rho, measured (ecbench)"),
        plt.Line2D([], [], color=INK2, linewidth=2, linestyle=(0, (2, 2)),
                   label="signed-Frobenius rho, measured (ecbench)"),
    ]
    ax.legend(handles=handles, fontsize=7, loc="upper left", ncol=2, frameon=False)
    ax.set_title(
        "Relation phase of each PDP variant, PROJECTED in S (pinned XOR ratio; marker size = g ∈ {0,2,4}),\n"
        "against measured rho on the same curve",
        fontsize=9,
        color=INK,
    )
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    analysis = json.loads(Path(sys.argv[1]).read_text())
    out = Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    per_g_panels(
        analysis,
        "solving_degree_per_target",
        "F4 solving degree",
        "Solving degree per target (4 targets per cell), by plan and guess count",
        out / "fig1-solving-degree.svg",
    )
    per_g_panels(
        analysis,
        "word_xors_per_target",
        "F4 64-bit word XORs (log)",
        "F4 work per target, all 2^g branches summed",
        out / "fig2-word-xors.svg",
        log=True,
    )
    projection_figure(analysis, out / "fig3-projection-vs-rho.svg")


if __name__ == "__main__":
    main()

"""Render the recorded frontier synthesis values; no research computation."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE.parent.parent / "frontier-data.json").read_text())

COLORS = {
    "factor base": "#167D9A",
    "F4/F5/F6": "#7649A8",
    "root finding": "#D27621",
}


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(HERE / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(HERE / f"{stem}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def cold_ops_gap() -> None:
    rows = DATA["cold_ops_gap"]
    fig, ax = plt.subplots(figsize=(10.8, 7.2))
    y = list(range(len(rows)))
    bars = ax.barh(
        y,
        [row["ratio_to_matched_rho"] for row in rows],
        color=[COLORS[row["group"]] for row in rows],
        edgecolor="#1B2430",
        linewidth=0.8,
    )
    for bar, row in zip(bars, rows):
        if row["lower_bound"]:
            bar.set_hatch("////")
        value = row["ratio_to_matched_rho"]
        ax.text(value * 1.08, bar.get_y() + bar.get_height() / 2, f"{value:.3g}x", va="center", fontsize=9)

    ax.axvline(1.0, color="#1B2430", linewidth=1.4, linestyle="--")
    ax.text(1.04, -0.85, "rho parity", fontsize=9, color="#1B2430")
    ax.set_xscale("log")
    ax.set_xlim(0.8, 900)
    ax.set_yticks(y, [row["label"] for row in rows])
    ax.invert_yaxis()
    ax.set_xlabel("Cold whole-pipeline S / matched rho S (log scale; lower is better)")
    ax.set_title("No measured IC candidate reaches the cold operations frontier")
    ax.grid(axis="x", which="both", alpha=0.22)
    legend = [Patch(facecolor=color, edgecolor="#1B2430", label=group) for group, color in COLORS.items()]
    legend.append(Patch(facecolor="white", edgecolor="#1B2430", hatch="////", label="lower bound: unpriced work remains"))
    ax.legend(handles=legend, loc="upper right", frameon=False, fontsize=9)
    fig.text(
        0.01,
        0.005,
        "Each bar uses its study's matched generic rho reference and reported estimator; curves, references, and estimators differ, so only parity is shared. "
        "Sources and derivations: frontier-data.json. Current audit: ECBHW1h7373da65e433.",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, "cold_ops_gap")


def online_boundary() -> None:
    rows = DATA["online_boundary"]
    labels = ["n61 compact orbit\n64 unseen targets", "n83 dense K1798\none public target"]
    x = [0, 1]
    width = 0.34
    online = [row["online_rho_over_ic"] for row in rows]
    cold = [row["setup_inclusive_rho_over_ic"] for row in rows]

    fig, ax = plt.subplots(figsize=(9.8, 6.2))
    online_bars = ax.bar([p - width / 2 for p in x], online, width, color="#2A9D6F", edgecolor="#1B2430", label="online rho / IC")
    cold_bars = ax.bar(
        [p + width / 2 for p in x],
        cold,
        width,
        color="#D1495B",
        edgecolor="#1B2430",
        label="cold / setup-inclusive diagnostic",
    )
    online_bars[1].set_hatch("\\\\")
    cold_bars[1].set_hatch("\\\\")

    ci = rows[0]["online_ci95"]
    ax.errorbar(
        x[0] - width / 2,
        online[0],
        yerr=[[online[0] - ci[0]], [ci[1] - online[0]]],
        fmt="none",
        ecolor="#1B2430",
        capsize=5,
        linewidth=1.2,
    )

    for bar, value in zip(online_bars, online):
        ax.text(bar.get_x() + bar.get_width() / 2, value * 1.16, f"{value:.3g}x", ha="center", va="bottom", fontsize=10)
    for bar, value, row in zip(cold_bars, cold, rows):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value / 1.45,
            f"1/{row['setup_inclusive_ic_over_rho']:.1f}",
            ha="center",
            va="top",
            fontsize=9,
        )

    ax.axhline(1.0, color="#1B2430", linewidth=1.4, linestyle="--")
    ax.text(1.45, 1.1, "parity", fontsize=9)
    ax.set_yscale("log")
    ax.set_ylim(0.001, 70)
    ax.set_xticks(x, labels)
    ax.set_ylabel("rho cost / IC cost (log scale; above 1 favours IC)")
    ax.set_title("Reusable setup reverses the apparent one-target boundary")
    ax.grid(axis="y", which="both", alpha=0.22)
    extra = Patch(facecolor="white", edgecolor="#1B2430", hatch="\\\\", label="exploratory shared-host N83")
    handles, legend_labels = ax.get_legend_handles_labels()
    ax.legend(
        handles + [extra],
        legend_labels + [extra.get_label()],
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        frameon=False,
    )
    fig.text(
        0.01,
        0.005,
        "n61 online: admitted isolated 95% interval [9.8, 19.1]; cold IC/rho median 83.4. "
        "n83: one-point exploratory online ratio; hatched setup diagnostic = reported rho online / complete IC process, not matched cold.",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.045, 0.79, 1))
    save(fig, "online_boundary")


if __name__ == "__main__":
    cold_ops_gap()
    online_boundary()

#!/usr/bin/env python3
"""Draw supplied native values; never estimate or generate research results.

Every series is supplied explicitly in figure_data.json, with native evidence
pointers, units, population and scope. This file handles axes, coordinates,
labels, colors and SVG/150 dpi PNG export only. It contains no aggregation,
quantile, ratio, selection, bootstrap or cryptographic operations.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
COLORS = ["#087F83", "#DA8749", "#63788F", "#A15D80", "#4B7EBD"]


def validate(spec):
    for key in ("name", "kind", "title", "source", "population"):
        if not spec.get(key):
            raise ValueError(f"Figure needs {key}: {spec.get('name')}")
    for series in spec.get("series", []):
        if len(series["x"]) != len(series["y"]):
            raise ValueError("x/y display arrays must have equal lengths")
        for bound in ("lower", "upper"):
            if bound in series and len(series[bound]) != len(series["y"]):
                raise ValueError("Every native display bound needs a matching point")


def draw_series(axis, spec):
    for index, series in enumerate(spec.get("series", [])):
        color = series.get("color", COLORS[index % len(COLORS)])
        x, y = series["x"], series["y"]
        label = series.get("label", "")
        kind = spec["kind"]
        if kind == "grouped_bars":
            # Offsets and widths are display geometry, not measured quantities.
            offset = series.get("display_offset", 0)
            axis.bar([v + offset for v in x], y,
                     width=series.get("display_width", 0.35), label=label,
                     color=color, alpha=series.get("alpha", 0.9))
        elif kind in ("points", "intervals"):
            axis.scatter(x, y, label=label, color=color,
                         s=series.get("marker_size", 28),
                         marker=series.get("marker", "o"),
                         alpha=series.get("alpha", 0.85), zorder=4)
        elif kind == "lines":
            axis.plot(x, y, label=label, color=color,
                      linestyle=series.get("linestyle", "-"),
                      marker=series.get("marker"),
                      linewidth=series.get("linewidth", 1.8))
        else:
            raise ValueError(f"Unknown presentation kind: {kind}")
        if "lower" in series and "upper" in series:
            # Endpoints are already recorded native intervals, not estimated.
            axis.vlines(x, series["lower"], series["upper"], color=color,
                        linewidth=1.6, alpha=0.85, zorder=3)


def plot(spec, output):
    import matplotlib.pyplot as plt
    validate(spec)
    fig, axis = plt.subplots(figsize=tuple(spec.get("figsize", [9.8, 5.4])))
    fig.patch.set_facecolor("white")
    axis.set_facecolor("#F8FAFC")
    draw_series(axis, spec)
    axis.set_title(spec["title"], loc="left", fontsize=13.5, fontweight="bold", pad=15)
    axis.set_xlabel(spec.get("xlabel", ""), labelpad=8)
    axis.set_ylabel(spec.get("ylabel", ""), labelpad=8)
    if "xticks" in spec:
        axis.set_xticks(spec["xticks"])
    if "xticklabels" in spec:
        axis.set_xticklabels(spec["xticklabels"], rotation=spec.get("label_rotation", 0))
    axis.set_xscale(spec.get("xscale", "linear"))
    axis.set_yscale(spec.get("yscale", "linear"))
    if "xlim" in spec:
        axis.set_xlim(spec["xlim"])
    if "ylim" in spec:
        axis.set_ylim(spec["ylim"])
    for marker in spec.get("vertical_markers", []):
        axis.axvline(marker["x"], color="#63788F", linestyle=":", linewidth=1.3)
        axis.text(marker["x"], marker.get("label_y", 0.96), marker["label"],
                  transform=axis.get_xaxis_transform(), rotation=90,
                  ha="right", va="top", fontsize=8.5, color="#63788F")
    for marker in spec.get("horizontal_markers", []):
        axis.axhline(marker["y"], color="#63788F", linestyle=":", linewidth=1.3)
        if marker.get("label"):
            axis.text(marker.get("label_x", 0.99), marker["y"], marker["label"],
                      transform=axis.get_yaxis_transform(), ha="right", va="bottom",
                      fontsize=8.5, color="#63788F")
    axis.grid(axis="y", linewidth=0.6, alpha=0.22)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)
    if spec.get("legend", True):
        axis.legend(loc=spec.get("legend_location", "best"), frameon=False,
                    fontsize=9, ncol=spec.get("legend_columns", 1))
    footer = spec["population"] + "\n" + spec["source"]
    if spec.get("scope_note"):
        footer = spec["scope_note"] + "\n" + footer
    fig.text(0.055, 0.017, footer, fontsize=8, color="#5B6A7A", va="bottom")
    fig.tight_layout(rect=[0, spec.get("footer_space", 0.16), 1, 0.98])
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{spec['name']}.svg", format="svg", facecolor="white")
    fig.savefig(output / f"{spec['name']}.png", dpi=150, facecolor="white")
    plt.close(fig)


def workflow(output):
    """Draw the stated accounting boundaries; no result is represented."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    fig, axis = plt.subplots(figsize=(10, 5.4))
    axis.set_xlim(0, 10)
    axis.set_ylim(0, 5.4)
    axis.axis("off")

    def box(x, y, width, height, text, color="#EFF6F7", border="#087F83", size=9):
        axis.add_patch(FancyBboxPatch((x, y), width, height,
                                     boxstyle="round,pad=0.08,rounding_size=0.10",
                                     facecolor=color, edgecolor=border, linewidth=1))
        axis.text(x + width/2, y + height/2, text, ha="center", va="center",
                  fontsize=size, color="#25354A", linespacing=1.35)

    def arrow(x1, y1, x2, y2):
        axis.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="->",
                                      mutation_scale=12, color="#63788F", linewidth=1.2))

    axis.text(0.2, 5.10, "Measured process boundaries", fontsize=15,
              fontweight="bold", color="#142B45")
    box(0.2, 4.03, 9.55, 0.65,
        "Inputs already available: exact public point(s) and retained K600 artifact\n"
        "Historical K600 selection / construction cost remains unavailable", size=10)
    box(0.2, 2.75, 2.05, 0.8, "Fresh singleton\n1 public target\n8 paired processes")
    box(2.70, 2.75, 3.02, 0.8, "Once-only setup\nParse, curve / basis, index,\nrank / logs or rho jumps")
    box(6.15, 2.75, 3.60, 0.8, "One complete target solve\nQ checks, attempts, recovery,\nbuilt-in verification, output / exit")
    arrow(2.29, 3.15, 2.58, 3.15)
    arrow(5.77, 3.15, 6.04, 3.15)
    box(0.2, 1.38, 2.05, 0.86, "Genuine shared cohort\n32 public targets\n3 paired cohorts")
    box(2.70, 1.38, 3.02, 0.86, "Once-only setup\nBoth arms initialize once;\nall preparation is charged")
    box(6.15, 1.38, 3.60, 0.86, "32 sequential target solves\nFresh target-dependent state / seed\nNo shared multi-target collision table")
    arrow(2.29, 1.81, 2.58, 1.81)
    arrow(5.77, 1.81, 6.04, 1.81)
    axis.text(0.2, 0.85, "Producer query timers are secondary: IC and rho Q-validation boundaries differ.",
              fontsize=10, color="#A15D80")
    axis.text(0.2, 0.35,
              "Outside both clocks: target generation, research / tuning, compilation, independent post-hoc audit.\n"
              "Cold denotes a fresh process, not an OS-cache flush. Complete process costs lead this report.",
              fontsize=9, color="#5B6A7A", linespacing=1.4)
    fig.tight_layout(pad=0.4)
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / "workload_boundaries.svg", format="svg", facecolor="white")
    fig.savefig(output / "workload_boundaries.png", dpi=150, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT/"figure_data.json")
    parser.add_argument("--output", type=Path, default=ROOT/"figures")
    args = parser.parse_args()
    payload = json.loads(args.data.read_text())
    if payload.get("evidence_status") != "NATIVE_RECORDS_VERIFIED_LAYOUT_ONLY":
        raise SystemExit("Native verified display data is not ready; plotting is withheld")
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                               "svg.fonttype": "none"})
    for figure in payload["figures"]:
        plot(figure, args.output)
    workflow(args.output)


if __name__ == "__main__":
    main()

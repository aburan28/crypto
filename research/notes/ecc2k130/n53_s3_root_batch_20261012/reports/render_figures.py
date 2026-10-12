"""Render native RESULT_DATA.json values only; no research computation."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPORTS = Path(__file__).resolve().parent
STUDY = REPORTS.parent
DATA = json.loads((STUDY / "RESULT_DATA.json").read_text())
RUN = REPORTS / "run-20261012T1737Z"
FAILURE = REPORTS / "run-20261012T1735Z"
PAIRS = [f"Q{q}/s{s}" for q in range(3) for s in range(3)]


def save(fig, directory, name):
    directory.mkdir(parents=True, exist_ok=True)
    for extension in ("svg", "png", "pdf"):
        kwargs = {"dpi": 150} if extension == "png" else {}
        destination = directory / f"{name}.{extension}"
        fig.savefig(destination, bbox_inches="tight", **kwargs)
        if extension == "svg":
            destination.write_text(
                "".join(line.rstrip(" \t\n") + "\n" for line in destination.read_text().splitlines())
            )
    plt.close(fig)


def ratio_panel(axis, names, keys, title, logarithmic=False):
    x = list(range(9))
    palette = ("#2070a0", "#b65723")
    for offset, (name, key) in enumerate(zip(names, keys)):
        entry = DATA[key]
        axis.plot(
            x,
            entry["values"],
            marker="o" if offset == 0 else "s",
            linewidth=1.2,
            markersize=4.5,
            color=palette[offset],
            label=f"{name} median {entry['median']:.3f}",
        )
        axis.axhline(
            entry["median"],
            color=palette[offset],
            linestyle=":",
            linewidth=1,
            alpha=0.8,
        )
    axis.axhline(1.0, color="#555555", linestyle="--", linewidth=0.8)
    axis.set_xticks(x, PAIRS, rotation=45, ha="right", fontsize=8)
    axis.set_ylabel("Ratio of in-process ms")
    axis.set_title(title)
    axis.grid(True, axis="y", alpha=0.22)
    axis.legend(frameon=False, fontsize=8, loc="best")
    if logarithmic:
        axis.set_yscale("log")


fig, axes = plt.subplots(2, 1, figsize=(11.2, 7.7), constrained_layout=True)
ratio_panel(
    axes[0],
    ["batch/control"],
    ["cold_candidate_over_control"],
    "Complete cold IC process (gate 0.90)",
)
axes[0].axhline(0.90, color="#945e9c", linestyle="-.", linewidth=0.9)
ratio_panel(
    axes[1],
    ["batch/control"],
    ["online_candidate_over_control"],
    "Target-online IC phase (gate 1.10)",
)
axes[1].axhline(1.10, color="#945e9c", linestyle="-.", linewidth=0.9)
fig.suptitle(
    "N53 K220 S3 root batching • 9 verified same-Q, same-rank-seed pairs\n"
    "Run 20261012T1737Z, shared Apple M4 Pro (L0); dotted = median, purple = frozen gate",
    fontsize=12,
)
save(fig, RUN, "paired_ratios")

fig, axes = plt.subplots(2, 1, figsize=(11.2, 7.7), constrained_layout=True)
ratio_panel(
    axes[0],
    ["control/rho", "batch/rho"],
    ["cold_control_over_rho", "cold_candidate_over_rho"],
    "Complete cold method against same-Q strong rho",
    logarithmic=True,
)
ratio_panel(
    axes[1],
    ["control/rho", "batch/rho"],
    ["online_control_over_rho", "online_candidate_over_rho"],
    "Target-online phase against same-Q strong rho",
)
fig.suptitle(
    "N53 K220 reference accounting • 9 verified public Q/seed pairs\n"
    "Run 20261012T1737Z, shared Apple M4 Pro (L0); dotted = medians, dashed = rho parity",
    fontsize=12,
)
save(fig, RUN, "rho_reference")

fig, ax = plt.subplots(figsize=(10.5, 3.0), constrained_layout=True)
ax.set_xlim(0, 10)
ax.set_ylim(0, 3)
ax.axis("off")
boxes = [
    (0.1, "Frozen source, 3 public Q,\n27-cell schedule"),
    (3.5, "Eight child outputs;\nninth interrupted"),
    (6.9, "time -l: sysctl denied;\nRSS absent, status incomplete"),
]
for x, label in boxes:
    ax.text(
        x + 1.45,
        1.55,
        label,
        ha="center",
        va="center",
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.5", "facecolor": "#edf3f7", "edgecolor": "#527492"},
    )
for start in (3.0, 6.4):
    ax.annotate(
        "",
        xy=(start + 0.4, 1.55),
        xytext=(start - 0.2, 1.55),
        arrowprops={"arrowstyle": "->", "color": "#527492", "lw": 1.4},
    )
ax.text(
    5.0,
    0.28,
    "Run 20261012T1735Z • measurement-wrapper infrastructure failure; "
    "raw attempt retained, excluded from timing comparison",
    ha="center",
    fontsize=9,
)
save(fig, FAILURE, "failure_path")

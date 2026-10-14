"""Render recorded native analyzer values; this script does no experiment work."""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: plot_recorded.py <analysis.json> <output-stem> <caption>")
    record = json.loads(Path(sys.argv[1]).read_text())
    arms = [arm for arm in record["arms"] if arm["status"] == "PASS_INDEPENDENT_REPLAY"]
    if not arms:
        raise SystemExit("no verified arms to plot")
    labels = [arm["arm"].replace("_", "\n") for arm in arms]
    colors = ["#5274a2" if not arm["root_table_dense_enabled"] else "#138f84" for arm in arms]
    series = [
        ("Mean target probes", [arm["mean_target_probes"] for arm in arms], "probes"),
        (
            "Allocated root slots",
            [arm["root_table_allocated_bytes"] / 2**30 for arm in arms],
            "GiB",
        ),
        (
            "Mean online time (exploratory)",
            [arm["mean_online_ms_exploratory"] for arm in arms],
            "ms",
        ),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.3), layout="constrained")
    for axis, (title, values, unit) in zip(axes, series):
        axis.bar(range(len(arms)), values, color=colors)
        axis.set_xticks(range(len(arms)), labels, fontsize=8)
        axis.set_ylabel(unit)
        axis.set_title(title, fontsize=11)
        axis.grid(axis="y", color="#dce3eb", linewidth=0.7)
        axis.set_axisbelow(True)
        axis.set_ylim(0, max(values) * 1.16)
        for index, value in enumerate(values):
            axis.annotate(f"{value:,.1f}", (index, value), ha="center", va="bottom", fontsize=7)
    n = arms[0]["targets"]
    fig.suptitle(f"{sys.argv[3]}  ·  {n} public target{'s' if n != 1 else ''} per arm", fontsize=13)
    stem = Path(sys.argv[2])
    fig.savefig(stem.with_suffix(".svg"))
    fig.savefig(stem.with_suffix(".png"), dpi=150)
    plt.close(fig)

    comparisons = [
        item for item in record.get("probe_comparisons", [])
        if item.get("status") == "PASS_PAIRED_PROBE_BOOTSTRAP"
    ]
    if comparisons:
        fig, axis = plt.subplots(figsize=(9.5, 4.4), layout="constrained")
        y_positions = list(range(len(comparisons)))
        means = [item["mean_difference_probes"] for item in comparisons]
        lower = [item["paired_bootstrap_95_percentile_probes"][0] for item in comparisons]
        upper = [item["paired_bootstrap_95_percentile_probes"][1] for item in comparisons]
        axis.errorbar(
            means,
            y_positions,
            xerr=[[mean - low for mean, low in zip(means, lower)],
                  [high - mean for mean, high in zip(means, upper)]],
            fmt="o",
            color="#138f84",
            capsize=4,
        )
        axis.axvline(0, color="#676f77", linewidth=1)
        axis.set_yticks(
            y_positions,
            [f"{item['left']} − {item['right']}" for item in comparisons],
            fontsize=9,
        )
        axis.invert_yaxis()
        axis.set_xlabel("Paired mean target probes (left − right); positive favors right")
        axis.set_title("64 public targets · 95% paired bootstrap intervals · 10,000 resamples")
        axis.grid(axis="x", color="#dce3eb", linewidth=0.7)
        axis.set_axisbelow(True)
        interval_stem = stem.with_name(stem.name + "_paired")
        fig.savefig(interval_stem.with_suffix(".svg"))
        fig.savefig(interval_stem.with_suffix(".png"), dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    main()

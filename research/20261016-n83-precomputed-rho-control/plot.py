#!/usr/bin/env python3
"""Lay out recorded online intervals from the Rust analyzer's analysis.json files.

Rendering only: every plotted number is read from records/<job>/analysis.json;
nothing is computed here beyond placing points.
usage: plot.py <run-records-dir> <out-dir> <file-stem> <run-label>
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ARMS = [
    ("ic", "IC (frozen producer)", "#7f3c8d"),
    ("rho_plain_fold", "rho, empty table (fold)", "#80868b"),
    ("bl_a_time_matched", "rho, precomputed table, time-matched (a)", "#11a579"),
    ("bl_a2_wall_equal", "rho, precomputed table, IC-wall-equal (a')", "#e68310"),
    ("bl_b_memory_matched", "rho, precomputed table, memory-matched (b)", "#3969ac"),
]
TICKS = {"ic": "IC", "rho_plain_fold": "rho\nempty", "bl_a_time_matched": "precomp\n(a) 650 s",
         "bl_a2_wall_equal": "precomp\n(a') = IC wall", "bl_b_memory_matched": "precomp\n(b) = IC RSS"}
JOBS = [("k600", "IC K600 job"), ("k3300", "IC K3300 job")]


def main():
    records, out, stem, run_label = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4]
    out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    for ax, (job, title) in zip(axes, JOBS):
        path = records / job / "analysis.json"
        if not path.exists():
            ax.set_title(f"{title}: no analysis.json")
            continue
        data = json.loads(path.read_text())
        present = [a for a in ARMS if a[0] in data["arms"]]
        for k, (key, label, color) in enumerate(present):
            arm = data["arms"][key]
            values = arm["online_ms"]
            xs = [k + (i - (len(values) - 1) / 2) * 0.025 for i in range(len(values))]
            ax.scatter(xs, values, s=14, color=color, alpha=0.75, label=label, zorder=3)
            lo, hi = arm["median_ms_ci95"]
            med = arm["median_ms"]
            ax.errorbar([k + 0.32], [med], yerr=[[med - lo], [hi - med]], fmt="D",
                        color=color, ms=6, capsize=4, zorder=4)
        ax.set_yscale("log")
        ax.set_xticks(range(len(present)))
        ax.set_xticklabels([TICKS[a[0]] for a in present], fontsize=8)
        lscpu = (records / job / "host" / "lscpu.txt").read_text().splitlines()
        cpu = next((l.split(":", 1)[1].strip() for l in lscpu if l.startswith("Model name")), "")
        ax.set_title(f"{title}, n = {data['targets']} points\n{cpu}", fontsize=9)
        ax.grid(axis="y", which="both", alpha=0.25)
    axes[0].set_ylabel("online interval per target (ms, log scale)")
    handles, labels = max((ax.get_legend_handles_labels() for ax in axes), key=lambda h: len(h[0]))
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, fontsize=8.5)
    fig.suptitle(f"N83 a=1 Koblitz, r = 8569786107849059, {run_label}: online time per target; "
                 "diamond = median, 95% bootstrap CI", fontsize=10)
    fig.tight_layout(rect=(0, 0.1, 1, 0.95))
    fig.savefig(out / f"{stem}.svg")
    fig.savefig(out / f"{stem}.png", dpi=150)


if __name__ == "__main__":
    main()

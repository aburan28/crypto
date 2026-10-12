#!/usr/bin/env python3
"""Plot native recorded N83 results only; no statistics or result computation."""
from pathlib import Path
import json
import os
import subprocess

os.environ.setdefault("MPLCONFIGDIR", "/tmp/n83-report-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
DATA = json.loads((ROOT / "data/report_data_secondary.json").read_text())
STAT = json.loads((ROOT / "data/statistics.json").read_text())
IC = DATA["methods"]["ic"]["rows"]
RHO = DATA["methods"]["rho"]["rows"]
NAVY, TEAL, BLUE, MUTED, GRID = "#142B45", "#087F83", "#527AAB", "#5B6A7A", "#DAE4EB"


def save(fig, name):
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", dpi=210, bbox_inches="tight")
    plt.close(fig)


def runtime_plot():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.0))
    fig.subplots_adjust(left=.07, right=.99, bottom=.24, top=.76, wspace=.26)
    positions = list(range(5))
    labels = [r["pair_id"] for r in IC]
    for ax, key, label, title in (
        (axes[0], "online_ms", "Complete prepared-online wall (ms)", "Prepared online"),
        (axes[1], "whole_process_wall_s", "Outer fresh-process wall (s; log axis)", "From retained base"),
    ):
        ax.bar([x-.17 for x in positions], [r[key] for r in RHO], .32, color=BLUE, label="Strong rho")
        ax.bar([x+.17 for x in positions], [r[key] for r in IC], .32, color=TEAL, label="IC / K600")
        ax.set_xticks(positions, labels)
        ax.set_xlabel("Paired seed round")
        ax.set_ylabel(label)
        ax.set_title(title, loc="left", fontsize=12.5, fontweight="bold", color=NAVY, pad=12)
        ax.grid(axis="y", color=GRID, linewidth=.7)
        ax.set_axisbelow(True)
        if key == "whole_process_wall_s":
            ax.set_yscale("log")
            ax.set_ylim(1, 2500)
            ticks = [1,10,100,1000]
            ax.yaxis.set_major_locator(FixedLocator(ticks))
            ax.yaxis.set_major_formatter(FixedFormatter([str(v) for v in ticks]))
            ax.minorticks_off()
        else:
            ax.set_ylim(0, 7800)
    fig.text(.025, .96, "One fixed Q, five paired seed rounds", fontsize=18,
             fontweight="bold", color=NAVY, va="top")
    fig.text(.025, .89, "Recorded complete intervals; the process metric includes measured preparation from retained input", fontsize=10.8, color=MUTED)
    axes[0].legend(loc="upper left", frameon=False, ncol=2, fontsize=9.5)
    fig.text(.025, .08, "n=5 pairs on one target and one host. Warm-ups excluded; all formal answers independently replayed.", fontsize=9.7, color=MUTED)
    fig.text(.025, .025, "Fresh-process/state cold: original base selection/construction omitted; OS page cache not flushed. run-20261012T0226Z", fontsize=9.2, color=MUTED)
    save(fig, "paired_runtime")


def interval_plot():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    fig.subplots_adjust(left=.16, right=.97, top=.75, bottom=.28, wspace=.55)
    online = DATA["primary_reciprocal_presentation"]["online"]
    process = DATA["primary_reciprocal_presentation"]["whole_process"]
    entries = [("IC / strong", online["ic_over_strong_primary_ratio"], online["ic_over_strong_primary_95_ci"], 1, TEAL),
               ("Strong A / B", STAT["aa"]["geometric_mean_paired_online_ratio"], STAT["aa"]["paired_percentile_bootstrap_95_ci"], 0, BLUE)]
    for label, center, interval, y, color in entries:
        # Differences are geometric coordinates for an errorbar, not new estimates.
        axes[0].errorbar(center,y,xerr=[[center-interval[0]],[interval[1]-center]],fmt="o",color=color,capsize=5,markersize=7,elinewidth=2)
    axes[0].set_yticks([1,0], ["IC / strong", "A/A control"])
    axes[0].set_ylim(-.6,1.6)
    axes[0].set_xlim(.85,2.05)
    axes[0].set_xticks([1,1.25,1.5,1.75,2])
    axes[0].set_xlabel("Prepared-online ratio")
    center = process["ic_over_strong_primary_ratio"]
    interval = process["ic_over_strong_primary_95_ci"]
    axes[1].errorbar(center,0,xerr=[[center-interval[0]],[interval[1]-center]],fmt="o",color=TEAL,capsize=5,markersize=7,elinewidth=2)
    axes[1].set_yticks([0],["IC / strong"])
    axes[1].set_ylim(-.6,.6)
    axes[1].set_xscale("log")
    axes[1].set_xlim(.85,400)
    ticks = [1,10,100,300]
    axes[1].xaxis.set_major_locator(FixedLocator(ticks))
    axes[1].xaxis.set_major_formatter(FixedFormatter([str(v) for v in ticks]))
    axes[1].minorticks_off()
    axes[1].set_xlabel("Fresh-process wall ratio (log axis)")
    for ax in axes:
        ax.axvline(1,color=NAVY,linestyle=(0,(4,3)),linewidth=1.3)
        ax.grid(axis="x",color=GRID,lw=.7)
        ax.tick_params(axis="y",length=0)
    fig.text(.025,.96,"Native paired geometric means and 95% intervals",fontsize=17.5,fontweight="bold",color=NAVY,va="top")
    fig.text(.025,.875,"IC/strong values are native reciprocals of the frozen strong/IC estimator; no new resampling",fontsize=10.6,color=MUTED)
    fig.text(.025,.105,"IC/strong > 1 means greater IC cost; parity is the dashed line. A/A compares two strong runs with the same seed.",fontsize=9.5,color=MUTED)
    fig.text(.025,.055,"Each interval is conditional on n=5 paired repeats of one Q and one host. It is not a fresh-target population interval.",fontsize=9.5,color=MUTED)
    fig.text(.025,.012,"Sources: statistics.json + native report_data_secondary.json · run-20261012T0226Z",fontsize=9.2,color=MUTED)
    save(fig,"paired_intervals")


def rasterize_diagrams():
    for name in ("three_spaces", "coefficient_transport"):
        subprocess.run(["inkscape", str(FIG/f"{name}.svg"), "--export-type=png",
                        f"--export-filename={FIG/name}.png", "--export-width=2200"], check=True)


if __name__ == "__main__":
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10.2,
                         "axes.spines.top":False,"axes.spines.right":False,
                         "svg.fonttype":"none","savefig.facecolor":"white"})
    FIG.mkdir(exist_ok=True)
    runtime_plot()
    interval_plot()
    rasterize_diagrams()
    print("Rendered recorded-value graphs and editable-diagram previews; no statistics computed.")

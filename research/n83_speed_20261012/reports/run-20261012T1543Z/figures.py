#!/usr/bin/env python3
"""Lay out native recorded values and the frozen campaign diagram only.

No algorithm, target generation, benchmark, selection, result calculation,
resampling, or independent answer replay occurs in this renderer.
Run only after the isolated native campaign has ended.
"""
from pathlib import Path
import json
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/n83-exact-scalar-report-mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/n83-exact-scalar-report-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures"
DATA = json.loads((ROOT / "data/statistics.json").read_text())
POLICY = json.loads((ROOT / "data/policy_selection.json").read_text())
ARITHMETIC = json.loads((ROOT / "data/arithmetic_selection.json").read_text())
AA = json.loads((ROOT / "data/noise_summary.json").read_text())
RUN = ROOT.name
NAVY, TEAL, BLUE, ORANGE, MUTED, GRID = "#142B45", "#087F83", "#527AAB", "#B66A32", "#5B6A7A", "#DAE4EB"


def save(fig, name):
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def caption(fig, text, source):
    fig.text(.025, .075, text, fontsize=9.2, color=MUTED)
    fig.text(.025, .025, f"Source: native {source} · {RUN} · recorded values only", fontsize=9.0, color=MUTED)


def campaign_flow():
    fig, ax = plt.subplots(figsize=(10, 4.45))
    ax.set_xlim(0, 10); ax.set_ylim(0, 4.2); ax.axis("off")
    boxes = [(0.15, 2.32, 2.65, 1.18, "Fresh native corpus", "3 training + 5 holdout points\nExclude all prior signed orbits", BLUE),
             (3.4, 2.32, 2.65, 1.18, "A/A: 10 processes", "Five identical pairs\nDiagnostic only; no selection", BLUE),
             (6.65, 2.32, 3.05, 1.18, "Arithmetic: 18 processes", "Reference / exact Barrett\n9 paired training cases", TEAL),
             (0.15, .46, 2.65, 1.18, "Policy: 81 processes", "Lanes 64 / 128 / 256\nDP 6 / 8 / 10; escape off", TEAL),
             (3.4, .46, 2.65, 1.18, "Frozen selection", "All cases / losses retained\nNo held-out target tuning", TEAL),
             (6.65, .46, 3.05, 1.18, "Holdout: 30 processes", "Current tuned / selected\n5 fresh points × 3 seeds × 2", BLUE)]
    for x, y, width, height, heading, body, color in boxes:
        ax.add_patch(FancyBboxPatch((x, y), width, height, boxstyle="round,pad=.08", edgecolor=color, facecolor="#F3F7FA", linewidth=1.4))
        ax.text(x + .13, y + height - .2, heading, va="top", fontsize=11.4, fontweight="bold", color=color)
        ax.text(x + .13, y + .65, body, va="top", fontsize=9.5, linespacing=1.55, color=NAVY)
    for start, end in [((2.94, 2.91), (3.25, 2.91)), ((6.17, 2.91), (6.51, 2.91)),
                       ((2.94, 1.04), (3.25, 1.04)), ((6.17, 1.04), (6.51, 1.04))]:
        ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle":"->", "color":MUTED, "lw":1.4})
    ax.plot([8.2, 8.2, 1.48], [2.1, 1.95, 1.95], color=MUTED, lw=1.2)
    ax.annotate("", xy=(1.48, 1.72), xytext=(1.48, 1.95), arrowprops={"arrowstyle":"->", "color":MUTED, "lw":1.2})
    fig.text(.025, .96, "Frozen 139-process follow-up: exact arithmetic, policy and fresh holdout", fontsize=16.1, fontweight="bold", color=NAVY, va="top")
    caption(fig, "Whole-process comparison against the previously tuned control. No new IC run or IC/rho ratio.", "PROTOCOL.md / stage manifests")
    save(fig, "campaign_flow")


def arithmetic_scores():
    rows = ARITHMETIC["candidate_scores"]
    fig, ax = plt.subplots(figsize=(10, 3.8)); fig.subplots_adjust(left=.32, right=.95, top=.72, bottom=.28)
    y = list(range(len(rows)-1, -1, -1))
    values = [r["primary_capped_arithmetic_mean_process_wall_s"] for r in rows]
    ax.barh(y, values, height=.57, color=[TEAL if r["candidate_id"] == ARITHMETIC["selected_candidate_id"] else BLUE for r in rows])
    ax.set_yticks(y, [r["candidate_id"] for r in rows]); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Native mean whole-process wall (s), 9 training cases per kernel")
    ax.grid(axis="x", color=GRID, lw=.7); ax.set_axisbelow(True)
    ax.margins(x=.18)
    for pos, value in zip(y, values):
        ax.text(value, pos, f"  {value:.6f}", va="center", fontsize=9, color=MUTED)
    fig.text(.025, .96, "Exact scalar arithmetic: all registered native training scores", fontsize=16.3, fontweight="bold", color=NAVY, va="top")
    fig.text(.025, .85, "Same Fold64/DP6/off policy, point and seed; recorded semantics must match in all nine pairs.", fontsize=10.3, color=MUTED)
    caption(fig, "Teal marks frozen native selection. Equality covers recorded counters and charges; no full step trace.", "arithmetic_selection.json")
    save(fig, "arithmetic_scores")


def aa_noise():
    values = AA["process_wall_a_over_b"]
    if not AA["all_complete_verified_and_counter_equal"] or values is None:
        raise ValueError("A/A diagnostic is incomplete; describe its failure instead of plotting a ratio")
    fig, ax = plt.subplots(figsize=(10, 4.0)); fig.subplots_adjust(left=.11, right=.96, top=.75, bottom=.27)
    ax.scatter(list(range(1, 6)), values["ratios"], s=52, color=BLUE, zorder=3)
    ax.axhline(1, color=NAVY, linestyle=(0, (4, 3)), lw=1.2)
    ax.axhline(values["median"], color=TEAL, lw=1.2, label=f"Native median {values['median']:.6f}")
    ax.axhspan(values["minimum"], values["maximum"], color=TEAL, alpha=.09)
    ax.set_xticks(list(range(1, 6))); ax.set_xlabel("Identical-run A/A pair"); ax.set_ylabel("Recorded A/B process-wall ratio")
    ax.grid(axis="y", color=GRID, lw=.7); ax.legend(frameon=False, fontsize=9)
    fig.text(.025, .96, "Five A/A pairs: fresh-process runtime variability", fontsize=16.3, fontweight="bold", color=NAVY, va="top")
    fig.text(.025, .86, f"Recorded min {values['minimum']:.6f}; max {values['maximum']:.6f}; range {values['range']:.6f}.", fontsize=10.5, color=MUTED)
    caption(fig, "One training point and seed; excluded from selection and holdout. Shading is observed range, not a CI.", "noise_summary.json")
    save(fig, "aa_noise")


def policy_scores():
    rows = POLICY["candidate_scores"]
    fig, ax = plt.subplots(figsize=(10, 5.0)); fig.subplots_adjust(left=.31, right=.94, top=.80, bottom=.23)
    y = list(range(len(rows)-1, -1, -1))
    values = [r["primary_capped_arithmetic_mean_process_wall_s"] for r in rows]
    ax.barh(y, values, color=[TEAL if r["candidate_id"] == POLICY["selected_candidate_id"] else BLUE for r in rows], height=.64)
    ax.set_yticks(y, [r["candidate_id"] for r in rows]); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Native arithmetic mean whole-process wall (s), 9 training cases per policy")
    ax.grid(axis="x", color=GRID, lw=.7); ax.set_axisbelow(True)
    ax.margins(x=.18)
    for pos, value in zip(y, values):
        ax.text(value, pos, f"  {value:.3f}", va="center", fontsize=9, color=MUTED)
    fig.text(.025, .96, "Expanded lane / distinguished-point grid: every registered score", fontsize=16.1, fontweight="bold", color=NAVY, va="top")
    fig.text(.025, .87, "Lane counts 64 / 128 / 256 × DP6 / DP8 / DP10; selected scalar kernel; escape fixed off.", fontsize=10.2, color=MUTED)
    caption(fig, "Teal marks the frozen native selection. Walk policies may change trajectories; total process cost is scored.", "policy_selection.json")
    save(fig, "policy_scores")


def holdout_cases():
    targets = DATA["per_target"]
    fig, ax = plt.subplots(figsize=(10, 4.8)); fig.subplots_adjust(left=.10, right=.97, top=.74, bottom=.26)
    x = list(range(len(targets)))
    for arm, shift, color, label in [("current", -.18, BLUE, "Current tuned: reference Fold64/DP6"), ("selected", .18, TEAL, "Frozen selected scalar / policy")]:
        ax.bar([a + shift for a in x], [t[arm]["mean_verified_actual_process_wall_s"] for t in targets], width=.32, color=color, label=label)
        for n, target in enumerate(targets):
            rows = [r for r in DATA["raw_rows"] if r["target_id"] == target["target_id"] and r["arm"] == arm]
            offsets = [-.08, 0, .08]
            ax.scatter([n + shift + offsets[j] for j in range(len(rows))], [r["elapsed_process_wall_s"] for r in rows], color=NAVY, s=18, zorder=3, edgecolor="white", linewidth=.4)
    ax.set_xticks(x, [f"Q{n+1}" for n in x]); ax.set_ylim(bottom=0)
    ax.set_ylabel("Whole-process wall (s)"); ax.set_xlabel("Fresh held-out target block")
    ax.grid(axis="y", color=GRID, lw=.7); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8.7, ncol=2, loc="lower left", bbox_to_anchor=(0,1.035), borderaxespad=0)
    fig.text(.025, .96, "Fresh holdout: all process times and native three-seed means", fontsize=16.5, fontweight="bold", color=NAVY, va="top")
    fig.text(.025, .87, "Five disjoint target blocks; dots retain all 30 individual solves. Bars are recorded native target means.", fontsize=10.2, color=MUTED)
    caption(fig, "Both arms share one reference-scalar monomorphization and executable; only lane count differs.", "statistics.json")
    save(fig, "holdout_cases")


def holdout_intervals():
    fields = [("Current/selected process: primary", DATA["primary_current_over_selected_process_wall"], TEAL),
              ("Current/selected online: secondary", DATA["secondary_current_over_selected_online"], BLUE),
              ("Arithmetic mean-time ratio: secondary", DATA["secondary_current_over_selected_arithmetic_process_wall"], ORANGE)]
    fig, ax = plt.subplots(figsize=(10, 4.1)); fig.subplots_adjust(left=.39, right=.96, top=.75, bottom=.27)
    for y, (_, value, color) in zip([2, 1, 0], fields):
        center = value["ratio"]; low, high = value["ci95_percentile"]
        # Endpoint differences place a native interval on the canvas only.
        ax.errorbar(center, y, xerr=[[center-low], [high-center]], fmt="o", color=color, capsize=5, markersize=7, elinewidth=2)
    ax.set_yticks([2, 1, 0], [a[0] for a in fields]); ax.tick_params(axis="y", length=0)
    ax.set_xscale("log"); ax.set_ylim(-.6, 2.6); ax.axvline(1, color=NAVY, linestyle=(0, (4, 3)), lw=1.3)
    ax.grid(axis="x", color=GRID, lw=.7); ax.set_xlabel("Recorded native current / selected cost ratio (log axis)")
    fig.text(.025, .96, "Native held-out ratios and 95% target-block intervals", fontsize=16.6, fontweight="bold", color=NAVY, va="top")
    fig.text(.025, .87, "Five paired target blocks; 20,000 native resamples preserve each target's three-seed means.", fontsize=10.3, color=MUTED)
    caption(fig, "Values above 1 favor selected. Geometric and arithmetic estimands retain separate definitions.", "statistics.json")
    save(fig, "holdout_intervals")


if __name__ == "__main__":
    if not DATA["all_verified"]:
        raise SystemExit("Incomplete holdout: use an explicit failure diagram and keep native null estimates.")
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10, "axes.spines.top":False, "axes.spines.right":False, "svg.fonttype":"none", "savefig.facecolor":"white"})
    OUT.mkdir(exist_ok=True)
    campaign_flow(); arithmetic_scores(); aa_noise(); policy_scores(); holdout_cases(); holdout_intervals()
    print("Rendered frozen diagram plus native recorded means/raw values/intervals only; no analysis or resampling.")

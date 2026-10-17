"""Render the lit2026 speedup journey from recorded entries.

Layout only: every ratio and interval is copied from an entry file, which
copies it from its run report. The board counts recorded steps by kind.
No experiment, measurement, fit, or cross-session ratio is computed here.

usage: python -I render_journey.py <program_dir>
writes <program_dir>/journey/journey.svg and <program_dir>/journey/board.svg
"""
import collections
import json
import pathlib
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def save_svg(figure, path):
    figure.savefig(path)
    # Matplotlib leaves spaces before newlines in path data. Normalising
    # that layout whitespace keeps generated SVG diffs clean.
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")


root = pathlib.Path(sys.argv[1])
ws = json.loads((root / "workstreams.json").read_text())
entries = []
for path in sorted((root / "journey" / "entries").glob("*.json")):
    entry = json.loads(path.read_text())
    entry["_file"] = path.name
    entries.append(entry)
superseded = {e.get("supersedes") for e in entries if e.get("supersedes")}
entries = [e for e in entries if e["_file"] not in superseded]
measured = [e for e in entries
            if (e.get("ratio_original_over_candidate") or {}).get("value") is not None]
measured.sort(key=lambda e: (int(e["workstream"][1:]), e["workload"]["id"],
                             e["metric"], e["date"], e["_file"]))

# A row per recorded comparison avoids overlapping labels for same-date
# steps. Units stay explicit beside every workload; no ratios are chained.
fig, ax = plt.subplots(figsize=(13, 7.2))
colors = {"counted ops": "#2673b8", "s": "#8b5ca6"}
markers = {"counted ops": "s", "s": "o"}
if measured:
    labels = []
    for row, entry in enumerate(measured):
        ratio = entry["ratio_original_over_candidate"]
        value = ratio["value"]
        unit = entry["unit"]
        color = colors.get(unit, "#555555")
        interval = ratio.get("interval") or [None, None]
        if interval[0] is not None and interval[1] is not None:
            ax.hlines(row, interval[0], interval[1], color=color, linewidth=1.5)
        point = ax.scatter([value], [row], color=color,
                           marker=markers.get(unit, "o"), s=24, zorder=3)
        point.set_gid(entry["_file"])
        labels.append(f'{entry["date"]}  {entry["workstream"]} {entry["slug"]} | '
                      f'{entry["workload"]["id"]} [{unit}]')
        if row % 2 == 0:
            ax.axhspan(row - 0.5, row + 0.5, color="#f2f4f7", zorder=0)
    ax.set_yticks(range(len(labels)), labels, fontsize=7.3)
    ax.set_ylim(len(labels) - 0.5, -0.5)
    ax.set_xscale("log")
    ax.axvline(1.0, color="#555555", linewidth=0.8, linestyle="--")
    units = sorted({e["unit"] for e in measured})
    ax.legend(handles=[Line2D([], [], color=colors.get(unit, "#555555"),
                             marker=markers.get(unit, "o"), linestyle="none",
                             label=unit) for unit in units],
              fontsize=8, frameon=False, loc="lower left", bbox_to_anchor=(0, 1.02),
              ncol=len(units), borderaxespad=0)
    ax.grid(axis="x", alpha=0.15)
    ax.set_xlabel("Recorded original / candidate ratio (>1 uses less of the recorded metric)",
                  fontsize=8)
    ax.tick_params(axis="y", length=0, pad=6)
    ax.tick_params(axis="x", labelsize=8)
else:
    ax.text(0.5, 0.5, "No measured comparisons recorded yet.",
            ha="center", va="center", transform=ax.transAxes)
    ax.set_xticks([])
    ax.set_yticks([])
fig.suptitle(f"lit2026 program baseline {ws['program_baseline_commit']}: recorded comparisons\n"
             f"{len(entries)} entries; {len(measured)} have baseline ratios. "
             "Each row retains its workload, unit and recorded date.", fontsize=10)
fig.text(0.01, 0.01, "Source: research/lit2026_speedup_program/journey/entries/*.json. "
         "Intervals are copied as recorded; these are separate stage or method comparisons.",
         fontsize=7, color="#555555")
fig.subplots_adjust(left=0.63, right=0.99, bottom=0.11, top=0.88)
save_svg(fig, root / "journey" / "journey.svg")
plt.close(fig)

# Separate axes keep the complete workstream descriptions outside the
# counting plot rather than letting long labels cover bars or get clipped.
kinds = ["speedup", "no-change", "regression", "diagnostic", "negative"]
kcol = {"speedup": "#3aa76d", "no-change": "#9aa0a6", "regression": "#d0453b",
        "diagnostic": "#2a6fdb", "negative": "#7a4cc2"}
counts = collections.Counter((e["workstream"], e["kind"]) for e in entries)
workstreams = ws["workstreams"]
fig, (ax, descriptions) = plt.subplots(1, 2, figsize=(13, 6),
                                       gridspec_kw={"width_ratios": [1, 2.5]})
left = [0] * len(workstreams)
for kind in kinds:
    values = [counts[(w["id"], kind)] for w in workstreams]
    ax.barh(range(len(workstreams)), values, left=left, color=kcol[kind],
            label=kind, edgecolor="white")
    left = [a + b for a, b in zip(left, values)]
ax.set_yticks(range(len(workstreams)), [w["id"] for w in workstreams], fontsize=9)
ax.set_ylim(len(workstreams) - 0.5, -0.5)
ax.set_xlim(0, max(left, default=0) + 1)
ax.set_xlabel("Recorded steps (count)", fontsize=9)
ax.xaxis.get_major_locator().set_params(integer=True)
ax.tick_params(axis="x", labelsize=8)
ax.grid(axis="x", alpha=0.15)
descriptions.set_xlim(0, 1)
descriptions.set_ylim(len(workstreams) - 0.5, -0.5)
descriptions.axis("off")
for row, workstream in enumerate(workstreams):
    heading = f'[{workstream["status"]}] ' + ", ".join(workstream["prs"])
    description = textwrap.fill(heading + " | " + workstream["title"], width=110)
    descriptions.text(0, row, description, fontsize=7.1, va="center", linespacing=1.05)
fig.suptitle(f"lit2026 workstreams: status and {len(entries)} recorded steps by kind", fontsize=11)
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, fontsize=8, frameon=False, loc="lower center",
           bbox_to_anchor=(0.5, 0.04), ncol=len(kinds))
fig.text(0.01, 0.01, "Source: workstreams.json and journey/entries/*.json. "
         "PR status does not close the research requirements recorded in the source reports.",
         fontsize=7, color="#555555")
fig.subplots_adjust(left=0.05, right=0.99, bottom=0.17, top=0.91, wspace=0.06)
save_svg(fig, root / "journey" / "board.svg")
plt.close(fig)
print(f"rendered {len(entries)} entries ({len(measured)} recorded baseline comparisons)")

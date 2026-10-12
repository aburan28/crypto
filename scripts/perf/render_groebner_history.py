#!/usr/bin/env python3
"""Render prepared native plot data; never execute benchmarks or derive speedups."""
import argparse
import json
import os
import tempfile
from pathlib import Path

# Rendering caches must also work in managed environments with read-only homes.
cache = Path(tempfile.mkdtemp(prefix="groebner-history-render-"))
os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(cache / "fontconfig"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import typst

matplotlib.rcParams.update({"svg.fonttype": "none", "svg.hashsalt": "groebner-history/v1"})


def typst_text(value, code=False):
    return ("text(" if code else "#text(") + json.dumps(str(value), ensure_ascii=False) + ")"


def draw_chart(group, chart, output):
    fig, ax = plt.subplots(figsize=(10.4, 3.35), layout="constrained")
    fig.supxlabel(chart["caption"], fontsize=7.1)
    ticks = {}
    boundaries = set()
    for series in chart["series"]:
        segment = []
        labelled = False

        def flush():
            nonlocal segment, labelled
            if segment:
                ax.plot(
                    [p["x"] for p in segment],
                    [p["value"] for p in segment],
                    color=series["color"], marker="o", markersize=5, linewidth=2,
                    label=series["label"] if not labelled else None,
                )
                labelled = True
                segment = []

        for point in series["points"]:
            ticks[point["x"]] = point["label"]
            if point["break_before"]:
                flush()
                if point["x"]:
                    boundaries.add(point["x"] - 0.5)
            if point["value"] is None:
                flush()
            else:
                segment.append(point)
        flush()
    for boundary in sorted(boundaries):
        ax.axvline(boundary, color="#a8b1bd", linewidth=0.8, linestyle="--", alpha=0.65)
    if not ax.lines:
        ax.text(0.5, 0.5, "No recorded measurements", transform=ax.transAxes, ha="center")
    ax.set_xticks(list(ticks), list(ticks.values()), fontsize=7.7)
    ax.set_ylabel(chart["unit"], fontsize=9)
    ax.set_xlabel("Measured revision order · dates and commits shown", fontsize=8)
    ax.set_title(chart["title"], fontsize=11, loc="left", pad=10)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}" if chart["id"].endswith(("instructions", "bytes", "requests")) else f"{value:,.1f}"))
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", color="#dbe1e8", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False, fontsize=9, loc="best")
    fig.savefig(output, metadata={"Date": None, "Title": group["title"], "Description": group["workload"]})
    plt.close(fig)


def table_rows(group):
    """Arrange recorded native labels; no new numerical result is computed."""
    charts = group["charts"]
    columns = ["UTC date / revision"]
    rows = {}
    for chart in charts:
        for series in chart["series"]:
            columns.append(series["label"] + (" Ir" if chart["id"].endswith("instructions") else
                           " ms" if chart["id"].endswith("runtime") else
                           " bytes" if chart["id"].endswith("bytes") else " requests"))
            for point in series["points"]:
                row = rows.setdefault(point["x"], [point["date"] + " / " + point["revision"][:8]])
                row.append(point["value_label"])
    return columns, [row for _, row in sorted(rows.items())]


def group_body(group, source_count):
    columns, rows = table_rows(group)
    body = ["= " + group["title"], typst_text(group["workload"]), ""]
    for chart in group["charts"]:
        body.append('#image("' + chart["id"] + '.svg", width: 100%)')
    widths = "1.7fr, " + ", ".join("1fr" for _ in columns[1:])
    cells = ",\n".join(typst_text(cell, code=True) for row in rows for cell in row)
    header = ", ".join(typst_text(cell, code=True) for cell in columns)
    body.append("#text(size: 7.5pt)[#table(columns: (" + widths + "), inset: 3pt, table.header(" + header + "),\n" + cells + ")]\n")
    body.extend(typst_text(note) + "\n" for note in group["notes"])
    if group["id"].startswith("boolean-"):
        body.append(typst_text(f"Evidence: {source_count} imported perf.history/v1 snapshots; historical_sources in plot-data.json retains raw filenames, SHA-256 hashes, hosts, harnesses and sample counts."))
    else:
        body.append(typst_text("Evidence: allocation_sources in plot-data.json; frozen native comparison and baseline/candidate JSONL records in research/f4_f5_row_build_20261011/records/."))
    body.append(typst_text("Source data and provenance: plot-data.json and REPORT.md. These plots do not establish an end-to-end ECDLP gain."))
    return "\n".join(body)


PREAMBLE = """#set document(date: none)
#set page(paper: "a4", margin: 15mm, numbering: "1")
#set text(font: "DejaVu Sans", size: 8.8pt)
#set par(leading: 0.5em, spacing: 0.8em)
#set heading(numbering: none)
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    data = json.loads(args.data.read_text())
    if data["schema"] != "groebner.plot-data/v1":
        raise ValueError("Expected native groebner.plot-data/v1 evidence")
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    (out / ".gitattributes").write_text("*.pdf binary\n*.svg -whitespace\n")
    encoded = json.dumps(data, indent=2, ensure_ascii=False) + "\n"
    (out / "plot-data.json").write_text(encoded)
    report = ["# Separate F4/F5 performance history", "",
              f"Frozen rendering of **{data['snapshot_count']} recorded main revisions**, plus the native allocation baseline/candidate audit. No new performance run or speedup claim.", "",
              "Plots retain each algorithm as its own series. CI runtime is diagnostic; incomplete historical profiler coverage is retained as recorded. Connecting lines stop when fingerprints or measurement conditions change. X spacing is revision order, not elapsed time.", "",
              "[Combined PDF](report.pdf) · [Editable source](report.typ) · [Native plot data and provenance](plot-data.json)", ""]
    combined = []
    for group in data["groups"]:
        for chart in group["charts"]:
            draw_chart(group, chart, out / (chart["id"] + ".svg"))
        body = group_body(group, data["snapshot_count"])
        source = out / (group["id"] + ".typ")
        source.write_text(PREAMBLE + body + "\n")
        typst.compile(str(source), output=str(out / (group["id"] + ".pdf")))
        combined.append(body)
        report.extend(["## " + group["title"], "", group["workload"], "",
                       f"[Group PDF]({group['id']}.pdf) · [Editable source]({group['id']}.typ)", ""])
        for chart in group["charts"]:
            report.extend([f"![{chart['title']}]({chart['id']}.svg)", ""])
        columns, rows = table_rows(group)
        report.extend(["| " + " | ".join(columns) + " |",
                       "| " + " | ".join("---" for _ in columns) + " |"])
        report.extend("| " + " | ".join(row) + " |" for row in rows)
        report.append("")
        report.extend("- " + note for note in group["notes"])
        report.append("")
    report.extend(["## Historical provenance", "",
                   "Raw snapshots: `docs/perf/history/` on the `perf-history` branch; frozen copies: `research/f4_f5_history_20261012/records/history/` on the implementation branch. Every SHA-256 below was calculated natively from the original file bytes.", "",
                   "| Snapshot | SHA-256 | Recorded samples |",
                   "| --- | --- | --- |"])
    report.extend(f"| `{s['file']}` | `{s['sha256']}` | {s['wall_samples']} |" for s in data["historical_sources"])
    report.extend(["", "Allocation provenance and hashes are in `plot-data.json`. Millisecond labels are the native conversion `wall_ns / 1,000,000`; no ratios, uncertainty estimates or fitted trends are calculated by the renderer.", "",
                   "Future history regeneration: native `tools/groebner-history`, followed by render-only `scripts/perf/render_groebner_history.py`; the workflow saves SVGs, PDFs, sources and JSON together.", ""])
    (out / "REPORT.md").write_text("\n".join(report))
    combined_source = out / "report.typ"
    combined_source.write_text(PREAMBLE + "\n#pagebreak()\n".join(combined) + "\n")
    typst.compile(str(combined_source), output=str(out / "report.pdf"))
    chart_count = sum(len(group["charts"]) for group in data["groups"])
    print(f"Rendered {len(data['groups'])} groups, {chart_count} SVGs and {len(data['groups']) + 1} PDFs into {out}")


if __name__ == "__main__":
    main()

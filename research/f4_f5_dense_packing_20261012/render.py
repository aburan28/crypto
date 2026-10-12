#!/usr/bin/env python3
"""Lay out recorded Rust results; execute no benchmark or result computation."""
import argparse
import json
import os
from pathlib import Path
import tempfile

with tempfile.TemporaryDirectory(prefix="f4-f5-packing-render-") as cache:
    os.environ.setdefault("MPLCONFIGDIR", cache)
    os.environ.setdefault("XDG_CACHE_HOME", cache)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import typst

    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--graphs-only", action="store_true")
    args = parser.parse_args()
    study = Path(__file__).resolve().parent
    data = json.loads((study / "records/comparison.json").read_text())
    candidate = (study / "records/candidate-revision.txt").read_text().strip()[:8]
    prior = json.loads((study.parent / "f4_f5_hash_20261012/records/comparison.json").read_text())
    fused = json.loads((study.parent / "f4_f5_fused_readback_20261012/records/attempt-2-fixed/comparison.json").read_text())
    args.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.hashsalt": "f4-f5-dense-packing-20261012",
                         "axes.spines.top": False, "axes.spines.right": False})
    colors = {"f4": "#326eab", "f5": "#d17b35", "criterion": "#7756a5"}
    run_id = "run-20261012T1503Z"
    arms = {row["arm"]: row for row in data["arms"]}

    def save(fig, name):
        fig.savefig(args.out / f"{name}.svg", metadata={"Date": None})
        fig.savefig(args.out / f"{name}.png", dpi=150, metadata={"Software": "matplotlib"})
        plt.close(fig)

    def history(name, selected, title):
        fig, ax = plt.subplots(figsize=(9.5, 4.6))
        for episode, positions in [(prior, [0, 1]), (fused, [3, 4]), (data, [6, 7])]:
            lookup = {row["arm"]: row for row in episode["arms"]}
            for arm in selected:
                row = lookup[arm]
                ax.plot(positions, [row["baseline_median_ir"], row["candidate_median_ir"]],
                        marker="o", linewidth=2, color=colors[arm], label=arm.upper() if episode is prior else None)
                for x, field in zip(positions, ["baseline_ir", "candidate_ir"]):
                    ax.scatter([x] * 5, row[field], s=17, color=colors[arm], zorder=3)
        ax.set_xticks([0, 1, 3, 4, 6, 7], ["baa62421\nbaseline", "0505fb45\nhashing",
                     "c2963dcf\nbaseline", "d0a1d2d8\nindex/fusion", "0677b1eb\nbaseline", candidate + "\npacking"])
        ax.set_xlim(-0.3, 7.3)
        ax.set_ylabel("Instruction references (Ir) / 11-system suite")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(6, 6))
        ax.set_title(title, loc="left", pad=14, fontweight="bold")
        ax.grid(axis="y", alpha=0.2)
        ax.legend()
        fig.text(0.06, 0.025, "Three independently matched episodes; gaps do not imply cumulative gains.\n"
                 "Each endpoint: 5 native paired samples and median; complete merged-worker dumps. " + run_id, fontsize=9)
        fig.subplots_adjust(left=0.12, bottom=0.25, right=0.97, top=0.86)
        save(fig, name)

    history("complete-instructions", ["f4", "f5"], "Complete Boolean F4/F5 calls: measured instruction history")
    history("criterion-instructions", ["criterion"], "Standalone F5 criterion: diagnostic instruction history")
    keys = [(case, seed) for case in ["n8_d3", "n8_d4", "n10_d4", "n10_d5"]
            for seed in [0x2545f4914f6cdd1d, 0x9e3779b97f4a7c15]]
    keys += [(case, 0xd1b54a32d192ed03) for case in ["holdout_n12_d3", "holdout_n12_d4", "holdout_shifted_d4"]]
    labels = ["8/d3/s1", "8/d3/s2", "8/d4/s1", "8/d4/s2", "10/d4/s1", "10/d4/s2",
              "10/d5/s1", "10/d5/s2", "H12/d3", "H12/d4", "Hshift/d4"]

    def panels(name, field, title, unit, cases=False):
        fig, axes = plt.subplots(3, 1, figsize=(9.5, 6.9), sharex=True)
        for ax, arm in zip(axes, ["f4", "f5", "criterion"]):
            lookup = {(row["case"], row["seed"]): row for row in data["case_verdicts" if cases else "memory_cases"] if row["arm"] == arm}
            rows = [lookup[key] for key in keys]
            if cases:
                ax.axhline(0, color="#777777", linewidth=0.7)
                ax.plot(range(11), [row["instruction_reduction_percent"] for row in rows],
                        marker="o", linewidth=1.2, color=colors[arm])
                for x, row in enumerate(rows):
                    ax.vlines(x, *row["instruction_reduction_ci95_percent"], color=colors[arm], linewidth=3)
            else:
                for role, marker, color in [("baseline", "o", "#326eab"), ("candidate", "x", "#d17b35")]:
                    ax.plot(range(11), [row[f"{role}_{field}"] for row in rows],
                            marker=marker, linewidth=1.4, color=color, label=role)
                ax.set_yscale("log")
            ax.set_ylabel(arm.upper() + "\n" + unit)
            ax.grid(axis="y", alpha=0.2)
        axes[0].set_title(title, loc="left", pad=14, fontweight="bold")
        if not cases:
            axes[0].legend()
        axes[-1].set_xticks(range(11), labels, rotation=35, ha="right")
        caption = "5 paired rounds; native reduction estimates and 95% bootstrap intervals. Positive = fewer Ir." if cases else "Median of 3 primary calls per case/variant; log axis. Native allocator units, not RSS."
        fig.text(0.06, 0.025, caption + "\nH = known frozen regression fixture; Hshift uses bits 52–63 and 6 multiplier variables.\n"
                 + run_id + " · Source: records/comparison.json.", fontsize=9)
        fig.subplots_adjust(left=0.13, bottom=0.23, right=0.97, top=0.9, hspace=0.18)
        save(fig, name)

    panels("case-instructions", None, "Per-case complete calls and criterion diagnostics", "% fewer Ir", cases=True)
    panels("allocation-requests", "allocations", "Allocation requests in each measured interval", "requests (log)")
    panels("requested-bytes", "allocated_bytes", "Cumulative requested allocation bytes", "bytes (log)")
    panels("peak-live", "peak_additional_live_bytes", "Peak live allocation above that call's entry", "bytes (log)")

    if args.graphs_only:
        raise SystemExit(0)
    validation = json.loads((study / "records/validation.json").read_text())

    url = "https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/"
    preamble = '#set document(date: none)\n#set page(paper: "a4", margin: 16mm)\n#set text(font: "DejaVu Sans", size: 9pt)\n#set par(leading: 0.5em, spacing: 0.65em)\n'
    result = " ".join(f'{arm.upper()}: {arms[arm]["baseline_median_ir"]:,.0f} to {arms[arm]["candidate_median_ir"]:,.0f} Ir; paired reduction {arms[arm]["instruction_reduction_percent"]:.4f}% (95% CI {arms[arm]["instruction_reduction_ci95_percent"][0]:.4f} to {arms[arm]["instruction_reduction_ci95_percent"][1]:.4f}%); {arms[arm]["verdict"].replace("_", " ")}.' for arm in ["f4", "f5", "criterion"])
    guard = f'Preregistered F5/complete-call guard: {data["f5_with_complete_call_guards_met"]}. Memory regression entries: {len(data["memory_regressions"])}. All case verdicts, including criterion regressions, remain in the native comparison.'
    sources = f'#text(size: 8pt)[{run_id}; baseline 0677b1eb; candidate {candidate}; AMD EPYC 9V74 Linux x86-64; Rust 1.99.0; Callgrind 3.24; one Rayon worker; 18 processes / 792 audited calls plus the 33-call hash control. No elapsed-time, process RSS or full-solver gain measured. Sources: #link("{url}records/comparison.json")[native results and raw hashes], #link("{url}PROTOCOL.md")[frozen protocol], #link("{url}REPORT.md")[report], #link("{url}records/validation.json")[local validation].]\n'
    sections = [
        ("complete-instructions", "Complete-call instruction history", result + " " + guard),
        ("criterion-instructions", "Criterion instruction diagnostics", "Criterion construction is measured separately from complete calls. One n8/d4/s2 case records 15.4494% more instructions; the cause is not established. Unchanged source does not establish unchanged generated instruction costs. " + result),
        ("case-instructions", "Per-case results and uncertainty", "All eleven cases and all three arms are shown, including neutral and more-instruction verdicts. The frozen gate requires no complete-call regression, with standalone criterion Ir reported separately. " + guard),
        ("allocation-requests", "Allocation requests", "Measured requests include packing, elimination and readback; setup, warmup, verification and JSON serialization are excluded. All primary and paired memory guards are checked, including criterion intervals. " + guard),
        ("requested-bytes", "Cumulative requested bytes", "These recorded native values include allocation and reallocation requests, not resident memory or allocation duration. Identical reduced outputs, ranks and counted word operations are required. " + guard),
        ("peak-live", "Additional peak live allocation", "The native comparator derives peak minus live-at-entry separately for each interval; inputs and retained warm outputs are outside this additional peak. Absolute peaks and entry bytes remain in raw JSONL. " + guard),
    ]
    pages = []
    for name, title, body in sections:
        page = f'= {title}\n#image("{name}.svg", width: 100%)\n\n{body}\n\n'
        page += f'Local validation: {validation["release_library"]} {validation["focused"]} {validation["examples"]}\n\n' + sources
        source = args.out / f"{name}.typ"
        source.write_text(preamble + page)
        typst.compile(str(source), output=str(args.out / f"{name}.pdf"))
        pages.append(page)
    source = args.out / "report.typ"
    source.write_text(preamble + "\n#pagebreak()\n".join(pages))
    typst.compile(str(source), output=str(args.out / "report.pdf"))

#!/usr/bin/env python3
"""Render recorded native values only; execute no benchmark or comparison."""
import argparse
import json
import os
from pathlib import Path
import tempfile

with tempfile.TemporaryDirectory(prefix="f4-f5-fused-render-") as cache:
    os.environ.setdefault("MPLCONFIGDIR", cache)
    os.environ.setdefault("XDG_CACHE_HOME", cache)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import typst

    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    study = Path(__file__).resolve().parent
    data = json.loads((study / "records/attempt-2-fixed/comparison.json").read_text())
    versions = json.loads((study / "records/attempt-2-fixed/production-revisions.json").read_text())
    prior = json.loads((study.parent / "f4_f5_hash_20261012/records/comparison.json").read_text())
    rejected = json.loads((study / "records/comparison.json").read_text())
    memory_rejection = json.loads((study / "records/attempt-2/comparison.json").read_text())
    validation = json.loads((study / "records/validation.json").read_text())
    args.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.hashsalt": "f4-f5-fused-20261012", "axes.spines.top": False,
                         "axes.spines.right": False})
    arms = {a["arm"]: a for a in data["arms"]}
    revisions = [versions["baseline"][:8] + "\nfresh baseline",
                 versions["candidate"][:8] + "\nindex + fusion"]
    colors = {"f4": "#326eab", "f5": "#d17b35", "criterion": "#7756a5"}

    def instruction_graph(name, selected, title):
        fig, ax = plt.subplots(figsize=(9.2, 4.4))
        for episode, positions in [(prior, [0, 1]), (data, [3, 4])]:
            lookup = {a["arm"]: a for a in episode["arms"]}
            for arm in selected:
                record = lookup[arm]
                ax.plot(positions, [record["baseline_median_ir"], record["candidate_median_ir"]],
                        marker="o", markersize=7, linewidth=2, color=colors[arm],
                        label=arm.upper() if episode is prior else None)
                for x, field in zip(positions, ["baseline_ir", "candidate_ir"]):
                    ax.scatter([x] * 5, record[field], color=colors[arm], s=18, zorder=3)
        ax.set_xticks([0, 1, 3, 4], ["baa62421\nbefore hashing", "0505fb45\nmixed F5 hashing"] + revisions)
        ax.set_xlim(-0.3, 4.3)
        ax.set_ylabel("Instruction references (Ir) / 11-system suite")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(6, 6))
        ax.set_title(title, loc="left", pad=14, fontweight="bold")
        ax.grid(axis="y", alpha=0.2)
        ax.legend(loc="best")
        fig.text(0.06, 0.035, "Each segment: 5 paired rounds, 11 systems, native medians and all samples; merged workers.\n"
                 "run-20261012T0543Z · Gap crosses unrelated main updates; both studies' comparison.json.", fontsize=9)
        fig.subplots_adjust(left=0.12, bottom=0.25, right=0.97, top=0.86)
        fig.savefig(args.out / f"{name}.svg", metadata={"Date": None})
        plt.close(fig)

    instruction_graph("complete-instructions", ["f4", "f5"], "Complete Boolean F4/F5 calls: native instruction history")
    instruction_graph("criterion-instructions", ["criterion"], "F5 criterion construction: native instruction history")

    seed1, seed2 = 0x2545f4914f6cdd1d, 0x9e3779b97f4a7c15
    holdout = 0xd1b54a32d192ed03
    order = [(case, seed) for case in ["n8_d3", "n8_d4", "n10_d4", "n10_d5"]
             for seed in [seed1, seed2]]
    order += [(case, holdout) for case in ["holdout_n12_d3", "holdout_n12_d4", "holdout_shifted_d4"]]
    labels = ["8/d3/s1", "8/d3/s2", "8/d4/s1", "8/d4/s2", "10/d4/s1", "10/d4/s2",
              "10/d5/s1", "10/d5/s2", "H12/d3", "H12/d4", "Hshift/d4"]
    def memory_graph(name, field, unit, title, scope):
        fig, axes = plt.subplots(3, 1, figsize=(9.5, 7), sharex=True)
        for ax, arm in zip(axes, ["f4", "f5", "criterion"]):
            lookup = {(r["case"], r["seed"]): r for r in data["memory_cases"] if r["arm"] == arm}
            records = [lookup[key] for key in order]
            for role, marker, color in [("baseline", "o", "#326eab"), ("candidate", "x", "#d17b35")]:
                ax.plot(range(11), [r[f"{role}_{field}"] for r in records],
                        marker=marker, linewidth=1.5, color=color, label=role)
            ax.set_yscale("log")
            ax.set_ylabel(f"{arm.upper()}\n{unit} (log scale)")
            ax.grid(axis="y", alpha=0.2)
        axes[0].set_title(title, loc="left", fontweight="bold", pad=14)
        axes[0].legend(loc="best")
        axes[-1].set_xticks(range(11), labels, rotation=35, ha="right")
        fig.text(0.06, 0.035, "run-20261012T0543Z · Median of 3 native primary calls per system and variant.\n"
                 "H = frozen holdout; Hshift uses bits 52–63 and 6 multiplier variables.\n"
                 + scope + " Source: records/attempt-2-fixed/comparison.json.", fontsize=9)
        fig.subplots_adjust(left=0.13, bottom=0.23, right=0.97, top=0.9, hspace=0.18)
        fig.savefig(args.out / f"{name}.svg", metadata={"Date": None})
        plt.close(fig)

    memory_graph("allocation-requests", "allocations", "requests",
                 "Allocation requests during complete calls and criterion construction",
                 "Allocator requests; not elapsed time.")
    memory_graph("peak-live", "peak_additional_live_bytes", "bytes",
                 "Peak live allocation above each call's entry",
                 "Same-interval peak minus entry; not RSS or cumulative traffic.")

    fig, ax = plt.subplots(figsize=(9.2, 4.4))
    for record in rejected["arms"]:
        arm = record["arm"]
        if arm == "criterion":
            continue
        ax.plot([0, 1], [record["baseline_median_ir"], record["candidate_median_ir"]],
                marker="o", linewidth=2, color=colors[arm], label=arm.upper())
        for x, field in [(0, "baseline_ir"), (1, "candidate_ir")]:
            ax.scatter([x] * 5, record[field], color=colors[arm], s=18, zorder=3)
    ax.set_xticks([0, 1], ["c2963dcf\nbaseline", "07151ff0\nfirst candidate: rejected"])
    ax.set_xlim(-0.18, 1.18)
    ax.set_ylabel("Instruction references (Ir) / 11-system suite")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(6, 6))
    ax.set_title("First attempt: F4 missed the frozen instruction gate", loc="left", fontweight="bold", pad=14)
    ax.grid(axis="y", alpha=0.2)
    ax.legend(loc="best")
    fig.text(0.06, 0.035, "run-20261012T0543Z · 5 paired rounds. F4 verdict neutral; two-arm acceptance failed.\n"
             "All first-attempt records retained. Sources: records/comparison.json and AMENDMENT-2.md.", fontsize=9)
    fig.subplots_adjust(left=0.12, bottom=0.25, right=0.97, top=0.86)
    fig.savefig(args.out / "rejected-attempt.svg", metadata={"Date": None})
    plt.close(fig)

    fallback = next(r for r in memory_rejection["memory_cases"]
                    if r["case"] == "holdout_shifted_d4" and r["arm"] == "f4")
    fig, ax = plt.subplots(figsize=(9.2, 4.4))
    ax.plot([0, 1], [fallback["baseline_peak_additional_live_bytes"],
                     fallback["candidate_peak_additional_live_bytes"]],
            marker="o", color="#a33d44", linewidth=2, markersize=8)
    ax.set_xticks([0, 1], ["c2963dcf\nbaseline", "07e568b1\nflat enum: rejected"])
    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(54088, 54128)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.set_ylabel("Peak additional live allocation bytes (zoomed axis)")
    ax.set_title("Second attempt: metadata enlarged the high-bit fallback", loc="left", fontweight="bold", pad=14)
    ax.grid(axis="y", alpha=0.2)
    fig.text(0.06, 0.035, "run-20261012T0543Z · Median of 3 primary high-bit F4 calls; confirmed in all 5 pairs.\n"
             "Memory gate failed. Sources: records/attempt-2/comparison.json and native layout-probe.json.", fontsize=9)
    fig.subplots_adjust(left=0.13, bottom=0.25, right=0.97, top=0.86)
    fig.savefig(args.out / "rejected-memory.svg", metadata={"Date": None})
    plt.close(fig)

    source_url = "https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/"
    preamble = '''#set document(date: none)
#set page(paper: "a4", margin: 16mm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(leading: 0.55em, spacing: 0.8em)
#set heading(numbering: none)
'''
    sources = f'''#text(size: 8pt)[Native sources: #link("{source_url}records/attempt-2-fixed/comparison.json")[final comparison and raw-record hashes], #link("{source_url}PROTOCOL.md")[frozen protocol], #link("{source_url}AMENDMENT-2.md")[first rejection], #link("{source_url}METADATA-FIX.md")[bounded correction], #link("{source_url}ACCEPTANCE-CLARIFICATION.md")[gate scope], #link("{source_url}TEST-FIX.md")[fixture repair], #link("{source_url}REPORT.md")[report and validation]. Raw Callgrind dumps and Rust JSONL remain beside the data. Earlier segment: #link("https://github.com/aburan28/crypto/blob/main/research/f4_f5_hash_20261012/records/comparison.json")[hashing study].]
#text(size: 8pt)[Validation: {validation["release_library"]} {validation["focused"]} {validation["history"]} Runtime and process RSS remain unmeasured.]
'''
    def verdict_text(arm):
        record = arms[arm]
        interval = record["instruction_reduction_ci95_percent"]
        return (f'{arm.upper()}: {record["baseline_median_ir"]:,.0f} → {record["candidate_median_ir"]:,.0f} Ir per suite. '
                f'Paired instruction reduction {record["instruction_reduction_percent"]:.2f}% '
                f'(95% bootstrap interval {interval[0]:.3f}% to {interval[1]:.3f}%). '
                f'Native verdict: {record["verdict"].replace("_", " ")}.')

    sections = [
        ("complete-instructions", "Complete Boolean F4/F5 calls",
         verdict_text("f4") + "\n\n" + verdict_text("f5") + "\n\n"
         "Five paired rounds of all eleven frozen systems, one Rayon worker. Both complete-call arms are optimized. "
         "The native verdict requires a 95% bootstrap interval excluding zero and an effect above twice the A/A noise estimate and 2%. "
         "Every reduced output, input hash, pruning report and word-XOR count matched. All complete-call cases have the fewer-instruction verdict. "
         "The clarified complete-call gate passes; the original gate over 33 cases fails on three standalone criterion cases. "
         "ACCEPTANCE-CLARIFICATION records this scope change before final measurement."),
        ("criterion-instructions", "F5 criterion construction",
         verdict_text("criterion") + "\n\n"
         "Criterion production source/rules are unchanged; the standalone aggregate is neutral. "
         "Three cases have more-instruction verdicts: the two n8/d3 seeds (reductions −2.31% and −5.16%) and holdout n12/d4 (−4.07%). "
         "Their cause is not established here, and the original all-arm gate fails. "
         "F5 counts full F4 rows and builds selected rows in one traversal; F4 uses a memory-bounded direct index and sorted readback terms. "
         "Pruning rules, signature order, row limits and counted XOR work are unchanged. "
         "Standalone criterion cost is a phase result; it is already included in complete F5 calls and is not added a second time."),
        ("allocation-requests", "Allocation requests",
         "The graph displays the native median allocation-request counts for each workload and variant, over three primary repetitions. "
         "F4 uses a bounded compact index and avoids a second term vector during readback. F5 avoids the separate full-row counting pass. "
         "Both complete-call arms include the unchanged packing and elimination work. "
         "Raw primary records retain requests and cumulative requested bytes separately. "
         "No memory metric regresses on any input in primary or paired controls."),
        ("peak-live", "Peak live allocation bytes",
         "Native primary records retain three repetitions per variant, all eleven inputs, live bytes at entry and the absolute live-byte peak. "
         "The graph shows native medians of peak minus the entry value over that same allocation interval. "
         "Entry records separate retained inputs and warm outputs from the call's additional peak. "
         "No allocation-request, requested-byte, absolute-peak or additional-live-byte regression occurred in primary or paired controls. "
         "The complete reduced outputs and retained warm inputs are identical, and the criterion allocation control is unchanged. "
         "Process RSS and elapsed runtime were not measured; the host lacks the pressure data needed for controlled timing."),
        ("rejected-attempt", "Rejected first attempt",
         "The first candidate combined ordered F4 readback with default fused F5 construction. "
         "The native comparison records 1.44% fewer F4 instructions, below the frozen 2% verdict threshold, so F4 is neutral. "
         "F5 recorded 11.13% fewer instructions; outputs/accounting matched and memory did not regress. "
         "The frozen two-arm acceptance failed. All profiles, the comparison and controls are retained unchanged. "
         "AMENDMENT-2 preregistered compact F4 column lookup before the second measurement, preserving inputs and gates. "
         "The final attempt reuses the two unchanged baseline profiles and runs fresh candidate primary and paired profiles."),
        ("rejected-memory", "Rejected metadata representation",
         "The second implementation recorded 5.51% fewer F4 instructions and 11.03% fewer F5 instructions, but failed the memory gate. "
         "The high-bit fallback's additional peak grew from 54,104 to 54,112 bytes across all primary and paired calls. "
         "The native probe measured the original index metadata at 48 bytes, the flat three-way enum at 56 bytes, and a nested representation at 48 bytes. "
         "The heap-resident column layout accounts for the observed eight-byte increase. "
         "METADATA-FIX explicitly extended the two-attempt stop rule by one bounded correction and rerun. "
         "It preserved inputs and gates, added a metadata-size regression assertion and retained all failed results. "
         "The superseded full-library run was interrupted, not counted as passed."),
    ]
    pages = []
    for name, title, body in sections:
        page = f'= {title}\n#text(size: 8pt, fill: rgb("#79542b"), weight: "bold")[Native instruction and allocation measurements · run-20261012T0543Z]\n\n'
        page += f'#image("{name}.svg", width: 100%)\n\n{body}\n\n{sources}'
        path = args.out / f"{name}.typ"
        path.write_text(preamble + page)
        typst.compile(str(path), output=str(args.out / f"{name}.pdf"))
        pages.append(page)
    combined = args.out / "report.typ"
    combined.write_text(preamble + "\n#pagebreak()\n".join(pages))
    typst.compile(str(combined), output=str(args.out / "report.pdf"))

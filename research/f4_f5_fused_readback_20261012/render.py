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
    data = json.loads((study / "records/comparison.json").read_text())
    versions = json.loads((study / "records/production-revisions.json").read_text())
    validation = json.loads((study / "records/validation.json").read_text())
    args.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.hashsalt": "f4-f5-fused-20261012", "axes.spines.top": False,
                         "axes.spines.right": False})
    arms = {a["arm"]: a for a in data["arms"]}
    revisions = [versions["baseline"][:8] + "\nproduction baseline",
                 versions["candidate"][:8] + "\nsorted readback + fused rows"]
    colors = {"f4": "#326eab", "f5": "#d17b35", "criterion": "#7756a5"}

    def instruction_graph(name, selected, title):
        fig, ax = plt.subplots(figsize=(9.2, 4.4))
        for arm in selected:
            record = arms[arm]
            ax.plot([0, 1], [record["baseline_median_ir"], record["candidate_median_ir"]],
                    marker="o", markersize=7, linewidth=2, color=colors[arm], label=arm.upper())
            for x, field in [(0, "baseline_ir"), (1, "candidate_ir")]:
                ax.scatter([x] * 5, record[field], color=colors[arm], s=18, zorder=3)
        ax.set_xticks([0, 1], revisions)
        ax.set_xlim(-0.18, 1.18)
        ax.set_ylabel("Instruction references (Ir) / 11-system suite")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(6, 6))
        ax.set_title(title, loc="left", pad=14, fontweight="bold")
        ax.grid(axis="y", alpha=0.2)
        ax.legend(loc="best")
        fig.text(0.06, 0.035, "run-20261012T0543Z · 5 paired A/B/A′ rounds; native medians and all samples.\n"
                 "Instrumented complete calls; no wall-time or RSS claim. Source: records/comparison.json.", fontsize=9)
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
                 + scope + " Sources: records/comparison.json and PROTOCOL.md.", fontsize=9)
        fig.subplots_adjust(left=0.13, bottom=0.23, right=0.97, top=0.9, hspace=0.18)
        fig.savefig(args.out / f"{name}.svg", metadata={"Date": None})
        plt.close(fig)

    memory_graph("allocation-requests", "allocations", "requests",
                 "Allocation requests during complete calls and criterion construction",
                 "Allocator requests; not elapsed time.")
    memory_graph("peak-live", "peak_additional_live_bytes", "bytes",
                 "Peak additional live allocation bytes above each call's entry",
                 "Same-interval peak minus entry; not RSS or cumulative traffic.")

    source_url = "https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/"
    preamble = '''#set document(date: none)
#set page(paper: "a4", margin: 16mm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(leading: 0.55em, spacing: 0.8em)
#set heading(numbering: none)
'''
    sources = f'''#text(size: 8pt)[Native sources: #link("{source_url}records/comparison.json")[comparison and raw-record hashes], #link("{source_url}PROTOCOL.md")[frozen protocol], #link("{source_url}REPORT.md")[result report and validation]. Raw Callgrind dumps and Rust JSONL remain beside the source data. This 11-system suite is separate from historical n=17/n=23 CI diagnostics.]
#text(size: 8pt)[Validation: {validation["release_library"]} {validation["focused"]} Runtime and process RSS remain unmeasured.]
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
         "Every reduced output, input hash, pruning report and word-XOR count matched. No per-workload instruction regression crossed the gate."),
        ("criterion-instructions", "F5 criterion construction",
         verdict_text("criterion") + "\n\n"
         "Criterion construction is the unchanged control arm for this iteration. "
         "F5 now counts all nonempty F4 rows and builds selected rows in one traversal; F4 adopts already sorted readback terms. "
         "Pruning rules, signature order, row limits and counted XOR work are unchanged. "
         "Standalone criterion cost is a phase result; it is already included in complete F5 calls and is not added a second time."),
        ("allocation-requests", "Allocation requests",
         "The graph displays the native median allocation-request counts for each workload and variant, over three primary repetitions. "
         "F4 avoids a second term vector and redundant sorting during readback. F5 avoids the separate full-row counting pass. "
         "Both complete-call arms include the unchanged packing and elimination work. "
         "Raw primary records retain requests and cumulative requested bytes separately. "
         "No memory metric regresses on any input in primary or paired controls."),
        ("peak-live", "Peak live allocation bytes",
         "Native primary records retain three repetitions per variant, all eleven inputs, live bytes at entry and the absolute live-byte peak. "
         "The graph shows native medians of peak minus the entry value over that same allocation interval. "
         "This separates a smaller retained warm criterion from F4's own allocation needs. "
         "No allocation-request, requested-byte, absolute-peak or additional-live-byte regression occurred in primary or paired controls. "
         "The complete reduced outputs and retained warm inputs are identical, and the criterion allocation control is unchanged. "
         "Process RSS and elapsed runtime were not measured; the host lacks the pressure data needed for controlled timing."),
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

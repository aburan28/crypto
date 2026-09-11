#!/usr/bin/env python3
"""Regenerate a receipt-checked Markdown report; optionally plot using matplotlib."""

import argparse
import hashlib
import json
from pathlib import Path
import statistics

import study


def load_verified(directory):
    manifest = json.loads((directory / "manifest.json").read_text())
    summary = json.loads((directory / "summary.json").read_text())
    raw = (directory / "trials.jsonl").read_bytes()
    if hashlib.sha256(raw).hexdigest() != summary["trials_sha256"]:
        raise ValueError("trial ledger hash mismatch")
    source = directory / "study_source.py"
    if not source.exists():
        source = Path(study.__file__)
    if hashlib.sha256(source.read_bytes()).hexdigest() != manifest["source_sha256"]:
        raise ValueError("study source does not match receipt")
    rows = [json.loads(line) for line in raw.splitlines()]
    regenerated = study.summarize(rows, manifest["coverage_budget"])
    if any(regenerated[k] != summary[k] for k in ("coverage", "arithmetic")):
        raise ValueError("summary does not match raw observations")
    return manifest, summary, rows


def ratio_text(value):
    if value is None:
        return "inconclusive"
    return f"{value['ratio']:.3f}× [{value['ci95'][0]:.3f}, {value['ci95'][1]:.3f}]"


def render(directory, plot=False):
    manifest, summary, rows = load_verified(directory)
    lines = ["# Synthetic performance study", "",
             "These are known numerical positive controls, not new cryptanalytic algorithms. "
             "There is no elliptic-curve arithmetic, key input, discrete-log solver, or GPU backend in this suite.", "",
             "## Operation-count distribution", "",
             "Each trial covers every state of a finite ring. The baseline takes nearest-neighbor random steps; "
             "the candidate refreshes to an independent uniform state. This changes the transition law. "
             "It is not a drop-in walk optimization for a collision solver. One transition counts as one operation; "
             "its machine cost need not be equal between variants. The initial state is free in both variants.", "",
             "Intervals are descriptive 95% paired percentile bootstrap intervals (1,000 resamples). "
             "Ratios are baseline/candidate. All capped trials contribute their full budget to the restricted mean "
             "E[min(T, budget)]. The ratio is not an estimate of uncensored mean completion-time speedup.", "",
             "| Split / states | Baseline completed | Candidate completed | Restricted mean ratio [95% CI] |",
             "|---|---:|---:|---:|"]
    for key, value in summary["coverage"].items():
        a, b = value["nearest"], value["refresh"]
        lines.append(f"| {key} | {a['completed']}/{a['trials']} | {b['completed']}/{b['trials']} | "
                     f"{ratio_text(value['restricted_mean_ratio'])} |")
    lines += ["", "## Exact arithmetic and full measured cost", "",
              "The candidate is ordinary Horner evaluation; the baseline sums modular powers. "
              "Both evaluate degree-24 polynomials modulo 257. Every sampled polynomial is checked "
              "at all 257 field elements against an independent unbounded-integer power-sum reference. "
              "Every timed output is also verified, and paired input/output digests must agree.", "",
              "Timing uses fresh child processes, identical fixtures within each pair, and randomized variant order. "
              "Full measured cost includes initialization, host input encoding/decoding, warmup, computation, "
              "output encoding, exhaustive verification, interpreter startup/teardown, and subprocess IPC. "
              "The final run also includes parent receipt parsing. Report generation and final artifact writes "
              "are outside per-trial timing. No GPU, PCIe, or network transfer has been measured.", "",
              "| Split | Verified and paired | Compute ratio [95% CI] | Full wall ratio [95% CI] | Decision |",
              "|---|---|---:|---:|---|"]
    for split, value in summary["arithmetic"].items():
        lines.append(f"| {split} | {value['all_verified_and_paired']} | {ratio_text(value['kernel_ratio'])} | "
                     f"{ratio_text(value['full_wall_ratio'])} | {value['decision']} |")
    lines += ["", "A decision requires every arithmetic timing and memory trial to complete and verify. "
              "Failed, incorrect, or timed-out trials remain in the ledger and block a speedup claim. "
              "Coverage caps are intentionally retained as censored observations. "
              "The holdout fixtures are disjoint from discovery, and variants are fixed before either split. "
              "Confidence intervals describe this process and machine sample, not all hardware or later workloads.", "",
              "## Memory (separate instrumented passes)", "",
              "| Split / variant | Python traced peak bytes | Process peak RSS bytes |",
              "|---|---:|---:|"]
    for split, value in summary["arithmetic"].items():
        for row in value["memory_passes"]:
            memory = row.get("memory", {})
            lines.append(f"| {split}/{row['variant']} | {memory.get('python_traced_peak_bytes')} | "
                         f"{memory.get('process_peak_rss_bytes')} |")
    lines += ["", "Memory timing is excluded from speed ratios. RSS is a process high-water mark including "
              "the interpreter; traced memory covers Python allocations and is not equivalent to RSS. "
              "One memory pass per variant and split is diagnostic, not a statistically established memory improvement.", "",
              "## Provenance", "", f"- Platform: `{manifest['platform']}`",
              f"- Python: `{manifest['python'].splitlines()[0]}`",
              f"- Study source SHA-256: `{manifest['source_sha256']}`",
              f"- Trial ledger SHA-256: `{summary['trials_sha256']}`",
              f"- Coverage trials per split/size/variant: {manifest['trials_per_split_and_size']}",
              f"- Coverage operation budget: {manifest['coverage_budget']}",
              f"- Arithmetic pairs per split: {manifest['arithmetic_pairs_per_split']}",
              f"- Timed evaluations per arithmetic trial: {manifest['arithmetic_count']}",
              f"- Child watchdog: {manifest['child_timeout_s']} seconds", "",
              "The host was not isolated from other workloads. No CPU affinity, frequency lock, or GPU "
              "measurement was used. Startup and scheduling variation may dominate small workloads. "
              "These results support only the stated finite numerical comparisons.", ""]
    if plot:
        plot_results(directory, rows)
        lines += ["![Holdout coverage and cost accounting](holdout.png)", ""]
    (directory / "REPORT.md").write_text("\n".join(lines))


def plot_results(directory, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, states in zip(axes[:2], (32, 128)):
        for variant, color in (("nearest", "#476885"), ("refresh", "#15917d")):
            rs = [r for r in rows if r["study"] == "coverage" and r["split"] == "holdout"
                  and r["states"] == states and r["variant"] == variant]
            times = sorted(r["operations"] for r in rs if r["status"] == "completed")
            budget = max(r["operations"] for r in rs)
            ax.step([0] + times + [budget], [0] + [(i + 1) / len(rs) for i in range(len(times))]
                    + [len(times) / len(rs)], where="post", label=variant, color=color)
        ax.set(title=f"{states}-state ring: holdout", xlabel="Transitions", ylabel="Fraction completed", ylim=(0, 1.02))
        ax.grid(alpha=.15)
        ax.legend()
    timing = {v: [r for r in rows if r["study"] == "arithmetic" and r["split"] == "holdout"
                  and not r["memory_pass"] and r["variant"] == v] for v in study.VARIANTS}
    if all(r["status"] == "completed" for rs in timing.values() for r in rs):
        categories = {"Compute": lambda r: r["phases_ns"]["compute"],
                      "Verification": lambda r: r["phases_ns"]["verify"],
                      "Setup + host encoding": lambda r: sum(v for k, v in r["phases_ns"].items()
                                                              if k not in ("compute", "verify")),
                      "Process + other": lambda r: r["process_overhead_ns"] + r["unattributed_worker_ns"]}
        bottoms = [0., 0.]
        for (label, fn), color in zip(categories.items(), ("#15917d", "#476885", "#e4b04e", "#bcc4cb")):
            heights = [statistics.mean(fn(r) for r in timing[v]) / 1e6 for v in study.VARIANTS]
            axes[2].bar(study.VARIANTS, heights, bottom=bottoms, label=label, color=color)
            bottoms = [a + b for a, b in zip(bottoms, heights)]
        axes[2].set(title="Arithmetic: mean full measured cost", ylabel="Milliseconds per child process")
        axes[2].legend(fontsize=8)
    fig.suptitle("Synthetic positive controls — no cryptanalytic speedup claim", fontsize=13)
    fig.tight_layout()
    fig.savefig(directory / "holdout.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    render(args.directory, args.plot)

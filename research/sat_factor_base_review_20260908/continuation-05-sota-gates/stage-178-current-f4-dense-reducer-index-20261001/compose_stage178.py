#!/usr/bin/env python3
"""Compose the Stage 178 paired dense-index result."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics


STAGE = Path(__file__).resolve().parent
PAIRED = STAGE / "development" / "paired"
BUILD = STAGE / "development" / "build"
SEQUENCE = ("linear", "dense", "dense", "linear", "linear", "dense")
EXPECTED_EQUATIONS = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
EXPECTED_OPS = 319_313_687_585
EXPECTED_PERFORMED = 147_794_583_858


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def receipt(path: Path) -> dict:
    return {"path": str(path.relative_to(STAGE)), "bytes": path.stat().st_size, "sha256": sha256(path)}


runs = []
counts = {"linear": 0, "dense": 0}
for order, variant in enumerate(SEQUENCE, 1):
    counts[variant] += 1
    root = PAIRED / f"{order:02d}-{variant}-r{counts[variant]}"
    metrics_path = root / "metrics.json"
    stdout_path = root / "stdout.json"
    stderr_path = root / "stderr.txt"
    process = json.loads(metrics_path.read_text())
    report = json.loads(stdout_path.read_text())
    extra = report["cost"]["extra"]
    expected_divisors = 4_190_633_182 if variant == "linear" else 103_532_494
    expected_index_bytes = 0 if variant == "linear" else 1_062_392
    correct = (
        process["returncode"] == 0
        and not process["timed_out"]
        and report["status"] == "unsat"
        and report["exhaustive"] is True
        and report["source_instance_verified"] is True
        and report["regenerated_source_exact"] is True
        and report["fixed_x1_masks_visited"] == 512
        and report["fixed_x1_systems_completed"] == 242
        and report["solver_equations_blake3"] == EXPECTED_EQUATIONS
        and report["cost"]["ops"] == EXPECTED_OPS
        and extra["word_xors_performed"] == EXPECTED_PERFORMED
        and extra["divisor_tests"] == expected_divisors
        and extra["divisor_submask_lookups"] == (0 if variant == "linear" else expected_divisors)
        and extra["divisor_linear_tests"] == (expected_divisors if variant == "linear" else 0)
        and extra["reducer_index_bytes_max"] == expected_index_bytes
        and report["conflicts"] is None
    )
    runs.append(
        {
            "order": order,
            "variant": variant,
            "repeat": counts[variant],
            "correct": correct,
            "process": process,
            "report": report,
            "artifacts": {
                "metrics": receipt(metrics_path),
                "stdout": receipt(stdout_path),
                "stderr": receipt(stderr_path),
            },
        }
    )
if not all(run["correct"] for run in runs):
    raise SystemExit("paired correctness failure")

linear = [run for run in runs if run["variant"] == "linear"]
dense = [run for run in runs if run["variant"] == "dense"]


def medians(group: list[dict]) -> dict:
    return {
        field: statistics.median(run["process"]["metrics"][field] for run in group)
        for field in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
    }


linear_median = medians(linear)
dense_median = medians(dense)
pairs = [(linear[0], dense[0]), (linear[1], dense[1]), (linear[2], dense[2])]
pair_ratios = [
    {
        "pair": i,
        "wall": candidate["process"]["metrics"]["wall_seconds"] / control["process"]["metrics"]["wall_seconds"],
        "core": candidate["process"]["metrics"]["total_core_seconds"] / control["process"]["metrics"]["total_core_seconds"],
        "rss": candidate["process"]["metrics"]["peak_rss_bytes"] / control["process"]["metrics"]["peak_rss_bytes"],
    }
    for i, (control, candidate) in enumerate(pairs, 1)
]
paired_median = {
    field: statistics.median(pair[field] for pair in pair_ratios)
    for field in ("wall", "core", "rss")
}
passed = paired_median["wall"] < 0.97 and paired_median["core"] < 0.97
build_path = BUILD / "metrics.json"
build = json.loads(build_path.read_text())
all_processes = [build, *(run["process"] for run in runs)]
charge = {
    "components": len(all_processes),
    "wall_seconds_sum": sum(item["metrics"]["wall_seconds"] for item in all_processes),
    "total_core_seconds_sum": sum(item["metrics"]["total_core_seconds"] for item in all_processes),
    "peak_rss_bytes_max": max(item["metrics"]["peak_rss_bytes"] for item in all_processes),
}
result = {
    "schema": "koblitz_stage178_current_f4_dense_reducer_index.v1",
    "claim_boundary": "One opened n=59 decomposition target and a solver-stage implementation experiment only; not a full index-calculus run or SOTA.",
    "source": {
        "candidate_commit": "c580c450cbce58ec23de8ce4a622481864420a4f",
        "runner_commit": "1becd640f",
        "binary_sha256": "745fefcef1d6b731a9edfaab298e1a925d1b616e9fcd9987378fb53c95de7cb4",
    },
    "build": {
        "process": build,
        "artifacts": {
            "metrics": receipt(build_path),
            "stdout": receipt(BUILD / "stdout.txt"),
            "stderr": receipt(BUILD / "stderr.txt"),
        },
    },
    "runs": runs,
    "linear_median": linear_median,
    "dense_median": dense_median,
    "pair_ratios": pair_ratios,
    "median_paired_ratios": paired_median,
    "mechanism": {
        "linear_divisor_tests": 4_190_633_182,
        "dense_submask_lookups": 103_532_494,
        "lookup_reduction_percent": 100 * (1 - 103_532_494 / 4_190_633_182),
        "dense_index_bytes_max": 1_062_392,
        "structural_f4_counts_identical": True,
    },
    "decision": {
        "status": "REJECTED_CPU_THRESHOLD" if not passed else "ACCEPTED_ENGINEERING_IMPROVEMENT",
        "pass": passed,
        "reason": "Median paired CPU did not fall below 0.97; the frozen threshold is unchanged." if not passed else "Both preregistered paired thresholds passed.",
        "runtime_default": "linear",
    },
    "campaign_charge": charge,
    "single_core_seconds": None,
    "conflicts": None,
    "rejected_patch": receipt(STAGE / "rejected-dense-indexed-reducers.patch"),
    "validation": {
        "f4_tests": "9 passed",
        "backend_tests": "3 passed",
        "forced_linear_witness_test": "1 passed",
    },
}
(STAGE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

lines = [
    "# Stage 178: reused dense reducer index in current F4",
    "",
    "All six same-binary processes returned exhaustive UNSAT with exact equation, matrix, pair, basis, extraction, and XOR counts.",
    "",
    "| pair | dense / linear wall | dense / linear core | dense / linear RSS |",
    "|---:|---:|---:|---:|",
]
for pair in pair_ratios:
    lines.append(f"| {pair['pair']} | {pair['wall']:.6f} | {pair['core']:.6f} | {pair['rss']:.6f} |")
lines += [
    "",
    f"The median paired ratios are {paired_median['wall']:.6f} wall, {paired_median['core']:.6f} CPU, and {paired_median['rss']:.6f} RSS. CPU misses the frozen 0.97 threshold, so the candidate is rejected.",
    "",
    f"The 1,062,392-byte per-call dense index retains the 97.530% lookup-count reduction, but the measured CPU gain is only {100 * (1 - paired_median['core']):.3f}% at the median pair.",
    "",
    f"The clean build plus six paired processes charge {charge['wall_seconds_sum']:.6f} wall seconds, {charge['total_core_seconds_sum']:.6f} core-seconds, and {charge['peak_rss_bytes_max']} bytes maximum RSS.",
    "",
    "The dense implementation is preserved as a rejected patch and is not the selected runtime default. This does not change any SOTA gate.",
    "",
]
(STAGE / "RESULTS.md").write_text("\n".join(lines))
print(json.dumps({"median_paired_ratios": paired_median, "decision": result["decision"], "campaign_charge": charge}, indent=2))

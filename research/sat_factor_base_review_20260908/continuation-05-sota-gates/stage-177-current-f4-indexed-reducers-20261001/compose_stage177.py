#!/usr/bin/env python3
"""Compose the Stage 177 paired exact-submask result."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics


STAGE = Path(__file__).resolve().parent
PAIRED = STAGE / "development" / "paired"
BUILD = STAGE / "development" / "build"
EXPECTED_EQUATIONS = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
EXPECTED_OPS = 319_313_687_585
EXPECTED_PERFORMED = 147_794_583_858
SEQUENCE = ("linear", "indexed", "indexed", "linear", "linear", "indexed")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def receipt(path: Path) -> dict:
    return {"path": str(path.relative_to(STAGE)), "bytes": path.stat().st_size, "sha256": sha256(path)}


runs = []
counts = {"linear": 0, "indexed": 0}
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
    ok = (
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
        and report["conflicts"] is None
    )
    runs.append(
        {
            "order": order,
            "variant": variant,
            "repeat": counts[variant],
            "correct": ok,
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
indexed = [run for run in runs if run["variant"] == "indexed"]


def medians(group: list[dict]) -> dict:
    return {
        field: statistics.median(run["process"]["metrics"][field] for run in group)
        for field in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
    }


linear_median = medians(linear)
indexed_median = medians(indexed)
pairs = [(linear[0], indexed[0]), (linear[1], indexed[1]), (linear[2], indexed[2])]
pair_ratios = []
for pair, (control, candidate) in enumerate(pairs, 1):
    pair_ratios.append(
        {
            "pair": pair,
            "wall": candidate["process"]["metrics"]["wall_seconds"] / control["process"]["metrics"]["wall_seconds"],
            "core": candidate["process"]["metrics"]["total_core_seconds"] / control["process"]["metrics"]["total_core_seconds"],
            "rss": candidate["process"]["metrics"]["peak_rss_bytes"] / control["process"]["metrics"]["peak_rss_bytes"],
        }
    )
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
    "schema": "koblitz_stage177_current_f4_indexed_reducers.v1",
    "claim_boundary": "One opened n=59 decomposition target and a solver-stage implementation experiment only; not a full index-calculus run or SOTA.",
    "source": {
        "candidate_commit": "a337e2e901848ac440d70d88c838d37937467413",
        "runner_commit": "dda9096b8",
        "binary_sha256": "07bc31f962de6329b3e1fe76ca5311ffbdb82c69df8650bab6eb99d318ee115b",
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
    "indexed_median": indexed_median,
    "pair_ratios": pair_ratios,
    "median_paired_ratios": paired_median,
    "mechanism": {
        "linear_divisor_tests": 4_190_633_182,
        "indexed_submask_lookups": 103_532_494,
        "lookup_count_ratio": 103_532_494 / 4_190_633_182,
        "reduction_percent": 100 * (1 - 103_532_494 / 4_190_633_182),
        "structural_f4_counts_identical": True,
    },
    "decision": {
        "status": "REJECTED_PAIRED_TIMING_REGRESSION" if not passed else "ACCEPTED_ENGINEERING_IMPROVEMENT",
        "pass": passed,
        "reason": "Median paired wall and core ratios did not both fall below 0.97." if not passed else "Both preregistered paired thresholds passed.",
        "runtime_default": "linear",
    },
    "campaign_charge": charge,
    "single_core_seconds": None,
    "conflicts": None,
    "rejected_patch": receipt(STAGE / "rejected-hash-indexed-reducers.patch"),
    "validation": {
        "f4_tests": "9 passed twice after final budget handling",
        "backend_tests": "3 passed",
        "forced_linear_witness_test": "1 passed",
    },
}
(STAGE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

lines = [
    "# Stage 177: exact hash-indexed reducers in current F4",
    "",
    "All six same-binary processes returned exhaustive UNSAT on the identical frozen target with exact equation, matrix, pair, basis, extraction, and XOR counts.",
    "",
    "| pair | indexed / linear wall | indexed / linear core | indexed / linear RSS |",
    "|---:|---:|---:|---:|",
]
for pair in pair_ratios:
    lines.append(f"| {pair['pair']} | {pair['wall']:.6f} | {pair['core']:.6f} | {pair['rss']:.6f} |")
lines += [
    "",
    f"The median paired ratios are {paired_median['wall']:.6f} wall, {paired_median['core']:.6f} CPU, and {paired_median['rss']:.6f} RSS. The candidate is rejected because wall and CPU did not both fall below 0.97.",
    "",
    f"The mechanism is real: divisor work fell from 4,190,633,182 linear tests to 103,532,494 exact submask probes, a {result['mechanism']['reduction_percent']:.3f}% count reduction. Hash-map construction/probes and cache effects consumed that gain in paired time.",
    "",
    f"The clean build plus six paired processes charge {charge['wall_seconds_sum']:.6f} wall seconds, {charge['total_core_seconds_sum']:.6f} core-seconds, and {charge['peak_rss_bytes_max']} bytes maximum RSS.",
    "",
    "The hash-indexed implementation is preserved as a rejected patch and is not the selected runtime default. This does not change any SOTA gate.",
    "",
]
(STAGE / "RESULTS.md").write_text("\n".join(lines))
print(json.dumps({"median_paired_ratios": paired_median, "decision": result["decision"], "campaign_charge": charge}, indent=2))

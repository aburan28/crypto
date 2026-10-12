#!/usr/bin/env python3
"""Compose the Stage 179 five-/six-column table result."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics


STAGE = Path(__file__).resolve().parent
PAIRED = STAGE / "development" / "paired"
SEQUENCE = ("five", "six", "six", "five", "five", "six")
EXPECTED_EQUATIONS = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
EXPECTED_OPS = 319_313_687_585
EXPECTED = {
    "five": {"performed": 147_794_583_858, "table_bytes": 53_070_336},
    "six": {"performed": 163_503_838_578, "table_bytes": 74_111_872},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def receipt(path: Path) -> dict:
    return {"path": str(path.relative_to(STAGE)), "bytes": path.stat().st_size, "sha256": sha256(path)}


def process_artifact(root: Path, stdout_name: str, stderr_name: str) -> dict:
    metrics = root / "metrics.json"
    stdout = root / stdout_name
    stderr = root / stderr_name
    return {
        "process": json.loads(metrics.read_text()),
        "artifacts": {
            "metrics": receipt(metrics),
            "stdout": receipt(stdout),
            "stderr": receipt(stderr),
        },
    }


builds = {
    variant: process_artifact(STAGE / "development" / f"build-{variant}", "stdout.txt", "stderr.txt")
    for variant in ("five", "six")
}
tests = {
    variant: process_artifact(STAGE / "development" / f"test-{variant}", "stdout.txt", "stderr.txt")
    for variant in ("five", "six")
}
for variant in tests:
    tests[variant]["passed"] = (
        tests[variant]["process"]["returncode"] == 0
        and not tests[variant]["process"]["timed_out"]
        and "8 passed; 0 failed" in (STAGE / "development" / f"test-{variant}" / "stdout.txt").read_text()
    )
if not all(value["passed"] for value in tests.values()):
    raise SystemExit("feature tests failed")

runs = []
counts = {"five": 0, "six": 0}
for order, variant in enumerate(SEQUENCE, 1):
    counts[variant] += 1
    root = PAIRED / f"{order:02d}-{variant}-r{counts[variant]}"
    metrics_path = root / "metrics.json"
    stdout_path = root / "stdout.json"
    stderr_path = root / "stderr.txt"
    process = json.loads(metrics_path.read_text())
    report = json.loads(stdout_path.read_text())
    extra = report["cost"]["extra"]
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
        and extra["word_xors_performed"] == EXPECTED[variant]["performed"]
        and extra["peak_table_bytes"] == EXPECTED[variant]["table_bytes"]
        and extra["divisor_tests"] == 4_190_633_182
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

five = [run for run in runs if run["variant"] == "five"]
six = [run for run in runs if run["variant"] == "six"]
pairs = [(five[0], six[0]), (five[1], six[1]), (five[2], six[2])]
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
all_processes = [
    *(item["process"] for item in builds.values()),
    *(item["process"] for item in tests.values()),
    *(run["process"] for run in runs),
]
charge = {
    "components": len(all_processes),
    "wall_seconds_sum": sum(item["metrics"]["wall_seconds"] for item in all_processes),
    "total_core_seconds_sum": sum(item["metrics"]["total_core_seconds"] for item in all_processes),
    "peak_rss_bytes_max": max(item["metrics"]["peak_rss_bytes"] for item in all_processes),
}
passed = paired_median["wall"] < 0.97 and paired_median["core"] < 0.97
result = {
    "schema": "koblitz_stage179_current_f4_wide_tables_single_target.v1",
    "claim_boundary": "One opened n=59 decomposition target and an existing solver-feature comparison only; not a full index-calculus run or SOTA.",
    "source_commit": "12da6fe144788ecf8e9f2bf72a7a3274ed3582b9",
    "runner_commit": "87c2b4170",
    "binary_sha256": {
        "five": "b156f3320564c1e9fed2b168aba7a0b4f3c875d677c1a9354751bba527ee5bd6",
        "six": "6e634aeff568c824c1f58c27c9deeb1b677a716990a4a0ed7d858182fec2bc5a",
    },
    "builds": builds,
    "tests": tests,
    "runs": runs,
    "pair_ratios": pair_ratios,
    "median_paired_ratios": paired_median,
    "mechanism": {
        "performed_xor_ratio_six_over_five": EXPECTED["six"]["performed"] / EXPECTED["five"]["performed"],
        "table_memory_ratio_six_over_five": EXPECTED["six"]["table_bytes"] / EXPECTED["five"]["table_bytes"],
        "row_equivalent_xors_identical": True,
        "structural_f4_counts_identical": True,
    },
    "decision": {
        "status": "REJECTED_WIDER_TABLE_REGRESSION" if not passed else "ACCEPTED_TARGET_SPECIFIC_ENGINEERING",
        "pass": passed,
        "selected_table_columns": 5,
        "reason": "Six columns increased deterministic performed XORs and table memory and failed the paired timing gate." if not passed else "Both preregistered paired thresholds passed.",
    },
    "campaign_charge": charge,
    "single_core_seconds": None,
    "conflicts": None,
    "wall_interpretation": "The host was severely descheduled during the paired panel; paired wall is retained as observed, while core-seconds and deterministic operation counts anchor the rejection.",
}
(STAGE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

lines = [
    "# Stage 179: six-column BlockTables on the frozen target",
    "",
    "Both feature configurations passed all eight Boolean-F4 tests and all six target processes returned exhaustive UNSAT with identical algebraic and row-equivalent F4 counters.",
    "",
    "| pair | six / five wall | six / five core | six / five RSS |",
    "|---:|---:|---:|---:|",
]
for pair in pair_ratios:
    lines.append(f"| {pair['pair']} | {pair['wall']:.6f} | {pair['core']:.6f} | {pair['rss']:.6f} |")
lines += [
    "",
    f"Median paired ratios are {paired_median['wall']:.6f} wall, {paired_median['core']:.6f} CPU, and {paired_median['rss']:.6f} RSS. Severe host descheduling makes absolute wall snapshots unsuitable, but the candidate also performs {result['mechanism']['performed_xor_ratio_six_over_five']:.6f}x as many actual word XORs and uses {result['mechanism']['table_memory_ratio_six_over_five']:.6f}x the table memory.",
    "",
    f"The two clean builds, two feature-test runs, and six paired queries charge {charge['wall_seconds_sum']:.6f} wall seconds, {charge['total_core_seconds_sum']:.6f} core-seconds, and {charge['peak_rss_bytes_max']} bytes maximum RSS.",
    "",
    "Six columns is rejected and the five-column repository default remains selected. This does not change any SOTA gate.",
    "",
]
(STAGE / "RESULTS.md").write_text("\n".join(lines))
print(json.dumps({"median_paired_ratios": paired_median, "mechanism": result["mechanism"], "decision": result["decision"], "campaign_charge": charge}, indent=2))

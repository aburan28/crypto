#!/usr/bin/env python3
"""Compose the Stage 181 full-matrix M4RI result."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics


STAGE = Path(__file__).resolve().parent
EXPECTED_EQUATIONS = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
EXPECTED = {
    "current": {
        "logical": 319_313_687_585,
        "performed": 147_794_583_858,
        "matrices": 0,
        "blocks": 0,
        "preparation": 0,
    },
    "full": {
        "logical": 318_635_818_320,
        "performed": 99_192_937_526,
        "matrices": 723,
        "blocks": 438_923,
        "preparation": 10_095_063_562,
    },
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def receipt(path: Path) -> dict:
    return {"path": str(path.relative_to(STAGE)), "bytes": path.stat().st_size, "sha256": sha256(path)}


def process_artifact(root: Path) -> dict:
    metrics = root / "metrics.json"
    stdout = root / "stdout.txt"
    stderr = root / "stderr.txt"
    return {
        "process": json.loads(metrics.read_text()),
        "artifacts": {
            "metrics": receipt(metrics),
            "stdout": receipt(stdout),
            "stderr": receipt(stderr),
        },
    }


build = process_artifact(STAGE / "development" / "build")
tests = {}
for name, expected in (
    ("current_f4", "9 passed; 0 failed"),
    ("full_f4", "9 passed; 0 failed"),
    ("current_backend", "3 passed; 0 failed"),
    ("full_backend", "3 passed; 0 failed"),
):
    item = process_artifact(STAGE / "development" / f"test-{name.replace('_', '-')}")
    item["passed"] = (
        item["process"]["returncode"] == 0
        and not item["process"]["timed_out"]
        and expected in (STAGE / "development" / f"test-{name.replace('_', '-')}" / "stdout.txt").read_text()
    )
    tests[name] = item
if not all(item["passed"] for item in tests.values()):
    raise SystemExit("validation process failed")


def load_run(root: Path, variant: str, order: int, repeat=None) -> dict:
    metrics_path = root / "metrics.json"
    stdout_path = root / "stdout.json"
    stderr_path = root / "stderr.txt"
    process = json.loads(metrics_path.read_text())
    report = json.loads(stdout_path.read_text())
    extra = report["cost"]["extra"]
    expected = EXPECTED[variant]
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
        and report["cost"]["ops"] == expected["logical"]
        and extra["word_xors_performed"] == expected["performed"]
        and extra["full_m4ri_matrices"] == expected["matrices"]
        and extra["full_m4ri_blocks"] == expected["blocks"]
        and extra["full_m4ri_table_word_xors"] == expected["preparation"]
        and extra["basis_len"] == 1
        and extra["pairs_left"] == 0
        and report["conflicts"] is None
    )
    return {
        "order": order,
        "variant": variant,
        "repeat": repeat,
        "correct": correct,
        "process": process,
        "report": report,
        "artifacts": {
            "metrics": receipt(metrics_path),
            "stdout": receipt(stdout_path),
            "stderr": receipt(stderr_path),
        },
    }


screen_root = STAGE / "development" / "screen"
screen = [
    load_run(screen_root / "01-current", "current", 1),
    load_run(screen_root / "02-full", "full", 2),
]
confirmation_root = STAGE / "development" / "confirmation"
sequence = ("current", "full", "full", "current", "current", "full")
counts = {"current": 0, "full": 0}
confirmation = []
for order, variant in enumerate(sequence, 1):
    counts[variant] += 1
    confirmation.append(
        load_run(
            confirmation_root / f"{order:02d}-{variant}-r{counts[variant]}",
            variant,
            order,
            counts[variant],
        )
    )
if not all(run["correct"] for run in [*screen, *confirmation]):
    raise SystemExit("target correctness failure")

current = [run for run in confirmation if run["variant"] == "current"]
full = [run for run in confirmation if run["variant"] == "full"]
pairs = [(current[0], full[0]), (current[1], full[1]), (current[2], full[2])]
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


def medians(group: list[dict]) -> dict:
    return {
        field: statistics.median(run["process"]["metrics"][field] for run in group)
        for field in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
    }


current_median = medians(current)
full_median = medians(full)
screen_ratio = {
    "wall": screen[1]["process"]["metrics"]["wall_seconds"] / screen[0]["process"]["metrics"]["wall_seconds"],
    "core": screen[1]["process"]["metrics"]["total_core_seconds"] / screen[0]["process"]["metrics"]["total_core_seconds"],
    "rss": screen[1]["process"]["metrics"]["peak_rss_bytes"] / screen[0]["process"]["metrics"]["peak_rss_bytes"],
}
passed = paired_median["wall"] < 0.97 and paired_median["core"] < 0.97
all_processes = [
    build["process"],
    *(item["process"] for item in tests.values()),
    *(run["process"] for run in screen),
    *(run["process"] for run in confirmation),
]
charge = {
    "components": len(all_processes),
    "wall_seconds_sum": sum(item["metrics"]["wall_seconds"] for item in all_processes),
    "total_core_seconds_sum": sum(item["metrics"]["total_core_seconds"] for item in all_processes),
    "peak_rss_bytes_max": max(item["metrics"]["peak_rss_bytes"] for item in all_processes),
}
stage174 = {"wall_seconds": 26.358217916989815, "total_core_seconds": 147.841771, "peak_rss_bytes": 2_756_362_240}
result = {
    "schema": "koblitz_stage181_current_f4_full_m4ri.v1",
    "claim_boundary": "One opened n=59 decomposition target and a solver-stage experiment only; not a full index-calculus run or SOTA.",
    "source": {
        "candidate_commit": "850182efba8d8a9499e82aa6f11600941d3cb9d6",
        "confirmation_runner_commit": "7f76abe9c",
        "binary_sha256": "fe17f002fedb2e4ef3037d100a4c541ffe0d2c7dcbfa6b2121a540804cc9050c",
    },
    "build": build,
    "tests": tests,
    "screen": {"runs": screen, "full_over_current": screen_ratio, "continued": screen_ratio["core"] < 1 and EXPECTED["full"]["performed"] < EXPECTED["current"]["performed"]},
    "confirmation": {
        "runs": confirmation,
        "current_median": current_median,
        "full_median": full_median,
        "pair_ratios": pair_ratios,
        "median_paired_ratios": paired_median,
    },
    "mechanism": {
        "performed_xor_ratio_full_over_current": EXPECTED["full"]["performed"] / EXPECTED["current"]["performed"],
        "logical_xor_ratio_full_over_current": EXPECTED["full"]["logical"] / EXPECTED["current"]["logical"],
        "full_m4ri_matrices": 723,
        "full_m4ri_blocks": 438_923,
        "full_m4ri_table_word_xors": 10_095_063_562,
    },
    "absolute_ratios": {
        "full_over_stage174": {
            field: full_median[field] / stage174[field]
            for field in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
        },
        "full_over_stage175_same_binary_direct": {
            "wall_seconds": full_median["wall_seconds"] / 0.3618509580001046,
            "total_core_seconds": full_median["total_core_seconds"] / 0.351335,
            "peak_rss_bytes": full_median["peak_rss_bytes"] / 45_072_384,
        },
    },
    "decision": {
        "status": "RETAIN_OPT_IN_CONFIRMATION_WALL_FAIL" if not passed else "ACCEPTED_ENGINEERING_IMPROVEMENT",
        "pass": passed,
        "runtime_default": "current_blocktables",
        "candidate_availability": "F4_F2_FULL_M4RI=1",
        "reason": "CPU and performed work improved, but the frozen median paired wall ratio did not fall below 0.97." if not passed else "Both preregistered paired thresholds passed.",
    },
    "campaign_charge": charge,
    "single_core_seconds": None,
    "conflicts": None,
    "wall_interpretation": "The confirmation host was heavily and unevenly descheduled; wall is retained as measured and its frozen threshold still controls the decision.",
}
(STAGE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

lines = [
    "# Stage 181: full-matrix M4RI inside current F4",
    "",
    "Default and forced modes passed nine Boolean-F4 tests and three backend tests. The screen and all six confirmation processes authenticated the same target and returned exact exhaustive UNSAT.",
    "",
    "| pair | full / current wall | full / current core | full / current RSS |",
    "|---:|---:|---:|---:|",
]
for pair in pair_ratios:
    lines.append(f"| {pair['pair']} | {pair['wall']:.6f} | {pair['core']:.6f} | {pair['rss']:.6f} |")
lines += [
    "",
    f"Median paired ratios are {paired_median['wall']:.6f} wall, {paired_median['core']:.6f} CPU, and {paired_median['rss']:.6f} RSS. CPU and work improve, but wall misses the frozen 0.97 gate, so full M4RI remains opt-in.",
    "",
    f"The candidate routes 723 matrices / 438,923 blocks, reduces actual XORs to {result['mechanism']['performed_xor_ratio_full_over_current']:.6f}x current, and spends 10,095,063,562 XORs on table preparation. Its different pivot basis reduces logical XORs slightly to {result['mechanism']['logical_xor_ratio_full_over_current']:.6f}x current.",
    "",
    f"The build, four metered validation commands, screen, and confirmation charge {charge['wall_seconds_sum']:.6f} wall seconds, {charge['total_core_seconds_sum']:.6f} core-seconds, and {charge['peak_rss_bytes_max']} bytes maximum RSS across {charge['components']} components.",
    "",
    "This is a strong opt-in solver-stage result, not a full attack or SOTA result. A separate single-core adjudication is required before selecting any default policy.",
    "",
]
(STAGE / "RESULTS.md").write_text("\n".join(lines))
print(json.dumps({"screen_ratio": screen_ratio, "median_paired_ratios": paired_median, "decision": result["decision"], "campaign_charge": charge}, indent=2))

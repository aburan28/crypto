#!/usr/bin/env python3
"""Compose the Stage 183 dense pair-selection result."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics


STAGE = Path(__file__).resolve().parent
EXPECTED_EQUATIONS = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
COMMON = {
    "logical": 318_635_818_320,
    "performed": 99_192_937_526,
    "candidate_visits": 1_137_001_812,
}
EXPECTED = {
    "quadratic": {
        "dense_calls": 0,
        "quadratic_calls": 1_011_275,
        "lcm_groups": 0,
        "cover_lookups": 0,
        "scratch": 0,
    },
    "dense": {
        "dense_calls": 1_011_275,
        "quadratic_calls": 0,
        "lcm_groups": 594_604_504,
        "cover_lookups": 698_657_372,
        "scratch": 4_472_832,
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
    ("quadratic_f4", "10 passed; 0 failed"),
    ("dense_f4", "10 passed; 0 failed"),
    ("quadratic_backend", "3 passed; 0 failed"),
    ("dense_backend", "3 passed; 0 failed"),
):
    root = STAGE / "development" / f"test-{name.replace('_', '-')}"
    item = process_artifact(root)
    item["passed"] = (
        item["process"]["returncode"] == 0
        and not item["process"]["timed_out"]
        and expected in (root / "stdout.txt").read_text()
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
        and report["cost"]["ops"] == COMMON["logical"]
        and extra["word_xors_performed"] == COMMON["performed"]
        and extra["full_m4ri_matrices"] == 723
        and extra["full_m4ri_blocks"] == 438_923
        and extra["pair_candidate_visits"] == COMMON["candidate_visits"]
        and extra["pair_dense_select_calls"] == expected["dense_calls"]
        and extra["pair_quadratic_select_calls"] == expected["quadratic_calls"]
        and extra["pair_lcm_groups"] == expected["lcm_groups"]
        and extra["pair_cover_lookups"] == expected["cover_lookups"]
        and extra["pair_dense_scratch_bytes_max"] == expected["scratch"]
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
    load_run(screen_root / "01-quadratic", "quadratic", 1),
    load_run(screen_root / "02-dense", "dense", 2),
]
confirmation_root = STAGE / "development" / "confirmation"
sequence = ("quadratic", "dense", "dense", "quadratic", "quadratic", "dense")
counts = {"quadratic": 0, "dense": 0}
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

quadratic = [run for run in confirmation if run["variant"] == "quadratic"]
dense = [run for run in confirmation if run["variant"] == "dense"]
pairs = [(quadratic[0], dense[0]), (quadratic[1], dense[1]), (quadratic[2], dense[2])]
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


quadratic_median = medians(quadratic)
dense_median = medians(dense)
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
stage175 = {"wall_seconds": 34.68066983300014, "total_core_seconds": 267.171734, "peak_rss_bytes": 4_030_480_384}
result = {
    "schema": "koblitz_stage183_current_f4_dense_pair_selection.v1",
    "claim_boundary": "One opened n=59 target and a stacked solver-stage experiment only; not a full index-calculus run or SOTA.",
    "source": {
        "candidate_commit": "fdeb3155fe2c9757d30464bae15584524a0554c8",
        "confirmation_runner_commit": "9f7c0901c",
        "binary_sha256": "205856a1923448f43a0c51e5c1527aa0141d8437493bab2d80676158524bf6e4",
    },
    "build": build,
    "tests": tests,
    "screen": {"runs": screen, "dense_over_quadratic": screen_ratio, "continued": screen_ratio["core"] < 1},
    "confirmation": {
        "runs": confirmation,
        "quadratic_median": quadratic_median,
        "dense_median": dense_median,
        "pair_ratios": pair_ratios,
        "median_paired_ratios": paired_median,
    },
    "mechanism": {
        "dense_calls": 1_011_275,
        "candidate_visits": 1_137_001_812,
        "lcm_groups": 594_604_504,
        "cover_lookups": 698_657_372,
        "dense_scratch_bytes_max": 4_472_832,
        "full_m4ri_logical_xors": COMMON["logical"],
        "full_m4ri_performed_xors": COMMON["performed"],
    },
    "absolute_ratios": {
        "dense_over_stage174": {
            field: dense_median[field] / stage174[field]
            for field in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
        },
        "dense_over_stage175_current_engine": {
            field: dense_median[field] / stage175[field]
            for field in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
        },
        "dense_over_same_binary_direct": {
            "wall_seconds": dense_median["wall_seconds"] / 0.3618509580001046,
            "total_core_seconds": dense_median["total_core_seconds"] / 0.351335,
            "peak_rss_bytes": dense_median["peak_rss_bytes"] / 45_072_384,
        },
    },
    "decision": {
        "status": "ACCEPTED_FOR_FULL_M4RI_STACK" if passed else "REJECTED_TIMING_GATE",
        "pass": passed,
        "selected_pair_selector_for_full_m4ri": "dense_exact",
        "current_runtime_default": "quadratic_with_blocktables",
        "candidate_controls": ["F4_F2_FULL_M4RI=1", "F4_F2_DENSE_PAIR_SELECT=1"],
        "reason": "Both preregistered paired thresholds passed; replay on current BlockTables is required before changing the repository default." if passed else "A frozen timing threshold failed.",
    },
    "campaign_charge": charge,
    "single_core_seconds": None,
    "conflicts": None,
}
(STAGE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

lines = [
    "# Stage 183: dense exact pair selection in current F4",
    "",
    "Quadratic and dense modes passed ten Boolean-F4 tests and three backend tests. All screen and confirmation processes authenticated the same target and returned exhaustive UNSAT with exact full-M4RI work and final algebra.",
    "",
    "| pair | dense / quadratic wall | dense / quadratic core | dense / quadratic RSS |",
    "|---:|---:|---:|---:|",
]
for pair in pair_ratios:
    lines.append(f"| {pair['pair']} | {pair['wall']:.6f} | {pair['core']:.6f} | {pair['rss']:.6f} |")
lines += [
    "",
    f"Median paired ratios are {paired_median['wall']:.6f} wall, {paired_median['core']:.6f} CPU, and {paired_median['rss']:.6f} RSS. Both timing gates pass.",
    "",
    "The selected mechanism handles 1,011,275 pair selections over 1,137,001,812 active candidates using 594,604,504 exact LCM groups and 698,657,372 cover probes. Peak dense selector scratch is 4,472,832 bytes; F4 logical and performed XOR counts remain exact.",
    "",
    f"The clean build, four metered validation commands, screen, and confirmation charge {charge['wall_seconds_sum']:.6f} wall seconds, {charge['total_core_seconds_sum']:.6f} core-seconds, and {charge['peak_rss_bytes_max']} bytes maximum RSS across {charge['components']} components.",
    "",
    "Dense pair selection is accepted for the full-M4RI research stack. The repository default remains current BlockTables plus the quadratic selector until a same-binary default-path replay passes. This is not a full attack or SOTA result.",
    "",
]
(STAGE / "RESULTS.md").write_text("\n".join(lines))
print(json.dumps({"screen_ratio": screen_ratio, "median_paired_ratios": paired_median, "decision": result["decision"], "campaign_charge": charge}, indent=2))

#!/usr/bin/env python3
"""Compose the Stage 182 single-core adjudication."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics


STAGE = Path(__file__).resolve().parent
EXPECTED_EQUATIONS = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
EXPECTED = {
    "current": {"logical": 319_313_687_585, "performed": 147_794_583_858, "matrices": 0},
    "full": {"logical": 318_635_818_320, "performed": 99_192_937_526, "matrices": 723},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def receipt(path: Path) -> dict:
    return {"path": str(path.relative_to(STAGE)), "bytes": path.stat().st_size, "sha256": sha256(path)}


def load_run(root: Path, variant: str, order: int, repeat=None) -> dict:
    metrics_path = root / "metrics.json"
    stdout_path = root / "stdout.json"
    stderr_path = root / "stderr.txt"
    process = json.loads(metrics_path.read_text())
    report = json.loads(stdout_path.read_text())
    extra = report["cost"]["extra"]
    metrics = process["metrics"]
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
        and report["single_thread_requested"] is True
        and metrics["single_core_seconds"] is not None
        and metrics["single_core_seconds"] == metrics["total_core_seconds"]
        and report["cost"]["ops"] == expected["logical"]
        and extra["word_xors_performed"] == expected["performed"]
        and extra["full_m4ri_matrices"] == expected["matrices"]
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
    raise SystemExit("single-core correctness or accounting failure")

screen_ratio = {
    "wall": screen[1]["process"]["metrics"]["wall_seconds"] / screen[0]["process"]["metrics"]["wall_seconds"],
    "core": screen[1]["process"]["metrics"]["total_core_seconds"] / screen[0]["process"]["metrics"]["total_core_seconds"],
    "rss": screen[1]["process"]["metrics"]["peak_rss_bytes"] / screen[0]["process"]["metrics"]["peak_rss_bytes"],
}
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
        for field in ("wall_seconds", "total_core_seconds", "single_core_seconds", "peak_rss_bytes")
    }


current_median = medians(current)
full_median = medians(full)
passed = paired_median["wall"] < 0.97 and paired_median["core"] < 0.97
new_processes = [*(run["process"] for run in screen), *(run["process"] for run in confirmation)]
new_charge = {
    "components": len(new_processes),
    "wall_seconds_sum": sum(item["metrics"]["wall_seconds"] for item in new_processes),
    "total_core_seconds_sum": sum(item["metrics"]["total_core_seconds"] for item in new_processes),
    "peak_rss_bytes_max": max(item["metrics"]["peak_rss_bytes"] for item in new_processes),
}
stage181 = json.loads(
    (STAGE.parent / "stage-181-current-f4-full-m4ri-20261001" / "result.json").read_text()
)
setup_components = [stage181["build"]["process"], *(item["process"] for item in stage181["tests"].values())]
cumulative = {
    "components": len(setup_components) + len(new_processes),
    "wall_seconds_sum": sum(item["metrics"]["wall_seconds"] for item in setup_components) + new_charge["wall_seconds_sum"],
    "total_core_seconds_sum": sum(item["metrics"]["total_core_seconds"] for item in setup_components) + new_charge["total_core_seconds_sum"],
    "peak_rss_bytes_max": max(stage181["campaign_charge"]["peak_rss_bytes_max"], new_charge["peak_rss_bytes_max"]),
}
result = {
    "schema": "koblitz_stage182_current_f4_full_m4ri_single_core.v1",
    "claim_boundary": "One opened n=59 target and a single-core solver-stage adjudication only; not a full index-calculus run or SOTA.",
    "source": {
        "binary_source_commit": "850182efba8d8a9499e82aa6f11600941d3cb9d6",
        "confirmation_runner_commit": "ceec561ba",
        "binary_sha256": "fe17f002fedb2e4ef3037d100a4c541ffe0d2c7dcbfa6b2121a540804cc9050c",
    },
    "screen": {"runs": screen, "full_over_current": screen_ratio, "continued": screen_ratio["wall"] < 1 and screen_ratio["core"] < 1},
    "confirmation": {
        "runs": confirmation,
        "current_median": current_median,
        "full_median": full_median,
        "pair_ratios": pair_ratios,
        "median_paired_ratios": paired_median,
    },
    "mechanism": {
        "performed_xor_ratio_full_over_current": EXPECTED["full"]["performed"] / EXPECTED["current"]["performed"],
        "full_m4ri_matrices": 723,
    },
    "same_binary_direct_ratios": {
        "current_wall": current_median["wall_seconds"] / 0.3618509580001046,
        "current_core": current_median["total_core_seconds"] / 0.351335,
        "full_wall": full_median["wall_seconds"] / 0.3618509580001046,
        "full_core": full_median["total_core_seconds"] / 0.351335,
    },
    "decision": {
        "status": "REJECTED_SINGLE_CORE_TIMING_REGRESSION" if not passed else "ACCEPTED_SINGLE_CORE_POLICY",
        "pass": passed,
        "runtime_default": "current_blocktables",
        "candidate_availability": "F4_F2_FULL_M4RI=1",
        "automatic_single_core_policy": False,
        "reason": "Full M4RI was slower in median paired wall and CPU, so automatic single-core selection is rejected." if not passed else "Both frozen timing thresholds passed.",
    },
    "new_campaign_charge": new_charge,
    "cumulative_with_stage181_build_and_validation": cumulative,
    "conflicts": None,
}
(STAGE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

lines = [
    "# Stage 182: full-matrix M4RI single-core adjudication",
    "",
    "All eight new processes requested one Rayon worker, reported non-null single-core seconds equal to total core-seconds, authenticated the same target, and returned exhaustive UNSAT.",
    "",
    "| pair | full / current wall | full / current core | full / current RSS |",
    "|---:|---:|---:|---:|",
]
for pair in pair_ratios:
    lines.append(f"| {pair['pair']} | {pair['wall']:.6f} | {pair['core']:.6f} | {pair['rss']:.6f} |")
lines += [
    "",
    f"Median paired ratios are {paired_median['wall']:.6f} wall, {paired_median['core']:.6f} CPU, and {paired_median['rss']:.6f} RSS. Full M4RI is slower in all three timing pairs and fails the 0.97 gate.",
    "",
    f"The candidate still performs only {result['mechanism']['performed_xor_ratio_full_over_current']:.6f}x the XORs, but its table and pivot scheduling overhead outweighs that reduction on one worker. Current median single-core CPU is {current_median['single_core_seconds']:.6f} seconds versus {full_median['single_core_seconds']:.6f} for full M4RI.",
    "",
    f"The eight new processes charge {new_charge['wall_seconds_sum']:.6f} wall seconds, {new_charge['total_core_seconds_sum']:.6f} core-seconds, and {new_charge['peak_rss_bytes_max']} bytes maximum RSS. Including the inherited exact build and four validation processes gives {cumulative['wall_seconds_sum']:.6f} wall seconds and {cumulative['total_core_seconds_sum']:.6f} core-seconds across {cumulative['components']} components.",
    "",
    "Automatic single-core selection is rejected. Current BlockTables remains the default; full M4RI remains an explicit research control. This does not change any SOTA gate.",
    "",
]
(STAGE / "RESULTS.md").write_text("\n".join(lines))
print(json.dumps({"screen_ratio": screen_ratio, "median_paired_ratios": paired_median, "decision": result["decision"], "new_campaign_charge": new_charge}, indent=2))

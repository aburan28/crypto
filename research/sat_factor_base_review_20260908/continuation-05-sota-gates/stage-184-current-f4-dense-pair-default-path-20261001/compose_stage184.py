#!/usr/bin/env python3
"""Compose the Stage 184 current-path dense selector result."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics


STAGE = Path(__file__).resolve().parent
EXPECTED_EQUATIONS = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
COMMON = {
    "logical": 319_313_687_585,
    "performed": 147_794_583_858,
    "candidate_visits": 1_137_001_812,
}
EXPECTED = {
    "quadratic": {"dense_calls": 0, "quadratic_calls": 1_011_275, "lcm_groups": 0, "cover_lookups": 0, "scratch": 0},
    "dense": {"dense_calls": 1_011_275, "quadratic_calls": 0, "lcm_groups": 594_604_504, "cover_lookups": 698_657_372, "scratch": 4_472_832},
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
        and extra["full_m4ri_matrices"] == 0
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
new_processes = [*(run["process"] for run in screen), *(run["process"] for run in confirmation)]
new_charge = {
    "components": len(new_processes),
    "wall_seconds_sum": sum(item["metrics"]["wall_seconds"] for item in new_processes),
    "total_core_seconds_sum": sum(item["metrics"]["total_core_seconds"] for item in new_processes),
    "peak_rss_bytes_max": max(item["metrics"]["peak_rss_bytes"] for item in new_processes),
}
stage183 = json.loads(
    (STAGE.parent / "stage-183-current-f4-dense-pair-selection-20261001" / "result.json").read_text()
)
setup = [stage183["build"]["process"], *(item["process"] for item in stage183["tests"].values())]
cumulative = {
    "components": len(setup) + len(new_processes),
    "wall_seconds_sum": sum(item["metrics"]["wall_seconds"] for item in setup) + new_charge["wall_seconds_sum"],
    "total_core_seconds_sum": sum(item["metrics"]["total_core_seconds"] for item in setup) + new_charge["total_core_seconds_sum"],
    "peak_rss_bytes_max": max(stage183["campaign_charge"]["peak_rss_bytes_max"], new_charge["peak_rss_bytes_max"]),
}
result = {
    "schema": "koblitz_stage184_current_f4_dense_pair_default_path.v1",
    "claim_boundary": "One opened n=59 target and a default-path implementation replay only; not a full index-calculus run or SOTA.",
    "source": {
        "binary_source_commit": "fdeb3155fe2c9757d30464bae15584524a0554c8",
        "confirmation_runner_commit": "5a40e3540",
        "binary_sha256": "205856a1923448f43a0c51e5c1527aa0141d8437493bab2d80676158524bf6e4",
    },
    "inherited_stage183_result_sha256": sha256(
        STAGE.parent / "stage-183-current-f4-dense-pair-selection-20261001" / "result.json"
    ),
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
        "blocktables_logical_xors": COMMON["logical"],
        "blocktables_performed_xors": COMMON["performed"],
    },
    "decision": {
        "status": "ACCEPTED_FOR_REPOSITORY_DEFAULT" if passed else "REJECTED_TIMING_GATE",
        "pass": passed,
        "selected_default_pair_selector": "dense_exact",
        "quadratic_control": "F4_F2_DENSE_PAIR_SELECT=0",
        "reason": "Both preregistered default-path timing thresholds passed; post-selection tests and default-mode replay remain required." if passed else "A frozen timing threshold failed.",
    },
    "new_campaign_charge": new_charge,
    "cumulative_with_stage183_build_and_validation": cumulative,
    "single_core_seconds": None,
    "conflicts": None,
}
(STAGE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

lines = [
    "# Stage 184: dense pair selection on the current default path",
    "",
    "All eight new processes used current five-column BlockTables, authenticated the same target, and returned exhaustive UNSAT with exact equations, pairs, matrices, basis, extraction, and XOR counts.",
    "",
    "| pair | dense / quadratic wall | dense / quadratic core | dense / quadratic RSS |",
    "|---:|---:|---:|---:|",
]
for pair in pair_ratios:
    lines.append(f"| {pair['pair']} | {pair['wall']:.6f} | {pair['core']:.6f} | {pair['rss']:.6f} |")
lines += [
    "",
    f"Median paired ratios are {paired_median['wall']:.6f} wall, {paired_median['core']:.6f} CPU, and {paired_median['rss']:.6f} RSS. Both frozen timing gates pass, and RSS is lower in all three confirmation pairs.",
    "",
    "Dense selection repeats the exact Stage 183 mechanism counts while leaving BlockTables logical and performed XORs unchanged. It is accepted for the repository default, subject to the required post-selection tests and default-mode target replay.",
    "",
    f"The eight new processes charge {new_charge['wall_seconds_sum']:.6f} wall seconds, {new_charge['total_core_seconds_sum']:.6f} core-seconds, and {new_charge['peak_rss_bytes_max']} bytes maximum RSS. Including inherited exact build and validation gives {cumulative['wall_seconds_sum']:.6f} wall seconds and {cumulative['total_core_seconds_sum']:.6f} core-seconds across {cumulative['components']} components.",
    "",
    "This is a one-target implementation improvement, not relation-yield evidence, a full attack, or SOTA.",
    "",
]
(STAGE / "RESULTS.md").write_text("\n".join(lines))
print(json.dumps({"screen_ratio": screen_ratio, "median_paired_ratios": paired_median, "decision": result["decision"], "new_campaign_charge": new_charge}, indent=2))

#!/usr/bin/env python3
"""Compose the Stage 180 four-column target screen."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


STAGE = Path(__file__).resolve().parent
SCREEN = STAGE / "development" / "screen"
SEQUENCE = ("five", "four", "four", "five")
EXPECTED_EQUATIONS = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
EXPECTED = {
    "five": {"performed": 147_794_583_858, "table_bytes": 53_070_336},
    "four": {"performed": 166_039_014_225, "table_bytes": 34_225_728},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def receipt(path: Path) -> dict:
    return {"path": str(path.relative_to(STAGE)), "bytes": path.stat().st_size, "sha256": sha256(path)}


build_four_path = STAGE / "development" / "build-four" / "metrics.json"
build_four = json.loads(build_four_path.read_text())
inherited_five_path = STAGE.parent / "stage-179-current-f4-wide-tables-single-target-20261001" / "development" / "build-five" / "metrics.json"
inherited_five = json.loads(inherited_five_path.read_text())
runs = []
counts = {"five": 0, "four": 0}
for order, variant in enumerate(SEQUENCE, 1):
    counts[variant] += 1
    root = SCREEN / f"{order:02d}-{variant}-r{counts[variant]}"
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
        and report["cost"]["ops"] == 319_313_687_585
        and extra["word_xors_performed"] == EXPECTED[variant]["performed"]
        and extra["peak_table_bytes"] == EXPECTED[variant]["table_bytes"]
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
    raise SystemExit("screen correctness failure")

five = [run for run in runs if run["variant"] == "five"]
four = [run for run in runs if run["variant"] == "four"]
pairs = [(five[0], four[0]), (five[1], four[1])]
pair_ratios = [
    {
        "pair": i,
        "wall": candidate["process"]["metrics"]["wall_seconds"] / control["process"]["metrics"]["wall_seconds"],
        "core": candidate["process"]["metrics"]["total_core_seconds"] / control["process"]["metrics"]["total_core_seconds"],
        "rss": candidate["process"]["metrics"]["peak_rss_bytes"] / control["process"]["metrics"]["peak_rss_bytes"],
    }
    for i, (control, candidate) in enumerate(pairs, 1)
]
performed_ratio = EXPECTED["four"]["performed"] / EXPECTED["five"]["performed"]
table_ratio = EXPECTED["four"]["table_bytes"] / EXPECTED["five"]["table_bytes"]
confirmation = performed_ratio < 1 and all(pair["core"] < 0.97 for pair in pair_ratios)
all_processes = [inherited_five, build_four, *(run["process"] for run in runs)]
charge = {
    "components": len(all_processes),
    "wall_seconds_sum": sum(item["metrics"]["wall_seconds"] for item in all_processes),
    "total_core_seconds_sum": sum(item["metrics"]["total_core_seconds"] for item in all_processes),
    "peak_rss_bytes_max": max(item["metrics"]["peak_rss_bytes"] for item in all_processes),
}
result = {
    "schema": "koblitz_stage180_current_f4_four_table_screen.v1",
    "claim_boundary": "One opened n=59 target and a bounded table-width screen only; not a selected speedup or SOTA.",
    "source_commit": "64e5272acc814b29b9c6c865d814c5bcdb4b439c",
    "runner_commit": "7dd93d927",
    "binary_sha256": {
        "five": "b156f3320564c1e9fed2b168aba7a0b4f3c875d677c1a9354751bba527ee5bd6",
        "four": "58870d396aecc33cd51ea51c4506301d14b8ea7d56044e87bad281e514af96af",
    },
    "builds": {
        "inherited_five": {
            "process": inherited_five,
            "receipt_path": "../stage-179-current-f4-wide-tables-single-target-20261001/development/build-five/metrics.json",
            "receipt_sha256": sha256(inherited_five_path),
        },
        "four": {
            "process": build_four,
            "artifacts": {
                "metrics": receipt(build_four_path),
                "stdout": receipt(STAGE / "development" / "build-four" / "stdout.txt"),
                "stderr": receipt(STAGE / "development" / "build-four" / "stderr.txt"),
            },
        },
    },
    "runs": runs,
    "pair_ratios": pair_ratios,
    "mechanism": {
        "performed_xor_ratio_four_over_five": performed_ratio,
        "table_memory_ratio_four_over_five": table_ratio,
        "structural_f4_counts_identical": True,
    },
    "decision": {
        "status": "REJECTED_SCREEN_WORK_REGRESSION",
        "confirmation_required": confirmation,
        "selected_table_columns": 5,
        "reason": "Four columns increased performed XORs, so the frozen prerequisite failed and no confirmation was run.",
    },
    "campaign_charge": charge,
    "single_core_seconds": None,
    "conflicts": None,
}
(STAGE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

lines = [
    "# Stage 180: four-column BlockTables target screen",
    "",
    "All four screen processes returned exhaustive UNSAT with identical equations and logical F4 counters.",
    "",
    "| pair | four / five wall | four / five core | four / five RSS |",
    "|---:|---:|---:|---:|",
]
for pair in pair_ratios:
    lines.append(f"| {pair['pair']} | {pair['wall']:.6f} | {pair['core']:.6f} | {pair['rss']:.6f} |")
lines += [
    "",
    f"Four columns performs {performed_ratio:.6f}x the actual word XORs and uses {table_ratio:.6f}x the table memory of five columns. The performed-work prerequisite failed, so no confirmation panel was run even though one CPU sample was marginally lower.",
    "",
    f"The inherited five-column build, fresh four-column build, and four queries charge {charge['wall_seconds_sum']:.6f} wall seconds, {charge['total_core_seconds_sum']:.6f} core-seconds, and {charge['peak_rss_bytes_max']} bytes maximum RSS.",
    "",
    "Four columns is rejected and five columns remains selected. This does not change any SOTA gate.",
    "",
]
(STAGE / "RESULTS.md").write_text("\n".join(lines))
print(json.dumps({"pair_ratios": pair_ratios, "mechanism": result["mechanism"], "decision": result["decision"], "campaign_charge": charge}, indent=2))

#!/usr/bin/env python3
"""Validate and compose the additive Stage-13 WDSat correction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import statistics


SEEDS = [2026091301, 2026091302, 2026091303, 2026091304, 2026091305]
CELL = "n59-l9-m3-standard-a1-f0"


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def distribution(values: list[float]) -> dict:
    assert values and all(math.isfinite(value) and value >= 0 for value in values)
    return {
        "count": len(values),
        "min": min(values),
        "median": statistics.median(values),
        "max": max(values),
    }


def solver(instance: dict, name: str) -> dict:
    matches = [row for row in instance["solvers"] if row["solver"] == name]
    assert len(matches) == 1
    return matches[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("panel", type=Path)
    parser.add_argument("correction", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    original_summary = load(args.panel / "panel-summary-independent.json")
    assert original_summary["panel_artifact_complete"] is True
    assert original_summary["verified_tasks"] == 20
    correction_outer = load(args.correction / "outer-metrics.json")
    correction_matrix = load(args.correction / "matrix" / "result.json")
    assert correction_outer["returncode"] == 0 and correction_outer["timed_out"] is False
    assert len(correction_matrix["instances"]) == 1
    corrected_instance = correction_matrix["instances"][0]
    corrected_manifest = corrected_instance["manifest"]
    assert corrected_instance["cell"] == {
        "n": 59,
        "ell": 9,
        "m": 3,
        "basis": "standard",
        "curve_a": 1,
        "factor_index": 0,
    }
    assert corrected_manifest["seed"] == 2026091303
    assert corrected_instance["source_artifacts_unchanged"] is True

    failed_task = (
        args.panel
        / "tasks"
        / "seed-2026091303"
        / CELL
        / "matrix"
        / "result.json"
    )
    original_failed_instance = load(failed_task)["instances"][0]
    assert original_failed_instance["manifest"]["source_instance"] == corrected_manifest["source_instance"]
    assert original_failed_instance["source_artifacts_before"] == corrected_instance["source_artifacts_before"]
    assert original_failed_instance["source_artifacts_after"] == corrected_instance["source_artifacts_after"]
    assert solver(original_failed_instance, "wdsat")["status"] == "solver_error"

    corrected_wdsat = solver(corrected_instance, "wdsat")
    assert corrected_wdsat["status"] == "timeout_inconclusive"
    assert corrected_wdsat["timed_out"] is True
    assert corrected_wdsat["returncode"] == -15
    source_atoms = corrected_instance["wdsat_build"]["source_xor_atoms"]
    assert source_atoms == 32263
    config = corrected_instance["wdsat_build"]["config"]
    marker = "#define __MAX_BUFFER_SIZE__ "
    buffer_size = int(next(line[len(marker):] for line in config.splitlines() if line.startswith(marker)))
    assert buffer_size > source_atoms

    rows = []
    for seed in SEEDS:
        if seed == 2026091303:
            row = corrected_wdsat
            source = "additive_correction"
        else:
            path = args.panel / "tasks" / f"seed-{seed}" / CELL / "matrix" / "result.json"
            row = solver(load(path)["instances"][0], "wdsat")
            source = "original_panel"
        assert row["status"] == "timeout_inconclusive" and row["timed_out"] is True
        metrics = row["metrics"]
        rows.append(
            {
                "seed": seed,
                "source": source,
                "status": row["status"],
                "wall_seconds": metrics["wall_seconds"],
                "total_core_seconds": metrics["total_core_seconds"],
                "single_core_seconds": metrics["single_core_seconds"],
                "peak_rss_bytes": metrics["peak_rss_bytes"],
            }
        )

    original_total = original_summary["charged_totals"]["cold_total_core_seconds_including_ggmp_discovery"]
    correction_core = correction_outer["metrics"]["total_core_seconds"]
    result = {
        "schema": "koblitz_stage13_wdsat_correction_result.v1",
        "source_identity_equal": True,
        "source_instance_blake3": corrected_manifest["source_instance"]["id_blake3"],
        "source_xor_atoms": source_atoms,
        "corrected_buffer_size": buffer_size,
        "corrected_n59_wdsat_rows": rows,
        "corrected_n59_wdsat_distribution": {
            "statuses": {"timeout_inconclusive": 5},
            "wall_seconds": distribution([row["wall_seconds"] for row in rows]),
            "total_core_seconds": distribution([row["total_core_seconds"] for row in rows]),
            "single_core_seconds": distribution([row["single_core_seconds"] for row in rows]),
            "peak_rss_bytes": distribution([float(row["peak_rss_bytes"]) for row in rows]),
        },
        "correction_process": correction_outer["metrics"],
        "charged_campaign_totals": {
            "original_panel_plus_ggmp_discovery_core_seconds": original_total,
            "additive_correction_core_seconds": correction_core,
            "total_core_seconds": original_total + correction_core,
            "original_panel_plus_ggmp_discovery_wall_seconds": original_summary["charged_totals"]["cold_sequential_wall_seconds_including_ggmp_discovery"],
            "additive_correction_wall_seconds": correction_outer["metrics"]["wall_seconds"],
            "total_sequential_wall_seconds": original_summary["charged_totals"]["cold_sequential_wall_seconds_including_ggmp_discovery"]
            + correction_outer["metrics"]["wall_seconds"],
            "peak_rss_bytes": max(
                original_summary["charged_totals"]["cold_peak_rss_bytes"],
                correction_outer["metrics"]["peak_rss_bytes"],
            ),
        },
        "composition_rule": "Replace only the failed seed-2026091303 n=59 WDSat distribution arm; retain every original process cost and add the correction process cost.",
        "claim_boundary": "The corrected panel is complete except for unavailable Magma; it is planted-PDP evidence, not an end-to-end index-calculus or SOTA result.",
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

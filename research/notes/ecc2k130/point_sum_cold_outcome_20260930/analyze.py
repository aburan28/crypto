#!/usr/bin/env python3
"""Recompute every six-cell complete-process CPU ratio from archived raw wait4 data."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path
import statistics
import tarfile

from archive import ARTIFACTS, CELLS, RUN_ID
from verify_archive import verify as verify_archive

CRITICAL_95 = {5: 2.7764451051977987, 20: 2.093024054408263}
ARM_ORDER = ("control_a", "point_sum", "rho", "control_b")


def paired(values: list[float]) -> dict:
    assert len(values) in CRITICAL_95 and all(value > 0 for value in values)
    logs = [math.log(value) for value in values]
    mean = statistics.mean(logs)
    half = CRITICAL_95[len(values)] * statistics.stdev(logs) / math.sqrt(len(values))
    return {"median": statistics.median(values),
            "interval_95pct": [math.exp(mean - half), math.exp(mean + half)]}


def member_json(archive: tarfile.TarFile, name: str) -> dict:
    item = archive.extractfile(name)
    assert item is not None, name
    return json.loads(item.read())


def analyze_cell(evidence: Path, cell: str) -> dict:
    path = evidence / "artifacts" / f"{ARTIFACTS[cell]}.tar.gz"
    with tarfile.open(path, "r:gz") as archive:
        run = member_json(archive, f"{cell}/cold_run.json")
        isolation = member_json(archive, f"{cell}/isolation.jsonl")
        hosted = member_json(archive, f"{cell}/receipt.json")
    replay = json.loads(gzip.decompress(
        (evidence / "replay" / f"{cell}.json.gz").read_bytes()))
    assert run["schema"] == "ecc2k130-point-sum-cold-run-v1"
    assert run["cell"] == cell and run["status"] == "PASS"
    assert run["mode"] == "measure"
    assert isolation["schema"] == "isolated-bench/1"
    assert isolation["label"] == f"point-sum-cold-{cell}"
    blocks = run["spec"]["blocks"]
    assert len(run["runs"]) == 4 * blocks
    by_block = {block: {} for block in range(blocks)}
    for item in run["runs"]:
        assert item["exit_code"] == 0 and item["stopped_for"] is None
        by_block[item["block"]][item["arm"]] = (
            item["child_user_cpu_seconds"] + item["child_system_cpu_seconds"])
    assert all(set(group) == set(ARM_ORDER) for group in by_block.values())
    aa = []
    point_over_control = []
    point_over_rho = []
    for group in by_block.values():
        a, b = group["control_a"], group["control_b"]
        aa.append(b / a)
        point_over_control.append(group["point_sum"] / math.sqrt(a * b))
        point_over_rho.append(group["point_sum"] / group["rho"])
    ratios = {"control_b_over_a": paired(aa),
              "point_over_control": paired(point_over_control),
              "point_over_rho": paired(point_over_rho)}
    aa_result = ratios["control_b_over_a"]
    aa_valid = (0.9 <= aa_result["median"] <= 1.1 and
                aa_result["interval_95pct"][0] <= 1 <= aa_result["interval_95pct"][1])
    zero_samples = isolation["contended_samples"] == 0
    remaining_user_threads = len(isolation["left_on_reserved"]["user_threads"])
    if cell == "n37_L1024":
        assert hosted["status"] == replay["status"] == "FAIL"
        assert hosted["error_type"] == replay["error_type"] == "AssertionError"
        frozen_replay = "FAIL"
        strict_eligible = False
    else:
        assert hosted["status"] == replay["status"] == "PASS"
        assert hosted["aa_valid"] == aa_valid
        assert hosted["uncontended"] == (zero_samples and remaining_user_threads == 0)
        strict_eligible = hosted["timing_eligible"]
        assert strict_eligible == (aa_valid and zero_samples and remaining_user_threads == 0)
        for key in ratios:
            for name in ("median",):
                assert math.isclose(ratios[key][name], hosted["paired"][key][name],
                                    rel_tol=1e-12, abs_tol=1e-12)
            for actual, expected in zip(ratios[key]["interval_95pct"],
                                        hosted["paired"][key]["interval_95pct"]):
                assert math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12)
        frozen_replay = "PASS"
    rho_ratio = ratios["point_over_rho"]["median"]
    return {"cell": cell, "n": run["spec"]["n"], "L": run["spec"]["L"],
            "K": run["spec"]["K"], "blocks": blocks, "children": len(run["runs"]),
            "producer_status": run["status"], "frozen_replay": frozen_replay,
            "aa_valid": aa_valid, "monitor_zero_contention_samples": zero_samples,
            "remaining_user_threads": remaining_user_threads,
            "strict_timing_eligible": strict_eligible,
            "paired_cpu_diagnostic": ratios,
            "required_point_cpu_reduction_to_match_rho_diagnostic": (
                max(0.0, 1 - 1 / rho_ratio)),
            "source_commit": run["source"]["commit"],
            "compact_binary_sha256": run["source"]["compact_binary_sha256"],
            "rho_binary_sha256": run["source"]["rho_binary_sha256"]}


def analyze(evidence: Path) -> dict:
    archive_receipt = verify_archive(evidence)
    assert archive_receipt["status"] == "PASS"
    rows = [analyze_cell(evidence, cell) for cell in CELLS]
    assert sum(row["children"] for row in rows) == 180
    assert all(row["aa_valid"] and row["monitor_zero_contention_samples"] for row in rows)
    assert not any(row["strict_timing_eligible"] for row in rows)
    assert all(row["source_commit"] == "d8d67d3b683334db105195b0fc74d8e660dc4b5e"
               for row in rows)
    return {"schema": "ecc2k130-point-sum-cold-analysis-v1",
            "run_id": RUN_ID, "status": "NO_ELIGIBLE_TIMING",
            "evidence_class": "complete_process_cpu_diagnostic_only",
            "fully_eligible_cells": 0, "rows": rows,
            "method_speedup": None, "rho_crossover": None,
            "n131_transfer": None, "common_operation_unit": None}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite analysis"
    result = analyze(args.evidence.resolve())
    args.out.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "run_id": RUN_ID,
                      "fully_eligible_cells": result["fully_eligible_cells"],
                      "cells": len(result["rows"])}, sort_keys=True))


if __name__ == "__main__":
    main()

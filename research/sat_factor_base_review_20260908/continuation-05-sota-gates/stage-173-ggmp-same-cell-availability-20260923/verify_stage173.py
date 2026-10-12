#!/usr/bin/env python3
"""Verify the Stage 173 GGMP same-cell availability evidence."""

from __future__ import annotations

import hashlib
import json
import math
import os
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEV = HERE / "development"
RESULT = HERE / "result.json"


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(value: bool, message: str) -> None:
    if not value:
        raise AssertionError(message)


def verify_artifacts(node: object) -> int:
    if isinstance(node, dict):
        if {"path", "bytes", "sha256"}.issubset(node):
            path = HERE / node["path"]
            require(path.is_file(), f"missing {path}")
            require(path.stat().st_size == node["bytes"], f"size {path}")
            require(sha256(path) == node["sha256"], f"hash {path}")
            return 1
        return sum(verify_artifacts(value) for value in node.values())
    if isinstance(node, list):
        return sum(verify_artifacts(value) for value in node)
    return 0


def metered() -> list[dict]:
    paths = sorted(set(DEV.rglob("metrics.json")) | set(DEV.rglob("*.metrics.json")))
    rows = []
    for path in paths:
        try:
            record = load(path)
        except (OSError, json.JSONDecodeError):
            continue
        metrics = record.get("metrics", record)
        if all(key in metrics for key in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")):
            rows.append(metrics)
    return rows


def main() -> None:
    result = load(RESULT)
    same = result["same_cell_availability"]
    cap = result["enumeration_cap_control"]
    positive = result["available_positive_control"]
    sr, cr, pr = same["report"], cap["report"], positive["report"]
    forbidden = {
        "target_constructed": False,
        "target_subgroup_enumerated": False,
        "discrete_log_labels_constructed": False,
        "relation_yield_used": False,
        "solver_timing_used": False,
    }
    checks = {
        "schema": result["schema"] == "koblitz_stage173_ggmp_same_cell_availability.v1",
        "status": result["status"] == "complete_ggmp_same_cell_structurally_unavailable",
        "same_cell": (
            sr["status"] == "unavailable_no_divisor_of_requested_dimension"
            and sr["complete_factor_degrees"] == [1, 58]
            and sr["factor_degrees"] == [1]
            and sr["requested_divisor_exists"] is False
            and sr["candidates"] == []
            and sr["selected"] is None
            and sr["forbidden_inputs"] == forbidden
        ),
        "cap_control": (
            cr["status"] == "unavailable_factor_enumeration_degree_cap"
            and cr["complete_factor_degrees"] == [1, 58]
            and cr["factor_degrees"] == [1]
            and cr["requested_divisor_exists"] is True
            and cr["requested_dimension"] == 58
            and cr["selected"] is None
            and cr["forbidden_inputs"] == forbidden
        ),
        "positive_control": (
            pr["status"] == "selected"
            and pr["complete_factor_degrees"] == [1, 5, 5, 5, 5, 5, 5]
            and pr["factor_degrees"] == pr["complete_factor_degrees"]
            and pr["requested_divisor_exists"] is True
            and len(pr["candidates"]) == 6
            and pr["selected"]["divisor_indices"] == [1]
            and pr["selected"]["m_cofactor_admissible"] is True
            and pr["forbidden_inputs"] == forbidden
        ),
        "clean_processes": all(
            item["process"]["returncode"] == 0
            and item["process"]["timed_out"] is False
            and item["process"]["orphan_group_terminated"] is False
            for item in (same, cap, positive)
        ),
        "claim_boundary": (
            result["ggmp_same_cell_benchmark_applicable"] is False
            and result["koblitz_index_calculus_sota"] is False
            and result["gates"]["all_seven_gates_passed"] is False
        ),
    }
    commit = result["source_revision"]["commit"]
    repo = HERE.parents[3]
    tree = subprocess.check_output(["git", "-C", str(repo), "rev-parse", f"{commit}^{{tree}}"], text=True).strip()
    checks["source_tree"] = tree == result["source_revision"]["tree"]
    archive = Path("/Volumes/SSD990/koblitz-ggmp-build25-c6205c56/source.tar")
    binary = Path(
        "/Volumes/SSD990/koblitz-ggmp-build25-c6205c56/bin/koblitz_public_factor_base_discovery"
    )
    external_source_or_binary_present = archive.is_file() or binary.is_file()
    checks["external_source_and_binary_if_present"] = (
        not external_source_or_binary_present
        or (
            archive.is_file()
            and binary.is_file()
            and sha256(archive) == result["source_revision"]["archive_sha256"]
            and sha256(binary) == result["source_revision"]["discovery_binary_sha256"]
        )
    )
    rows = metered()
    accounting = result["stage173_development_accounting"]
    checks["development_accounting"] = (
        len(rows) == accounting["resource_components"]
        and math.isclose(sum(row["wall_seconds"] for row in rows), accounting["summed_wall_seconds"], abs_tol=1e-12)
        and math.isclose(sum(row["total_core_seconds"] for row in rows), accounting["total_core_seconds"], abs_tol=1e-12)
        and max(row["peak_rss_bytes"] for row in rows) == accounting["peak_rss_bytes"]
    )
    tests = result["validation"]["consumer_tests"]
    checks["consumer_tests"] = (
        tests["unknown_scalar_consumers"]["process"]["returncode"] == 0
        and tests["production_consumers"]["process"]["returncode"] == 0
        and "Ran 12 tests" in (HERE / tests["unknown_scalar_consumers"]["artifacts"]["stderr"]["path"]).read_text()
        and "Ran 11 tests" in (HERE / tests["production_consumers"]["artifacts"]["stderr"]["path"]).read_text()
        and '"checks": 28' in (HERE / tests["production_consumers"]["artifacts"]["stdout"]["path"]).read_text()
    )
    artifact_count = verify_artifacts(result)
    checks["artifacts"] = artifact_count > 0
    require(all(checks.values()), json.dumps({key: value for key, value in checks.items() if not value}))
    verification = {
        "schema": "koblitz_stage173_ggmp_same_cell_availability_verification.v1",
        "status": "pass",
        "result_sha256": sha256(RESULT),
        "artifact_receipts_verified": artifact_count,
        "external_source_and_binary_present": external_source_or_binary_present,
        "checks": checks,
    }
    if os.environ.get("KIC_VERIFY_NO_WRITE") != "1":
        (HERE / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
    print(json.dumps(verification, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

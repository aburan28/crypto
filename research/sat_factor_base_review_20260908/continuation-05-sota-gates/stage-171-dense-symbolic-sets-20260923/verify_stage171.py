#!/usr/bin/env python3
"""Fail-closed verification of the Stage 171 native-F4 result."""

from __future__ import annotations

import hashlib
import json
import math
import statistics
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEV = HERE / "development"
RESULT = HERE / "result.json"


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=0.0, abs_tol=1e-12)


def stable_report(report: dict) -> dict:
    report = json.loads(json.dumps(report))
    report.pop("timing_ns", None)
    report.pop("solver_description", None)
    report.pop("solver_symbolic_monomial_sets", None)
    cost = report["cost"]
    cost.pop("wall_ns", None)
    for key in (
        "build_ns",
        "eliminate_ns",
        "pair_update_ns",
        "dense_symbolic_set_steps",
        "dense_symbolic_set_bytes_max",
        "m4ri_scratch_bytes_max",
    ):
        cost["extra"].pop(key, None)
    return report


def metric_receipts() -> list[dict]:
    records = []
    paths = sorted(set(DEV.rglob("metrics.json")) | set(DEV.rglob("*.metrics.json")))
    for path in paths:
        if "development/clean-t14-b64/clean-t14-b64/" in path.relative_to(HERE).as_posix():
            continue
        try:
            record = load(path)
        except (OSError, json.JSONDecodeError):
            continue
        metrics = record.get("metrics", record)
        if all(key in metrics for key in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")):
            records.append(metrics)
    return records


def verify_artifacts(node: object) -> int:
    checked = 0
    if isinstance(node, dict):
        if set(("path", "bytes", "sha256")).issubset(node):
            path = HERE / node["path"]
            require(path.is_file(), f"missing artifact {path}")
            require(path.stat().st_size == node["bytes"], f"size mismatch {path}")
            require(sha256(path) == node["sha256"], f"hash mismatch {path}")
            checked += 1
        else:
            checked += sum(verify_artifacts(value) for value in node.values())
    elif isinstance(node, list):
        checked += sum(verify_artifacts(value) for value in node)
    return checked


def main() -> None:
    result = load(RESULT)
    checks: dict[str, bool] = {}
    checks["schema"] = result["schema"] == "koblitz_stage171_dense_symbolic_sets.v1"
    checks["status"] = result["status"] == "complete_dense_symbolic_sets_true_negative"
    checks["claim_boundary"] = (
        result["koblitz_index_calculus_sota"] is False
        and result["full_cost_gate_passed"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["gates"]["all_seven_gates_passed"] is False
    )

    commit = result["source_revision"]["commit"]
    repo = HERE.parents[3]
    tree = subprocess.check_output(
        ["git", "-C", str(repo), "rev-parse", f"{commit}^{{tree}}"], text=True
    ).strip()
    checks["source_commit_and_tree"] = tree == result["source_revision"]["tree"]
    external_archive = Path("/Volumes/SSD990/koblitz-native-f4-build21-6de26d8c/source.tar")
    external_backend = Path(
        "/Volumes/SSD990/koblitz-native-f4-build21-6de26d8c/bin/koblitz_pdp_backend"
    )
    checks["external_clean_archive"] = (
        external_archive.is_file()
        and sha256(external_archive) == result["source_revision"]["archive_sha256"]
    )
    checks["external_clean_backend"] = (
        external_backend.is_file()
        and sha256(external_backend) == result["source_revision"]["backend_sha256"]
    )

    selected = result["selected_native_f4"]
    control = result["same_binary_hash_control"]
    fingerprint = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
    clean = result["clean_paired_repeats"]["candidates"] + result["clean_paired_repeats"]["controls"]
    checks["all_clean_terminals"] = all(
        item["process"]["returncode"] == 0
        and item["process"]["timed_out"] is False
        and item["report"]["status"] == "unsat"
        and item["report"]["exhaustive"] is True
        and item["report"]["solver_equations_blake3"] == fingerprint
        for item in clean
    )
    reference = stable_report(clean[0]["report"])
    checks["same_scientific_report"] = all(stable_report(item["report"]) == reference for item in clean)
    checks["selected_schedule"] = (
        selected["report"]["solver_x1_batch_size"] == 512
        and selected["report"]["solver_x1_batches_completed"] == 1
        and selected["report"]["solver_rayon_threads_requested"] == 12
        and selected["report"]["fixed_x1_systems_completed"] == 242
    )
    checks["same_binary_control"] = (
        selected["process"]["command"] == control["process"]["command"]
        and selected["report"]["cost"]["extra"]["dense_symbolic_set_steps"] == 1204
        and control["report"]["cost"]["extra"]["dense_symbolic_set_steps"] == 0
    )

    pairs = result["clean_paired_repeats"]["pair_ratios"]
    recomputed = {
        key: statistics.median(pair[key] for pair in pairs)
        for key in pairs[0]
        if key != "pair"
    }
    checks["paired_medians"] = all(
        close(recomputed[key], result["clean_paired_repeats"]["median_pair_ratios"][key])
        for key in recomputed
    )
    checks["measured_improvement"] = (
        recomputed["candidate_over_control_wall"] < 1
        and recomputed["candidate_over_control_core"] < 1
        and recomputed["candidate_over_control_rss"] < 1
        and recomputed["candidate_over_control_build_ns"] < 1
    )

    tests = result["validation"]
    checks["tests"] = (
        tests["boolean_f4_tests"]["passed"] == 13
        and tests["fixed_x1_backend_tests"]["passed"] == 3
        and "13 passed; 0 failed" in (HERE / tests["boolean_f4_tests"]["artifacts"]["stdout"]["path"]).read_text()
        and "3 passed; 0 failed" in (HERE / tests["fixed_x1_backend_tests"]["artifacts"]["stdout"]["path"]).read_text()
    )

    receipts = metric_receipts()
    accounting = result["development_accounting"]
    checks["development_accounting"] = (
        len(receipts) == accounting["metered_components"]
        and close(sum(m["wall_seconds"] for m in receipts), accounting["summed_wall_seconds"])
        and close(sum(m["total_core_seconds"] for m in receipts), accounting["total_core_seconds"])
        and max(m["peak_rss_bytes"] for m in receipts) == accounting["maximum_peak_rss_bytes"]
        and accounting["complete_total"] is None
    )
    campaign = result["campaign_accounting"]
    checks["campaign_accounting"] = (
        campaign["complete_campaign_cost"] is None
        and campaign["measured_lower_bound"]["resource_components"]
        == campaign["stage170_measured_lower_bound"]["resource_components"]
        + accounting["metered_components"]
        and close(
            campaign["measured_lower_bound"]["summed_wall_seconds"],
            campaign["stage170_measured_lower_bound"]["summed_wall_seconds"]
            + accounting["summed_wall_seconds"],
        )
    )

    artifact_count = verify_artifacts(result)
    checks["artifact_receipts"] = artifact_count > 0
    require(all(checks.values()), json.dumps({k: v for k, v in checks.items() if not v}))
    verification = {
        "schema": "koblitz_stage171_dense_symbolic_sets_verification.v1",
        "status": "pass",
        "result_sha256": sha256(RESULT),
        "artifact_receipts_verified": artifact_count,
        "checks": checks,
    }
    (HERE / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
    print(json.dumps(verification, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fail-closed verification of the Stage 172 single-core result."""

from __future__ import annotations

import hashlib
import json
import math
import os
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULT = HERE / "result.json"


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(value: bool, message: str) -> None:
    if not value:
        raise AssertionError(message)


def stable(report: dict) -> dict:
    report = json.loads(json.dumps(report))
    report.pop("timing_ns", None)
    report.pop("solver_description", None)
    report.pop("solver_symbolic_monomial_sets", None)
    report["cost"].pop("wall_ns", None)
    for key in (
        "build_ns",
        "eliminate_ns",
        "pair_update_ns",
        "dense_symbolic_set_steps",
        "dense_symbolic_set_bytes_max",
    ):
        report["cost"]["extra"].pop(key, None)
    return report


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


def main() -> None:
    result = load(RESULT)
    candidate = result["selected_single_core_native_f4"]
    control = result["same_binary_single_core_hash_control"]
    direct = result["same_binary_single_core_direct_mitm"]
    cr, hr, dr = candidate["report"], control["report"], direct["report"]
    cm, hm, dm = candidate["process"]["metrics"], control["process"]["metrics"], direct["process"]["metrics"]
    checks = {
        "schema": result["schema"] == "koblitz_stage172_native_f4_single_core.v1",
        "status": result["status"] == "complete_native_f4_single_core_true_negative",
        "single_thread_contract": all(
            report["single_thread_requested"] is True
            and report["solver_rayon_threads_requested"] == 1
            and report["solver_x1_batch_size"] == 1
            for report in (cr, hr)
        ),
        "valid_single_core_metrics": all(
            metrics["single_core_seconds"] == metrics["total_core_seconds"]
            and metrics["single_core_seconds"] > 0
            for metrics in (cm, hm, dm)
        ),
        "exact_f4_terminals": all(
            report["status"] == "unsat"
            and report["exhaustive"] is True
            and report["fixed_x1_systems_completed"] == 242
            and report["solver_equations_blake3"]
            == "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
            for report in (cr, hr)
        ),
        "same_scientific_report": stable(cr) == stable(hr),
        "direct_terminal": dr["status"] == "unsat" and dr["exhaustive"] is True,
        "same_binary": (
            candidate["process"]["command"][0]
            == control["process"]["command"][0]
            == direct["process"]["command"][0]
        ),
        "dense_mode_exercised": (
            cr["cost"]["extra"]["dense_symbolic_set_steps"] == 1204
            and hr["cost"]["extra"]["dense_symbolic_set_steps"] == 0
        ),
        "measured_single_core_improvement": (
            cm["wall_seconds"] < hm["wall_seconds"]
            and cm["total_core_seconds"] < hm["total_core_seconds"]
            and cr["cost"]["extra"]["build_ns"] < hr["cost"]["extra"]["build_ns"]
        ),
        "claim_boundary": (
            result["koblitz_index_calculus_sota"] is False
            and result["full_cost_gate_passed"] is False
            and result["gates"]["all_seven_gates_passed"] is False
        ),
    }
    commit = result["source_revision"]["commit"]
    repo = HERE.parents[3]
    tree = subprocess.check_output(["git", "-C", str(repo), "rev-parse", f"{commit}^{{tree}}"], text=True).strip()
    checks["source_tree"] = tree == result["source_revision"]["tree"]
    accounting = result["stage172_execution_accounting"]
    checks["accounting"] = (
        accounting["resource_components"] == 3
        and math.isclose(
            accounting["summed_wall_seconds"],
            cm["wall_seconds"] + hm["wall_seconds"] + dm["wall_seconds"],
            abs_tol=1e-12,
        )
        and math.isclose(
            accounting["total_core_seconds"],
            cm["total_core_seconds"] + hm["total_core_seconds"] + dm["total_core_seconds"],
            abs_tol=1e-12,
        )
    )
    artifact_count = verify_artifacts(result)
    checks["artifacts"] = artifact_count == 9
    require(all(checks.values()), json.dumps({key: value for key, value in checks.items() if not value}))
    verification = {
        "schema": "koblitz_stage172_native_f4_single_core_verification.v1",
        "status": "pass",
        "result_sha256": sha256(RESULT),
        "artifact_receipts_verified": artifact_count,
        "checks": checks,
    }
    if os.environ.get("KIC_VERIFY_NO_WRITE") != "1":
        (HERE / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
    print(json.dumps(verification, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

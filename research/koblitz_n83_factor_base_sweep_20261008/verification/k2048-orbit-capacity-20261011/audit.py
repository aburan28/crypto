#!/usr/bin/env python3
"""Independently cross-check the retained K2048 construction receipts."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / "structured-orbit-search-gate-20261010/evidence-01/independent-audit.json"
EXPECTED = {
    "required_max_variables": 111478,
    "required_domain_clauses": 892470,
    "sat_variables": 111478,
    "sat_clauses": 1140557,
    "sat_xor_rows": 4980,
    "sat_s3_nodes": 4,
    "domain_representatives": 2048,
    "legal_x_coordinates": 169984,
    "memory_cgroup_limit_bytes": 4 * 1024**3,
    "memory_cgroup_swap_limit_bytes": 0,
    "source_attestation_mode": "supervised_snapshot",
    "status": "PASS_model_construction_only",
    "solver_search_executed": False,
    "relation_stage_executed": False,
    "rank_stage_executed": False,
    "total_index_calculus_runtime_ms": None,
    "selected_best_total_runtime": None,
}


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def main() -> int:
    errors: list[str] = []
    build, identity, budget = (read(HERE / name) for name in
        ("build.json", "input-identity.json", "budget.json"))
    first, corrected = (HERE / name for name in ("run", "run-corrected"))
    first_outer, outer = (read(path / "outer.json") for path in (first, corrected))
    first_config, config = (read(path / "config.json") for path in (first, corrected))
    worker = read(corrected / "worker.json")
    manifest = read(HERE / "input/manifest.json")
    prior = read(PRIOR)
    if build["status"] != "PASS_build" or build["exit_code"] != 0 or not build["cleanup_verified"]:
        errors.append("build or build cleanup")
    if build["wall_elapsed_seconds_L0"] > 600:
        errors.append("build wall cap")
    for label, path, bound_outer, bound_config in (
        ("first", first, first_outer, first_config),
        ("corrected", corrected, outer, config),
    ):
        if bound_config["binary_sha256"] != build["binary_sha256"] or bound_config["source_commit"] != build["source_commit"]:
            errors.append(f"{label} binary or source binding")
        if sha(path / "config.json") != bound_outer["config_sha256"]:
            errors.append(f"{label} config digest")
        for name in ("stdout", "stderr"):
            if sha(path / f"{name}.log") != bound_outer[f"{name}_sha256"]:
                errors.append(f"{label} {name} digest")
        for relative, expected_hash in bound_config["source_sha256"].items():
            if sha(path / "source" / relative) != expected_hash:
                errors.append(f"{label} source {relative}")
        if bound_outer["cleanup_error"] is not None or bound_outer["source_changed"] or bound_outer["binary_changed"] or bound_outer["config_changed"]:
            errors.append(f"{label} cleanup or snapshot changed")
        if bound_outer["process_wall_ms"] > 180000:
            errors.append(f"{label} wall cap")
    if first_outer["status"] != "PRODUCER_FAILURE_worker_exit" or first_outer["worker_exit_code"] != 1 or (first / "worker.json").exists():
        errors.append("first failed attempt classification")
    if b"No such file or directory" not in (first / "stderr.log").read_bytes():
        errors.append("first failed attempt stderr")
    if outer["status"] != "PASS_model_construction_only" or outer["worker_exit_code"] != 0 or outer["worker_error"] is not None:
        errors.append("corrected outer status")
    if sha(corrected / "worker.json") != outer["worker_sha256"]:
        errors.append("worker digest")
    for key, value in EXPECTED.items():
        if worker.get(key) != value:
            errors.append(f"worker {key}")
    if worker["source_commit"] != build["source_commit"] or worker["point_set_blake3"] != manifest["bases"][0]["point_set_blake3"]:
        errors.append("worker source or point-set binding")
    if not 0 < worker["memory_cgroup_peak_bytes"] <= 4 * 1024**3 or worker["memory_cgroup_peak_bytes"] != outer["worker_cgroup_peak_bytes"]:
        errors.append("memory peak")
    for name, receipt in identity["inputs"].items():
        if name != "s3_object" and (sha(HERE / "input" / name) != receipt["sha256"] or
                                    (HERE / "input" / name).stat().st_size != receipt["bytes"]):
            errors.append(f"input {name}")
    object_receipt = identity["inputs"]["s3_object"]
    if object_receipt["s3_uri"] != manifest["bases"][0]["s3_uri"] or object_receipt["compressed_blake3_claim"] != manifest["bases"][0]["compressed_blake3"] or object_receipt["bytes"] != manifest["bases"][0]["bytes"]:
        errors.append("S3 object identity")
    if sha(PRIOR) != budget["prior_audit_sha256"] or prior["cumulative_charge_seconds"] != budget["prior_charge_seconds"]:
        errors.append("prior budget binding")
    if budget["prior_charge_seconds"] != 9000 or budget["round_conservative_charge_seconds"] != 900 or budget["cumulative_conservative_charge_seconds"] != 9900 or budget["remaining_authorized_seconds"] != 0:
        errors.append("conservative budget")
    if budget["build_phase_observed_seconds_L0"] != build["wall_elapsed_seconds_L0"] or budget["first_container_observed_seconds_L0"] != first_outer["process_wall_ms"] / 1000 or budget["corrected_container_observed_seconds_L0"] != outer["process_wall_ms"] / 1000:
        errors.append("observed budget intervals")
    result = {"schema": "n83.k2048-orbit-capacity-audit/v1", "status": "PASS" if not errors else "FAIL",
              "errors": errors, "source_commit": build["source_commit"],
              "worker_sha256": outer["worker_sha256"], "build_binary_sha256": build["binary_sha256"],
              "cumulative_conservative_charge_seconds": budget["cumulative_conservative_charge_seconds"],
              "selected_best_total_runtime": None}
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())

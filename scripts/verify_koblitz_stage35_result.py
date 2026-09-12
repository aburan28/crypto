#!/usr/bin/env python3
"""Verify and freeze the hosted Stage-35 algebraic walked-descent result."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path, PurePosixPath
import tarfile
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage35_algebraic_walk as stage35


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
ARCHIVE = STAGE / "koblitz-stage35-algebraic-walk-evidence-34693802792.tar.gz"
ARCHIVE_SHA = STAGE / "koblitz-stage35-algebraic-walk-evidence-34693802792.tar.gz.sha256"
WORKFLOW = STAGE / "stage-35-workflow-34693802792.json"
ARTIFACTS = STAGE / "stage-35-artifacts-34693802792.json"
DEFAULT_OUTPUT = STAGE / "stage-35-algebraic-walk-result-20260912"

RUN_ID = 34693802792
RUN_COMMIT = "11dc73b8f1aed58d0d6f29d1d95bbc6a156a1387"
ARTIFACT_ID = 10298541006
ARTIFACT_DIGEST = "sha256:d67480bbb4f4ac7ee805fedff050df95eacb3456a487e4dc2a2a7c7b4a1e33af"
ARCHIVE_ROOT = "koblitz-stage35-algebraic-walk-34693802792"
SCHEMA = "koblitz_stage35_hosted_result.v1"
SEAL_SCHEMA = "koblitz_stage35_hosted_result_seal.v1"


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def timestamp(value: str) -> datetime:
    require(isinstance(value, str) and value.endswith("Z"), "workflow timestamp is invalid")
    return datetime.fromisoformat(value[:-1] + "+00:00")


def archive_identity() -> dict[str, Any]:
    fields = ARCHIVE_SHA.read_text().split()
    relative = str(ARCHIVE.relative_to(REPO))
    require(len(fields) == 2 and fields[1] == relative, "Stage-35 archive sidecar is invalid")
    digest = custody.sha256_file(ARCHIVE, "Stage-35 archive")
    require(fields[0] == digest, "Stage-35 archive hash changed")
    return {"path": relative, "bytes": ARCHIVE.stat().st_size, "sha256": digest}


def extract(destination: Path) -> Path:
    with tarfile.open(ARCHIVE, "r:gz") as source:
        members = source.getmembers()
        require(bool(members), "Stage-35 archive is empty")
        roots: set[str] = set()
        for member in members:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and ".." not in path.parts, f"unsafe archive path: {member.name}")
            require(member.isfile() or member.isdir(), f"unsupported archive member: {member.name}")
            if path.parts:
                roots.add(path.parts[0])
        require(roots == {ARCHIVE_ROOT}, "Stage-35 archive root changed")
        source.extractall(destination, filter="data")
    return destination / ARCHIVE_ROOT


def workflow_accounting() -> dict[str, Any]:
    value = load(WORKFLOW, "Stage-35 workflow metadata")
    require(value.get("databaseId") == RUN_ID and value.get("headSha") == RUN_COMMIT, "Stage-35 workflow identity changed")
    require(value.get("status") == "completed" and value.get("conclusion") == "success", "Stage-35 workflow is incomplete")
    jobs = value.get("jobs")
    require(isinstance(jobs, list) and {row.get("name") for row in jobs} == {"validate-control", "algebraic-walk"}, "Stage-35 workflow jobs changed")
    for row in jobs:
        require(row.get("status") == "completed" and row.get("conclusion") == "success", f"Stage-35 job failed: {row.get('name')}")
    production = next(row for row in jobs if row["name"] == "algebraic-walk")
    return {
        "run_id": RUN_ID,
        "head_sha": RUN_COMMIT,
        "url": value.get("url"),
        "workflow_wall_seconds": (timestamp(value["updatedAt"]) - timestamp(value["createdAt"])).total_seconds(),
        "production_job_id": production.get("databaseId"),
        "production_job_url": production.get("url"),
        "production_job_wall_seconds": (timestamp(production["completedAt"]) - timestamp(production["startedAt"])).total_seconds(),
    }


def artifact_accounting() -> dict[str, Any]:
    value = load(ARTIFACTS, "Stage-35 artifact metadata")
    rows = value.get("artifacts")
    require(value.get("total_count") == 1 and isinstance(rows, list) and len(rows) == 1, "Stage-35 artifact count changed")
    row = rows[0]
    require(row.get("id") == ARTIFACT_ID and row.get("name") == f"koblitz-stage35-algebraic-walk-{RUN_ID}", "Stage-35 artifact identity changed")
    require(row.get("digest") == ARTIFACT_DIGEST and row.get("expired") is False, "Stage-35 artifact digest or state changed")
    source = row.get("workflow_run", {})
    require(source.get("id") == RUN_ID and source.get("head_sha") == RUN_COMMIT, "Stage-35 artifact source changed")
    return {key: row[key] for key in ("id", "name", "size_in_bytes", "digest", "expires_at")}


def result() -> dict[str, Any]:
    archive = archive_identity()
    workflow = workflow_accounting()
    artifact = artifact_accounting()
    with tempfile.TemporaryDirectory() as directory:
        root = extract(Path(directory))
        verification = stage35.verify(root / "stage35-build", root / "stage35-run")
        embedded = load(root / "stage35-verification.json", "embedded Stage-35 verification")
        require(embedded == verification, "embedded and recomputed Stage-35 verifications differ")
        build = load(root / "stage35-build/result.json", "Stage-35 build result")
        run = load(root / "stage35-run/result.json", "Stage-35 run result")
        build_seal_sha256 = custody.sha256_file(root / "stage35-build/result-seal.json", "Stage-35 build seal")
        run_seal_sha256 = custody.sha256_file(root / "stage35-run/result-seal.json", "Stage-35 run seal")
    workflow_result = run["workflow_result"]
    baseline = next(row["vs_rho"] for row in workflow_result["stages"] if row["stage"] == "baseline")
    predecessor = run["predecessor_stage33"]
    old_online = predecessor["ic_over_rho_charged_wall_ratio"]
    new_online = verification["ic_over_rho_online_wall_ratio"]
    require(new_online < 1.0 < old_online, "Stage-35 did not cross the finite online boundary")
    require(baseline["verdict"]["charged_crossover"] is True, "Stage-35 online crossover flag changed")
    require(baseline["verdict"]["amortised_crossover"] is False and baseline["verdict"]["whole_process_crossover"] is False, "Stage-35 widened the amortised or whole-process claim")
    return {
        "schema": SCHEMA,
        "status": "hosted_n41_algebraic_walk_verified",
        "archive": archive,
        "artifact": artifact,
        "workflow": workflow,
        "verification": verification,
        "build_seal_sha256": build_seal_sha256,
        "run_seal_sha256": run_seal_sha256,
        "source_commit": verification["source_commit"],
        "factor_base": {
            "algebraically_defined": True,
            "target_subgroup_enumerated": False,
            "target_or_log_labels_used_for_selection": False,
            "selected": verification["factor_base"],
        },
        "collection_summands": 3,
        "descent_summands": 2,
        "walked_probe_cap": 100_000_000,
        "relations": verification["relations"],
        "targets_verified": verification["targets_verified"],
        "target_scalars_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
        "descent_trials": verification["descent_trials"],
        "rho_iterations": verification["rho_iterations"],
        "comparison": {
            "stage33_online_ic_over_rho_wall_ratio": old_online,
            "stage35_online_ic_over_rho_wall_ratio": new_online,
            "online_ratio_improvement_factor": old_online / new_online,
            "stage35_online_crossover": True,
            "stage35_amortised_ic_over_rho_wall_ratio": verification["ic_over_rho_amortised_wall_ratio"],
            "stage35_amortised_crossover": False,
            "stage35_whole_process_crossover": False,
            "stage35_full_available_wall_over_same_targets_rho": verification["full_available_wall_over_same_targets_rho"],
            "rho_seconds_total": baseline["rho"]["seconds_total"],
            "descent_seconds_total": baseline["ic"]["descent_seconds_total"],
            "precompute_seconds": baseline["ic"]["precompute_seconds"],
            "claim_boundary": "finite_post_precomputation_online_crossover_only",
        },
        "resources": {
            "scientific_single_core_elapsed_seconds": verification["single_core_elapsed_seconds"],
            "scientific_outer_core_seconds": verification["scientific_outer_core_seconds"],
            "scientific_process_core_seconds": run["workflow_process"]["metrics"]["total_core_seconds"],
            "scientific_process_wall_seconds": run["workflow_process"]["metrics"]["wall_seconds"],
            "scientific_process_peak_rss_bytes": run["workflow_process"]["metrics"]["peak_rss_bytes"],
            "build_process_total_core_seconds": build["process_totals"]["total_core_seconds"],
            "build_outer_total_core_seconds": build["outer_resources"]["total_core_seconds"],
            "charged_total_core_seconds_available": verification["charged_total_core_seconds_available"],
            "charged_sequential_wall_seconds_available": verification["charged_sequential_wall_seconds_available"],
            "maximum_sampled_process_tree_rss_bytes": verification["maximum_sampled_process_tree_rss_bytes"],
            "workflow_wall_seconds": workflow["workflow_wall_seconds"],
            "production_job_wall_seconds": workflow["production_job_wall_seconds"],
        },
        "sat_conflicts": None,
        "sat_conflict_semantics": "not applicable to the pair-table individual-log arm; Phase B retains SAT conflicts",
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
        "amortised_crossover": False,
        "whole_process_crossover": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "remaining_blockers": [
            "licensed Magma F4 has not executed the same 160-input packet",
            "fresh build plus scientific wall remains more than 200 times the rho wall for these five targets",
            "unaffiliated reproduction and source-pinned novelty review remain absent",
            "the crossover changes a finite constant, not the asymptotic exponent",
        ],
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-35 result output must be new")
    output.mkdir(parents=False)
    custody.write_json_new(output / "verification.json", result())
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    payload = {
        "schema": SEAL_SCHEMA,
        "status": "hosted_result_frozen",
        "workflow_run_id": RUN_ID,
        "workflow_commit": RUN_COMMIT,
        "inventory": inventory,
        "inventory_sha256": custody.canonical_sha256(inventory),
    }
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "result-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "verification.json", "committed Stage-35 verification")
    require(committed == result(), "committed Stage-35 verification differs from source evidence")
    seal = load(output / "result-seal.json", "Stage-35 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-35 result seal is invalid")
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    require(inventory == seal.get("inventory") and custody.canonical_sha256(inventory) == seal.get("inventory_sha256"), "Stage-35 result inventory changed")
    return committed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    create = sub.add_parser("compose")
    create.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    check = sub.add_parser("verify")
    check.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    try:
        value = freeze(args.output.resolve()) if args.command == "compose" else verify(args.output.resolve())
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, VerificationError, stage35.Stage35Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage35-result: {error}")


if __name__ == "__main__":
    main()

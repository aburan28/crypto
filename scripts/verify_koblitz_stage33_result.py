#!/usr/bin/env python3
"""Verify and freeze the hosted Stage-33 n=41 public-target result."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path, PurePosixPath
import tarfile
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage33_n41_unknown_scalar as stage33


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
ARCHIVE = STAGE / "koblitz-stage33-n41-public-evidence-34653116265.tar.gz"
ARCHIVE_SHA = STAGE / "koblitz-stage33-n41-public-evidence-34653116265.tar.gz.sha256"
WORKFLOW = STAGE / "stage-33-n41-public-workflow-34653116265.json"
ARTIFACTS = STAGE / "stage-33-n41-public-artifacts-34653116265.json"
DEFAULT_OUTPUT = STAGE / "stage-33-n41-public-result-20260911"

RUN_ID = 34653116265
RUN_COMMIT = "18ba75e940cba564e6a08213ff29f6bdba93b577"
ARTIFACT_ID = 10284816498
ARTIFACT_DIGEST = "sha256:9096d4747fc168a5c246af552d9399bc410e6ca61e3b634f866457c4e813b7a7"
ARCHIVE_ROOT = "koblitz-stage33-n41-public-34653116265"
SCHEMA = "koblitz_stage33_hosted_result.v1"
SEAL_SCHEMA = "koblitz_stage33_hosted_result_seal.v1"


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
    require(len(fields) == 2 and fields[1] == relative, "Stage-33 archive sidecar is invalid")
    digest = custody.sha256_file(ARCHIVE, "Stage-33 archive")
    require(fields[0] == digest, "Stage-33 archive hash changed")
    return {"path": relative, "bytes": ARCHIVE.stat().st_size, "sha256": digest}


def extract(destination: Path) -> Path:
    with tarfile.open(ARCHIVE, "r:gz") as source:
        members = source.getmembers()
        require(bool(members), "Stage-33 archive is empty")
        roots: set[str] = set()
        for member in members:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and ".." not in path.parts, f"unsafe archive path: {member.name}")
            require(member.isfile() or member.isdir(), f"unsupported archive member: {member.name}")
            if path.parts:
                roots.add(path.parts[0])
        require(roots == {ARCHIVE_ROOT}, "Stage-33 archive root changed")
        source.extractall(destination, filter="data")
    return destination / ARCHIVE_ROOT


def workflow_accounting() -> dict[str, Any]:
    value = load(WORKFLOW, "Stage-33 workflow metadata")
    require(value.get("databaseId") == RUN_ID and value.get("headSha") == RUN_COMMIT, "Stage-33 workflow identity changed")
    require(value.get("status") == "completed" and value.get("conclusion") == "success", "Stage-33 workflow is incomplete")
    jobs = value.get("jobs")
    require(isinstance(jobs, list) and {row.get("name") for row in jobs} == {"validate-control", "n41-public-unknown-scalar"}, "Stage-33 workflow jobs changed")
    for row in jobs:
        require(row.get("status") == "completed" and row.get("conclusion") == "success", f"Stage-33 job failed: {row.get('name')}")
    production = next(row for row in jobs if row["name"] == "n41-public-unknown-scalar")
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
    value = load(ARTIFACTS, "Stage-33 artifact metadata")
    rows = value.get("artifacts")
    require(value.get("total_count") == 1 and isinstance(rows, list) and len(rows) == 1, "Stage-33 artifact count changed")
    row = rows[0]
    require(row.get("id") == ARTIFACT_ID, "Stage-33 artifact ID changed")
    require(row.get("name") == f"koblitz-stage33-n41-public-{RUN_ID}", "Stage-33 artifact name changed")
    require(row.get("digest") == ARTIFACT_DIGEST and row.get("expired") is False, "Stage-33 artifact digest or state changed")
    source = row.get("workflow_run", {})
    require(source.get("id") == RUN_ID and source.get("head_sha") == RUN_COMMIT, "Stage-33 artifact source changed")
    return {key: row[key] for key in ("id", "name", "size_in_bytes", "digest", "expires_at")}


def result() -> dict[str, Any]:
    archive = archive_identity()
    workflow = workflow_accounting()
    artifact = artifact_accounting()
    with tempfile.TemporaryDirectory() as directory:
        root = extract(Path(directory))
        verification = stage33.verify(root / "stage33-build", root / "stage33-run")
        embedded = load(root / "stage33-verification.json", "embedded Stage-33 verification")
        require(embedded == verification, "embedded and recomputed Stage-33 verifications differ")
        build = load(root / "stage33-build/result.json", "Stage-33 build result")
        run = load(root / "stage33-run/result.json", "Stage-33 run result")
        run_seal_sha256 = custody.sha256_file(root / "stage33-run/result-seal.json", "Stage-33 run seal")
        build_seal_sha256 = custody.sha256_file(root / "stage33-build/result-seal.json", "Stage-33 build seal")
    charged_ratio = verification["rho_over_ic_charged_ratio"]
    amortised_ratio = verification["rho_over_ic_amortised_ratio"]
    require(charged_ratio > 0 and amortised_ratio > 0, "Stage-33 ratios are invalid")
    return {
        "schema": SCHEMA,
        "status": "hosted_n41_public_unknown_scalar_verified",
        "archive": archive,
        "artifact": artifact,
        "workflow": workflow,
        "verification": verification,
        "build_seal_sha256": build_seal_sha256,
        "run_seal_sha256": run_seal_sha256,
        "source_commit": verification["source_commit"],
        "factor_base_discovery": {
            "algebraic_recipe": verification["factor_base"]["spec"],
            "target_subgroup_enumerated": False,
            "factor_base_logs_known_by_construction": False,
            "future_target_scalars_available_to_search": False,
            "sampled_targets": verification["factor_base_search"]["targets"],
            "exhaustive_targets": verification["factor_base_search"]["exhaustive_targets"],
            "elapsed_seconds": verification["factor_base_search"]["elapsed_ms"] / 1000.0,
            "selected": verification["factor_base"],
        },
        "relations": verification["relations"],
        "linear_algebra": verification["linear_algebra"],
        "targets_verified": verification["targets_verified"],
        "target_scalars_constructed_or_supplied": False,
        "descent_trials": verification["descent_trials"],
        "rho_iterations": verification["rho_iterations"],
        "ic_over_rho_charged_wall_ratio": 1.0 / charged_ratio,
        "ic_over_rho_amortised_wall_ratio": 1.0 / amortised_ratio,
        "resources": {
            "scientific_single_core_elapsed_seconds": verification["single_core_elapsed_seconds"],
            "scientific_outer_core_seconds": verification["workflow_outer_core_seconds"],
            "build_process_total_core_seconds": build["process_totals"]["total_core_seconds"],
            "build_outer_total_core_seconds": build["outer_resources"]["total_core_seconds"],
            "charged_total_core_seconds_available": verification["charged_total_core_seconds_available"],
            "charged_sequential_wall_seconds_available": verification["charged_sequential_wall_seconds_available"],
            "maximum_sampled_process_tree_rss_bytes": verification["maximum_sampled_process_tree_rss_bytes"],
            "scientific_process_peak_rss_bytes": run["workflow_process"]["metrics"]["peak_rss_bytes"],
            "scientific_process_core_seconds": run["workflow_process"]["metrics"]["total_core_seconds"],
            "scientific_process_wall_seconds": run["workflow_process"]["metrics"]["wall_seconds"],
            "workflow_wall_seconds": workflow["workflow_wall_seconds"],
            "production_job_wall_seconds": workflow["production_job_wall_seconds"],
        },
        "sat_conflicts": None,
        "sat_conflict_semantics": "not applicable to the direct pair-table end-to-end arm; the same-instance SAT conflict matrix is retained in Phase B",
        "preinstalled_operating_system_and_rust_toolchain_charged": False,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "remaining_blockers": [
            "licensed Magma F4 has not executed the same 160-input packet",
            "the larger n=59 regime remains a PDP benchmark rather than an end-to-end unknown-scalar run",
            "preinstalled operating-system and Rust-toolchain acquisition remain outside the Stage-33 charge",
            "unaffiliated reproduction and source-pinned novelty review remain absent",
        ],
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-33 result output must be new")
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
    expected = result()
    committed = load(output / "verification.json", "committed Stage-33 verification")
    require(committed == expected, "committed Stage-33 verification differs from source evidence")
    seal = load(output / "result-seal.json", "committed Stage-33 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "committed Stage-33 result seal is invalid")
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    require(inventory == seal.get("inventory") and custody.canonical_sha256(inventory) == seal.get("inventory_sha256"), "committed Stage-33 result inventory changed")
    return committed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    compose = sub.add_parser("compose")
    compose.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    check = sub.add_parser("verify")
    check.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    try:
        value = freeze(args.output.resolve()) if args.command == "compose" else verify(args.output.resolve())
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, VerificationError, stage33.Stage33Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage33-result: {error}")


if __name__ == "__main__":
    main()

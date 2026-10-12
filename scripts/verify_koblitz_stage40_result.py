#!/usr/bin/env python3
"""Verify the committed hosted Stage-40 parallel precompute control."""

from __future__ import annotations

from datetime import datetime
import json
from pathlib import Path, PurePosixPath
import tarfile
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage40_parallel_precompute as stage40


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
ARCHIVE = STAGE / "koblitz-stage40-parallel-precompute-evidence-34700889826.tar.gz"
ARCHIVE_SHA = STAGE / "koblitz-stage40-parallel-precompute-evidence-34700889826.tar.gz.sha256"
WORKFLOW = STAGE / "stage-40-workflow-34700889826.json"
ARTIFACTS = STAGE / "stage-40-artifacts-34700889826.json"
RESULT = STAGE / "stage-40-parallel-precompute-result-20260912"

RUN_ID = 34700889826
RUN_COMMIT = "62eb3e1fc008cf33bb9e910d7a61993cd8b6921f"
ARTIFACT_ID = 10300312749
ARTIFACT_DIGEST = "sha256:83473ff7c3dcef85789624fcfc33b144bc4de7886e88b49d42b34425e4e38bc4"
ARCHIVE_ROOT = "koblitz-stage40-parallel-precompute-34700889826"
SEAL_SCHEMA = "koblitz_stage40_hosted_result_seal.v1"


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
    require(len(fields) == 2 and fields[1] == relative, "Stage-40 archive sidecar is invalid")
    digest = custody.sha256_file(ARCHIVE, "Stage-40 archive")
    require(fields[0] == digest, "Stage-40 archive hash changed")
    return {"path": relative, "bytes": ARCHIVE.stat().st_size, "sha256": digest}


def extract(destination: Path) -> Path:
    with tarfile.open(ARCHIVE, "r:gz") as source:
        members = source.getmembers()
        require(bool(members), "Stage-40 archive is empty")
        roots: set[str] = set()
        for member in members:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and ".." not in path.parts, f"unsafe archive path: {member.name}")
            require(member.isfile() or member.isdir(), f"unsupported archive member: {member.name}")
            if path.parts:
                roots.add(path.parts[0])
        require(roots == {ARCHIVE_ROOT}, "Stage-40 archive root changed")
        source.extractall(destination, filter="data")
    return destination / ARCHIVE_ROOT


def verify_workflow(committed: dict[str, Any]) -> None:
    value = load(WORKFLOW, "Stage-40 workflow metadata")
    require(value.get("databaseId") == RUN_ID and value.get("headSha") == RUN_COMMIT, "Stage-40 workflow identity changed")
    require(value.get("status") == "completed" and value.get("conclusion") == "success", "Stage-40 workflow is incomplete")
    jobs = value.get("jobs")
    require(isinstance(jobs, list) and {row.get("name") for row in jobs} == {"validate-control", "parallel-precompute"}, "Stage-40 workflow jobs changed")
    require(all(row.get("status") == "completed" and row.get("conclusion") == "success" for row in jobs), "a Stage-40 job failed")
    production = next(row for row in jobs if row["name"] == "parallel-precompute")
    expected = committed["workflow"]
    require(expected["run_id"] == RUN_ID and expected["head_sha"] == RUN_COMMIT and expected["url"] == value.get("url"), "Stage-40 workflow summary changed")
    require(expected["production_job_id"] == production.get("databaseId") and expected["production_job_url"] == production.get("url"), "Stage-40 production job changed")
    require(expected["workflow_wall_seconds"] == (timestamp(value["updatedAt"]) - timestamp(value["createdAt"])).total_seconds(), "Stage-40 workflow wall changed")
    require(expected["production_job_wall_seconds"] == (timestamp(production["completedAt"]) - timestamp(production["startedAt"])).total_seconds(), "Stage-40 production wall changed")


def verify_artifact(committed: dict[str, Any]) -> None:
    value = load(ARTIFACTS, "Stage-40 artifact metadata")
    rows = value.get("artifacts")
    require(value.get("total_count") == 1 and isinstance(rows, list) and len(rows) == 1, "Stage-40 artifact count changed")
    row = rows[0]
    require(row.get("id") == ARTIFACT_ID and row.get("name") == f"koblitz-stage40-parallel-precompute-{RUN_ID}", "Stage-40 artifact identity changed")
    require(row.get("digest") == ARTIFACT_DIGEST and row.get("expired") is False, "Stage-40 artifact digest or state changed")
    source = row.get("workflow_run", {})
    require(source.get("id") == RUN_ID and source.get("head_sha") == RUN_COMMIT, "Stage-40 artifact source changed")
    require(committed["artifact"] == {key: row[key] for key in ("id", "name", "size_in_bytes", "digest", "expires_at")}, "Stage-40 artifact summary changed")


def verify_result() -> dict[str, Any]:
    committed = load(RESULT / "verification.json", "committed Stage-40 result")
    require(committed.get("schema") == "koblitz_stage40_hosted_result.v1" and committed.get("status") == "hosted_abba_parallel_precompute_verified", "Stage-40 result label changed")
    require(committed.get("source_commit") == RUN_COMMIT, "Stage-40 source commit changed")
    require(committed.get("archive") == archive_identity(), "Stage-40 archive summary changed")
    verify_workflow(committed)
    verify_artifact(committed)
    with tempfile.TemporaryDirectory() as directory:
        root = extract(Path(directory))
        verification = stage40.verify(root / "stage40-build", root / "stage40-run")
        embedded = load(root / "stage40-verification.json", "embedded Stage-40 verification")
        require(embedded == verification == committed.get("verification"), "Stage-40 recomputed verification differs")
        run = load(root / "stage40-run/result.json", "Stage-40 run result")
    summary = run["summary"]
    require(committed.get("summary") == summary, "Stage-40 committed summary changed")
    comparison = committed["comparison"]
    require(comparison["precompute_wall_speedup"] == summary["precompute_wall_speedup"] and comparison["five_target_ic_wall_speedup"] == summary["five_target_ic_wall_speedup"], "Stage-40 wall speedup changed")
    require(comparison["process_core_seconds_increase_factor"] == summary["four_core_process_core_seconds_median"] / summary["single_core_process_core_seconds_median"], "Stage-40 CPU cost changed")
    new_amortised = summary["four_core_amortised_ic_over_rho_wall_ratio_median"]
    require(comparison["improvement_factor_over_stage39"] == comparison["stage39_ic_over_rho_amortised_wall_ratio"] / new_amortised, "Stage-40 improvement over Stage 39 changed")
    require(comparison["improvement_factor_over_stage35"] == comparison["stage35_ic_over_rho_amortised_wall_ratio"] / new_amortised, "Stage-40 improvement over Stage 35 changed")
    require(committed.get("same_instance_in_every_cell") is True and committed.get("mathematical_outputs_identical") is True, "Stage-40 instance equality changed")
    require(committed.get("factor_base_algebraically_defined") is True and committed.get("target_subgroup_enumerated_for_factor_base") is False, "Stage-40 factor-base boundary changed")
    require(committed.get("target_scalars_constructed_or_supplied") is False and committed.get("factor_base_logs_known_by_construction") is False, "Stage-40 scalar-label boundary changed")
    require(committed.get("full_cost_gate_passed") is False and committed.get("koblitz_index_calculus_sota") is False, "Stage-40 claim widened")
    seal = load(RESULT / "result-seal.json", "Stage-40 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-40 result seal is invalid")
    inventory = custody.all_regular_inventory(RESULT, {"result-seal.json"})
    require(inventory == seal.get("inventory") and custody.canonical_sha256(inventory) == seal.get("inventory_sha256"), "Stage-40 result inventory changed")
    return committed


if __name__ == "__main__":
    try:
        print(json.dumps(verify_result(), indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, VerificationError, stage40.Stage40Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage40-result: {error}")

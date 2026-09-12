#!/usr/bin/env python3
"""Verify the committed hosted Stage-39 fixed algebraic holdout."""

from __future__ import annotations

from datetime import datetime
import json
from pathlib import Path, PurePosixPath
import tarfile
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage39_fixed_algebraic_holdout as stage39


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
ARCHIVE = STAGE / "koblitz-stage39-fixed-algebraic-evidence-34699822015.tar.gz"
ARCHIVE_SHA = STAGE / "koblitz-stage39-fixed-algebraic-evidence-34699822015.tar.gz.sha256"
WORKFLOW = STAGE / "stage-39-workflow-34699822015.json"
ARTIFACTS = STAGE / "stage-39-artifacts-34699822015.json"
RESULT = STAGE / "stage-39-fixed-algebraic-result-20260912"

RUN_ID = 34699822015
RUN_COMMIT = "6f7ff404d615e11df73c0a0e0607bf504b2e09af"
ARTIFACT_ID = 10300470572
ARTIFACT_DIGEST = "sha256:c97eb111c5c4904ec1c42093a8bd2f673367e60cc5370a4dd23123a47ec2e098"
ARCHIVE_ROOT = "koblitz-stage39-fixed-algebraic-34699822015"
SEAL_SCHEMA = "koblitz_stage39_hosted_result_seal.v1"


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
    require(len(fields) == 2 and fields[1] == relative, "Stage-39 archive sidecar is invalid")
    digest = custody.sha256_file(ARCHIVE, "Stage-39 archive")
    require(fields[0] == digest, "Stage-39 archive hash changed")
    return {"path": relative, "bytes": ARCHIVE.stat().st_size, "sha256": digest}


def extract(destination: Path) -> Path:
    with tarfile.open(ARCHIVE, "r:gz") as source:
        members = source.getmembers()
        require(bool(members), "Stage-39 archive is empty")
        roots: set[str] = set()
        for member in members:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and ".." not in path.parts, f"unsafe archive path: {member.name}")
            require(member.isfile() or member.isdir(), f"unsupported archive member: {member.name}")
            if path.parts:
                roots.add(path.parts[0])
        require(roots == {ARCHIVE_ROOT}, "Stage-39 archive root changed")
        source.extractall(destination, filter="data")
    return destination / ARCHIVE_ROOT


def verify_workflow(committed: dict[str, Any]) -> None:
    value = load(WORKFLOW, "Stage-39 workflow metadata")
    require(value.get("databaseId") == RUN_ID and value.get("headSha") == RUN_COMMIT, "Stage-39 workflow identity changed")
    require(value.get("status") == "completed" and value.get("conclusion") == "success", "Stage-39 workflow is incomplete")
    jobs = value.get("jobs")
    require(isinstance(jobs, list) and {row.get("name") for row in jobs} == {"validate-control", "fixed-algebraic-holdout"}, "Stage-39 workflow jobs changed")
    require(all(row.get("status") == "completed" and row.get("conclusion") == "success" for row in jobs), "a Stage-39 job failed")
    production = next(row for row in jobs if row["name"] == "fixed-algebraic-holdout")
    expected = committed["workflow"]
    require(expected["run_id"] == RUN_ID and expected["head_sha"] == RUN_COMMIT and expected["url"] == value.get("url"), "Stage-39 workflow summary changed")
    require(expected["production_job_id"] == production.get("databaseId") and expected["production_job_url"] == production.get("url"), "Stage-39 production job changed")
    require(expected["workflow_wall_seconds"] == (timestamp(value["updatedAt"]) - timestamp(value["createdAt"])).total_seconds(), "Stage-39 workflow wall changed")
    require(expected["production_job_wall_seconds"] == (timestamp(production["completedAt"]) - timestamp(production["startedAt"])).total_seconds(), "Stage-39 production wall changed")


def verify_artifact(committed: dict[str, Any]) -> None:
    value = load(ARTIFACTS, "Stage-39 artifact metadata")
    rows = value.get("artifacts")
    require(value.get("total_count") == 1 and isinstance(rows, list) and len(rows) == 1, "Stage-39 artifact count changed")
    row = rows[0]
    require(row.get("id") == ARTIFACT_ID and row.get("name") == f"koblitz-stage39-fixed-algebraic-{RUN_ID}", "Stage-39 artifact identity changed")
    require(row.get("digest") == ARTIFACT_DIGEST and row.get("expired") is False, "Stage-39 artifact digest or state changed")
    source = row.get("workflow_run", {})
    require(source.get("id") == RUN_ID and source.get("head_sha") == RUN_COMMIT, "Stage-39 artifact source changed")
    require(committed["artifact"] == {key: row[key] for key in ("id", "name", "size_in_bytes", "digest", "expires_at")}, "Stage-39 artifact summary changed")


def verify_result() -> dict[str, Any]:
    committed = load(RESULT / "verification.json", "committed Stage-39 result")
    require(committed.get("schema") == "koblitz_stage39_hosted_result.v1" and committed.get("status") == "hosted_fixed_algebraic_holdout_verified", "Stage-39 result label changed")
    require(committed.get("source_commit") == RUN_COMMIT, "Stage-39 source commit changed")
    require(committed.get("archive") == archive_identity(), "Stage-39 archive summary changed")
    verify_workflow(committed)
    verify_artifact(committed)
    with tempfile.TemporaryDirectory() as directory:
        root = extract(Path(directory))
        verification = stage39.verify(root / "stage39-build", root / "stage39-run")
        embedded = load(root / "stage39-verification.json", "embedded Stage-39 verification")
        require(embedded == verification == committed.get("verification"), "Stage-39 recomputed verification differs")
        run = load(root / "stage39-run/result.json", "Stage-39 run result")
    workflow = run["workflow_result"]
    baseline = next(row["vs_rho"] for row in workflow["stages"] if row["stage"] == "baseline")
    logs = next(row for row in workflow["stages"] if row["stage"] == "logs")
    relations = committed["relations"]
    require(relations["accepted"] == logs["relations"] and relations["trials"] == logs["trials"], "Stage-39 relation totals changed")
    require(relations["summands_scanned"] == logs["summands_scanned"] and relations["columns"] == logs["columns"], "Stage-39 relation charge changed")
    comparison = committed["comparison"]
    new_amortised = 1.0 / baseline["ratio"]["amortised"]
    require(comparison["stage39_ic_over_rho_amortised_wall_ratio"] == new_amortised, "Stage-39 amortized ratio changed")
    require(comparison["improvement_factor_over_stage38"] == comparison["stage38_ic_over_rho_amortised_wall_ratio"] / new_amortised, "Stage-39 improvement over Stage 38 changed")
    require(comparison["improvement_factor_over_stage35"] == comparison["stage35_ic_over_rho_amortised_wall_ratio"] / new_amortised, "Stage-39 improvement over Stage 35 changed")
    require(committed.get("factor_base_discovery") == verification["factor_base_discovery"], "Stage-39 factor-base discovery changed")
    invalid = committed.get("invalid_predecessor", {})
    require(invalid.get("run_id") == 34699355146 and invalid.get("status") == "completed_invalid" and invalid.get("scientific_result_admitted") is False, "Stage-39 invalid predecessor boundary changed")
    require(committed.get("factor_base_algebraically_defined") is True and committed.get("target_subgroup_enumerated_for_factor_base") is False, "Stage-39 factor-base boundary changed")
    require(committed.get("target_scalars_constructed_or_supplied") is False and committed.get("factor_base_logs_known_by_construction") is False, "Stage-39 scalar-label boundary changed")
    require(committed.get("full_cost_gate_passed") is False and committed.get("koblitz_index_calculus_sota") is False, "Stage-39 claim widened")
    seal = load(RESULT / "result-seal.json", "Stage-39 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-39 result seal is invalid")
    inventory = custody.all_regular_inventory(RESULT, {"result-seal.json"})
    require(inventory == seal.get("inventory") and custody.canonical_sha256(inventory) == seal.get("inventory_sha256"), "Stage-39 result inventory changed")
    return committed


if __name__ == "__main__":
    try:
        print(json.dumps(verify_result(), indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, VerificationError, stage39.Stage39Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage39-result: {error}")

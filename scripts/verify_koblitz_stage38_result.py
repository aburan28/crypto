#!/usr/bin/env python3
"""Verify the committed hosted Stage-38 algebraic window sweep."""

from __future__ import annotations

from datetime import datetime
import json
from pathlib import Path, PurePosixPath
import tarfile
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage38_collection_window as stage38


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
ARCHIVE = STAGE / "koblitz-stage38-algebraic-collection-evidence-34698066872.tar.gz"
ARCHIVE_SHA = STAGE / "koblitz-stage38-algebraic-collection-evidence-34698066872.tar.gz.sha256"
WORKFLOW = STAGE / "stage-38-workflow-34698066872.json"
ARTIFACTS = STAGE / "stage-38-artifacts-34698066872.json"
RESULT = STAGE / "stage-38-algebraic-collection-result-20260912"

RUN_ID = 34698066872
RUN_COMMIT = "647481caed3c703d270800922d1e262a5b683834"
ARTIFACT_ID = 10299351870
ARTIFACT_DIGEST = "sha256:0676d18bdf9fb3bb8995b9df570d2e4500e1bafe8785813b5ce7e425ef444839"
ARCHIVE_ROOT = "koblitz-stage38-algebraic-collection-34698066872"
SEAL_SCHEMA = "koblitz_stage38_hosted_result_seal.v1"


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
    require(len(fields) == 2 and fields[1] == relative, "Stage-38 archive sidecar is invalid")
    digest = custody.sha256_file(ARCHIVE, "Stage-38 archive")
    require(fields[0] == digest, "Stage-38 archive hash changed")
    return {"path": relative, "bytes": ARCHIVE.stat().st_size, "sha256": digest}


def extract(destination: Path) -> Path:
    with tarfile.open(ARCHIVE, "r:gz") as source:
        members = source.getmembers()
        require(bool(members), "Stage-38 archive is empty")
        roots: set[str] = set()
        for member in members:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and ".." not in path.parts, f"unsafe archive path: {member.name}")
            require(member.isfile() or member.isdir(), f"unsupported archive member: {member.name}")
            if path.parts:
                roots.add(path.parts[0])
        require(roots == {ARCHIVE_ROOT}, "Stage-38 archive root changed")
        source.extractall(destination, filter="data")
    return destination / ARCHIVE_ROOT


def verify_workflow(committed: dict[str, Any]) -> None:
    value = load(WORKFLOW, "Stage-38 workflow metadata")
    require(value.get("databaseId") == RUN_ID and value.get("headSha") == RUN_COMMIT, "Stage-38 workflow identity changed")
    require(value.get("status") == "completed" and value.get("conclusion") == "success", "Stage-38 workflow is incomplete")
    jobs = value.get("jobs")
    require(isinstance(jobs, list) and {row.get("name") for row in jobs} == {"validate-control", "collection-window-sweep"}, "Stage-38 workflow jobs changed")
    require(all(row.get("status") == "completed" and row.get("conclusion") == "success" for row in jobs), "a Stage-38 job failed")
    production = next(row for row in jobs if row["name"] == "collection-window-sweep")
    expected = committed["workflow"]
    require(expected["run_id"] == RUN_ID and expected["head_sha"] == RUN_COMMIT and expected["url"] == value.get("url"), "Stage-38 workflow summary changed")
    require(expected["production_job_id"] == production.get("databaseId") and expected["production_job_url"] == production.get("url"), "Stage-38 production job changed")
    require(expected["workflow_wall_seconds"] == (timestamp(value["updatedAt"]) - timestamp(value["createdAt"])).total_seconds(), "Stage-38 workflow wall changed")
    require(expected["production_job_wall_seconds"] == (timestamp(production["completedAt"]) - timestamp(production["startedAt"])).total_seconds(), "Stage-38 production wall changed")


def verify_artifact(committed: dict[str, Any]) -> None:
    value = load(ARTIFACTS, "Stage-38 artifact metadata")
    rows = value.get("artifacts")
    require(value.get("total_count") == 1 and isinstance(rows, list) and len(rows) == 1, "Stage-38 artifact count changed")
    row = rows[0]
    require(row.get("id") == ARTIFACT_ID and row.get("name") == f"koblitz-stage38-algebraic-collection-{RUN_ID}", "Stage-38 artifact identity changed")
    require(row.get("digest") == ARTIFACT_DIGEST and row.get("expired") is False, "Stage-38 artifact digest or state changed")
    source = row.get("workflow_run", {})
    require(source.get("id") == RUN_ID and source.get("head_sha") == RUN_COMMIT, "Stage-38 artifact source changed")
    require(committed["artifact"] == {key: row[key] for key in ("id", "name", "size_in_bytes", "digest", "expires_at")}, "Stage-38 artifact summary changed")


def verify_result() -> dict[str, Any]:
    committed = load(RESULT / "verification.json", "committed Stage-38 result")
    require(committed.get("schema") == "koblitz_stage38_hosted_result.v1" and committed.get("status") == "hosted_six_window_algebraic_sweep_verified", "Stage-38 result label changed")
    require(committed.get("source_commit") == RUN_COMMIT, "Stage-38 source commit changed")
    require(committed.get("archive") == archive_identity(), "Stage-38 archive summary changed")
    verify_workflow(committed)
    verify_artifact(committed)
    with tempfile.TemporaryDirectory() as directory:
        root = extract(Path(directory))
        verification = stage38.verify(root / "stage38-build", root / "stage38-run")
        embedded = load(root / "stage38-verification.json", "embedded Stage-38 verification")
        require(embedded == verification == committed.get("verification"), "Stage-38 recomputed verification differs")
        run = load(root / "stage38-run/result.json", "Stage-38 run result")
    full = next(row for row in run["rows"] if row["collection_window"] is None)
    winner = next(row for row in run["rows"] if row["label"] == run["winner"]["label"])
    require(winner["collection_window"] == 149 and run["winner"]["selection_metric"] == "collection_and_logs_seconds", "Stage-38 selected window changed")
    selection = committed["selection"]
    require(selection["winner"] == winner and selection["full_scan"] == full, "Stage-38 selected rows changed")
    require(selection["collection_and_logs_speedup"] == full["collection_and_logs_seconds"] / winner["collection_and_logs_seconds"], "Stage-38 collection speedup changed")
    require(selection["five_target_ic_speedup"] == full["five_target_ic_seconds"] / winner["five_target_ic_seconds"], "Stage-38 five-target speedup changed")
    require(selection["lookup_reduction"] == full["summands_scanned"] / winner["summands_scanned"], "Stage-38 lookup reduction changed")
    require(committed.get("factor_base_algebraically_defined") is True and committed.get("target_subgroup_enumerated_for_factor_base") is False, "Stage-38 factor-base boundary changed")
    require(committed.get("target_scalars_constructed_or_supplied") is False and committed.get("factor_base_logs_known_by_construction") is False, "Stage-38 scalar-label boundary changed")
    require(committed.get("full_cost_gate_passed") is False and committed.get("koblitz_index_calculus_sota") is False, "Stage-38 claim widened")
    seal = load(RESULT / "result-seal.json", "Stage-38 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-38 result seal is invalid")
    inventory = custody.all_regular_inventory(RESULT, {"result-seal.json"})
    require(inventory == seal.get("inventory") and custody.canonical_sha256(inventory) == seal.get("inventory_sha256"), "Stage-38 result inventory changed")
    return committed


if __name__ == "__main__":
    try:
        print(json.dumps(verify_result(), indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, VerificationError, stage38.Stage38Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage38-result: {error}")

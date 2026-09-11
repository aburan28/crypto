#!/usr/bin/env python3
"""Verify the committed Stage-26 score and five-artifact archive."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import tarfile
import tempfile

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage26_affinity_cell as cell_tool
import score_koblitz_stage26_affinity as score_tool


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
SCORE_ROOT = STAGE / "stage-26-affinity-matrix-result-20260911"
WORKFLOW = STAGE / "stage-26-affinity-workflow-34632018379.json"
ARTIFACTS = STAGE / "stage-26-affinity-artifacts-34632018379.json"
ARCHIVE = STAGE / "koblitz-stage26-terminal-evidence-successor-01-20260911.tar.gz"
ARCHIVE_SHA = STAGE / "koblitz-stage26-terminal-evidence-successor-01-20260911.tar.gz.sha256"


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load(path: Path) -> dict:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"JSON object required: {path}")
    return value


def verify_archive_hash() -> str:
    fields = ARCHIVE_SHA.read_text().split()
    require(len(fields) == 2 and fields[1] == str(ARCHIVE.relative_to(REPO)), "Stage-26 checksum record changed")
    custody.require_hex64(fields[0], "Stage-26 archive SHA-256")
    actual = hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
    require(actual == fields[0], "Stage-26 archive SHA-256 changed")
    return actual


def extract(destination: Path) -> Path:
    with tarfile.open(ARCHIVE, "r:gz") as source:
        for member in source.getmembers():
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and ".." not in path.parts, "unsafe Stage-26 archive path")
            require(member.isfile() or member.isdir(), "Stage-26 archive contains a link or special file")
        source.extractall(destination, filter="data")
    root = destination / "koblitz-stage26-run-34632018379-download-20260911"
    require(root.is_dir(), "Stage-26 archive root changed")
    return root


def verify() -> dict:
    archive_sha256 = verify_archive_hash()
    seal = load(SCORE_ROOT / "score-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == score_tool.SEAL_SCHEMA
        and seal.get("status") == "score_frozen"
        and seal.get("workflow_run_id") == score_tool.EXPECTED_RUN
        and seal.get("workflow_commit") == score_tool.EXPECTED_COMMIT
        and claimed == custody.canonical_sha256(payload),
        "Stage-26 score seal is invalid",
    )
    inventory = custody.all_regular_inventory(SCORE_ROOT, {"score-seal.json"})
    require(inventory == seal.get("inventory") and custody.canonical_sha256(inventory) == seal.get("inventory_sha256"), "Stage-26 score inventory changed")
    score = load(SCORE_ROOT / "score.json")
    require(score.get("schema") == score_tool.SCHEMA and score.get("status") == "complete_verified_truth_scored_four_cell_panel", "Stage-26 score status changed")
    require(score.get("instances") == 160 and score.get("backend_rows") == 480, "Stage-26 scored counts changed")
    require(score.get("classification_counts") == {"inconclusive": 219, "true_negative": 120, "true_positive": 141}, "Stage-26 classifications changed")
    require(score.get("full_cost_gate_passed") is False and score.get("koblitz_index_calculus_sota") is False, "Stage-26 score widened its conclusion")
    cells = set(score.get("cells", {}))
    require(cells == {
        "n31-l5-m3-standard-a1-f0", "n31-l5-m3-ggmp-a0-f0",
        "n41-l5-m3-standard-a1-f0", "n59-l9-m3-standard-a1-f0",
    }, "Stage-26 score cell inventory changed")
    workflow = score_tool.workflow_accounting(WORKFLOW, cells)
    artifacts = score_tool.artifact_accounting(ARTIFACTS, cells)
    require(score.get("workflow") == workflow and score.get("artifacts") == artifacts, "Stage-26 workflow or artifact binding changed")
    require(workflow.get("workflow_wall_seconds") == 11088.0 and workflow.get("parallel_cell_job_span_seconds") == 10750.0, "Stage-26 campaign wall changed")
    resources = score.get("cell_outer_resources", {})
    require(resources.get("summed_total_core_seconds") == 19414.190243999998 and resources.get("summed_single_core_elapsed_seconds") == 19418.302146035, "Stage-26 cell resource totals changed")
    require(score.get("charged_total_core_seconds_available") == 20611.878126, "Stage-26 charged core total changed")
    with tempfile.TemporaryDirectory() as directory:
        root = extract(Path(directory))
        paths = {
            "n31-l5-m3-standard-a1-f0": root / "n31-standard/stage26-n31-l5-m3-standard-a1-f0",
            "n31-l5-m3-ggmp-a0-f0": root / "n31-ggmp/stage26-n31-l5-m3-ggmp-a0-f0",
            "n41-l5-m3-standard-a1-f0": root / "n41/stage26-n41-l5-m3-standard-a1-f0",
            "n59-l9-m3-standard-a1-f0": root / "n59/stage26-n59-l9-m3-standard-a1-f0",
        }
        require((root / "tools/koblitz_pdp_backend").is_file() and (root / "tools/wdsat_solver").is_file() and (root / "tools/cryptominisat5").is_file(), "Stage-26 tool archive is incomplete")
        for cell, path in paths.items():
            verification = cell_tool.verify(path)
            require(score["cells"][cell]["verification"] == verification, f"Stage-26 {cell} verification changed")
            require(score["cells"][cell]["result_seal_sha256"] == custody.sha256_file(path / "result-seal.json", "Stage-26 cell seal"), f"Stage-26 {cell} seal changed")
    require(score.get("matched_direct_mitm", {}).get("verification", {}).get("status") == "committed_score_and_four_cell_archive_verified", "matched direct-MITM binding changed")
    require(score.get("natural_relation_yield", {}).get("independent_internal_payload_replay", {}).get("status") == "independent_internal_factor_base_pair_oracle_and_witness_replay_verified", "relation-yield replay binding changed")
    require(score.get("n31_unknown_scalar_index_calculus_and_rho", {}).get("verification", {}).get("status") == "validated_factor_base_and_five_public_unknown_scalars_verified", "n=31 unknown-scalar binding changed")
    return {
        "schema": "koblitz_stage26_affinity_result_verification.v1",
        "status": "score_artifacts_and_four_cell_archive_verified",
        "archive_sha256": archive_sha256,
        "score_sha256": custody.sha256_file(SCORE_ROOT / "score.json", "Stage-26 score"),
        "score_inventory_sha256": seal["inventory_sha256"],
        "workflow_run_id": score_tool.EXPECTED_RUN,
        "workflow_commit": score_tool.EXPECTED_COMMIT,
        "cells": 4,
        "instances": 160,
        "backend_rows": 480,
        "classification_counts": score["classification_counts"],
        "workflow_wall_seconds": workflow["workflow_wall_seconds"],
        "parallel_cell_job_span_seconds": workflow["parallel_cell_job_span_seconds"],
        "charged_total_core_seconds_available": score["charged_total_core_seconds_available"],
        "truth_rescored_without_oracle": False,
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


if __name__ == "__main__":
    try:
        print(json.dumps(verify(), indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, VerificationError, cell_tool.Stage26Error, score_tool.Stage26ScoreError, custody.PhaseBError) as error:
        raise SystemExit(f"stage26-result: {error}")

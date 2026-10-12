#!/usr/bin/env python3
"""Verify the committed Stage-27 score metadata and terminal evidence archive."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import tarfile
import tempfile

import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_stage27_direct_mitm as stage27
import score_koblitz_stage27_direct_mitm as score_tool


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
SCORE_ROOT = STAGE / "stage-27-direct-mitm-result-20260911"
WORKFLOW = STAGE / "stage-27-direct-mitm-workflow-34633920325.json"
ARTIFACTS = STAGE / "stage-27-direct-mitm-artifacts-34633920325.json"
ARCHIVE = STAGE / "koblitz-stage27-terminal-evidence-successor-01-20260911.tar.gz"
ARCHIVE_SHA = STAGE / "koblitz-stage27-terminal-evidence-successor-01-20260911.tar.gz.sha256"


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def read_json(path: Path) -> dict:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"JSON object required: {path}")
    return value


def archive_sha256() -> str:
    fields = ARCHIVE_SHA.read_text().split()
    require(len(fields) == 2 and fields[1] == str(ARCHIVE.relative_to(REPO)), "Stage-27 archive checksum file changed")
    expected = fields[0]
    phase_b.require_hex64(expected, "Stage-27 archive SHA-256")
    actual = hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
    require(actual == expected, "Stage-27 archive SHA-256 changed")
    return actual


def extract(destination: Path) -> Path:
    with tarfile.open(ARCHIVE, "r:gz") as source:
        for member in source.getmembers():
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and ".." not in path.parts, "unsafe Stage-27 archive path")
            require(member.isfile() or member.isdir(), "Stage-27 archive contains a link or special file")
        source.extractall(destination, filter="data")
    root = destination / "koblitz-stage27-run-34633920325-download-20260911"
    require(root.is_dir(), "Stage-27 archive root changed")
    return root


def verify() -> dict:
    archive_digest = archive_sha256()
    score_seal = read_json(SCORE_ROOT / "score-seal.json")
    payload = dict(score_seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(score_seal.get("schema") == score_tool.SEAL_SCHEMA and claimed == phase_b.canonical_sha256(payload), "Stage-27 score seal is invalid")
    inventory = phase_b.all_regular_inventory(SCORE_ROOT, {"score-seal.json"})
    require(inventory == score_seal.get("inventory") and phase_b.canonical_sha256(inventory) == score_seal.get("inventory_sha256"), "Stage-27 score inventory changed")
    score = read_json(SCORE_ROOT / "score.json")
    require(score.get("schema") == score_tool.SCHEMA and score.get("status") == "complete_verified_truth_scored_direct_mitm_panel", "Stage-27 score status changed")
    require(score.get("instances") == 160 and score.get("outcomes") == 160, "Stage-27 score counts changed")
    require(score.get("classification_counts") == {"true_negative": 80, "true_positive": 80}, "Stage-27 classification counts changed")
    require(score.get("full_cost_gate_passed") is False and score.get("koblitz_index_calculus_sota") is False, "Stage-27 score widened its claim")
    cells = set(score.get("cells", {}))
    require(cells == packet_cells(), "Stage-27 score cell inventory changed")
    workflow = score_tool.workflow_receipt(WORKFLOW, cells)
    artifacts = score_tool.artifact_receipt(ARTIFACTS, cells)
    require(score.get("workflow") == workflow, "Stage-27 workflow summary changed")
    with tempfile.TemporaryDirectory() as directory:
        root = extract(Path(directory))
        result_paths = sorted(root.glob("*/stage27-*/result.json"))
        require(len(result_paths) == 4, "Stage-27 archive does not contain four cell results")
        seen = set()
        for result_path in result_paths:
            cell_root = result_path.parent
            result = read_json(result_path)
            cell = result.get("cell_id")
            require(cell in cells and cell not in seen, "Stage-27 archive cell identity changed")
            seen.add(cell)
            verification = stage27.verify(cell_root)
            require(score["cells"][cell]["verification"] == verification, "Stage-27 embedded cell verification changed")
            require(score["cells"][cell]["result_seal_sha256"] == phase_b.sha256_file(cell_root / "result-seal.json", "Stage-27 archived cell seal"), "Stage-27 archived cell seal differs from score")
            require(score["cells"][cell]["artifact"] == artifacts[cell], "Stage-27 artifact binding changed")
        require(seen == cells, "Stage-27 archive cell set changed")
    return {
        "schema": "koblitz_stage27_direct_mitm_result_verification.v1",
        "status": "committed_score_and_four_cell_archive_verified",
        "archive_sha256": archive_digest,
        "score_inventory_sha256": score_seal["inventory_sha256"],
        "workflow_run_id": score_tool.RUN_ID,
        "workflow_commit": score_tool.RUN_COMMIT,
        "cells": 4,
        "instances": 160,
        "outcomes": 160,
        "truth_rescored_without_oracle": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


def packet_cells() -> set[str]:
    return {
        "n31-l5-m3-standard-a1-f0",
        "n31-l5-m3-ggmp-a0-f0",
        "n41-l5-m3-standard-a1-f0",
        "n59-l9-m3-standard-a1-f0",
    }


if __name__ == "__main__":
    try:
        print(json.dumps(verify(), indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, VerificationError, stage27.Stage27Error, score_tool.Stage27ScoreError, phase_b.PhaseBError) as error:
        raise SystemExit(f"stage27-result: {error}")

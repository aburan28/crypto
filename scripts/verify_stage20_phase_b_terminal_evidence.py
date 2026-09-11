#!/usr/bin/env python3
"""Verify the additive, self-contained Stage-20 Phase-B terminal evidence bundle."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import stat
from typing import Any


REPO = Path(__file__).resolve().parents[1]
HERE = REPO / "research" / "sat_factor_base_review_20260908" / "continuation-05-sota-gates"
DEFAULT_BUNDLE = HERE / "stage-20-phase-b-terminal-evidence-successor-04-20260910"
SEAL_SCHEMA = "koblitz_pdp_phase_b_terminal_evidence_seal.v1"


class EvidenceError(RuntimeError):
    pass


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def regular_bytes(path: Path, context: str) -> bytes:
    metadata = path.lstat()
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
        raise EvidenceError(f"{context} must be a regular non-symlink file: {path}")
    if metadata.st_nlink != 1:
        raise EvidenceError(f"{context} must not be hard-linked: {path}")
    return path.read_bytes()


def identity(path: Path, context: str = "file") -> dict[str, Any]:
    data = regular_bytes(path, context)
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def bundle_identity(root: Path, relative: str, context: str) -> dict[str, Any]:
    value = identity(root / relative, context)
    value["path"] = relative
    return value


def read_json(path: Path, context: str) -> tuple[dict[str, Any], bytes]:
    data = regular_bytes(path, context)

    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, child in pairs:
            if key in value:
                raise EvidenceError(f"duplicate JSON key {key!r} in {context}")
            value[key] = child
        return value

    try:
        value = json.loads(data, object_pairs_hook=no_duplicates)
    except json.JSONDecodeError as error:
        raise EvidenceError(f"invalid JSON in {context}: {error}") from error
    if not isinstance(value, dict):
        raise EvidenceError(f"{context} must contain an object")
    return value, data


def inventory(root: Path) -> list[dict[str, Any]]:
    records = []
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root).as_posix()
        metadata = path.lstat()
        if stat.S_ISLNK(metadata.st_mode):
            raise EvidenceError(f"bundle contains symlink {relative}")
        if stat.S_ISDIR(metadata.st_mode):
            continue
        if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
            raise EvidenceError(f"bundle contains inadmissible file {relative}")
        if relative == "bundle-seal.json":
            continue
        data = path.read_bytes()
        records.append({"path": relative, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
    return records


def validate_self_hash(value: dict[str, Any], field: str, context: str) -> None:
    payload = {key: child for key, child in value.items() if key != field}
    if value.get(field) != canonical_sha256(payload):
        raise EvidenceError(f"{context} self-hash is invalid")


def validate_process_metrics(value: dict[str, Any], context: str) -> None:
    if value.get("returncode") != 0 or value.get("timed_out") is not False or value.get(
        "orphan_group_terminated"
    ) is not False:
        raise EvidenceError(f"{context} did not terminate cleanly")
    metrics = value.get("metrics", {})
    user = metrics.get("user_seconds")
    system = metrics.get("system_seconds")
    core = metrics.get("total_core_seconds")
    alias = metrics.get("single_core_seconds")
    if any(isinstance(item, bool) or not isinstance(item, (int, float)) for item in (user, system, core, alias)):
        raise EvidenceError(f"{context} CPU fields are invalid")
    if not math.isclose(core, user + system, rel_tol=0, abs_tol=1e-9) or alias != core:
        raise EvidenceError(f"{context} CPU fields do not balance")
    if metrics.get("meter") != "fresh-process getrusage(RUSAGE_CHILDREN)":
        raise EvidenceError(f"{context} meter identity changed")


def verify(bundle: Path) -> dict[str, Any]:
    bundle = bundle.resolve(strict=True)
    if bundle.is_symlink() or not bundle.is_dir():
        raise EvidenceError("bundle root must be a real directory")
    seal, _ = read_json(bundle / "bundle-seal.json", "bundle seal")
    if seal.get("schema") != SEAL_SCHEMA or seal.get("status") != "terminal_scored_evidence_frozen":
        raise EvidenceError("bundle seal is not terminal")
    validate_self_hash(seal, "seal_payload_sha256", "bundle seal")
    actual_inventory = inventory(bundle)
    if actual_inventory != seal.get("inventory") or canonical_sha256(actual_inventory) != seal.get(
        "inventory_sha256"
    ):
        raise EvidenceError("bundle inventory changed after sealing")

    archived_sources = seal.get("archived_sources", {})
    expected_archived_sources = {
        "report": bundle_identity(
            bundle, "source/STAGE20_PHASE_B_RESULTS_20260910.md", "archived result report"
        ),
        "summary": bundle_identity(
            bundle,
            "source/stage-20-phase-b-result-summary-20260910.json",
            "archived result summary",
        ),
        "scorer_source": bundle_identity(
            bundle, "source/score_koblitz_blind_pdp_phase_b.py", "archived scorer source"
        ),
        "process_meter_source": bundle_identity(
            bundle, "source/process_meter.py", "archived process meter source"
        ),
        "verifier_source": bundle_identity(
            bundle,
            "source/verify_stage20_phase_b_terminal_evidence.py",
            "archived evidence verifier source",
        ),
    }
    if archived_sources != expected_archived_sources:
        raise EvidenceError("archived result or source identity changed")

    score, score_bytes = read_json(bundle / "score/score.json", "score")
    score_seal, _ = read_json(bundle / "score/score-seal.json", "score seal")
    if (
        score_seal.get("schema") != "koblitz_pdp_phase_b_score_seal.v1"
        or score_seal.get("status") != "post_run_truth_scoring_complete"
        or score_seal.get("score_path") != "score.json"
        or score_seal.get("score_bytes") != len(score_bytes)
        or score_seal.get("score_sha256") != hashlib.sha256(score_bytes).hexdigest()
    ):
        raise EvidenceError("score seal does not bind the score")
    rows = score.get("rows")
    if (
        score.get("schema") != "koblitz_pdp_phase_b_score.v1"
        or score.get("selected_instances") != 160
        or score.get("backend_rows") != 480
        or score.get("full_panel_scored") is not True
        or score.get("full_cost_gate_passed") is not False
        or score.get("independent_external_reproduction_satisfied") is not False
        or not isinstance(rows, list)
        or len(rows) != 480
    ):
        raise EvidenceError("score uses an unexpected terminal schema or claim boundary")
    classifications = Counter(row.get("classification") for row in rows)
    expected_outcomes = {
        "true_positive": 141,
        "true_negative": 120,
        "inconclusive": 219,
    }
    if classifications != expected_outcomes:
        raise EvidenceError(f"scored outcome totals changed: {dict(classifications)}")
    if any(row.get("classification") in {"false_positive", "false_negative"} for row in rows):
        raise EvidenceError("score contains a false terminal classification")

    run_seal, _ = read_json(bundle / "run/run-seal.json", "run seal")
    validate_self_hash(run_seal, "seal_payload_sha256", "run seal")
    if (
        run_seal.get("schema") != "koblitz_pdp_phase_b_run_seal.v1"
        or run_seal.get("status") != "solver_outputs_frozen"
        or run_seal.get("selected_instance_count") != 160
        or run_seal.get("backend_outcomes") != 480
        or run_seal.get("full_panel_complete") is not True
        or score.get("phase_b_run_inventory_sha256") != run_seal.get("inventory_sha256")
    ):
        raise EvidenceError("run seal and score are not consistently bound")

    phase_a_seal, phase_a_seal_bytes = read_json(bundle / "inputs/phase-a-seal.json", "Phase-A seal")
    phase_a_protocol = regular_bytes(bundle / "inputs/phase-a-protocol.json", "Phase-A protocol")
    phase_b_protocol = regular_bytes(bundle / "inputs/phase-b-protocol.json", "Phase-B protocol")
    blind_bundle = regular_bytes(bundle / "inputs/blind-bundle.json", "blind bundle")
    solver_binding, _ = read_json(bundle / "inputs/solver-binding.json", "solver binding")
    if (
        hashlib.sha256(phase_a_seal_bytes).hexdigest() != score.get("phase_a_seal_sha256")
        or hashlib.sha256(phase_a_protocol).hexdigest() != phase_a_seal.get("protocol_sha256")
        or hashlib.sha256(blind_bundle).hexdigest() != phase_a_seal.get("blind_bundle_sha256")
        or hashlib.sha256(phase_b_protocol).hexdigest() != score.get("protocol_sha256")
        or solver_binding.get("blind_bundle_sha256") != phase_a_seal.get("blind_bundle_sha256")
    ):
        raise EvidenceError("protocol, Phase-A seal, or blind-input binding changed")

    scorer_metrics, _ = read_json(bundle / "accounting/scorer.metrics.json", "scorer metrics")
    validate_process_metrics(scorer_metrics, "scorer")
    scorer_command = scorer_metrics.get("command", [])
    if (
        not isinstance(scorer_command, list)
        or not any(str(item).endswith("score_koblitz_blind_pdp_phase_b.py") for item in scorer_command)
        or "--protocol" not in scorer_command
        or "--run-root" not in scorer_command
        or "--oracle-ledger" not in scorer_command
        or "--outer-metrics" not in scorer_command
    ):
        raise EvidenceError("scorer command receipt is incomplete")
    scorer_stdout, _ = read_json(bundle / "accounting/scorer.stdout", "scorer stdout")
    if scorer_stdout != score_seal or regular_bytes(bundle / "accounting/scorer.stderr", "scorer stderr"):
        raise EvidenceError("scorer output streams do not match the terminal score seal")

    expected_builds = {
        "rust": ("builds/rust-receipt.json", "koblitz_pdp_phase_b_tool_build_receipt.v1"),
        "cryptominisat": (
            "builds/cryptominisat-receipt.json",
            "koblitz_pdp_phase_b_tool_build_receipt.v1",
        ),
        "wdsat": ("builds/wdsat-receipt.json", "koblitz_pdp_phase_b_wdsat_build_receipt.v2"),
    }
    for name, (relative, schema) in expected_builds.items():
        receipt, _ = read_json(bundle / relative, f"{name} build receipt")
        if receipt.get("schema") != schema or receipt.get("status") != "completed":
            raise EvidenceError(f"{name} build receipt is not terminal")
        validate_self_hash(receipt, "receipt_payload_sha256", f"{name} build receipt")
        metrics, _ = read_json(bundle / f"builds/{name}-outer.metrics.json", f"{name} outer build")
        validate_process_metrics(metrics, f"{name} outer build")
    wdsat_seal, _ = read_json(bundle / "builds/wdsat-build-seal.json", "WDSat build seal")
    validate_self_hash(wdsat_seal, "seal_payload_sha256", "WDSat build seal")
    if wdsat_seal.get("status") != "build_frozen":
        raise EvidenceError("WDSat build seal is not terminal")

    summary, _ = read_json(
        bundle / "source/stage-20-phase-b-result-summary-20260910.json",
        "archived result summary",
    )
    catalog = summary.get("receipt_catalog")
    if not isinstance(catalog, list) or len(catalog) != 21:
        raise EvidenceError("campaign receipt catalog must contain 21 entries")
    totals = {"all_core": 0.0, "all_wall": 0.0, "success_core": 0.0, "success_wall": 0.0}
    for entry in catalog:
        label = entry.get("label")
        receipt_path = bundle / f"campaign-receipts/{label}.metrics.json"
        receipt, receipt_bytes = read_json(receipt_path, f"campaign receipt {label}")
        if hashlib.sha256(receipt_bytes).hexdigest() != entry.get("sha256"):
            raise EvidenceError(f"campaign receipt {label} hash changed")
        metrics = receipt.get("metrics", {})
        if (
            receipt.get("returncode") != entry.get("returncode")
            or not math.isclose(
                metrics.get("total_core_seconds"),
                entry.get("total_core_seconds"),
                rel_tol=0,
                abs_tol=1e-12,
            )
            or not math.isclose(
                metrics.get("wall_seconds"),
                entry.get("wall_seconds"),
                rel_tol=0,
                abs_tol=1e-12,
            )
            or metrics.get("peak_rss_bytes") != entry.get("peak_rss_bytes")
        ):
            raise EvidenceError(f"campaign receipt {label} differs from its catalog")
        totals["all_core"] += metrics["total_core_seconds"]
        totals["all_wall"] += metrics["wall_seconds"]
        if entry.get("included_in_successful_terminal_path") is True:
            totals["success_core"] += metrics["total_core_seconds"]
            totals["success_wall"] += metrics["wall_seconds"]
    expected_campaign = summary["campaign_accounting"]
    if (
        not math.isclose(totals["all_core"], expected_campaign["all_listed_receipts"]["total_core_seconds"], abs_tol=1e-9)
        or not math.isclose(totals["all_wall"], expected_campaign["all_listed_receipts"]["summed_receipt_wall_seconds"], abs_tol=1e-9)
        or not math.isclose(totals["success_core"], expected_campaign["successful_terminal_path"]["total_core_seconds"], abs_tol=1e-9)
        or not math.isclose(totals["success_wall"], expected_campaign["successful_terminal_path"]["summed_receipt_wall_seconds"], abs_tol=1e-9)
    ):
        raise EvidenceError("campaign accounting does not reconstruct from retained receipts")

    if seal.get("cross_bindings") != {
        "phase_b_implementation_commit": "62be8cf6dfdaf1aec54e0a8f290ed55b822229a0",
        "phase_a_implementation_commit": "47235e51b74a6fa8f3c8dc85d68bf886e20a1e88",
        "run_seal_file_sha256": identity(bundle / "run/run-seal.json")["sha256"],
        "run_inventory_sha256": run_seal["inventory_sha256"],
        "score_sha256": score_seal["score_sha256"],
        "score_seal_file_sha256": identity(bundle / "score/score-seal.json")["sha256"],
        "scorer_metrics_sha256": identity(bundle / "accounting/scorer.metrics.json")["sha256"],
        "wdsat_build_receipt_sha256": identity(bundle / "builds/wdsat-receipt.json")["sha256"],
    }:
        raise EvidenceError("terminal cross-bindings changed")
    if seal.get("claim_boundary") != {
        "full_cost_gate_passed": False,
        "independent_external_reproduction_satisfied": False,
        "koblitz_index_calculus_sota": False,
    }:
        raise EvidenceError("terminal evidence widened its claim boundary")
    return {
        "schema": "koblitz_pdp_phase_b_terminal_evidence_verification.v1",
        "status": "pass",
        "bundle_inventory_sha256": seal["inventory_sha256"],
        "score_sha256": score_seal["score_sha256"],
        "run_seal_file_sha256": identity(bundle / "run/run-seal.json")["sha256"],
        "outcomes": dict(classifications),
        "campaign_totals": totals,
        "claim_boundary": seal["claim_boundary"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    args = parser.parse_args()
    try:
        print(json.dumps(verify(args.bundle), indent=2, sort_keys=True))
    except (OSError, EvidenceError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()

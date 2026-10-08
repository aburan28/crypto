#!/usr/bin/env python3
"""Verify the committed Stage-165 fixed-X1 Hamming-weight schedule result."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import build_koblitz_phase_b_native_f4 as build_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4
import verify_koblitz_stage164_dense_f4 as stage164_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE = GATES / "stage-165-x1-weight-order-20260923"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
COMMIT = "341e46d02d49da94ca24837f7116c812ff2299d3"


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"JSON object required: {path}")
    return value


def verify_score() -> None:
    root = STAGE / "selected-score"
    seal = load(root / "score-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(root, {"score-seal.json"})
    require(
        seal.get("schema") == native_f4.SCORE_SEAL_SCHEMA
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and load(root / "score.json").get("classification_counts") == {"true_positive": 1},
        "Stage-165 score changed",
    )


def verify() -> dict[str, Any]:
    require(stage164_verify.verify().get("status") == "verified", "Stage-164 failed verification")
    seal = load(STAGE / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    require(
        seal.get("schema") == "koblitz_stage165_x1_weight_order_seal.v1"
        and seal.get("status") == "stage165_x1_weight_order_result_frozen"
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and phase_b.sha256_file(STAGE / "result.json", "Stage-165 result")
        == seal.get("result_sha256"),
        "Stage-165 seal changed",
    )
    result = load(STAGE / "result.json")
    require(
        result.get("schema") == "koblitz_stage165_x1_weight_order.v1"
        and result.get("status") == "complete_weight_order_true_positive"
        and result["instance"]["blind_instance_id"] == BLIND_ID
        and result["factor_base"]["factor_base_logs_known_by_construction"] is False
        and result["factor_base"]["target_subgroup_enumerated"] is False
        and result["mechanism"]["target_independent"] is True
        and result["mechanism"]["uses_discrete_log_labels"] is False
        and result["mechanism"]["uses_known_witness"] is False
        and "post hoc" in result["mechanism"]["selection_caveat"],
        "Stage-165 identity or schedule boundary changed",
    )
    run_seal, plan = native_f4.validate_run_seal(STAGE / "selected-run")
    require(
        run_seal["status"] == "native_f4_run_frozen_before_truth_scoring"
        and plan["selected_blind_instance_ids"] == [BLIND_ID]
        and plan["source_revision"]
        == {"commit": COMMIT, "dirty": False, "porcelain": []},
        "Stage-165 blind run changed",
    )
    tasks = sorted((STAGE / "selected-run/tasks").glob("*/task-result.json"))
    require(len(tasks) == 1, "Stage-165 task count changed")
    selected = load(tasks[0])["result"]
    report = selected["backend_report"]
    metrics = selected["metrics"]
    extra = report["cost"]["extra"]
    require(
        selected["status"] == "sat"
        and selected["source_model_valid"] is True
        and selected["source_witness_valid"] is True
        and report["solver_x1_order"] == "hamming_weight_then_mask"
        and report["solver_x1_order_target_independent"] is True
        and report["solver_x1_order_control"] == "PQ_F4_X1_ORDER=ascending"
        and report["fixed_x1_masks_visited"] == 85
        and report["fixed_x1_nonrational_skipped"] == 46
        and report["fixed_x1_systems_constructed"] == 39
        and report["fixed_x1_systems_completed"] == 39
        and report["solver_equations_total"] == 2_301
        and report["solver_terms_total"] == 2_295_166
        and report["solver_equations_blake3"]
        == "bd95f10ffd0ab09a2384a2a283fbb85f886df79cd33776b3a3a5a852ccf3bd13"
        and report["cost"]["ops"] == 16_821_055_616
        and extra["f4_calls"] == 39,
        "Stage-165 exact weighted-order result changed",
    )
    require(
        metrics["wall_seconds"] == 21.512089291005395
        and metrics["total_core_seconds"] == 21.484441
        and metrics["peak_rss_bytes"] == 948_502_528,
        "Stage-165 selected resources changed",
    )
    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path)
    build_tool.validate_receipt(
        receipt_path, receipt["binaries"], receipt["rust_source_objects"]
    )
    require(
        receipt["source_commit"] == COMMIT
        and receipt["source_clean"] is True
        and receipt["resources"]["summed_process_wall_seconds"] == 177.882565666012
        and receipt["resources"]["total_core_seconds"] == 174.176185,
        "Stage-165 clean build changed",
    )
    verify_score()
    comparisons = result["comparisons"]
    require(
        comparisons["same_binary_ascending_control"]["selected_wall_speedup"] > 1.60
        and comparisons["three_run_medians"]["wall_speedup"] > 1.60
        and comparisons["selected_over_stage164"]["wall_speedup"] > 1.69
        and comparisons["selected_over_same_host_direct_mitm"]["wall_ratio"] > 7
        and comparisons["selected_cold_build_plus_run"]["wall_seconds"]
        == 200.07609191601665,
        "Stage-165 comparisons changed",
    )
    campaign = result["campaign_accounting"]
    require(
        campaign["measured_lower_bound"]["summed_wall_seconds"] == 6902.64765482998
        and campaign["measured_lower_bound"]["total_core_seconds"] == 6817.038349
        and campaign["measured_lower_bound"]["peak_rss_bytes"] == 6_310_576_128
        and campaign["measured_lower_bound"]["resource_components"] == 150
        and campaign["complete_campaign_cost"] is None,
        "Stage-165 campaign accounting changed",
    )
    development = load(STAGE / "development/manifest.json")
    require(
        development["arms"]["clean-ascending"]["decision"]
        == "same_binary_ascending_control"
        and development["arms"]["weight-r1"]["decision"]
        == "rejected_stale_binary_environment_override_not_present",
        "Stage-165 development boundary changed",
    )
    for arm in development["arms"].values():
        for artifact in arm["artifacts"].values():
            path = STAGE / "development" / artifact["path"]
            require(
                path.stat().st_size == artifact["bytes"]
                and phase_b.sha256_file(path, "Stage-165 arm artifact")
                == artifact["sha256"],
                "Stage-165 arm artifact changed",
            )
    require(
        result["gates"]["all_seven_gates_passed"] is False
        and result["licensed_magma_f4_complete"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["full_cost_gate_passed"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "Stage-165 claim boundary widened",
    )
    gate_status = (GATES / "GATE_STATUS.md").read_text()
    scoreboard = (REPO / "docs/index-calculus-scoreboard.html").read_text()
    require(
        "Stage 165 target-independent sparse fixed-X1 order" in gate_status
        and "6,902.647655 sequential wall-seconds" in gate_status,
        "gate status lacks Stage 165",
    )
    require(
        "still 7.17&times; direct MITM" in scoreboard and "21.512089" in scoreboard,
        "scoreboard lacks Stage 165",
    )
    return {
        "schema": "koblitz_stage165_x1_weight_order_verification.v1",
        "status": "verified",
        "result_sha256": seal["result_sha256"],
        "inventory_sha256": seal["inventory_sha256"],
        "blind_instance_id": BLIND_ID,
        "source_commit": receipt["source_commit"],
        "selected_wall_seconds": metrics["wall_seconds"],
        "selected_core_seconds": metrics["total_core_seconds"],
        "selected_peak_rss_bytes": metrics["peak_rss_bytes"],
        "selected_word_xors": report["cost"]["ops"],
        "same_binary_ascending_speedup": comparisons["same_binary_ascending_control"][
            "selected_wall_speedup"
        ],
        "same_host_f4_over_mitm_wall": comparisons[
            "selected_over_same_host_direct_mitm"
        ]["wall_ratio"],
        "all_seven_gates_passed": False,
        "koblitz_index_calculus_sota": False,
    }


if __name__ == "__main__":
    try:
        print(json.dumps(verify(), indent=2, sort_keys=True))
    except (
        OSError,
        ValueError,
        KeyError,
        IndexError,
        VerificationError,
        phase_b.PhaseBError,
    ) as error:
        raise SystemExit(f"stage165-x1-order: {error}")

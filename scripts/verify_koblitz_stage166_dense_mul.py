#!/usr/bin/env python3
"""Verify the committed Stage-166 dense F4 monomial-product result."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import build_koblitz_phase_b_native_f4 as build_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4
import verify_koblitz_stage165_x1_order as stage165_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE = GATES / "stage-166-dense-product-cancellation-20260923"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
COMMIT = "ac3e385cec9c90b489017f6012ebe0b50ca06405"
BACKEND_SHA256 = "013d3525557a7e3c811f11b1b189d52aaafc782512c35e6c3a5f98947652ba00"


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
        and load(root / "score.json").get("classification_counts")
        == {"true_positive": 1},
        "Stage-166 score changed",
    )


def verify() -> dict[str, Any]:
    require(stage165_verify.verify().get("status") == "verified", "Stage-165 failed verification")
    seal = load(STAGE / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    require(
        seal.get("schema") == "koblitz_stage166_dense_product_cancellation_seal.v1"
        and seal.get("status") == "stage166_dense_product_result_frozen"
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and phase_b.sha256_file(STAGE / "result.json", "Stage-166 result")
        == seal.get("result_sha256"),
        "Stage-166 seal changed",
    )
    result = load(STAGE / "result.json")
    require(
        result.get("schema") == "koblitz_stage166_dense_product_cancellation.v1"
        and result.get("status") == "complete_dense_product_true_positive"
        and result["instance"]["blind_instance_id"] == BLIND_ID
        and result["factor_base"]["factor_base_logs_known_by_construction"] is False
        and result["factor_base"]["target_subgroup_enumerated"] is False
        and result["mechanism"]["uses_target_labels"] is False
        and result["mechanism"]["uses_discrete_log_labels"] is False
        and result["mechanism"]["uses_known_witness"] is False
        and result["mechanism"]["same_binary_control"] == "PQ_F4_DISABLE_DENSE_MUL=1",
        "Stage-166 identity or mechanism boundary changed",
    )
    run_seal, plan = native_f4.validate_run_seal(STAGE / "selected-run")
    require(
        run_seal["status"] == "native_f4_run_frozen_before_truth_scoring"
        and plan["selected_blind_instance_ids"] == [BLIND_ID]
        and plan["source_revision"]
        == {"commit": COMMIT, "dirty": False, "porcelain": []},
        "Stage-166 blind run changed",
    )
    tasks = sorted((STAGE / "selected-run/tasks").glob("*/task-result.json"))
    require(len(tasks) == 1, "Stage-166 task count changed")
    selected = load(tasks[0])["result"]
    report = selected["backend_report"]
    metrics = selected["metrics"]
    extra = report["cost"]["extra"]
    require(
        selected["status"] == "sat"
        and selected["source_model_valid"] is True
        and selected["source_witness_valid"] is True
        and report["solver_x1_order"] == "hamming_weight_then_mask"
        and report["solver_monomial_multiply"]
        == "epoch_dense_parity_cancel_then_encoded_u32_order_sort"
        and report["solver_monomial_multiply_control"] == "PQ_F4_DISABLE_DENSE_MUL=1"
        and report["fixed_x1_masks_visited"] == 85
        and report["fixed_x1_nonrational_skipped"] == 46
        and report["fixed_x1_systems_constructed"] == 39
        and report["fixed_x1_systems_completed"] == 39
        and report["solver_equations_total"] == 2_301
        and report["solver_terms_total"] == 2_295_166
        and report["solver_equations_blake3"]
        == "bd95f10ffd0ab09a2384a2a283fbb85f886df79cd33776b3a3a5a852ccf3bd13"
        and report["cost"]["ops"] == 16_821_055_616
        and extra["f4_calls"] == 39
        and extra["dense_mul_calls"] == 818_171
        and extra["dense_mul_input_terms"] == 651_502_431
        and extra["dense_mul_output_terms"] == 549_909_017
        and extra["dense_mul_cancelled_terms"] == 101_593_414
        and extra["dense_mul_scratch_bytes_max"] == 1_067_584,
        "Stage-166 exact dense-product result changed",
    )
    require(
        metrics["wall_seconds"] == 19.88775295900996
        and metrics["total_core_seconds"] == 19.858235999999998
        and metrics["peak_rss_bytes"] == 900_169_728,
        "Stage-166 selected resources changed",
    )
    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path)
    build_tool.validate_receipt(
        receipt_path, receipt["binaries"], receipt["rust_source_objects"]
    )
    require(
        receipt["source_commit"] == COMMIT
        and receipt["source_clean"] is True
        and receipt["binaries"]["backend"]["sha256"] == BACKEND_SHA256
        and receipt["resources"]["summed_process_wall_seconds"] == 172.554674500992
        and receipt["resources"]["total_core_seconds"] == 167.177895,
        "Stage-166 clean build changed",
    )
    verify_score()
    comparisons = result["comparisons"]
    require(
        comparisons["same_binary_sorted_product_control"]["selected_wall_speedup"]
        > 1.08
        and comparisons["same_binary_sorted_product_control"][
            "same_scientific_report_outside_declared_product_and_timing_fields"
        ]
        is True
        and comparisons["three_run_medians"]["wall_speedup"] > 1.07
        and comparisons["selected_over_stage165"]["wall_speedup"] > 1.08
        and 6.62
        < comparisons["selected_over_same_host_direct_mitm"]["wall_ratio"]
        < 6.64
        and comparisons["selected_cold_build_plus_run"]["wall_seconds"]
        == 193.16252025098908,
        "Stage-166 comparisons changed",
    )
    campaign = result["campaign_accounting"]
    require(
        campaign["measured_lower_bound"]["summed_wall_seconds"]
        == 7400.453979371927
        and campaign["measured_lower_bound"]["total_core_seconds"]
        == 7307.628306000001
        and campaign["measured_lower_bound"]["peak_rss_bytes"] == 6_310_576_128
        and campaign["measured_lower_bound"]["resource_components"] == 180
        and campaign["complete_campaign_cost"] is None,
        "Stage-166 campaign accounting changed",
    )
    development = load(STAGE / "development/manifest.json")
    require(
        development["arms"]["clean-control"]["decision"]
        == "same_binary_sorted_product_control"
        and development["arms"]["dense-v1-candidate-r1"]["decision"]
        == "rejected_two_array_dense_parity_wall_regression"
        and development["arms"]["m4ri-block8"]["decision"]
        == "retained_block_width_sweep_no_new_default",
        "Stage-166 development boundary changed",
    )
    for arm in development["arms"].values():
        for artifact in arm["artifacts"].values():
            path = STAGE / "development" / artifact["path"]
            require(
                path.stat().st_size == artifact["bytes"]
                and phase_b.sha256_file(path, "Stage-166 arm artifact")
                == artifact["sha256"],
                "Stage-166 arm artifact changed",
            )
    rejected = load(STAGE / "rejected-build-default-cache/manifest.json")
    require(
        rejected["status"] == "rejected_missing_dependency_in_default_offline_cache"
        and rejected["source_commit"] == COMMIT
        and rejected["resource_summary"]["resource_components"] == 6
        and load(STAGE / "rejected-build-default-cache/vendor.metrics.json")["returncode"]
        == 101,
        "Stage-166 failed-build accounting changed",
    )
    require(
        result["gates"]["all_seven_gates_passed"] is False
        and result["licensed_magma_f4_complete"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["full_cost_gate_passed"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "Stage-166 claim boundary widened",
    )
    gate_status = (GATES / "GATE_STATUS.md").read_text()
    scoreboard = (REPO / "docs/index-calculus-scoreboard.html").read_text()
    require(
        "Stage 166 dense F4 product cancellation" in gate_status
        and "7,400.453979 sequential wall-seconds" in gate_status,
        "gate status lacks Stage 166",
    )
    require(
        "still 6.63&times; direct MITM" in scoreboard and "19.887753" in scoreboard,
        "scoreboard lacks Stage 166",
    )
    return {
        "schema": "koblitz_stage166_dense_product_verification.v1",
        "status": "verified",
        "result_sha256": seal["result_sha256"],
        "inventory_sha256": seal["inventory_sha256"],
        "blind_instance_id": BLIND_ID,
        "source_commit": receipt["source_commit"],
        "selected_wall_seconds": metrics["wall_seconds"],
        "selected_core_seconds": metrics["total_core_seconds"],
        "selected_peak_rss_bytes": metrics["peak_rss_bytes"],
        "selected_word_xors": report["cost"]["ops"],
        "same_binary_control_speedup": comparisons[
            "same_binary_sorted_product_control"
        ]["selected_wall_speedup"],
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
        raise SystemExit(f"stage166-dense-product: {error}")

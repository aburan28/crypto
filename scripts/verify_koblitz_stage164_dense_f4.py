#!/usr/bin/env python3
"""Verify the committed Stage-164 dense native-F4 result."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import build_koblitz_phase_b_native_f4 as build_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4
import verify_koblitz_stage163_m4ri as stage163_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE = GATES / "stage-164-native-f4-dense-20260923"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
SOURCE_ID = "fb0ea6bb30be3dea841c609a2668939128b4b799d4c0e051b0b52dc6da05f0fa"
COMMIT = "9845f2df87c1aa759254c42739f11eef632ab2f2"


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
        "Stage-164 score changed",
    )


def verify() -> dict[str, Any]:
    require(stage163_verify.verify().get("status") == "verified", "Stage-163 failed verification")
    seal = load(STAGE / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    require(
        seal.get("schema") == "koblitz_stage164_native_f4_dense_seal.v1"
        and seal.get("status") == "stage164_native_f4_dense_result_frozen"
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and phase_b.sha256_file(STAGE / "result.json", "Stage-164 result")
        == seal.get("result_sha256"),
        "Stage-164 seal changed",
    )
    result = load(STAGE / "result.json")
    require(
        result.get("schema") == "koblitz_stage164_native_f4_dense.v1"
        and result.get("status") == "complete_dense_f4_true_positive"
        and result["instance"]["blind_instance_id"] == BLIND_ID
        and result["instance"]["source_instance_id"] == SOURCE_ID
        and result["factor_base"]["factor_base_logs_known_by_construction"] is False
        and result["factor_base"]["target_subgroup_enumerated"] is False,
        "Stage-164 identity or algebraic-base boundary changed",
    )
    run_seal, plan = native_f4.validate_run_seal(STAGE / "selected-run")
    require(
        run_seal["status"] == "native_f4_run_frozen_before_truth_scoring"
        and plan["selected_blind_instance_ids"] == [BLIND_ID]
        and plan["source_revision"]
        == {"commit": COMMIT, "dirty": False, "porcelain": []},
        "Stage-164 blind run changed",
    )
    tasks = sorted((STAGE / "selected-run/tasks").glob("*/task-result.json"))
    require(len(tasks) == 1, "Stage-164 task count changed")
    selected = load(tasks[0])["result"]
    report = selected["backend_report"]
    metrics = selected["metrics"]
    extra = report["cost"]["extra"]
    require(
        selected["status"] == "sat"
        and selected["source_model_valid"] is True
        and selected["source_witness_valid"] is True
        and report["solver_echelon_policy"].startswith("shape_selected_m4ri_default_block8")
        and report["solver_equations_blake3"]
        == "d6969c5f89b5be6a7eacbdac45799fc734a273c9034949a482723b4422504ae7"
        and report["solver_terms_total"] == 3_864_601
        and report["cost"]["op_unit"]
        == "word XORs (elimination, including M4RI tables)"
        and report["cost"]["ops"] == 27_264_366_281
        and extra["m4ri_matrices"] == 196
        and extra["m4ri_block_width_max"] == 8
        and extra["m4ri_table_word_xors"] == 2_509_320_857
        and extra["m4ri_trimmed_word_xors_avoided"] == 2_237_851_750
        and extra["pair_dense_select_calls"] == 269_235
        and extra["pair_cover_lookups"] == 186_105_200
        and extra["dense_column_matrices"] == 391
        and extra["divisor_tests"] == 29_134_673,
        "Stage-164 exact dense-F4 result changed",
    )
    require(
        metrics["wall_seconds"] == 36.45048054200015
        and metrics["total_core_seconds"] == 36.319424999999995
        and metrics["peak_rss_bytes"] == 946_388_992,
        "Stage-164 selected resources changed",
    )
    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path)
    build_tool.validate_receipt(
        receipt_path, receipt["binaries"], receipt["rust_source_objects"]
    )
    require(
        receipt["source_commit"] == COMMIT
        and receipt["source_clean"] is True
        and receipt["resources"]["summed_process_wall_seconds"] == 176.585437122994
        and receipt["resources"]["total_core_seconds"] == 172.196192,
        "Stage-164 clean build changed",
    )
    verify_score()
    comparisons = result["comparisons"]
    require(
        comparisons["same_binary_streaming_control"]["selected_wall_speedup"] > 1.32
        and comparisons["same_binary_batch_disabled"]["selected_wall_speedup"] > 1.06
        and comparisons["selected_over_stage163"]["wall_speedup"] > 1.50
        and comparisons["selected_over_same_host_direct_mitm"]["wall_ratio"] > 12
        and comparisons["selected_cold_dependency_build_plus_run"]["wall_seconds"]
        == 214.04431670698943,
        "Stage-164 comparisons changed",
    )
    campaign = result["campaign_accounting"]
    require(
        campaign["measured_lower_bound"]["summed_wall_seconds"] == 6464.947902204968
        and campaign["measured_lower_bound"]["total_core_seconds"] == 6383.990738
        and campaign["measured_lower_bound"]["peak_rss_bytes"] == 6_310_576_128
        and campaign["measured_lower_bound"]["resource_components"] == 132
        and campaign["complete_campaign_cost"] is None,
        "Stage-164 campaign accounting changed",
    )
    development = load(STAGE / "development/manifest.json")
    require(
        development["arms"]["clean-streaming"]["decision"]
        == "same_clean_binary_control"
        and development["arms"]["clean-batch-disabled"]["decision"]
        == "same_clean_binary_control"
        and development["arms"]["dependency-fetch"]["decision"]
        == "separately_charged_dependency_acquisition"
        and len(development["failed_clean_builds"]) == 2,
        "Stage-164 development boundary changed",
    )
    for arm in development["arms"].values():
        for artifact in arm["artifacts"].values():
            path = STAGE / "development" / artifact["path"]
            require(
                path.stat().st_size == artifact["bytes"]
                and phase_b.sha256_file(path, "Stage-164 arm artifact")
                == artifact["sha256"],
                "Stage-164 arm artifact changed",
            )
    for failed in development["failed_clean_builds"].values():
        root = STAGE / "development" / failed["path"]
        failed_inventory = phase_b.all_regular_inventory(root)
        require(
            failed_inventory == failed["inventory"]
            and phase_b.canonical_sha256(failed_inventory)
            == failed["inventory_sha256"],
            "Stage-164 failed-build inventory changed",
        )
    require(
        result["gates"]["all_seven_gates_passed"] is False
        and result["licensed_magma_f4_complete"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["full_cost_gate_passed"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "Stage-164 claim boundary widened",
    )
    gate_status = (GATES / "GATE_STATUS.md").read_text()
    scoreboard = (REPO / "docs/index-calculus-scoreboard.html").read_text()
    require(
        "Stage 164 dense native-F4 pipeline" in gate_status
        and "6,464.947902 sequential wall-seconds" in gate_status,
        "gate status lacks Stage 164",
    )
    require(
        "still 12.15&times; direct MITM" in scoreboard and "36.450481" in scoreboard,
        "scoreboard lacks Stage 164",
    )
    return {
        "schema": "koblitz_stage164_native_f4_dense_verification.v1",
        "status": "verified",
        "result_sha256": seal["result_sha256"],
        "inventory_sha256": seal["inventory_sha256"],
        "blind_instance_id": BLIND_ID,
        "source_commit": receipt["source_commit"],
        "selected_wall_seconds": metrics["wall_seconds"],
        "selected_core_seconds": metrics["total_core_seconds"],
        "selected_peak_rss_bytes": metrics["peak_rss_bytes"],
        "selected_word_xors": report["cost"]["ops"],
        "same_binary_streaming_speedup": comparisons["same_binary_streaming_control"][
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
        raise SystemExit(f"stage164-dense-f4: {error}")

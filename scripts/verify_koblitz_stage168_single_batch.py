#!/usr/bin/env python3
"""Verify the committed Stage-168 one-batch fixed-X1 result."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import build_koblitz_phase_b_native_f4 as build_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4
import verify_koblitz_stage167_parallel_x1 as stage167_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE = GATES / "stage-168-single-batch-39-20260923"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
COMMIT = "509e43c79bcf36425f4329a09c34fd20c75949ed"
BACKEND_SHA256 = "58767973f4bc1fcf30d5bbca231d1d5647e6a79f6282c55f1c3874a962820931"
FINGERPRINT = "bd95f10ffd0ab09a2384a2a283fbb85f886df79cd33776b3a3a5a852ccf3bd13"


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
        "Stage-168 score changed",
    )


def verify() -> dict[str, Any]:
    require(stage167_verify.verify().get("status") == "verified", "Stage-167 failed verification")
    seal = load(STAGE / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    require(
        seal.get("schema") == "koblitz_stage168_single_batch_39_seal.v1"
        and seal.get("status") == "stage168_single_batch_result_frozen"
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and phase_b.sha256_file(STAGE / "result.json", "Stage-168 result")
        == seal.get("result_sha256"),
        "Stage-168 seal changed",
    )
    result = load(STAGE / "result.json")
    require(
        result.get("schema") == "koblitz_stage168_single_batch_39.v1"
        and result.get("status") == "complete_one_batch_true_positive"
        and result["instance"]["blind_instance_id"] == BLIND_ID
        and result["factor_base"]["factor_base_logs_known_by_construction"] is False
        and result["factor_base"]["target_subgroup_enumerated"] is False
        and result["mechanism"]["selected_batch_size"] == 39
        and result["mechanism"]["selected_rayon_threads"] == 12
        and result["mechanism"]["uses_target_labels"] is False
        and result["mechanism"]["uses_discrete_log_labels"] is False
        and result["mechanism"]["uses_known_witness"] is False
        and "post hoc" in result["mechanism"]["selection_caveat"],
        "Stage-168 identity or claim boundary changed",
    )
    run_seal, plan = native_f4.validate_run_seal(STAGE / "selected-run")
    require(
        run_seal["status"] == "native_f4_run_frozen_before_truth_scoring"
        and plan["source_revision"]
        == {"commit": COMMIT, "dirty": False, "porcelain": []}
        and plan["selected_blind_instance_ids"] == [BLIND_ID]
        and plan["x1_batch_size"] == 39
        and plan["rayon_threads_requested"] == 12
        and plan["single_thread_requested"] is False,
        "Stage-168 run plan changed",
    )
    tasks = sorted((STAGE / "selected-run/tasks").glob("*/task-result.json"))
    require(len(tasks) == 1, "Stage-168 task count changed")
    selected = load(tasks[0])["result"]
    report = selected["backend_report"]
    metrics = result["selected_parallel_native_f4"]["metrics"]
    require(
        selected["status"] == "sat"
        and selected["source_model_valid"] is True
        and selected["source_witness_valid"] is True
        and report["solver_x1_batch_size"] == 39
        and report["solver_rayon_threads_requested"] == 12
        and report["solver_x1_batches_completed"] == 1
        and report["solver_x1_speculative_systems_completed"] == 0
        and report["fixed_x1_masks_visited"] == 85
        and report["fixed_x1_nonrational_skipped"] == 46
        and report["fixed_x1_systems_constructed"] == 39
        and report["fixed_x1_systems_completed"] == 39
        and report["solver_equations_total"] == 2_301
        and report["solver_terms_total"] == 2_295_166
        and report["solver_equations_blake3"] == FINGERPRINT
        and report["cost"]["ops"] == 16_821_055_616
        and report["cost"]["extra"]["f4_calls"] == 39,
        "Stage-168 exact result changed",
    )
    require(
        metrics["wall_seconds"] == 4.987009749995195
        and metrics["total_core_seconds"] == 28.179809
        and metrics["single_core_seconds"] is None
        and metrics["peak_rss_bytes"] == 2_763_997_184,
        "Stage-168 selected resources changed",
    )
    summary = load(STAGE / "selected-run/run-summary.json")
    require(
        summary["f4_process_resources"]["single_core_seconds"] is None
        and summary["f4_process_resources"]["summed_process_wall_seconds"]
        == 4.987009749995
        and summary["f4_process_resources"]["total_core_seconds"] == 28.179809,
        "Stage-168 process accounting changed",
    )
    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path)
    build_tool.validate_receipt(
        receipt_path, receipt["binaries"], receipt["rust_source_objects"]
    )
    require(
        receipt["source_commit"] == COMMIT
        and receipt["source_clean"] is True
        and receipt["binaries"]["backend"]["sha256"] == BACKEND_SHA256,
        "Stage-168 build changed",
    )
    verify_score()
    comparisons = result["comparisons"]
    require(
        comparisons["selected_over_stage167_parallel"]["wall_speedup"] > 1.12
        and comparisons["selected_over_stage167_parallel"]["core_ratio"] < 0.78
        and comparisons["selected_over_stage167_parallel"]["rss_ratio"] < 0.92
        and comparisons["selected_over_single"]["wall_speedup"] > 3.97
        and 1.66
        < comparisons["selected_over_same_host_direct_mitm"]["wall_ratio"]
        < 1.67
        and comparisons["three_run_median"]["wall_seconds"]
        == 4.744984999997541
        and comparisons["cold_build_plus_selected_run"]["wall_seconds"]
        == 172.14679829198545
        and comparisons["cold_build_plus_selected_run"][
            "build_already_charged_in_stage167_campaign"
        ]
        is True,
        "Stage-168 comparisons changed",
    )
    campaign = result["campaign_accounting"]
    require(
        campaign["reused_build_already_charged"] is True
        and campaign["measured_lower_bound"]["summed_wall_seconds"]
        == 7988.745489327874
        and campaign["measured_lower_bound"]["total_core_seconds"]
        == 8381.199829000003
        and campaign["measured_lower_bound"]["peak_rss_bytes"] == 6_310_576_128
        and campaign["measured_lower_bound"]["resource_components"] == 223
        and campaign["complete_campaign_cost"] is None,
        "Stage-168 campaign accounting changed",
    )
    development = load(STAGE / "development/manifest.json")
    require(
        development["arms"]["batch39-t12-r1"]["decision"]
        == "selected_thread_count_repeat"
        and development["arms"]["batch39-t10-r1"]["decision"]
        == "thread_count_sweep_control",
        "Stage-168 development boundary changed",
    )
    for arm in development["arms"].values():
        require(arm["valid_single_core_seconds"] is None, "parallel arm invented single-core time")
        for artifact in arm["artifacts"].values():
            path = STAGE / "development" / artifact["path"]
            require(
                path.stat().st_size == artifact["bytes"]
                and phase_b.sha256_file(path, "Stage-168 arm artifact")
                == artifact["sha256"],
                "Stage-168 arm artifact changed",
            )
    require(
        result["gates"]["all_seven_gates_passed"] is False
        and result["licensed_magma_f4_complete"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["full_cost_gate_passed"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "Stage-168 claim boundary widened",
    )
    gate_status = (GATES / "GATE_STATUS.md").read_text()
    scoreboard = (REPO / "docs/index-calculus-scoreboard.html").read_text()
    require(
        "Stage 168 one deterministic batch of 39 systems" in gate_status
        and "7,988.745489 sequential wall-seconds" in gate_status,
        "gate status lacks Stage 168",
    )
    require(
        "parallel F4 still 1.66&times; direct MITM" in scoreboard
        and "4.987010" in scoreboard,
        "scoreboard lacks Stage 168",
    )
    return {
        "schema": "koblitz_stage168_single_batch_verification.v1",
        "status": "verified",
        "result_sha256": seal["result_sha256"],
        "inventory_sha256": seal["inventory_sha256"],
        "blind_instance_id": BLIND_ID,
        "source_commit": receipt["source_commit"],
        "parallel_wall_seconds": metrics["wall_seconds"],
        "parallel_total_core_seconds": metrics["total_core_seconds"],
        "parallel_single_core_seconds": None,
        "parallel_peak_rss_bytes": metrics["peak_rss_bytes"],
        "stage167_wall_speedup": comparisons["selected_over_stage167_parallel"][
            "wall_speedup"
        ],
        "same_host_parallel_f4_over_mitm_wall": comparisons[
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
        raise SystemExit(f"stage168-single-batch: {error}")

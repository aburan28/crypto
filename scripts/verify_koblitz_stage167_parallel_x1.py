#!/usr/bin/env python3
"""Verify the committed Stage-167 deterministic parallel fixed-X1 result."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import build_koblitz_phase_b_native_f4 as build_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4
import verify_koblitz_stage166_dense_mul as stage166_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE = GATES / "stage-167-parallel-x1-batches-20260923"
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


def verify_score(root: Path, context: str) -> None:
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
        f"Stage-167 {context} score changed",
    )


def only_result(root: Path) -> dict[str, Any]:
    tasks = sorted((root / "tasks").glob("*/task-result.json"))
    require(len(tasks) == 1, "Stage-167 task count changed")
    task = load(tasks[0])
    require(
        task.get("blind_instance_id") == BLIND_ID,
        "Stage-167 task identity changed",
    )
    return task["result"]


def verify() -> dict[str, Any]:
    require(stage166_verify.verify().get("status") == "verified", "Stage-166 failed verification")
    seal = load(STAGE / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    require(
        seal.get("schema") == "koblitz_stage167_parallel_x1_batches_seal.v1"
        and seal.get("status") == "stage167_parallel_x1_result_frozen"
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and phase_b.sha256_file(STAGE / "result.json", "Stage-167 result")
        == seal.get("result_sha256"),
        "Stage-167 seal changed",
    )
    result = load(STAGE / "result.json")
    require(
        result.get("schema") == "koblitz_stage167_parallel_x1_batches.v1"
        and result.get("status") == "complete_parallel_and_single_true_positive"
        and result["instance"]["blind_instance_id"] == BLIND_ID
        and result["factor_base"]["factor_base_logs_known_by_construction"] is False
        and result["factor_base"]["target_subgroup_enumerated"] is False
        and result["mechanism"]["selected_batch_size"] == 13
        and result["mechanism"]["selected_rayon_threads"] == 13
        and result["mechanism"]["uses_target_labels"] is False
        and result["mechanism"]["uses_discrete_log_labels"] is False
        and result["mechanism"]["uses_known_witness"] is False
        and "post hoc" in result["mechanism"]["selection_caveat"],
        "Stage-167 identity or claim boundary changed",
    )

    parallel_seal, parallel_plan = native_f4.validate_run_seal(
        STAGE / "selected-parallel-run"
    )
    single_seal, single_plan = native_f4.validate_run_seal(STAGE / "single-thread-run")
    require(
        parallel_seal["status"] == single_seal["status"]
        == "native_f4_run_frozen_before_truth_scoring"
        and parallel_plan["source_revision"]
        == single_plan["source_revision"]
        == {"commit": COMMIT, "dirty": False, "porcelain": []}
        and parallel_plan["x1_batch_size"] == 13
        and parallel_plan["rayon_threads_requested"] == 13
        and parallel_plan["single_thread_requested"] is False
        and single_plan["x1_batch_size"] == 1
        and single_plan["rayon_threads_requested"] == 1
        and single_plan["single_thread_requested"] is True,
        "Stage-167 run plans changed",
    )
    parallel = only_result(STAGE / "selected-parallel-run")
    single = only_result(STAGE / "single-thread-run")
    parallel_report = parallel["backend_report"]
    single_report = single["backend_report"]
    for report, batch, threads, batches in (
        (parallel_report, 13, 13, 3),
        (single_report, 1, 1, 39),
    ):
        require(
            report["solver_x1_batch_size"] == batch
            and report["solver_rayon_threads_requested"] == threads
            and report["solver_x1_batches_completed"] == batches
            and report["solver_x1_speculative_systems_completed"] == 0
            and report["fixed_x1_masks_visited"] == 85
            and report["fixed_x1_systems_completed"] == 39
            and report["solver_equations_blake3"] == FINGERPRINT
            and report["cost"]["ops"] == 16_821_055_616
            and report["cost"]["extra"]["f4_calls"] == 39
            and report["witness_points"] == parallel_report["witness_points"],
            "Stage-167 exact report changed",
        )
    require(
        parallel["status"] == single["status"] == "sat"
        and parallel["source_model_valid"] is single["source_model_valid"] is True
        and parallel["source_witness_valid"] is single["source_witness_valid"] is True,
        "Stage-167 validation terminal changed",
    )
    parallel_metrics = result["selected_parallel_native_f4"]["metrics"]
    single_metrics = result["single_thread_control"]["metrics"]
    require(
        parallel_metrics["wall_seconds"] == 5.589951040994492
        and parallel_metrics["total_core_seconds"] == 36.302718
        and parallel_metrics["single_core_seconds"] is None
        and parallel_metrics["peak_rss_bytes"] == 3_033_563_136
        and single_metrics["wall_seconds"] == 19.82845924999856
        and single_metrics["total_core_seconds"] == 19.805719000000003
        and single_metrics["single_core_seconds"] == 19.805719000000003
        and single_metrics["peak_rss_bytes"] == 904_085_504,
        "Stage-167 resources changed",
    )
    parallel_summary = load(STAGE / "selected-parallel-run/run-summary.json")
    single_summary = load(STAGE / "single-thread-run/run-summary.json")
    require(
        parallel_summary["f4_process_resources"]["single_core_seconds"] is None
        and single_summary["f4_process_resources"]["single_core_seconds"]
        == 19.805719,
        "Stage-167 single-core semantics changed",
    )

    receipt_path = STAGE / "selected-parallel-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path)
    build_tool.validate_receipt(
        receipt_path, receipt["binaries"], receipt["rust_source_objects"]
    )
    require(
        receipt["source_commit"] == COMMIT
        and receipt["source_clean"] is True
        and receipt["binaries"]["backend"]["sha256"] == BACKEND_SHA256
        and receipt["resources"]["summed_process_wall_seconds"]
        == 166.910593999986
        and receipt["resources"]["total_core_seconds"] == 163.261594,
        "Stage-167 clean build changed",
    )
    verify_score(STAGE / "selected-parallel-score", "parallel")
    verify_score(STAGE / "single-thread-score", "single")

    comparisons = result["comparisons"]
    require(
        comparisons["parallel_over_single"]["wall_speedup"] > 3.54
        and comparisons["parallel_over_single"]["core_ratio"] > 1.83
        and comparisons["three_run_medians"]["wall_speedup"] > 3.56
        and 1.86
        < comparisons["parallel_over_same_host_direct_mitm"]["wall_ratio"]
        < 1.87
        and comparisons["cold_build_plus_parallel_run"]["wall_seconds"]
        == 173.19033541699142
        and comparisons["cold_build_plus_parallel_run"]["total_core_seconds"]
        == 199.81424900000002,
        "Stage-167 comparisons changed",
    )
    campaign = result["campaign_accounting"]
    require(
        campaign["measured_lower_bound"]["summed_wall_seconds"]
        == 7948.233269577871
        and campaign["measured_lower_bound"]["total_core_seconds"]
        == 8144.184517000002
        and campaign["measured_lower_bound"]["peak_rss_bytes"] == 6_310_576_128
        and campaign["measured_lower_bound"]["resource_components"] == 215
        and campaign["complete_campaign_cost"] is None,
        "Stage-167 campaign accounting changed",
    )
    development = load(STAGE / "development/manifest.json")
    require(
        development["arms"]["pattern-tail-candidate-r1"]["decision"]
        == "rejected_pattern_tail_variable_loop_regression"
        and development["arms"]["x1-stale-batch2-r1"]["decision"]
        == "rejected_stale_binary_batch_override_absent"
        and development["arms"]["x1-batch13-r1"]["decision"]
        == "parallel_batch13_development_repeat",
        "Stage-167 development boundary changed",
    )
    for arm in development["arms"].values():
        for artifact in arm["artifacts"].values():
            path = STAGE / "development" / artifact["path"]
            require(
                path.stat().st_size == artifact["bytes"]
                and phase_b.sha256_file(path, "Stage-167 arm artifact")
                == artifact["sha256"],
                "Stage-167 arm artifact changed",
            )
    require(
        result["gates"]["all_seven_gates_passed"] is False
        and result["licensed_magma_f4_complete"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["full_cost_gate_passed"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "Stage-167 claim boundary widened",
    )
    gate_status = (GATES / "GATE_STATUS.md").read_text()
    scoreboard = (REPO / "docs/index-calculus-scoreboard.html").read_text()
    require(
        "Stage 167 deterministic parallel fixed-X1 batches" in gate_status
        and "7,948.233270 sequential wall-seconds" in gate_status,
        "gate status lacks Stage 167",
    )
    require(
        "parallel F4 still 1.86&times; direct MITM" in scoreboard
        and "5.589951" in scoreboard,
        "scoreboard lacks Stage 167",
    )
    return {
        "schema": "koblitz_stage167_parallel_x1_verification.v1",
        "status": "verified",
        "result_sha256": seal["result_sha256"],
        "inventory_sha256": seal["inventory_sha256"],
        "blind_instance_id": BLIND_ID,
        "source_commit": receipt["source_commit"],
        "parallel_wall_seconds": parallel_metrics["wall_seconds"],
        "parallel_total_core_seconds": parallel_metrics["total_core_seconds"],
        "parallel_single_core_seconds": None,
        "parallel_peak_rss_bytes": parallel_metrics["peak_rss_bytes"],
        "single_wall_seconds": single_metrics["wall_seconds"],
        "single_core_seconds": single_metrics["single_core_seconds"],
        "single_peak_rss_bytes": single_metrics["peak_rss_bytes"],
        "same_binary_wall_speedup": comparisons["parallel_over_single"][
            "wall_speedup"
        ],
        "same_host_parallel_f4_over_mitm_wall": comparisons[
            "parallel_over_same_host_direct_mitm"
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
        raise SystemExit(f"stage167-parallel-x1: {error}")

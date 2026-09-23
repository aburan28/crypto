#!/usr/bin/env python3
"""Verify the committed Stage-170 parallel fixed-X1 construction result."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import build_koblitz_phase_b_native_f4 as build_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4
import verify_koblitz_stage169_fresh_holdout as stage169_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE = GATES / "stage-170-parallel-fixed-x1-construction-20260923"
BLIND_ID = "b-421e22a9c1c3b9d56396c8bbd0e46185bbebc32de0306f2bee6e9585703c2be4"
SOURCE_SYSTEM_ID = "x-af85283c5728becd64890c6eaa654ccb164cb499fcc4cde7b5cfd1484b4e0efa"
COMMIT = "080d3c2bfe9fefaf86bd10c8279a7fb9ebe000a0"
BACKEND_SHA256 = "5515eb70f24c195172707dfaa0dce99bd80bb5e21958747af6b8221597ec8b80"
FINGERPRINT = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
WORD_XORS = 99_199_976_264


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"JSON object required: {path}")
    return value


def validate_report(report: dict[str, Any], context: str) -> None:
    require(
        report.get("status") == "unsat"
        and report.get("exhaustive") is True
        and report.get("fixed_x1_masks_visited") == 512
        and report.get("fixed_x1_nonrational_skipped") == 270
        and report.get("fixed_x1_systems_constructed") == 242
        and report.get("fixed_x1_systems_completed") == 242
        and report.get("solver_equations_total") == 14_278
        and report.get("solver_terms_total") == 14_515_915
        and report.get("solver_equations_blake3") == FINGERPRINT
        and report.get("algebraic_roots") == 0
        and report.get("witness_points") is None
        and report["cost"]["ops"] == WORD_XORS
        and report["cost"]["extra"]["f4_calls"] == 242,
        f"{context} scientific result changed",
    )


def verify_score() -> None:
    root = STAGE / "selected-score"
    seal = load(root / "score-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(root, {"score-seal.json"})
    score = load(root / "score.json")
    require(
        seal.get("schema") == native_f4.SCORE_SEAL_SCHEMA
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and score.get("classification_counts") == {"true_negative": 1}
        and score["rows"][0]["blind_instance_id"] == BLIND_ID
        and score["rows"][0]["source_system_id"] == SOURCE_SYSTEM_ID
        and score["rows"][0]["target_class"] == "nondecomposable",
        "Stage-170 score changed",
    )


def verify() -> dict[str, Any]:
    require(stage169_verify.verify().get("status") == "verified", "Stage-169 failed verification")
    seal = load(STAGE / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    require(
        seal.get("schema")
        == "koblitz_stage170_parallel_fixed_x1_construction_seal.v1"
        and seal.get("status") == "stage170_parallel_construction_result_frozen"
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and phase_b.sha256_file(STAGE / "result.json", "Stage-170 result")
        == seal.get("result_sha256"),
        "Stage-170 seal changed",
    )
    result = load(STAGE / "result.json")
    require(
        result.get("schema") == "koblitz_stage170_parallel_fixed_x1_construction.v1"
        and result.get("status") == "complete_parallel_construction_true_negative"
        and result["instance"]["blind_instance_id"] == BLIND_ID
        and result["factor_base"]["factor_base_logs_known_by_construction"] is False
        and result["factor_base"]["target_subgroup_enumerated"] is False
        and result["mechanism"]["target_independent"] is True
        and result["mechanism"]["uses_truth_or_witness"] is False
        and result["mechanism"]["preserves_batch_and_mask_order"] is True,
        "Stage-170 identity or mechanism changed",
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
        "Stage-170 selected run plan changed",
    )
    tasks = sorted((STAGE / "selected-run/tasks").glob("*/task-result.json"))
    require(len(tasks) == 1, "Stage-170 task count changed")
    task = load(tasks[0])
    require(
        task["blind_instance_id"] == BLIND_ID
        and task["source_system_id"] == SOURCE_SYSTEM_ID,
        "Stage-170 task identity changed",
    )
    report = task["result"]["backend_report"]
    validate_report(report, "selected run")
    selected = result["selected_native_f4"]
    metrics = selected["metrics"]
    require(
        report["solver_parallel_construction"] is True
        and report["solver_x1_batch_size"] == 39
        and report["solver_rayon_threads_requested"] == 12
        and report["solver_x1_batches_completed"] == 7
        and metrics["wall_seconds"] == 17.131543916999362
        and metrics["total_core_seconds"] == 176.171665
        and metrics["single_core_seconds"] is None
        and metrics["peak_rss_bytes"] == 2_904_489_984
        and selected["construction_wall_seconds"] == 0.528025834,
        "Stage-170 selected resources changed",
    )
    summary = load(STAGE / "selected-run/run-summary.json")
    require(
        summary["f4_process_resources"]["single_core_seconds"] is None,
        "parallel selected run invented single-core time",
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
        "Stage-170 build changed",
    )
    verify_score()

    control = result["same_binary_serial_control"]
    posthoc = result["posthoc_single_target_ceiling"]
    require(
        control["binary_sha256"] == BACKEND_SHA256
        and control["parallel_construction"] is False
        and control["batch_size"] == 39
        and control["rayon_threads"] == 12
        and control["metrics"]["wall_seconds"] == 20.796178000004147
        and control["metrics"]["total_core_seconds"] == 172.00142200000002
        and control["metrics"]["single_core_seconds"] is None
        and control["construction_wall_seconds"] == 4.62954329
        and posthoc["parallel_construction"] is True
        and posthoc["batch_size"] == 64
        and posthoc["rayon_threads"] == 14
        and posthoc["metrics"]["wall_seconds"] == 14.661977374998969
        and posthoc["metrics"]["single_core_seconds"] is None,
        "Stage-170 clean controls changed",
    )
    for arm in (control, posthoc):
        for artifact in arm["artifacts"].values():
            path = STAGE / "clean-controls" / artifact["path"]
            require(
                path.stat().st_size == artifact["bytes"]
                and phase_b.sha256_file(path, "Stage-170 clean artifact")
                == artifact["sha256"],
                "Stage-170 clean artifact changed",
            )

    comparisons = result["comparisons"]
    require(
        1.213 < comparisons["clean_same_binary"]["wall_speedup"] < 1.215
        and comparisons["clean_same_binary"]["construction_speedup"] > 8.7
        and comparisons["three_run_medians"]["wall_speedup"] > 1.2
        and comparisons["selected_over_stage169_preregistered_baseline"][
            "wall_speedup"
        ]
        > 1.19
        and 5.80 < comparisons["selected_over_direct_mitm"]["wall_ratio"] < 5.81
        and comparisons["posthoc_over_direct_mitm"]["wall_ratio"] > 4.9
        and comparisons["cold_build_plus_selected"]["wall_seconds"]
        == 186.05967537601788,
        "Stage-170 comparisons changed",
    )
    campaign = result["campaign_accounting"]
    require(
        campaign["stage170_incremental_measured"]["resource_components"] == 35
        and campaign["stage170_incremental_measured"]["summed_wall_seconds"]
        == 546.884947297047
        and campaign["stage170_incremental_measured"]["total_core_seconds"]
        == 3621.256367
        and campaign["measured_lower_bound"]["resource_components"] == 264
        and campaign["measured_lower_bound"]["summed_wall_seconds"]
        == 8808.374611916948
        and campaign["measured_lower_bound"]["total_core_seconds"]
        == 12423.234481000003
        and campaign["measured_lower_bound"]["peak_rss_bytes"] == 6_310_576_128
        and campaign["complete_campaign_cost"] is None,
        "Stage-170 campaign accounting changed",
    )
    rejected = result["rejected_build"]
    require(
        rejected["status"] == "rejected_missing_dependency_in_default_offline_cache"
        and rejected["resource_summary"]["resource_components"] == 6,
        "Stage-170 rejected build changed",
    )
    require(
        result["gates"]["all_seven_gates_passed"] is False
        and result["licensed_magma_f4_complete"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["full_cost_gate_passed"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "Stage-170 claim boundary widened",
    )
    gate_status = (GATES / "GATE_STATUS.md").read_text()
    scoreboard = (REPO / "docs/index-calculus-scoreboard.html").read_text()
    require(
        "Stage 170 parallel fixed-X1 construction" in gate_status
        and "12,423.234481 core-seconds" in gate_status,
        "gate status lacks Stage 170",
    )
    require(
        "parallel fixed-X1 construction" in scoreboard
        and "17.131544" in scoreboard,
        "scoreboard lacks Stage 170",
    )
    return {
        "schema": "koblitz_stage170_parallel_construction_verification.v1",
        "status": "verified",
        "result_sha256": seal["result_sha256"],
        "inventory_sha256": seal["inventory_sha256"],
        "blind_instance_id": BLIND_ID,
        "source_commit": receipt["source_commit"],
        "classification": "true_negative",
        "selected_wall_seconds": metrics["wall_seconds"],
        "selected_total_core_seconds": metrics["total_core_seconds"],
        "selected_single_core_seconds": None,
        "selected_peak_rss_bytes": metrics["peak_rss_bytes"],
        "serial_control_wall_seconds": control["metrics"]["wall_seconds"],
        "clean_wall_speedup": comparisons["clean_same_binary"]["wall_speedup"],
        "posthoc_wall_seconds": posthoc["metrics"]["wall_seconds"],
        "direct_mitm_wall_ratio": comparisons["selected_over_direct_mitm"][
            "wall_ratio"
        ],
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
        raise SystemExit(f"stage170-parallel-construction: {error}")

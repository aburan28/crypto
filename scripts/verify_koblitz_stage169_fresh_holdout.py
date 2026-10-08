#!/usr/bin/env python3
"""Verify the committed Stage-169 preregistered fresh-target holdout."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import build_koblitz_phase_b_native_f4 as build_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4
import verify_koblitz_stage168_single_batch as stage168_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE = GATES / "stage-169-fresh-target-holdout-20260923"
BLIND_ID = "b-421e22a9c1c3b9d56396c8bbd0e46185bbebc32de0306f2bee6e9585703c2be4"
SOURCE_SYSTEM_ID = "x-af85283c5728becd64890c6eaa654ccb164cb499fcc4cde7b5cfd1484b4e0efa"
COMMIT = "509e43c79bcf36425f4329a09c34fd20c75949ed"
BACKEND_SHA256 = "58767973f4bc1fcf30d5bbca231d1d5647e6a79f6282c55f1c3874a962820931"


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
    root = STAGE / "holdout-score"
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
        and score["classification_counts"] == {"true_negative": 1}
        and score["rows"][0]["blind_instance_id"] == BLIND_ID
        and score["rows"][0]["source_system_id"] == SOURCE_SYSTEM_ID
        and score["rows"][0]["target_class"] == "nondecomposable",
        "Stage-169 score changed",
    )


def verify() -> dict[str, Any]:
    require(stage168_verify.verify().get("status") == "verified", "Stage-168 failed verification")
    seal = load(STAGE / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    require(
        seal.get("schema") == "koblitz_stage169_fresh_target_holdout_seal.v1"
        and seal.get("status") == "stage169_fresh_holdout_result_frozen"
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and phase_b.sha256_file(STAGE / "result.json", "Stage-169 result")
        == seal.get("result_sha256"),
        "Stage-169 seal changed",
    )
    result = load(STAGE / "result.json")
    prereg = load(STAGE / "preregistration.json")
    require(
        result.get("schema") == "koblitz_stage169_fresh_target_holdout.v1"
        and result.get("status") == "complete_preregistered_true_negative"
        and result["instance"]["blind_instance_id"] == BLIND_ID
        and result["instance"]["source_system_id"] == SOURCE_SYSTEM_ID
        and result["instance"]["target_class_opened_after_run_seal"]
        == "nondecomposable"
        and result["factor_base"]["factor_base_logs_known_by_construction"] is False
        and result["factor_base"]["target_subgroup_enumerated"] is False
        and prereg["status"] == "frozen_before_execution_and_truth_scoring"
        and prereg["selection"]["selected_blind_instance_id"] == BLIND_ID
        and prereg["selection"]["selected_bundle_ordinal_zero_based"] == 13
        and prereg["selection"]["selected_instance_truth_consulted_for_selection"]
        is False
        and prereg["execution"]["x1_batch_size"] == 39
        and prereg["execution"]["rayon_threads"] == 12,
        "Stage-169 identity or preregistration changed",
    )
    run_seal, plan = native_f4.validate_run_seal(STAGE / "holdout-run")
    require(
        run_seal["status"] == "native_f4_run_frozen_before_truth_scoring"
        and plan["source_revision"]
        == {"commit": COMMIT, "dirty": False, "porcelain": []}
        and plan["selected_blind_instance_ids"] == [BLIND_ID]
        and plan["x1_batch_size"] == 39
        and plan["rayon_threads_requested"] == 12
        and plan["single_thread_requested"] is False,
        "Stage-169 run plan changed",
    )
    tasks = sorted((STAGE / "holdout-run/tasks").glob("*/task-result.json"))
    require(len(tasks) == 1, "Stage-169 task count changed")
    selected = load(tasks[0])["result"]
    report = selected["backend_report"]
    metrics = result["fresh_holdout_native_f4"]["metrics"]
    require(
        selected["status"] == "unsat"
        and report["status"] == "unsat"
        and report["exhaustive"] is True
        and report["fixed_x1_values"] == 512
        and report["fixed_x1_masks_visited"] == 512
        and report["fixed_x1_nonrational_skipped"] == 270
        and report["fixed_x1_systems_constructed"] == 242
        and report["fixed_x1_systems_completed"] == 242
        and report["solver_x1_batches_completed"] == 7
        and report["solver_x1_speculative_systems_completed"] == 0
        and report["solver_equations_total"] == 14_278
        and report["solver_terms_total"] == 14_515_915
        and report["solver_equations_blake3"]
        == "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
        and report["algebraic_roots"] == 0
        and report["witness_points"] is None
        and report["cost"]["ops"] == 99_199_976_264
        and report["cost"]["extra"]["f4_calls"] == 242,
        "Stage-169 exhaustive F4 terminal changed",
    )
    require(
        metrics["wall_seconds"] == 20.53981879200728
        and metrics["total_core_seconds"] == 169.0941
        and metrics["single_core_seconds"] is None
        and metrics["peak_rss_bytes"] == 2_872_541_184,
        "Stage-169 F4 resources changed",
    )
    summary = load(STAGE / "holdout-run/run-summary.json")
    require(summary["f4_process_resources"]["single_core_seconds"] is None, "parallel holdout invented single-core time")
    receipt_path = STAGE / "holdout-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path)
    build_tool.validate_receipt(
        receipt_path, receipt["binaries"], receipt["rust_source_objects"]
    )
    require(
        receipt["source_commit"] == COMMIT
        and receipt["binaries"]["backend"]["sha256"] == BACKEND_SHA256,
        "Stage-169 build changed",
    )
    verify_score()
    table = {row["backend"]: row for row in result["same_target_table"]}
    require(
        table["direct-mitm"]["status"] == "unsat"
        and table["direct-mitm"]["wall_seconds"] == 2.9531290830054786
        and table["native-xor"]["status"] == "unknown_inconclusive"
        and table["native-xor"]["conflicts"] == 100_000
        and table["wdsat"]["wall_seconds"] == 120.00246820801112
        and table["wdsat"]["status"] == "timeout_inconclusive"
        and table["cryptominisat"]["wall_seconds"] == 120.00753962500312
        and table["cryptominisat"]["status"] == "timeout_inconclusive"
        and table["licensed-magma-f4"]["status"] == "not_run"
        and table["ggmp"]["status"] == "not_same_instance",
        "Stage-169 same-instance table changed",
    )
    controls = load(STAGE / "controls/manifest.json")
    require(
        controls["wdsat"]["decision"]
        == "rejected_long_absolute_path_fopen_failure"
        and controls["wdsat-valid"]["decision"]
        == "accepted_exact_input_timeout"
        and controls["cryptominisat"]["decision"]
        == "accepted_exact_input_timeout",
        "Stage-169 control decisions changed",
    )
    for name in ("direct-mitm", "native-xor", "wdsat", "wdsat-valid", "cryptominisat"):
        for artifact in controls[name]["artifacts"].values():
            path = STAGE / "controls" / artifact["path"]
            require(
                path.stat().st_size == artifact["bytes"]
                and phase_b.sha256_file(path, "Stage-169 control artifact")
                == artifact["sha256"],
                "Stage-169 control artifact changed",
            )
    comparisons = result["comparisons"]
    require(
        6.95 < comparisons["f4_over_direct_mitm"]["wall_ratio"] < 6.96
        and comparisons["f4_over_direct_mitm"]["core_ratio"] > 57
        and comparisons["selected_positive_vs_fresh_negative"]["wall_ratio"]
        > 4.11
        and comparisons["selected_positive_vs_fresh_negative"][
            "system_count_ratio"
        ]
        > 6.2
        and comparisons["cold_build_plus_holdout"]["wall_seconds"]
        == 187.6906990839922,
        "Stage-169 comparisons changed",
    )
    campaign = result["campaign_accounting"]
    require(
        campaign["reused_build_already_charged"] is True
        and campaign["measured_lower_bound"]["summed_wall_seconds"]
        == 8261.489664619901
        and campaign["measured_lower_bound"]["total_core_seconds"]
        == 8801.978114000003
        and campaign["measured_lower_bound"]["peak_rss_bytes"] == 6_310_576_128
        and campaign["measured_lower_bound"]["resource_components"] == 229
        and campaign["complete_campaign_cost"] is None,
        "Stage-169 campaign accounting changed",
    )
    require(
        result["fresh_target_holdout_complete"] is True
        and result["gates"]["2_same_instance_solver_matrix"].startswith(
            "partial_two_native_f4_n59_targets"
        )
        and result["gates"]["all_seven_gates_passed"] is False
        and result["licensed_magma_f4_complete"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["full_cost_gate_passed"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "Stage-169 claim boundary widened",
    )
    gate_status = (GATES / "GATE_STATUS.md").read_text()
    scoreboard = (REPO / "docs/index-calculus-scoreboard.html").read_text()
    require(
        "Stage 169 preregistered fresh-target holdout" in gate_status
        and "8,261.489665 sequential wall-seconds" in gate_status,
        "gate status lacks Stage 169",
    )
    require(
        "preregistered fresh n = 59 holdout" in scoreboard
        and "20.539819" in scoreboard,
        "scoreboard lacks Stage 169",
    )
    return {
        "schema": "koblitz_stage169_fresh_holdout_verification.v1",
        "status": "verified",
        "result_sha256": seal["result_sha256"],
        "inventory_sha256": seal["inventory_sha256"],
        "blind_instance_id": BLIND_ID,
        "source_commit": receipt["source_commit"],
        "classification": "true_negative",
        "f4_wall_seconds": metrics["wall_seconds"],
        "f4_total_core_seconds": metrics["total_core_seconds"],
        "f4_single_core_seconds": None,
        "f4_peak_rss_bytes": metrics["peak_rss_bytes"],
        "f4_systems_completed": report["fixed_x1_systems_completed"],
        "f4_word_xors": report["cost"]["ops"],
        "direct_mitm_wall_seconds": table["direct-mitm"]["wall_seconds"],
        "f4_over_direct_mitm_wall": comparisons["f4_over_direct_mitm"][
            "wall_ratio"
        ],
        "fresh_target_holdout_complete": True,
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
        raise SystemExit(f"stage169-fresh-holdout: {error}")

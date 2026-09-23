#!/usr/bin/env python3
"""Verify the committed Stage-162 grouped critical-pair result."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import build_koblitz_phase_b_native_f4 as build_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4
import verify_koblitz_stage161_fast_hash as stage161_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE = GATES / "stage-162-grouped-pair-selection-20260923"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
SOURCE_ID = "fb0ea6bb30be3dea841c609a2668939128b4b799d4c0e051b0b52dc6da05f0fa"


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
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256"),
        "Stage-162 score seal changed",
    )
    require(
        load(root / "score.json").get("classification_counts") == {"true_positive": 1},
        "Stage-162 score changed",
    )


def verify() -> dict[str, Any]:
    require(stage161_verify.verify().get("status") == "verified", "Stage-161 failed verification")
    seal = load(STAGE / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    require(
        seal.get("schema") == "koblitz_stage162_grouped_pair_selection_seal.v1"
        and seal.get("status") == "stage162_grouped_pair_result_frozen"
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and phase_b.sha256_file(STAGE / "result.json", "Stage-162 result")
        == seal.get("result_sha256"),
        "Stage-162 result seal changed",
    )
    result = load(STAGE / "result.json")
    require(
        result.get("schema") == "koblitz_stage162_grouped_pair_selection.v1"
        and result.get("status") == "complete_grouped_pair_selection_true_positive"
        and result["instance"]["blind_instance_id"] == BLIND_ID
        and result["instance"]["source_instance_id"] == SOURCE_ID
        and result["factor_base"]["factor_base_logs_known_by_construction"] is False
        and result["factor_base"]["target_subgroup_enumerated"] is False,
        "Stage-162 identity or algebraic-base boundary changed",
    )
    run_seal, plan = native_f4.validate_run_seal(STAGE / "selected-run")
    require(
        run_seal["status"] == "native_f4_run_frozen_before_truth_scoring"
        and plan["selected_blind_instance_ids"] == [BLIND_ID]
        and plan["source_revision"]
        == {"commit": "29524b34ef5bfcbb726a927c362e825c3203c390", "dirty": False, "porcelain": []},
        "Stage-162 blind run changed",
    )
    tasks = sorted((STAGE / "selected-run/tasks").glob("*/task-result.json"))
    require(len(tasks) == 1, "Stage-162 task count changed")
    task = load(tasks[0])
    selected = task["result"]
    report = selected["backend_report"]
    metrics = selected["metrics"]
    require(
        task["blind_instance_id"] == BLIND_ID
        and selected["status"] == "sat"
        and selected["source_model_valid"] is True
        and selected["source_witness_valid"] is True
        and report["solver_pair_selector"]
        == "grouped_lcm_exact_submask_equivalent_to_quadratic_update"
        and report["solver_equations_blake3"]
        == "d6969c5f89b5be6a7eacbdac45799fc734a273c9034949a482723b4422504ae7"
        and report["solver_terms_total"] == 3_864_601
        and report["cost"]["extra"]["f4_calls"] == 65
        and report["cost"]["extra"]["pairs_chain_skipped"] == 229_620_889
        and report["cost"]["extra"]["pairs_product_skipped"] == 70_978_814
        and report["cost"]["ops"] == 87_513_949_370,
        "Stage-162 exact F4 result changed",
    )
    require(
        metrics["wall_seconds"] == 57.86289354200562
        and metrics["total_core_seconds"] == 57.736423
        and metrics["peak_rss_bytes"] == 968_146_944,
        "Stage-162 selected resources changed",
    )
    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path)
    build_tool.validate_receipt(receipt_path, receipt["binaries"], receipt["rust_source_objects"])
    require(
        receipt["source_commit"] == "29524b34ef5bfcbb726a927c362e825c3203c390"
        and receipt["source_clean"] is True,
        "Stage-162 clean build changed",
    )
    verify_score()
    comparisons = result["comparisons"]
    require(
        comparisons["selected_over_stage161"]["wall_ratio"] < 0.955
        and comparisons["selected_over_stage161"]["same_f4_structure"] is True
        and comparisons["selected_over_same_host_direct_mitm"]["wall_ratio"] > 19
        and comparisons["selected_cold_build_plus_run"]["wall_seconds"]
        == 225.3537250400185,
        "Stage-162 comparisons changed",
    )
    campaign = result["campaign_accounting"]
    require(
        campaign["measured_lower_bound"]["summed_wall_seconds"] == 3029.348390705014
        and campaign["measured_lower_bound"]["total_core_seconds"] == 2981.63296
        and campaign["measured_lower_bound"]["peak_rss_bytes"] == 6_310_576_128
        and campaign["complete_campaign_cost"] is None,
        "Stage-162 campaign accounting changed",
    )
    development = load(STAGE / "development/manifest.json")
    require(
        development["profile"]["samples"] == 3327
        and development["profile"]["top_of_stack"]["state_insert"] == 695
        and development["candidate"]["decision"] == "promoted_after_clean_rebuild",
        "Stage-162 development record changed",
    )
    for section in (development["profile"], development["candidate"]):
        for artifact in section["artifacts"].values():
            path = STAGE / "development" / artifact["path"]
            require(
                path.stat().st_size == artifact["bytes"]
                and phase_b.sha256_file(path, "Stage-162 development artifact")
                == artifact["sha256"],
                "Stage-162 development artifact changed",
            )
    require(
        result["gates"]["all_seven_gates_passed"] is False
        and result["licensed_magma_f4_complete"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["full_cost_gate_passed"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "Stage-162 claim boundary widened",
    )
    gate_status = (GATES / "GATE_STATUS.md").read_text()
    scoreboard = (REPO / "docs/index-calculus-scoreboard.html").read_text()
    require(
        "Stage 162 grouped critical-pair selection" in gate_status
        and "3,029.348391 sequential wall-seconds" in gate_status,
        "gate status lacks Stage 162",
    )
    require(
        "still 19.29&times; direct MITM" in scoreboard and "57.862894" in scoreboard,
        "scoreboard lacks Stage 162",
    )
    return {
        "schema": "koblitz_stage162_pair_selection_verification.v1",
        "status": "verified",
        "result_sha256": seal["result_sha256"],
        "inventory_sha256": seal["inventory_sha256"],
        "blind_instance_id": BLIND_ID,
        "source_commit": receipt["source_commit"],
        "selected_wall_seconds": metrics["wall_seconds"],
        "selected_core_seconds": metrics["total_core_seconds"],
        "selected_peak_rss_bytes": metrics["peak_rss_bytes"],
        "whole_target_speedup": comparisons["selected_over_stage161"]["wall_speedup"],
        "same_host_f4_over_mitm_wall": comparisons["selected_over_same_host_direct_mitm"]["wall_ratio"],
        "all_seven_gates_passed": False,
        "koblitz_index_calculus_sota": False,
    }


if __name__ == "__main__":
    try:
        print(json.dumps(verify(), indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, IndexError, VerificationError, phase_b.PhaseBError) as error:
        raise SystemExit(f"stage162-pair-selection: {error}")

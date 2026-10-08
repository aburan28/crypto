#!/usr/bin/env python3
"""Verify the committed Stage-163 block-4 M4RI result."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import build_koblitz_phase_b_native_f4 as build_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4
import verify_koblitz_stage162_pair_selection as stage162_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE = GATES / "stage-163-m4ri-echelons-20260923"
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
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and load(root / "score.json").get("classification_counts") == {"true_positive": 1},
        "Stage-163 score changed",
    )


def verify() -> dict[str, Any]:
    require(stage162_verify.verify().get("status") == "verified", "Stage-162 failed verification")
    seal = load(STAGE / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    require(
        seal.get("schema") == "koblitz_stage163_m4ri_echelons_seal.v1"
        and seal.get("status") == "stage163_m4ri_result_frozen"
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and phase_b.sha256_file(STAGE / "result.json", "Stage-163 result")
        == seal.get("result_sha256"),
        "Stage-163 seal changed",
    )
    result = load(STAGE / "result.json")
    require(
        result.get("schema") == "koblitz_stage163_m4ri_echelons.v1"
        and result.get("status") == "complete_m4ri_true_positive"
        and result["instance"]["blind_instance_id"] == BLIND_ID
        and result["instance"]["source_instance_id"] == SOURCE_ID
        and result["factor_base"]["factor_base_logs_known_by_construction"] is False
        and result["factor_base"]["target_subgroup_enumerated"] is False,
        "Stage-163 identity or algebraic-base boundary changed",
    )
    run_seal, plan = native_f4.validate_run_seal(STAGE / "selected-run")
    require(
        run_seal["status"] == "native_f4_run_frozen_before_truth_scoring"
        and plan["selected_blind_instance_ids"] == [BLIND_ID]
        and plan["source_revision"]
        == {"commit": "04e49c930caccea31bad7172aa8f567d96744b17", "dirty": False, "porcelain": []},
        "Stage-163 blind run changed",
    )
    tasks = sorted((STAGE / "selected-run/tasks").glob("*/task-result.json"))
    require(len(tasks) == 1, "Stage-163 task count changed")
    selected = load(tasks[0])["result"]
    report = selected["backend_report"]
    metrics = selected["metrics"]
    require(
        selected["status"] == "sat"
        and selected["source_model_valid"] is True
        and selected["source_witness_valid"] is True
        and report["solver_echelon_policy"].startswith("block4_m4ri")
        and report["solver_equations_blake3"]
        == "d6969c5f89b5be6a7eacbdac45799fc734a273c9034949a482723b4422504ae7"
        and report["solver_terms_total"] == 3_864_601
        and report["cost"]["op_unit"]
        == "word XORs (elimination, including M4RI tables)"
        and report["cost"]["extra"]["m4ri_matrices"] == 196
        and report["cost"]["extra"]["m4ri_table_word_xors"] == 316_714_298
        and report["cost"]["ops"] == 45_879_309_338,
        "Stage-163 exact M4RI result changed",
    )
    require(
        metrics["wall_seconds"] == 54.77377895799873
        and metrics["total_core_seconds"] == 54.662604
        and metrics["peak_rss_bytes"] == 915_668_992,
        "Stage-163 resources changed",
    )
    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path)
    build_tool.validate_receipt(receipt_path, receipt["binaries"], receipt["rust_source_objects"])
    require(
        receipt["source_commit"] == "04e49c930caccea31bad7172aa8f567d96744b17"
        and receipt["source_clean"] is True,
        "Stage-163 clean build changed",
    )
    verify_score()
    comparisons = result["comparisons"]
    require(
        comparisons["same_binary_m4ri_over_streaming"]["wall_ratio"] < 0.935
        and comparisons["same_binary_m4ri_over_streaming"]["word_xor_ratio"] < 0.525
        and comparisons["selected_over_stage162"]["wall_ratio"] < 0.947
        and comparisons["selected_over_same_host_direct_mitm"]["wall_ratio"] > 18
        and comparisons["selected_cold_build_plus_run"]["wall_seconds"]
        == 223.24849729199894,
        "Stage-163 comparisons changed",
    )
    campaign = result["campaign_accounting"]
    require(
        campaign["measured_lower_bound"]["summed_wall_seconds"] == 3589.954839372015
        and campaign["measured_lower_bound"]["total_core_seconds"] == 3532.297699
        and campaign["measured_lower_bound"]["peak_rss_bytes"] == 6_310_576_128
        and campaign["complete_campaign_cost"] is None,
        "Stage-163 campaign accounting changed",
    )
    development = load(STAGE / "development/manifest.json")
    require(
        development["arms"]["streaming-control"]["decision"] == "same_binary_control"
        and development["arms"]["m4ri-provisional"]["decision"]
        == "promoted_after_clean_rebuild"
        and development["parser_rejected_clean_run"]["scientific_terminal_promoted"] is False,
        "Stage-163 development boundary changed",
    )
    for arm in development["arms"].values():
        for artifact in arm["artifacts"].values():
            path = STAGE / "development" / artifact["path"]
            require(
                path.stat().st_size == artifact["bytes"]
                and phase_b.sha256_file(path, "Stage-163 arm artifact") == artifact["sha256"],
                "Stage-163 arm artifact changed",
            )
    failed_root = STAGE / "development" / development["parser_rejected_clean_run"]["path"]
    failed_inventory = phase_b.all_regular_inventory(failed_root)
    require(
        phase_b.canonical_sha256(failed_inventory)
        == development["parser_rejected_clean_run"]["inventory_sha256"],
        "Stage-163 failed-run inventory changed",
    )
    require(
        result["gates"]["all_seven_gates_passed"] is False
        and result["licensed_magma_f4_complete"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["full_cost_gate_passed"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "Stage-163 claim boundary widened",
    )
    gate_status = (GATES / "GATE_STATUS.md").read_text()
    scoreboard = (REPO / "docs/index-calculus-scoreboard.html").read_text()
    require(
        "Stage 163 shape-selected block-4 M4RI" in gate_status
        and "3,589.954839 sequential wall-seconds" in gate_status,
        "gate status lacks Stage 163",
    )
    require(
        "still 18.26&times; direct MITM" in scoreboard and "54.773779" in scoreboard,
        "scoreboard lacks Stage 163",
    )
    return {
        "schema": "koblitz_stage163_m4ri_verification.v1",
        "status": "verified",
        "result_sha256": seal["result_sha256"],
        "inventory_sha256": seal["inventory_sha256"],
        "blind_instance_id": BLIND_ID,
        "source_commit": receipt["source_commit"],
        "selected_wall_seconds": metrics["wall_seconds"],
        "selected_core_seconds": metrics["total_core_seconds"],
        "selected_peak_rss_bytes": metrics["peak_rss_bytes"],
        "selected_word_xors": report["cost"]["ops"],
        "same_binary_wall_speedup": comparisons["same_binary_m4ri_over_streaming"]["wall_speedup"],
        "same_host_f4_over_mitm_wall": comparisons["selected_over_same_host_direct_mitm"]["wall_ratio"],
        "all_seven_gates_passed": False,
        "koblitz_index_calculus_sota": False,
    }


if __name__ == "__main__":
    try:
        print(json.dumps(verify(), indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, IndexError, VerificationError, phase_b.PhaseBError) as error:
        raise SystemExit(f"stage163-m4ri: {error}")

#!/usr/bin/env python3
"""Verify the committed Stage-160 fixed-X1 construction result."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import build_koblitz_phase_b_native_f4 as build_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4
import verify_koblitz_stage159_native_f4 as stage159_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE = GATES / "stage-160-fixed-x1-constant-specialisation-20260922"
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
    require(
        seal.get("schema") == native_f4.SCORE_SEAL_SCHEMA
        and claimed == phase_b.canonical_sha256(payload),
        "Stage-160 score seal changed",
    )
    inventory = phase_b.all_regular_inventory(root, {"score-seal.json"})
    require(
        inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256"),
        "Stage-160 score inventory changed",
    )
    score = load(root / "score.json")
    require(
        score.get("classification_counts") == {"true_positive": 1}
        and score.get("rows", [{}])[0].get("blind_instance_id") == BLIND_ID,
        "Stage-160 score changed",
    )


def verify() -> dict[str, Any]:
    stage159 = stage159_verify.verify()
    require(stage159.get("status") == "verified", "Stage-159 baseline failed verification")

    seal = load(STAGE / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == "koblitz_stage160_fixed_x1_constant_specialisation_seal.v1"
        and seal.get("status") == "stage160_fixed_x1_result_frozen"
        and claimed == phase_b.canonical_sha256(payload),
        "Stage-160 result seal changed",
    )
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    require(
        inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256"),
        "Stage-160 inventory changed",
    )
    require(
        phase_b.sha256_file(STAGE / "result.json", "Stage-160 result")
        == seal.get("result_sha256"),
        "Stage-160 result hash changed",
    )

    result = load(STAGE / "result.json")
    require(
        result.get("schema") == "koblitz_stage160_fixed_x1_constant_specialisation.v1"
        and result.get("status") == "complete_fixed_x1_constant_specialisation_true_positive",
        "Stage-160 status changed",
    )
    require(
        result["instance"]["blind_instance_id"] == BLIND_ID
        and result["instance"]["source_instance_id"] == SOURCE_ID
        and result["factor_base"]["factor_base_logs_known_by_construction"] is False
        and result["factor_base"]["target_subgroup_enumerated"] is False
        and result["factor_base"]["discrete_log_labels_used"] is False,
        "Stage-160 public-algebra boundary changed",
    )

    run_seal, plan = native_f4.validate_run_seal(STAGE / "selected-run")
    require(
        run_seal["status"] == "native_f4_run_frozen_before_truth_scoring"
        and plan["selected_blind_instance_ids"] == [BLIND_ID]
        and plan["source_revision"]
        == {"commit": "3cb24b82e8a059a293f6371178bd0d3647272a63", "dirty": False, "porcelain": []},
        "Stage-160 blind run changed",
    )
    task_paths = sorted((STAGE / "selected-run/tasks").glob("*/task-result.json"))
    require(len(task_paths) == 1, "Stage-160 task count changed")
    task = load(task_paths[0])
    selected = task["result"]
    report = selected["backend_report"]
    metrics = selected["metrics"]
    require(
        task["blind_instance_id"] == BLIND_ID
        and task["source_instance_id"] == SOURCE_ID
        and selected["status"] == "sat"
        and selected["source_model_valid"] is True
        and selected["source_witness_valid"] is True
        and report["solver_constructor"] == "fixed_x1_constant_linear_specialisation"
        and report["solver_equations_blake3"]
        == "d6969c5f89b5be6a7eacbdac45799fc734a273c9034949a482723b4422504ae7"
        and report["solver_terms_total"] == 3_864_601
        and report["fixed_x1_systems_completed"] == 65
        and report["cost"]["extra"]["f4_calls"] == 65
        and report["cost"]["ops"] == 87_513_949_370
        and report["timing_ns"]["fixed_x1_s4_construction"] == 973_558_384,
        "Stage-160 exact F4 result changed",
    )
    require(
        metrics["wall_seconds"] == 73.75799025000015
        and metrics["total_core_seconds"] == 73.585054
        and metrics["peak_rss_bytes"] == 910_934_016,
        "Stage-160 selected resources changed",
    )
    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path)
    build_tool.validate_receipt(receipt_path, receipt["binaries"], receipt["rust_source_objects"])
    require(
        receipt["source_commit"] == "3cb24b82e8a059a293f6371178bd0d3647272a63"
        and receipt["source_clean"] is True,
        "Stage-160 clean build changed",
    )
    verify_score()

    comparisons = result["comparisons"]
    require(
        comparisons["selected_over_stage159"]["wall_ratio"] < 0.737
        and comparisons["selected_over_stage159"]["construction_ratio"] < 0.037
        and comparisons["selected_over_stage159"]["same_equation_fingerprint"] is True
        and comparisons["selected_over_stage159"]["same_word_xors"] is True
        and comparisons["selected_over_same_host_direct_mitm"]["wall_ratio"] > 24
        and comparisons["selected_cold_build_plus_run"]["wall_seconds"]
        == 241.26225399900434,
        "Stage-160 comparisons changed",
    )
    campaign = result["campaign_accounting"]
    require(
        campaign["measured_lower_bound"]["summed_wall_seconds"] == 2319.762065749015
        and campaign["measured_lower_bound"]["total_core_seconds"] == 2281.647854
        and campaign["measured_lower_bound"]["peak_rss_bytes"] == 6_310_576_128
        and campaign["complete_campaign_cost"] is None
        and campaign["unmetered_development_commands_present"] is True,
        "Stage-160 campaign accounting changed",
    )
    development = load(STAGE / "development/manifest.json")
    require(
        development["status"] == "metered_development_arms_retained"
        and [arm["decision"] for arm in development["arms"]]
        == ["rejected", "rejected", "rejected", "rejected", "promoted_after_clean_rebuild"],
        "Stage-160 development decisions changed",
    )
    for arm in development["arms"]:
        for artifact in arm["artifacts"].values():
            path = STAGE / "development" / artifact["path"]
            require(
                path.stat().st_size == artifact["bytes"]
                and phase_b.sha256_file(path, "Stage-160 development artifact")
                == artifact["sha256"],
                "Stage-160 development artifact changed",
            )

    require(
        result["gates"]["all_seven_gates_passed"] is False
        and result["licensed_magma_f4_complete"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["full_cost_gate_passed"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "Stage-160 claim boundary widened",
    )
    gate_status = (GATES / "GATE_STATUS.md").read_text()
    scoreboard = (REPO / "docs/index-calculus-scoreboard.html").read_text()
    require(
        "Stage 160 fixed-X1 constant-linear construction" in gate_status
        and "2,319.762066 sequential wall-seconds" in gate_status,
        "gate status lacks Stage 160",
    )
    require(
        'id="phase-b-native-f4-single-target"' in scoreboard
        and "73.757990" in scoreboard,
        "canonical scoreboard no longer retains Stage 160",
    )
    return {
        "schema": "koblitz_stage160_fixed_x1_verification.v1",
        "status": "verified",
        "result_sha256": seal["result_sha256"],
        "inventory_sha256": seal["inventory_sha256"],
        "blind_instance_id": BLIND_ID,
        "source_commit": receipt["source_commit"],
        "selected_wall_seconds": metrics["wall_seconds"],
        "selected_core_seconds": metrics["total_core_seconds"],
        "selected_peak_rss_bytes": metrics["peak_rss_bytes"],
        "construction_speedup": comparisons["selected_over_stage159"]["construction_speedup"],
        "whole_target_speedup": comparisons["selected_over_stage159"]["wall_speedup"],
        "same_host_f4_over_mitm_wall": comparisons["selected_over_same_host_direct_mitm"][
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
        raise SystemExit(f"stage160-fixed-x1: {error}")

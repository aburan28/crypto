#!/usr/bin/env python3
"""Verify the committed Stage-159 native-F4 single-target result."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import build_koblitz_phase_b_native_f4 as build_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE = GATES / "stage-159-native-f4-single-target-20260922"
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


def verify_score(name: str, classification: str) -> None:
    root = STAGE / name
    seal = load(root / "score-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == native_f4.SCORE_SEAL_SCHEMA
        and claimed == phase_b.canonical_sha256(payload),
        f"{name} score seal changed",
    )
    inventory = phase_b.all_regular_inventory(root, {"score-seal.json"})
    require(
        inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256"),
        f"{name} inventory changed",
    )
    score = load(root / "score.json")
    require(
        score.get("classification_counts") == {classification: 1}
        and score.get("rows", [{}])[0].get("blind_instance_id") == BLIND_ID,
        f"{name} classification changed",
    )


def verify() -> dict[str, Any]:
    seal = load(STAGE / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == "koblitz_stage159_native_f4_single_target_seal.v1"
        and seal.get("status") == "stage159_native_f4_result_frozen"
        and claimed == phase_b.canonical_sha256(payload),
        "Stage-159 result seal changed",
    )
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    require(
        inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256"),
        "Stage-159 inventory changed",
    )
    require(
        phase_b.sha256_file(STAGE / "result.json", "Stage-159 result")
        == seal.get("result_sha256"),
        "Stage-159 result hash changed",
    )

    result = load(STAGE / "result.json")
    require(
        result.get("schema") == "koblitz_stage159_native_f4_single_target.v1"
        and result.get("status") == "complete_single_target_native_f4_true_positive",
        "Stage-159 result status changed",
    )
    require(
        result["instance"]["blind_instance_id"] == BLIND_ID
        and result["instance"]["source_instance_id"] == SOURCE_ID
        and result["instance"]["n"] == 59
        and result["instance"]["ell"] == 9
        and result["instance"]["m"] == 3,
        "Stage-159 instance identity changed",
    )
    factor_base = result["factor_base"]
    require(
        factor_base["factor_base_logs_known_by_construction"] is False
        and factor_base["target_subgroup_enumerated"] is False
        and factor_base["discrete_log_labels_used"] is False,
        "Stage-159 factor-base knowledge boundary changed",
    )
    selected = result["selected_native_f4"]
    report = selected["report"]
    metrics = selected["metrics"]
    require(
        selected["status"] == "sat"
        and selected["classification"] == "true_positive"
        and selected["source_model_valid"] is True
        and selected["source_witness_valid"] is True
        and report["fixed_x1_nonrational_skipped"] == 73
        and report["fixed_x1_systems_completed"] == 65
        and report["cost"]["extra"]["f4_calls"] == 65
        and report["cost"]["ops"] == 87_513_949_370,
        "Stage-159 native-F4 terminal changed",
    )
    require(
        metrics["wall_seconds"] == 100.17938004100142
        and metrics["total_core_seconds"] == 99.920352
        and metrics["peak_rss_bytes"] == 921_681_920,
        "Stage-159 selected resources changed",
    )
    selected_build_path = STAGE / "attempt-05-rational-x1/tool-builds/rust/receipt.json"
    selected_build = load(selected_build_path)
    build_tool.validate_receipt(
        selected_build_path,
        selected_build["binaries"],
        selected_build["rust_source_objects"],
    )
    rows = {row["backend"]: row for row in result["same_target_table"]}
    require(
        rows["direct-mitm"]["status"] == "sat"
        and rows["direct-mitm"]["wall_seconds"] == 2.999660749999748
        and rows["native-xor"]["status"] == "unknown_inconclusive"
        and rows["wdsat"]["status"] == "timeout_inconclusive"
        and rows["cryptominisat"]["status"] == "timeout_inconclusive"
        and rows["magma-f4"]["wall_seconds"] is None
        and rows["ggmp"]["status"] == "not_same_instance",
        "Stage-159 same-target table changed",
    )
    comparisons = result["comparisons"]
    require(
        comparisons["selected_over_full_fixed_x1"]["wall_ratio"] < 0.488
        and comparisons["selected_over_full_fixed_x1"]["core_ratio"] < 0.488
        and comparisons["selected_over_full_fixed_x1"]["word_xor_ratio"] < 0.48
        and comparisons["selected_over_full_fixed_x1"]["same_verified_witness_x"] is True,
        "Stage-159 selected improvement changed",
    )
    require(
        comparisons["selected_over_same_host_direct_mitm"]["wall_ratio"] > 33
        and comparisons["selected_over_same_host_direct_mitm"]["core_ratio"] > 33
        and comparisons["rejected_trace_over_selected"]["wall_ratio"] > 1.18
        and comparisons["rejected_trace_over_selected"]["word_xor_ratio"] > 1.58,
        "Stage-159 reference or rejected-control boundary changed",
    )
    campaign = result["campaign_accounting"]
    require(
        campaign["measured_lower_bound"]["summed_wall_seconds"] == 1596.064122959007
        and campaign["measured_lower_bound"]["total_core_seconds"] == 1565.839519
        and campaign["measured_lower_bound"]["peak_rss_bytes"] == 6_310_576_128
        and campaign["complete_campaign_cost"] is None
        and campaign["unmetered_development_commands_present"] is True,
        "Stage-159 campaign accounting changed",
    )
    require(
        result["gates"]["all_seven_gates_passed"] is False
        and result["licensed_magma_f4_complete"] is False
        and result["full_cost_gate_passed"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "Stage-159 claim boundary widened",
    )

    controls = STAGE / "controls"
    provenance = load(controls / "provenance.json")
    require(
        provenance["blind_instance_id"] == BLIND_ID
        and provenance["source_instance_id"] == SOURCE_ID
        and phase_b.sha256_file(controls / "stage26-task-result.json", "Stage-26 compact task")
        == provenance["stage26"]["task_sha256"]
        and phase_b.sha256_file(controls / "stage27-task-result.json", "Stage-27 compact task")
        == provenance["stage27"]["task_sha256"]
        and provenance["stage26"]["verification"]["status"] == "verified"
        and provenance["stage27"]["verification"]["status"] == "verified",
        "Stage-159 imported control provenance changed",
    )
    control26 = load(controls / "stage26-task-result.json")
    control27 = load(controls / "stage27-task-result.json")
    require(
        control26["source_instance_id"] == control27["source_instance_id"] == SOURCE_ID
        and control26["source_artifacts_before"] == control26["source_artifacts_after"]
        and control27["source_before"] == control27["source_after"],
        "Stage-159 compact control source custody changed",
    )

    for name in (
        "attempt-03-fixed-x1",
        "attempt-04-fixed-x1-long",
        "attempt-05-rational-x1",
        "attempt-06-trace",
    ):
        run_seal, plan = native_f4.validate_run_seal(STAGE / name)
        require(
            run_seal["status"] == "native_f4_run_frozen_before_truth_scoring"
            and plan["selected_blind_instance_ids"] == [BLIND_ID],
            f"{name} changed",
        )
    verify_score("attempt-03-score", "inconclusive")
    verify_score("attempt-04-score", "true_positive")
    verify_score("attempt-05-score", "true_positive")
    verify_score("attempt-06-score", "true_positive")

    gate_status = (GATES / "GATE_STATUS.md").read_text()
    scoreboard = (REPO / "docs/index-calculus-scoreboard.html").read_text()
    require("Stage 159 native-F4 single-target supplement" in gate_status, "gate status not updated")
    require('id="phase-b-native-f4-single-target"' in scoreboard, "scoreboard not updated")
    return {
        "schema": "koblitz_stage159_native_f4_verification.v1",
        "status": "verified",
        "result_sha256": seal["result_sha256"],
        "inventory_sha256": seal["inventory_sha256"],
        "blind_instance_id": BLIND_ID,
        "selected_wall_seconds": metrics["wall_seconds"],
        "selected_core_seconds": metrics["total_core_seconds"],
        "selected_peak_rss_bytes": metrics["peak_rss_bytes"],
        "selected_word_xors": report["cost"]["ops"],
        "same_host_f4_over_mitm_wall": comparisons["selected_over_same_host_direct_mitm"][
            "wall_ratio"
        ],
        "all_seven_gates_passed": False,
        "koblitz_index_calculus_sota": False,
    }


if __name__ == "__main__":
    try:
        print(json.dumps(verify(), indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, IndexError, VerificationError, phase_b.PhaseBError) as error:
        raise SystemExit(f"stage159-native-f4: {error}")

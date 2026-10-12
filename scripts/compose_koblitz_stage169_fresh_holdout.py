#!/usr/bin/env python3
"""Compose the sealed Stage-169 preregistered fresh-target F4 holdout."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
from typing import Any

import build_koblitz_phase_b_native_f4 as build_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4
import verify_koblitz_stage168_single_batch as stage168_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE168 = GATES / "stage-168-single-batch-39-20260923"
STAGE = GATES / "stage-169-fresh-target-holdout-20260923"
PREREG = GATES / "stage-169-fresh-target-holdout-preregistration-20260923.json"
STAGE20_SCORE = (
    GATES
    / "stage-20-phase-b-terminal-evidence-successor-04-20260910/score/score.json"
)
STAGE27_SCORE = GATES / "stage-27-direct-mitm-result-20260911/score.json"
BLIND_ID = "b-421e22a9c1c3b9d56396c8bbd0e46185bbebc32de0306f2bee6e9585703c2be4"
SOURCE_SYSTEM_ID = "x-af85283c5728becd64890c6eaa654ccb164cb499fcc4cde7b5cfd1484b4e0efa"
SOURCE_INSTANCE_ID = "954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7"
COMMIT = "509e43c79bcf36425f4329a09c34fd20c75949ed"
BACKEND_SHA256 = "58767973f4bc1fcf30d5bbca231d1d5647e6a79f6282c55f1c3874a962820931"
WDSAT_SHA256 = "1de1e3328af3f00963cd5339c07cc4ab18efd55d1749c953ecaafcefc54fd185"
CMS_SHA256 = "6c509f09622f103d8a3ad90afc151e1c4031275c052d7f8af481465b4f27f2af"
FINGERPRINT = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"


class Stage169Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage169Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be an object")
    return value


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def only_task(root: Path) -> dict[str, Any]:
    paths = sorted((root / "tasks").glob("*/task-result.json"))
    require(len(paths) == 1, "holdout run must contain one task")
    task = load(paths[0], "holdout task")
    require(
        task.get("blind_instance_id") == BLIND_ID
        and task.get("source_system_id") == SOURCE_SYSTEM_ID
        and task.get("source_instance_id") == SOURCE_INSTANCE_ID,
        "holdout task identity changed",
    )
    return task


def verify_score(root: Path) -> dict[str, Any]:
    seal = load(root / "score-seal.json", "holdout score seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(root, {"score-seal.json"})
    score = load(root / "score.json", "holdout score")
    require(
        seal.get("schema") == native_f4.SCORE_SEAL_SCHEMA
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and score.get("classification_counts") == {"true_negative": 1}
        and score["rows"][0]["blind_instance_id"] == BLIND_ID
        and score["rows"][0]["target_class"] == "nondecomposable"
        and score["rows"][0]["solver_status"] == "unsat",
        "holdout score changed",
    )
    return score


def copy_control(root: Path, output: Path, name: str) -> tuple[dict[str, Any], dict[str, Any]]:
    artifacts = {}
    for suffix in ("metrics", "stdout", "stderr"):
        source = root / f"{name}.{suffix}"
        require(source.is_file(), f"missing control artifact: {source}")
        destination = output / f"{name}.{suffix}"
        shutil.copyfile(source, destination)
        artifacts[suffix] = {
            "path": destination.name,
            "bytes": destination.stat().st_size,
            "sha256": sha256(destination),
        }
    return artifacts, load(output / f"{name}.metrics", f"{name} metrics")


def sum_resources(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "resource_components": len(rows),
        "summed_wall_seconds": round(math.fsum(row["wall_seconds"] for row in rows), 12),
        "total_core_seconds": round(
            math.fsum(row["total_core_seconds"] for row in rows), 12
        ),
        "peak_rss_bytes": max(row["peak_rss_bytes"] for row in rows),
    }


def historical_rows() -> dict[str, Any]:
    stage20 = load(STAGE20_SCORE, "Stage-20 score")
    stage27 = load(STAGE27_SCORE, "Stage-27 direct score")
    solver_rows = [
        row for row in stage20["rows"] if row.get("blind_instance_id") == BLIND_ID
    ]
    direct_rows = [
        row for row in stage27["rows"] if row.get("blind_instance_id") == BLIND_ID
    ]
    require(
        {row["backend"] for row in solver_rows}
        == {"native-xor", "wdsat", "cryptominisat"}
        and len(direct_rows) == 1
        and all(row["target_class"] == "nondecomposable" for row in solver_rows)
        and direct_rows[0]["target_class"] == "nondecomposable"
        and direct_rows[0]["classification"] == "true_negative",
        "historical same-instance rows changed",
    )
    return {
        "stage20_score": {
            "path": str(STAGE20_SCORE.relative_to(REPO)),
            "sha256": sha256(STAGE20_SCORE),
            "rows": solver_rows,
        },
        "stage27_score": {
            "path": str(STAGE27_SCORE.relative_to(REPO)),
            "sha256": sha256(STAGE27_SCORE),
            "row": direct_rows[0],
        },
    }


def write_markdown(result: dict[str, Any]) -> None:
    f4 = result["fresh_holdout_native_f4"]
    table = result["same_target_table"]
    direct = next(row for row in table if row["backend"] == "direct-mitm")
    lines = [
        "# Stage 169: preregistered fresh n=59 target holdout",
        "",
        result["claim_boundary"],
        "",
        "The holdout was selected from the blind bundle by presentation order and committed before execution. Batch size 39, twelve Rayon threads, a 300-second F4 budget, and no retry after truth were frozen in the preregistration.",
        "",
        f"Native F4 exhaustively visits all {f4['report']['fixed_x1_masks_visited']} masks, skips {f4['report']['fixed_x1_nonrational_skipped']} non-rational values, and completes {f4['report']['fixed_x1_systems_completed']} F4 systems. It returns UNSAT in {f4['metrics']['wall_seconds']:.6f} wall seconds / {f4['metrics']['total_core_seconds']:.6f} core-seconds / {f4['metrics']['peak_rss_bytes']} bytes peak RSS. Post-run scoring classifies it as a true negative.",
        "",
        f"Same-host direct MITM also returns UNSAT in {direct['wall_seconds']:.6f} wall seconds. Native XOR is inconclusive at 100,000 conflicts, while WDSat and CryptoMiniSat reach their 120-second watchdogs. Licensed Magma remains unexecuted.",
        "",
        "This holdout proves that the fixed policy can complete a fresh nondecomposable target and produce a correct exhaustive UNSAT terminal. Its 20.54-second wall cost is materially above the post-hoc 4.99-second positive target and direct MITM, so it rejects any inference that Stage 168 represents generic fresh-target time.",
        "",
        "One fresh holdout is finite evidence, not an expected-time distribution or full panel. Licensed Magma, the full native-F4 panel, end-to-end IC/rho crossover, and unaffiliated reproduction and novelty review remain open. This is not a Koblitz index-calculus SOTA result.",
        "",
    ]
    phase_b.write_new(STAGE / "RESULTS.md", ("\n".join(lines)).encode())


def compose(run_source: Path, score_source: Path, controls_root: Path) -> dict[str, Any]:
    require(not STAGE.exists(), f"refusing to overwrite {STAGE}")
    require(stage168_verify.verify().get("status") == "verified", "Stage-168 failed verification")
    baseline = load(STAGE168 / "result.json", "Stage-168 result")
    prereg = load(PREREG, "Stage-169 preregistration")
    require(
        prereg["status"] == "frozen_before_execution_and_truth_scoring"
        and prereg["preregistration_parent_commit"]
        == "02338cb7c55b9657a5612a407247bf3bab2e71b3"
        and prereg["selection"]["selected_blind_instance_id"] == BLIND_ID
        and prereg["selection"]["selected_bundle_ordinal_zero_based"] == 13
        and prereg["selection"]["selected_instance_truth_consulted_for_selection"]
        is False
        and prereg["execution"]["x1_batch_size"] == 39
        and prereg["execution"]["rayon_threads"] == 12,
        "preregistration changed",
    )
    STAGE.mkdir()
    shutil.copyfile(PREREG, STAGE / "preregistration.json")
    shutil.copytree(run_source, STAGE / "holdout-run")
    shutil.copytree(score_source, STAGE / "holdout-score")
    controls = STAGE / "controls"
    controls.mkdir()
    control_docs: dict[str, Any] = {}
    resource_rows: list[dict[str, Any]] = []
    for name in (
        "direct-mitm",
        "native-xor",
        "wdsat",
        "wdsat-valid",
        "cryptominisat",
    ):
        artifacts, process = copy_control(controls_root, controls, name)
        resource_rows.append(process["metrics"])
        control_docs[name] = {
            "artifacts": artifacts,
            "process": process,
        }
    direct_report = load(controls / "direct-mitm.stdout", "direct report")
    native_report = load(controls / "native-xor.stdout", "native report")
    require(
        control_docs["direct-mitm"]["process"]["returncode"] == 0
        and direct_report["status"] == "unsat"
        and direct_report["exhaustive"] is True
        and direct_report["source_instance_id"] == SOURCE_INSTANCE_ID
        and control_docs["native-xor"]["process"]["returncode"] == 0
        and native_report["status"] == "unknown_inconclusive"
        and native_report["stats"]["conflicts"] == 100_000
        and native_report["source_instance_id"] == SOURCE_INSTANCE_ID,
        "local direct or native-XOR control changed",
    )
    require(
        control_docs["wdsat"]["process"]["returncode"] == 1
        and control_docs["wdsat"]["process"]["timed_out"] is False
        and "unkown filename" in (controls / "wdsat.stdout").read_text()
        and control_docs["wdsat-valid"]["process"]["timed_out"] is True
        and control_docs["wdsat-valid"]["process"]["returncode"] == -15
        and control_docs["cryptominisat"]["process"]["timed_out"] is True
        and control_docs["cryptominisat"]["process"]["returncode"] == -15,
        "external SAT control terminal changed",
    )
    control_docs["wdsat"]["decision"] = "rejected_long_absolute_path_fopen_failure"
    control_docs["wdsat-valid"]["decision"] = "accepted_exact_input_timeout"
    control_docs["cryptominisat"]["decision"] = "accepted_exact_input_timeout"
    control_docs["native-xor"]["decision"] = "accepted_conflict_cap_inconclusive"
    control_docs["direct-mitm"]["decision"] = "accepted_exhaustive_true_negative"
    control_docs["binary_identities"] = {
        "backend_sha256": BACKEND_SHA256,
        "wdsat_sha256": WDSAT_SHA256,
        "cryptominisat_sha256": CMS_SHA256,
    }
    control_docs["input_bindings"] = {
        name: {
            "path": artifact["path"],
            "bytes": artifact["bytes"],
            "blake3": artifact["manifest_blake3"],
            "sha256": artifact["sha256"],
        }
        for name, artifact in only_task(STAGE / "holdout-run")["source_artifacts_before"][
            "exports"
        ].items()
    }
    phase_b.write_json_new(controls / "manifest.json", control_docs)
    history = historical_rows()
    phase_b.write_json_new(STAGE / "historical-same-instance.json", history)

    run_seal, plan = native_f4.validate_run_seal(STAGE / "holdout-run")
    task = only_task(STAGE / "holdout-run")
    selected = task["result"]
    report = selected["backend_report"]
    score = verify_score(STAGE / "holdout-score")
    raw_metrics = selected["metrics"]
    metrics = copy.deepcopy(raw_metrics)
    metrics["single_core_seconds"] = None
    require(
        plan["source_revision"] == {"commit": COMMIT, "dirty": False, "porcelain": []}
        and plan["selected_blind_instance_ids"] == [BLIND_ID]
        and plan["x1_batch_size"] == 39
        and plan["rayon_threads_requested"] == 12
        and selected["status"] == "unsat"
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
        and report["solver_equations_blake3"] == FINGERPRINT
        and report["algebraic_roots"] == 0
        and report["witness_points"] is None
        and report["cost"]["ops"] == 99_199_976_264
        and report["cost"]["extra"]["f4_calls"] == 242
        and score["classification_counts"] == {"true_negative": 1},
        "holdout result changed",
    )
    summary = load(STAGE / "holdout-run/run-summary.json", "holdout summary")
    require(
        summary["f4_process_resources"]["single_core_seconds"] is None,
        "parallel holdout invented single-core time",
    )
    receipt_path = STAGE / "holdout-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path, "holdout build receipt")
    build_tool.validate_receipt(
        receipt_path, receipt["binaries"], receipt["rust_source_objects"]
    )
    require(
        receipt["source_commit"] == COMMIT
        and receipt["binaries"]["backend"]["sha256"] == BACKEND_SHA256,
        "holdout build changed",
    )
    same_target_table = [
        {
            "backend": "native-f4-batch39-t12",
            "status": "unsat",
            "classification": "true_negative",
            "wall_seconds": metrics["wall_seconds"],
            "total_core_seconds": metrics["total_core_seconds"],
            "single_core_seconds": None,
            "peak_rss_bytes": metrics["peak_rss_bytes"],
            "conflicts": None,
            "operations": {"word_xors": report["cost"]["ops"], "f4_calls": 242},
        },
        {
            "backend": "direct-mitm",
            "status": "unsat",
            "classification": "true_negative",
            "wall_seconds": control_docs["direct-mitm"]["process"]["metrics"]["wall_seconds"],
            "total_core_seconds": control_docs["direct-mitm"]["process"]["metrics"]["total_core_seconds"],
            "single_core_seconds": control_docs["direct-mitm"]["process"]["metrics"]["single_core_seconds"],
            "peak_rss_bytes": control_docs["direct-mitm"]["process"]["metrics"]["peak_rss_bytes"],
            "conflicts": None,
            "operations": {
                "group_additions": direct_report["group_additions"],
                "pair_entries": direct_report["pair_entries"],
            },
        },
        {
            "backend": "native-xor",
            "status": "unknown_inconclusive",
            "classification": "inconclusive",
            "wall_seconds": control_docs["native-xor"]["process"]["metrics"]["wall_seconds"],
            "total_core_seconds": control_docs["native-xor"]["process"]["metrics"]["total_core_seconds"],
            "single_core_seconds": control_docs["native-xor"]["process"]["metrics"]["single_core_seconds"],
            "peak_rss_bytes": control_docs["native-xor"]["process"]["metrics"]["peak_rss_bytes"],
            "conflicts": native_report["stats"]["conflicts"],
        },
        {
            "backend": "wdsat",
            "status": "timeout_inconclusive",
            "classification": "inconclusive",
            "wall_seconds": control_docs["wdsat-valid"]["process"]["metrics"]["wall_seconds"],
            "total_core_seconds": control_docs["wdsat-valid"]["process"]["metrics"]["total_core_seconds"],
            "single_core_seconds": control_docs["wdsat-valid"]["process"]["metrics"]["single_core_seconds"],
            "peak_rss_bytes": control_docs["wdsat-valid"]["process"]["metrics"]["peak_rss_bytes"],
            "conflicts": None,
        },
        {
            "backend": "cryptominisat",
            "status": "timeout_inconclusive",
            "classification": "inconclusive",
            "wall_seconds": control_docs["cryptominisat"]["process"]["metrics"]["wall_seconds"],
            "total_core_seconds": control_docs["cryptominisat"]["process"]["metrics"]["total_core_seconds"],
            "single_core_seconds": control_docs["cryptominisat"]["process"]["metrics"]["single_core_seconds"],
            "peak_rss_bytes": control_docs["cryptominisat"]["process"]["metrics"]["peak_rss_bytes"],
            "conflicts": None,
        },
        {
            "backend": "licensed-magma-f4",
            "status": "not_run",
            "classification": None,
            "wall_seconds": None,
            "total_core_seconds": None,
            "single_core_seconds": None,
            "peak_rss_bytes": None,
            "conflicts": None,
        },
        {
            "backend": "ggmp",
            "status": "not_same_instance",
            "classification": None,
            "wall_seconds": None,
            "total_core_seconds": None,
            "single_core_seconds": None,
            "peak_rss_bytes": None,
            "conflicts": None,
        },
    ]
    direct_metrics = same_target_table[1]
    comparisons = {
        "f4_over_direct_mitm": {
            "wall_ratio": metrics["wall_seconds"] / direct_metrics["wall_seconds"],
            "core_ratio": metrics["total_core_seconds"] / direct_metrics["total_core_seconds"],
            "rss_ratio": metrics["peak_rss_bytes"] / direct_metrics["peak_rss_bytes"],
        },
        "f4_over_native_xor": {
            "wall_ratio": metrics["wall_seconds"] / same_target_table[2]["wall_seconds"],
            "core_ratio": metrics["total_core_seconds"] / same_target_table[2]["total_core_seconds"],
            "native_xor_inconclusive": True,
            "f4_conclusive": True,
        },
        "selected_positive_vs_fresh_negative": {
            "wall_ratio": metrics["wall_seconds"]
            / baseline["selected_parallel_native_f4"]["metrics"]["wall_seconds"],
            "core_ratio": metrics["total_core_seconds"]
            / baseline["selected_parallel_native_f4"]["metrics"]["total_core_seconds"],
            "system_count_ratio": report["fixed_x1_systems_completed"]
            / baseline["selected_parallel_native_f4"]["report"][
                "fixed_x1_systems_completed"
            ],
        },
        "cold_build_plus_holdout": {
            "wall_seconds": receipt["resources"]["summed_process_wall_seconds"]
            + summary["outer_resources"]["wall_seconds"],
            "total_core_seconds": receipt["resources"]["total_core_seconds"]
            + summary["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": max(
                receipt["resources"]["largest_child_peak_rss_bytes"],
                metrics["peak_rss_bytes"],
            ),
            "build_already_charged_in_stage167_campaign": True,
        },
    }
    prior = baseline["campaign_accounting"]["measured_lower_bound"]
    incremental_rows = list(resource_rows)
    incremental_rows.append(
        {
            "wall_seconds": summary["outer_resources"]["wall_seconds"],
            "total_core_seconds": summary["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": metrics["peak_rss_bytes"],
        }
    )
    incremental = sum_resources(incremental_rows)
    campaign = {
        "stage168_measured_lower_bound": prior,
        "stage169_incremental_measured": incremental,
        "reused_build_already_charged": True,
        "measured_lower_bound": {
            "resource_components": prior["resource_components"] + incremental["resource_components"],
            "summed_wall_seconds": round(
                prior["summed_wall_seconds"] + incremental["summed_wall_seconds"], 12
            ),
            "total_core_seconds": round(
                prior["total_core_seconds"] + incremental["total_core_seconds"], 12
            ),
            "peak_rss_bytes": max(prior["peak_rss_bytes"], incremental["peak_rss_bytes"]),
        },
        "complete_campaign_cost": None,
        "unmetered_scope": [
            "preregistration composition and Git commit",
            "post-run score composition and documentation",
        ],
        "reason_complete_is_null": "some orchestration commands lack complete outer receipts",
    }
    result = {
        "schema": "koblitz_stage169_fresh_target_holdout.v1",
        "status": "complete_preregistered_true_negative",
        "date": "2026-09-23",
        "claim_boundary": prereg["claim_boundary"],
        "preregistration": {
            "path": "preregistration.json",
            "sha256": sha256(STAGE / "preregistration.json"),
            "parent_commit": prereg["preregistration_parent_commit"],
            "truth_consulted_for_selection": False,
        },
        "instance": {
            "blind_instance_id": BLIND_ID,
            "source_system_id": SOURCE_SYSTEM_ID,
            "source_instance_id": SOURCE_INSTANCE_ID,
            "cell_id": "n59-l9-m3-standard-a1-f0",
            "n": 59,
            "ell": 9,
            "m": 3,
            "target_class_opened_after_run_seal": "nondecomposable",
        },
        "factor_base": baseline["factor_base"],
        "fresh_holdout_native_f4": {
            "source_revision": plan["source_revision"],
            "run_inventory_sha256": run_seal["inventory_sha256"],
            "status": selected["status"],
            "classification": "true_negative",
            "metrics": metrics,
            "report": report,
        },
        "same_target_table": same_target_table,
        "comparisons": comparisons,
        "controls": control_docs,
        "historical_same_instance": history,
        "predecessor": {
            "stage168_result_sha256": phase_b.sha256_file(
                STAGE168 / "result.json", "Stage-168 result"
            ),
            "stage168_seal_sha256": phase_b.sha256_file(
                STAGE168 / "result-seal.json", "Stage-168 seal"
            ),
        },
        "campaign_accounting": campaign,
        "gates": {
            **baseline["gates"],
            "2_same_instance_solver_matrix": (
                "partial_two_native_f4_n59_targets_one_preregistered_holdout_"
                "licensed_magma_and_full_panel_missing"
            ),
        },
        "fresh_target_holdout_complete": True,
        "licensed_magma_f4_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    phase_b.write_json_new(STAGE / "result.json", result)
    write_markdown(result)
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    seal = {
        "schema": "koblitz_stage169_fresh_target_holdout_seal.v1",
        "status": "stage169_fresh_holdout_result_frozen",
        "result_sha256": phase_b.sha256_file(STAGE / "result.json", "Stage-169 result"),
        "inventory": inventory,
        "inventory_sha256": phase_b.canonical_sha256(inventory),
    }
    seal["seal_payload_sha256"] = phase_b.canonical_sha256(seal)
    phase_b.write_json_new(STAGE / "result-seal.json", seal)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--score", type=Path, required=True)
    parser.add_argument("--controls-root", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = compose(
            args.run.resolve(strict=True),
            args.score.resolve(strict=True),
            args.controls_root.resolve(strict=True),
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (
        OSError,
        ValueError,
        KeyError,
        IndexError,
        Stage169Error,
        phase_b.PhaseBError,
    ) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Compose the sealed Stage-165 fixed-X1 Hamming-weight schedule result."""

from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
import shutil
from typing import Any

import build_koblitz_phase_b_native_f4 as build_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4
import verify_koblitz_stage164_dense_f4 as stage164_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE164 = GATES / "stage-164-native-f4-dense-20260923"
STAGE = GATES / "stage-165-x1-weight-order-20260923"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
SOURCE_ID = "fb0ea6bb30be3dea841c609a2668939128b4b799d4c0e051b0b52dc6da05f0fa"
COMMIT = "341e46d02d49da94ca24837f7116c812ff2299d3"
WEIGHT_FINGERPRINT = "bd95f10ffd0ab09a2384a2a283fbb85f886df79cd33776b3a3a5a852ccf3bd13"
ASCENDING_FINGERPRINT = "d6969c5f89b5be6a7eacbdac45799fc734a273c9034949a482723b4422504ae7"


class Stage165Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage165Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be an object")
    return value


def only_task(root: Path) -> dict[str, Any]:
    paths = sorted((root / "tasks").glob("*/task-result.json"))
    require(len(paths) == 1, "selected run must contain one task")
    task = load(paths[0], "selected task")
    require(
        task.get("blind_instance_id") == BLIND_ID
        and task.get("source_instance_id") == SOURCE_ID,
        "selected task identity changed",
    )
    return task


def verify_score(root: Path) -> None:
    seal = load(root / "score-seal.json", "score seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(root, {"score-seal.json"})
    require(
        seal.get("schema") == native_f4.SCORE_SEAL_SCHEMA
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and load(root / "score.json", "score").get("classification_counts")
        == {"true_positive": 1},
        "selected score changed",
    )


def copy_file(source: Path, destination: Path, context: str) -> dict[str, Any]:
    require(source.is_file(), f"missing {context}: {source}")
    shutil.copyfile(source, destination)
    return {
        "path": destination.name,
        "bytes": destination.stat().st_size,
        "sha256": phase_b.sha256_file(destination, context),
    }


def sum_resources(rows: list[dict[str, Any]]) -> dict[str, Any]:
    require(rows, "resource rows are empty")
    return {
        "resource_components": len(rows),
        "summed_wall_seconds": round(math.fsum(row["wall_seconds"] for row in rows), 12),
        "total_core_seconds": round(math.fsum(row["total_core_seconds"] for row in rows), 12),
        "peak_rss_bytes": max(row["peak_rss_bytes"] for row in rows),
    }


def copy_development(root: Path, output: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    output.mkdir()
    arms: dict[str, Any] = {}
    resource_rows: list[dict[str, Any]] = []
    for metrics_source in sorted(root.glob("stage165-*.metrics")):
        prefix = metrics_source.name.removesuffix(".metrics")
        name = prefix.removeprefix("stage165-")
        artifacts = {
            suffix: copy_file(
                root / f"{prefix}.{suffix}",
                output / f"{name}.{suffix}",
                f"Stage-165 {name} {suffix}",
            )
            for suffix in ("metrics", "stdout", "stderr")
        }
        metrics_doc = load(output / artifacts["metrics"]["path"], f"{name} metrics")
        metrics = metrics_doc["metrics"]
        require(
            metrics_doc.get("returncode") == 0 and metrics_doc.get("timed_out") is False,
            f"development arm failed: {name}",
        )
        report = load(output / artifacts["stdout"]["path"], f"{name} stdout")
        require(
            report.get("status") == "sat" and report.get("source_witness_valid") is True,
            f"development terminal changed: {name}",
        )
        if name == "weight-r1":
            arm_decision = "rejected_stale_binary_environment_override_not_present"
        elif name == "clean-ascending" or name.startswith("ascending-"):
            arm_decision = "same_binary_ascending_control"
            require(report.get("solver_x1_order") == "ascending_bitmask", f"bad ascending arm: {name}")
        else:
            arm_decision = "weight_order_development_repeat"
            require(
                report.get("solver_x1_order") == "hamming_weight_then_mask",
                f"bad weight arm: {name}",
            )
        resource_rows.append(metrics)
        arms[name] = {
            "decision": arm_decision,
            "metrics": metrics,
            "cost": report["cost"],
            "fixed_x1_masks_visited": report.get("fixed_x1_masks_visited"),
            "fixed_x1_systems_constructed": report["fixed_x1_systems_constructed"],
            "equation_fingerprint": report["solver_equations_blake3"],
            "artifacts": artifacts,
        }
    require(
        {"weight-valid-r1", "weight-r2", "weight-r3", "ascending-r1", "ascending-r2", "ascending-r3", "clean-ascending", "weight-r1"}
        == set(arms),
        "Stage-165 development arm set changed",
    )
    manifest = {
        "schema": "koblitz_stage165_development.v1",
        "status": "weight_order_repeats_and_ascending_controls_retained",
        "arms": arms,
        "complete_development_cost": None,
        "reason_complete_is_null": (
            "incremental compilation, example tests, composition, documentation and Git "
            "operations lack complete outer receipts"
        ),
    }
    phase_b.write_json_new(output / "manifest.json", manifest)
    return manifest, resource_rows


def median(values: list[float]) -> float:
    ordered = sorted(values)
    require(len(ordered) == 3, "three-run median requires exactly three rows")
    return ordered[1]


def write_markdown(result: dict[str, Any]) -> None:
    selected = result["selected_native_f4"]
    metrics = selected["metrics"]
    comparisons = result["comparisons"]
    lines = [
        "# Stage 165: target-independent sparse fixed-X1 order",
        "",
        result["claim_boundary"],
        "",
        "The fixed-X1 outer loop now orders algebraic-subspace coefficient masks by Hamming weight and then numeric mask. The order depends only on the published polynomial-basis coefficients, not on the target, truth class, witness, subgroup enumeration, factor-base logs, or a known scalar. `PQ_F4_X1_ORDER=ascending` retains the prior same-binary schedule.",
        "",
        f"The clean blind F4 process takes {metrics['wall_seconds']:.6f} wall seconds, {metrics['total_core_seconds']:.6f} core-seconds, and {metrics['peak_rss_bytes']} bytes peak RSS. It visits 85 masks, constructs 39 rational fixed-X1 systems, and charges {selected['report']['cost']['ops']:,} elimination-and-table word XORs before returning the same exact curve-verified witness.",
        "",
        f"The clean ascending control takes {comparisons['same_binary_ascending_control']['wall_seconds']:.6f} wall seconds and {comparisons['same_binary_ascending_control']['total_core_seconds']:.6f} core-seconds. The selected clean wall speedup is {comparisons['same_binary_ascending_control']['selected_wall_speedup']:.3f}x; the interleaved three-run development medians give {comparisons['three_run_medians']['wall_speedup']:.3f}x.",
        "",
        f"Direct MITM still wins by {comparisons['selected_over_same_host_direct_mitm']['wall_ratio']:.2f}x wall and {comparisons['selected_over_same_host_direct_mitm']['core_ratio']:.2f}x CPU. The order was selected post hoc after inspecting this target, so the result is a target-specific engineering optimization despite the order itself being target-independent; it is not evidence of generic expected-time scaling.",
        "",
        "Licensed Magma F4, the complete same-instance solver panel, full index-calculus/rho crossover, and unaffiliated reproduction and novelty review remain open. This is not a Koblitz index-calculus SOTA result.",
        "",
    ]
    phase_b.write_new(STAGE / "RESULTS.md", ("\n".join(lines)).encode())


def compose(run_source: Path, score_source: Path, development_root: Path) -> dict[str, Any]:
    require(not STAGE.exists(), f"refusing to overwrite {STAGE}")
    require(stage164_verify.verify().get("status") == "verified", "Stage-164 failed verification")
    baseline = load(STAGE164 / "result.json", "Stage-164 result")
    STAGE.mkdir()
    shutil.copytree(run_source, STAGE / "selected-run")
    shutil.copytree(score_source, STAGE / "selected-score")
    development, development_rows = copy_development(development_root, STAGE / "development")
    run_seal, plan = native_f4.validate_run_seal(STAGE / "selected-run")
    selected = only_task(STAGE / "selected-run")["result"]
    verify_score(STAGE / "selected-score")
    report = selected["backend_report"]
    metrics = selected["metrics"]
    extra = report["cost"]["extra"]
    require(
        selected["status"] == "sat"
        and selected["source_model_valid"] is True
        and selected["source_witness_valid"] is True
        and report["solver_x1_order"] == "hamming_weight_then_mask"
        and report["solver_x1_order_target_independent"] is True
        and report["fixed_x1_masks_visited"] == 85
        and report["fixed_x1_nonrational_skipped"] == 46
        and report["fixed_x1_systems_constructed"] == 39
        and report["fixed_x1_systems_completed"] == 39
        and report["solver_equations_total"] == 2_301
        and report["solver_terms_total"] == 2_295_166
        and report["solver_equations_blake3"] == WEIGHT_FINGERPRINT
        and report["cost"]["ops"] == 16_821_055_616
        and extra["f4_calls"] == 39,
        "selected weighted-order result changed",
    )
    old = baseline["selected_native_f4"]
    require(
        report["witness_points"] == old["report"]["witness_points"],
        "Stage-165 changed the exact witness",
    )
    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path, "selected receipt")
    build_tool.validate_receipt(
        receipt_path, receipt["binaries"], receipt["rust_source_objects"]
    )
    require(
        receipt["source_commit"] == COMMIT and receipt["source_clean"] is True,
        "selected source changed",
    )
    run_summary = load(STAGE / "selected-run/run-summary.json", "run summary")
    ascending = development["arms"]["clean-ascending"]
    require(
        ascending["equation_fingerprint"] == ASCENDING_FINGERPRINT
        and ascending["fixed_x1_masks_visited"] == 138
        and ascending["fixed_x1_systems_constructed"] == 65
        and ascending["cost"]["ops"] == 27_264_366_281,
        "clean ascending control changed",
    )
    weight_development = [development["arms"][name] for name in ("weight-valid-r1", "weight-r2", "weight-r3")]
    ascending_development = [development["arms"][name] for name in ("ascending-r1", "ascending-r2", "ascending-r3")]
    weight_wall_median = median([arm["metrics"]["wall_seconds"] for arm in weight_development])
    ascending_wall_median = median([arm["metrics"]["wall_seconds"] for arm in ascending_development])
    weight_core_median = median([arm["metrics"]["total_core_seconds"] for arm in weight_development])
    ascending_core_median = median([arm["metrics"]["total_core_seconds"] for arm in ascending_development])
    table = copy.deepcopy(baseline["same_target_table"])
    f4 = next(row for row in table if row["backend"] == "native-f4-selected")
    f4.update(
        {
            "formulation": "fixed-X1 direct S4; Hamming-weight coefficient order; dense F4 pipeline",
            "wall_seconds": metrics["wall_seconds"],
            "total_core_seconds": metrics["total_core_seconds"],
            "peak_rss_bytes": metrics["peak_rss_bytes"],
            "operations": {
                "word_xors_including_m4ri_tables": report["cost"]["ops"],
                "f4_calls": extra["f4_calls"],
                "fixed_x1_masks_visited": report["fixed_x1_masks_visited"],
            },
        }
    )
    mitm = next(row for row in table if row["backend"] == "direct-mitm")
    comparisons = {
        "same_binary_ascending_control": {
            "wall_seconds": ascending["metrics"]["wall_seconds"],
            "total_core_seconds": ascending["metrics"]["total_core_seconds"],
            "peak_rss_bytes": ascending["metrics"]["peak_rss_bytes"],
            "selected_wall_ratio": metrics["wall_seconds"] / ascending["metrics"]["wall_seconds"],
            "selected_wall_speedup": ascending["metrics"]["wall_seconds"] / metrics["wall_seconds"],
            "selected_core_ratio": metrics["total_core_seconds"] / ascending["metrics"]["total_core_seconds"],
            "selected_word_xor_ratio": report["cost"]["ops"] / ascending["cost"]["ops"],
            "same_verified_witness": True,
            "equation_prefixes_differ_by_schedule": True,
        },
        "three_run_medians": {
            "weight_wall_seconds": weight_wall_median,
            "ascending_wall_seconds": ascending_wall_median,
            "wall_speedup": ascending_wall_median / weight_wall_median,
            "weight_core_seconds": weight_core_median,
            "ascending_core_seconds": ascending_core_median,
            "core_speedup": ascending_core_median / weight_core_median,
        },
        "selected_over_stage164": {
            "wall_ratio": metrics["wall_seconds"] / old["metrics"]["wall_seconds"],
            "wall_speedup": old["metrics"]["wall_seconds"] / metrics["wall_seconds"],
            "core_ratio": metrics["total_core_seconds"] / old["metrics"]["total_core_seconds"],
            "rss_ratio": metrics["peak_rss_bytes"] / old["metrics"]["peak_rss_bytes"],
            "word_xor_ratio": report["cost"]["ops"] / old["report"]["cost"]["ops"],
            "same_verified_witness": True,
        },
        "selected_over_same_host_direct_mitm": {
            "wall_ratio": metrics["wall_seconds"] / mitm["wall_seconds"],
            "core_ratio": metrics["total_core_seconds"] / mitm["total_core_seconds"],
            "rss_ratio": metrics["peak_rss_bytes"] / mitm["peak_rss_bytes"],
        },
        "selected_cold_build_plus_run": {
            "wall_seconds": receipt["resources"]["summed_process_wall_seconds"]
            + run_summary["outer_resources"]["wall_seconds"],
            "total_core_seconds": receipt["resources"]["total_core_seconds"]
            + run_summary["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": max(
                receipt["resources"]["largest_child_peak_rss_bytes"],
                metrics["peak_rss_bytes"],
            ),
        },
    }
    prior = baseline["campaign_accounting"]["measured_lower_bound"]
    incremental_rows = list(development_rows)
    incremental_rows.extend(row["metrics"] for row in receipt["build_processes"])
    incremental_rows.append(
        {
            "wall_seconds": run_summary["outer_resources"]["wall_seconds"],
            "total_core_seconds": run_summary["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": metrics["peak_rss_bytes"],
        }
    )
    incremental = sum_resources(incremental_rows)
    campaign = {
        "stage164_measured_lower_bound": prior,
        "stage165_incremental_measured": incremental,
        "measured_lower_bound": {
            "resource_components": prior["resource_components"] + incremental["resource_components"],
            "summed_wall_seconds": round(prior["summed_wall_seconds"] + incremental["summed_wall_seconds"], 12),
            "total_core_seconds": round(prior["total_core_seconds"] + incremental["total_core_seconds"], 12),
            "peak_rss_bytes": max(prior["peak_rss_bytes"], incremental["peak_rss_bytes"]),
        },
        "unmetered_development_commands_present": True,
        "unmetered_scope": [
            "incremental compilation and example tests",
            "composition, documentation and Git operations",
        ],
        "complete_campaign_cost": None,
        "reason_complete_is_null": "some development commands lack complete outer receipts",
    }
    result = {
        "schema": "koblitz_stage165_x1_weight_order.v1",
        "status": "complete_weight_order_true_positive",
        "date": "2026-09-23",
        "claim_boundary": (
            "Strong target-specific internal engineering on one finite public toy PDP. The "
            "order is target-independent, but it was selected post hoc after inspecting this "
            "target. This is not generic expected-time evidence, a full solver panel, a full "
            "index-calculus/rho crossover, independent reproduction, novelty review, or SOTA."
        ),
        "instance": baseline["instance"],
        "factor_base": baseline["factor_base"],
        "mechanism": {
            "name": "fixed_x1_hamming_weight_then_mask",
            "definition": "sort coefficient masks by (popcount(mask), mask)",
            "target_independent": True,
            "uses_discrete_log_labels": False,
            "uses_known_witness": False,
            "same_binary_control": "PQ_F4_X1_ORDER=ascending",
            "selection_caveat": "chosen post hoc after inspecting this target; generalization untested",
            "permutation_test": "fixed_x1_weight_order_is_an_exact_target_independent_permutation",
        },
        "selected_native_f4": {
            "source_revision": plan["source_revision"],
            "run_inventory_sha256": run_seal["inventory_sha256"],
            "status": selected["status"],
            "classification": "true_positive",
            "source_model_valid": selected["source_model_valid"],
            "source_witness_valid": selected["source_witness_valid"],
            "metrics": metrics,
            "report": report,
        },
        "same_target_table": table,
        "comparisons": comparisons,
        "development": development,
        "predecessor": {
            "stage164_result_sha256": phase_b.sha256_file(STAGE164 / "result.json", "Stage-164 result"),
            "stage164_seal_sha256": phase_b.sha256_file(STAGE164 / "result-seal.json", "Stage-164 seal"),
        },
        "campaign_accounting": campaign,
        "gates": baseline["gates"],
        "licensed_magma_f4_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    phase_b.write_json_new(STAGE / "result.json", result)
    write_markdown(result)
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    seal = {
        "schema": "koblitz_stage165_x1_weight_order_seal.v1",
        "status": "stage165_x1_weight_order_result_frozen",
        "result_sha256": phase_b.sha256_file(STAGE / "result.json", "Stage-165 result"),
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
    parser.add_argument("--development-root", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = compose(
            args.run.resolve(strict=True),
            args.score.resolve(strict=True),
            args.development_root.resolve(strict=True),
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, IndexError, Stage165Error, phase_b.PhaseBError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()

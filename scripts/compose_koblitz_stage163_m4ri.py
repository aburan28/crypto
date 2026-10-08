#!/usr/bin/env python3
"""Compose the sealed Stage-163 block-4 M4RI result."""

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
import verify_koblitz_stage162_pair_selection as stage162_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE162 = GATES / "stage-162-grouped-pair-selection-20260923"
STAGE = GATES / "stage-163-m4ri-echelons-20260923"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
SOURCE_ID = "fb0ea6bb30be3dea841c609a2668939128b4b799d4c0e051b0b52dc6da05f0fa"
CELL = "n59-l9-m3-standard-a1-f0"


class Stage163Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage163Error(message)


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
        and task.get("source_instance_id") == SOURCE_ID
        and task.get("cell_id") == CELL,
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


def metric_rows(root: Path) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(root.rglob("*.metrics.json")):
        value = load(path, str(path))
        if isinstance(value.get("metrics"), dict):
            rows.append(value["metrics"])
    return rows


def sum_resources(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "resource_components": len(rows),
        "summed_wall_seconds": round(math.fsum(row["wall_seconds"] for row in rows), 12),
        "total_core_seconds": round(math.fsum(row["total_core_seconds"] for row in rows), 12),
        "peak_rss_bytes": max(row["peak_rss_bytes"] for row in rows),
    }


def copy_development(
    root: Path,
    failed_run: Path,
    output: Path,
) -> dict[str, Any]:
    output.mkdir()
    arms = {}
    for name, prefix in (
        ("streaming-control", "stage163-streaming-control-f4"),
        ("m4ri-provisional", "stage163-m4ri-f4"),
    ):
        artifacts = {
            suffix: copy_file(
                root / f"{prefix}.{suffix}",
                output / f"{name}.{suffix}",
                f"{name} {suffix}",
            )
            for suffix in ("metrics", "stdout", "stderr")
        }
        metrics = load(output / artifacts["metrics"]["path"], f"{name} metrics")["metrics"]
        report = load(output / artifacts["stdout"]["path"], f"{name} stdout")
        require(
            report.get("status") == "sat"
            and report.get("source_witness_valid") is True
            and report.get("solver_equations_blake3")
            == "d6969c5f89b5be6a7eacbdac45799fc734a273c9034949a482723b4422504ae7",
            f"{name} terminal changed",
        )
        arms[name] = {
            "metrics": metrics,
            "cost": report["cost"],
            "artifacts": artifacts,
        }
    arms["streaming-control"]["decision"] = "same_binary_control"
    arms["m4ri-provisional"]["decision"] = "promoted_after_clean_rebuild"

    failed_destination = output / "parser-rejected-clean-run"
    shutil.copytree(failed_run, failed_destination)
    require(
        not (failed_destination / "run-seal.json").exists()
        and not list((failed_destination / "tasks").glob("*/task-result.json")),
        "parser-rejected run unexpectedly has a scientific seal",
    )
    failed_task = next((failed_destination / "tasks").iterdir())
    failed_report = load(failed_task / "native-f4.stdout", "failed run backend output")
    require(
        failed_report.get("status") == "sat"
        and failed_report.get("source_witness_valid") is True
        and failed_report["cost"]["op_unit"]
        == "word XORs (elimination, including M4RI tables)",
        "failed run backend output changed",
    )
    failed_receipt = load(
        failed_destination / "tool-builds/rust/receipt.json", "failed run build receipt"
    )
    failed_tasks = sum_resources(metric_rows(failed_task))
    failed_inventory = phase_b.all_regular_inventory(failed_destination)
    manifest = {
        "schema": "koblitz_stage163_development.v1",
        "status": "same_binary_pair_and_parser_rejection_retained",
        "arms": arms,
        "parser_rejected_clean_run": {
            "status": "backend_contract_error_before_parser_compatibility_fix",
            "backend_status": "sat",
            "scientific_terminal_promoted": False,
            "source_commit": failed_receipt["source_commit"],
            "build_resources": failed_receipt["resources"],
            "task_process_resources": failed_tasks,
            "inventory_sha256": phase_b.canonical_sha256(failed_inventory),
            "inventory": failed_inventory,
            "path": failed_destination.name,
        },
        "complete_development_cost": None,
        "reason_complete_is_null": "incremental compile and focused tests lack outer receipts",
    }
    phase_b.write_json_new(output / "manifest.json", manifest)
    return manifest


def write_markdown(result: dict[str, Any]) -> None:
    s = result["selected_native_f4"]
    c = result["comparisons"]
    lines = [
        "# Stage 163: shape-selected block-4 M4RI F4 elimination",
        "",
        result["claim_boundary"],
        "",
        "Large dense matrices now use block-4 Method of Four Russians elimination when rows >= 128, columns >= 256, and columns <= 4*rows. Smaller or wider matrices retain streaming elimination. `PQ_F4_DISABLE_M4RI=1` is the same-binary control. Four large synthetic shapes pass exact pivot-column and canonical-row-space equality against streaming elimination.",
        "",
        f"The same-binary M4RI arm takes {c['same_binary_m4ri_over_streaming']['m4ri_wall_seconds']:.6f} seconds versus {c['same_binary_m4ri_over_streaming']['streaming_wall_seconds']:.6f} for streaming, a {c['same_binary_m4ri_over_streaming']['wall_speedup']:.3f}x gain. The clean selected process takes {s['metrics']['wall_seconds']:.6f} wall seconds, {s['metrics']['total_core_seconds']:.6f} core-seconds, and {s['metrics']['peak_rss_bytes']} bytes RSS.",
        "",
        f"Selected M4RI routes {s['report']['cost']['extra']['m4ri_matrices']} matrices and reduces counted table-inclusive word XORs from 87,513,949,370 to {s['report']['cost']['ops']:,}. Direct MITM still wins by {c['selected_over_same_host_direct_mitm']['wall_ratio']:.2f}x wall and {c['selected_over_same_host_direct_mitm']['core_ratio']:.2f}x CPU. Clean build plus run costs {c['selected_cold_build_plus_run']['wall_seconds']:.3f} wall seconds.",
        "",
        "Licensed Magma, the full native-F4 panel, end-to-end IC/rho cost, and independent external reproduction remain open. This does not establish a SOTA.",
        "",
    ]
    phase_b.write_new(STAGE / "RESULTS.md", ("\n".join(lines)).encode())


def compose(
    run_source: Path,
    score_source: Path,
    development_root: Path,
    failed_run: Path,
) -> dict[str, Any]:
    require(not STAGE.exists(), f"refusing to overwrite {STAGE}")
    require(stage162_verify.verify().get("status") == "verified", "Stage-162 failed verification")
    baseline = load(STAGE162 / "result.json", "Stage-162 result")
    STAGE.mkdir()
    shutil.copytree(run_source, STAGE / "selected-run")
    shutil.copytree(score_source, STAGE / "selected-score")
    development = copy_development(development_root, failed_run, STAGE / "development")
    run_seal, plan = native_f4.validate_run_seal(STAGE / "selected-run")
    selected = only_task(STAGE / "selected-run")["result"]
    verify_score(STAGE / "selected-score")
    report = selected["backend_report"]
    metrics = selected["metrics"]
    require(
        selected["status"] == "sat"
        and selected["source_model_valid"] is True
        and selected["source_witness_valid"] is True
        and report["solver_echelon_policy"].startswith("block4_m4ri")
        and report["cost"]["extra"]["m4ri_matrices"] == 196
        and report["cost"]["extra"]["m4ri_table_word_xors"] == 316_714_298
        and report["cost"]["ops"] == 45_879_309_338,
        "selected M4RI result changed",
    )
    old = baseline["selected_native_f4"]
    require(
        report["solver_equations_blake3"] == old["report"]["solver_equations_blake3"]
        and report["solver_terms_total"] == old["report"]["solver_terms_total"]
        and report["witness_points"] == old["report"]["witness_points"],
        "M4RI changed equations or witness",
    )
    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path, "selected receipt")
    build_tool.validate_receipt(receipt_path, receipt["binaries"], receipt["rust_source_objects"])
    require(
        receipt["source_commit"] == "04e49c930caccea31bad7172aa8f567d96744b17"
        and receipt["source_clean"] is True,
        "selected source changed",
    )
    run_summary = load(STAGE / "selected-run/run-summary.json", "run summary")
    table = copy.deepcopy(baseline["same_target_table"])
    f4 = next(row for row in table if row["backend"] == "native-f4-selected")
    f4.update(
        {
            "formulation": "fixed-X1 direct S4; fast masks; grouped pairs; shape-selected block-4 M4RI",
            "wall_seconds": metrics["wall_seconds"],
            "total_core_seconds": metrics["total_core_seconds"],
            "peak_rss_bytes": metrics["peak_rss_bytes"],
            "operations": {
                "word_xors_including_m4ri_tables": report["cost"]["ops"],
                "f4_calls": report["cost"]["extra"]["f4_calls"],
                "m4ri_matrices": report["cost"]["extra"]["m4ri_matrices"],
            },
        }
    )
    mitm = next(row for row in table if row["backend"] == "direct-mitm")
    control = development["arms"]["streaming-control"]
    provisional = development["arms"]["m4ri-provisional"]
    comparisons = {
        "same_binary_m4ri_over_streaming": {
            "streaming_wall_seconds": control["metrics"]["wall_seconds"],
            "m4ri_wall_seconds": provisional["metrics"]["wall_seconds"],
            "wall_ratio": provisional["metrics"]["wall_seconds"] / control["metrics"]["wall_seconds"],
            "wall_speedup": control["metrics"]["wall_seconds"] / provisional["metrics"]["wall_seconds"],
            "core_ratio": provisional["metrics"]["total_core_seconds"] / control["metrics"]["total_core_seconds"],
            "word_xor_ratio": provisional["cost"]["ops"] / control["cost"]["ops"],
            "same_equation_fingerprint": True,
            "same_verified_witness": True,
        },
        "selected_over_stage162": {
            "wall_ratio": metrics["wall_seconds"] / old["metrics"]["wall_seconds"],
            "wall_speedup": old["metrics"]["wall_seconds"] / metrics["wall_seconds"],
            "core_ratio": metrics["total_core_seconds"] / old["metrics"]["total_core_seconds"],
            "rss_ratio": metrics["peak_rss_bytes"] / old["metrics"]["peak_rss_bytes"],
            "word_xor_ratio": report["cost"]["ops"] / old["report"]["cost"]["ops"],
            "same_equation_fingerprint": True,
            "same_verified_witness": True,
        },
        "selected_over_same_host_direct_mitm": {
            "wall_ratio": metrics["wall_seconds"] / mitm["wall_seconds"],
            "core_ratio": metrics["total_core_seconds"] / mitm["total_core_seconds"],
            "rss_ratio": metrics["peak_rss_bytes"] / mitm["peak_rss_bytes"],
        },
        "selected_cold_build_plus_run": {
            "wall_seconds": receipt["resources"]["summed_process_wall_seconds"] + run_summary["outer_resources"]["wall_seconds"],
            "total_core_seconds": receipt["resources"]["total_core_seconds"] + run_summary["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": max(receipt["resources"]["largest_child_peak_rss_bytes"], metrics["peak_rss_bytes"]),
        },
    }
    prior = baseline["campaign_accounting"]["measured_lower_bound"]
    failed = development["parser_rejected_clean_run"]
    incremental = [
        control["metrics"], provisional["metrics"],
        {"wall_seconds": failed["build_resources"]["summed_process_wall_seconds"], "total_core_seconds": failed["build_resources"]["total_core_seconds"], "peak_rss_bytes": failed["build_resources"]["largest_child_peak_rss_bytes"]},
        {"wall_seconds": failed["task_process_resources"]["summed_wall_seconds"], "total_core_seconds": failed["task_process_resources"]["total_core_seconds"], "peak_rss_bytes": failed["task_process_resources"]["peak_rss_bytes"]},
        {"wall_seconds": receipt["resources"]["summed_process_wall_seconds"], "total_core_seconds": receipt["resources"]["total_core_seconds"], "peak_rss_bytes": receipt["resources"]["largest_child_peak_rss_bytes"]},
        {"wall_seconds": run_summary["outer_resources"]["wall_seconds"], "total_core_seconds": run_summary["outer_resources"]["total_core_seconds"], "peak_rss_bytes": metrics["peak_rss_bytes"]},
    ]
    inc = sum_resources(incremental)
    campaign = {
        "stage162_measured_lower_bound": prior,
        "stage163_incremental_measured": inc,
        "measured_lower_bound": {
            "resource_components": prior["resource_components"] + inc["resource_components"],
            "summed_wall_seconds": round(prior["summed_wall_seconds"] + inc["summed_wall_seconds"], 12),
            "total_core_seconds": round(prior["total_core_seconds"] + inc["total_core_seconds"], 12),
            "peak_rss_bytes": max(prior["peak_rss_bytes"], inc["peak_rss_bytes"]),
        },
        "unmetered_development_commands_present": True,
        "unmetered_scope": ["incremental compiles and focused tests", "failed runner outer-driver overhead", "score/evidence composition and Git operations"],
        "complete_campaign_cost": None,
        "reason_complete_is_null": "some development commands lack complete outer receipts",
    }
    result = {
        "schema": "koblitz_stage163_m4ri_echelons.v1",
        "status": "complete_m4ri_true_positive",
        "date": "2026-09-23",
        "claim_boundary": "Strong internal engineering and one finite public toy-PDP native-F4 M4RI improvement. It is not a full solver panel, direct-MITM improvement, full IC/rho crossover, independent reproduction, or Koblitz index-calculus SOTA.",
        "instance": baseline["instance"],
        "factor_base": baseline["factor_base"],
        "mechanism": {
            "name": "shape_selected_block4_m4ri",
            "selection": "rows >= 128, columns >= 256, columns <= 4*rows",
            "same_binary_disable": "PQ_F4_DISABLE_M4RI=1",
            "row_space_test": "m4ri_and_streaming_have_the_same_row_space",
            "large_matrix_test_shapes": 4,
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
            "stage162_result_sha256": phase_b.sha256_file(STAGE162 / "result.json", "Stage-162 result"),
            "stage162_seal_sha256": phase_b.sha256_file(STAGE162 / "result-seal.json", "Stage-162 seal"),
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
        "schema": "koblitz_stage163_m4ri_echelons_seal.v1",
        "status": "stage163_m4ri_result_frozen",
        "result_sha256": phase_b.sha256_file(STAGE / "result.json", "Stage-163 result"),
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
    parser.add_argument("--failed-run", type=Path, required=True)
    args = parser.parse_args()
    try:
        value = compose(
            args.run.resolve(strict=True), args.score.resolve(strict=True),
            args.development_root.resolve(strict=True), args.failed_run.resolve(strict=True),
        )
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, IndexError, Stage163Error, phase_b.PhaseBError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()

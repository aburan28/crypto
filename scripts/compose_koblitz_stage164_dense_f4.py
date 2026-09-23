#!/usr/bin/env python3
"""Compose the sealed Stage-164 dense native-F4 optimization result."""

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
import verify_koblitz_stage163_m4ri as stage163_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE163 = GATES / "stage-163-m4ri-echelons-20260923"
STAGE = GATES / "stage-164-native-f4-dense-20260923"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
SOURCE_ID = "fb0ea6bb30be3dea841c609a2668939128b4b799d4c0e051b0b52dc6da05f0fa"
CELL = "n59-l9-m3-standard-a1-f0"
COMMIT = "9845f2df87c1aa759254c42739f11eef632ab2f2"
EQUATIONS = "d6969c5f89b5be6a7eacbdac45799fc734a273c9034949a482723b4422504ae7"


class Stage164Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage164Error(message)


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
    require(rows, "resource rows are empty")
    return {
        "resource_components": len(rows),
        "summed_wall_seconds": round(math.fsum(row["wall_seconds"] for row in rows), 12),
        "total_core_seconds": round(math.fsum(row["total_core_seconds"] for row in rows), 12),
        "peak_rss_bytes": max(row["peak_rss_bytes"] for row in rows),
    }


def decision(name: str) -> str:
    if name in {"fused-active-f4", "hash-grouped-pairs-f4", "pivot-ends-f4", "unique-lcms-f4"}:
        return "rejected_wall_regression"
    if name in {"clean-streaming", "clean-batch-disabled", "clean-b7"}:
        return "same_clean_binary_control"
    if name.startswith("finaltrim-b") or name.startswith("m4ri-b") or name.startswith("final-b"):
        return "block_width_ablation"
    if name == "dependency-fetch":
        return "separately_charged_dependency_acquisition"
    if name == "trimmed-tables-f4":
        return "promoted_after_clean_rebuild"
    return "retained_development_checkpoint"


def copy_development(
    development_root: Path,
    failed_builds: list[Path],
    output: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    output.mkdir()
    arms: dict[str, Any] = {}
    resource_rows: list[dict[str, Any]] = []
    for metrics_source in sorted(development_root.glob("stage164-*.metrics")):
        prefix = metrics_source.name.removesuffix(".metrics")
        name = prefix.removeprefix("stage164-")
        require(name not in arms, f"duplicate development arm: {name}")
        artifacts = {}
        for suffix in ("metrics", "stdout", "stderr"):
            source = development_root / f"{prefix}.{suffix}"
            artifacts[suffix] = copy_file(
                source,
                output / f"{name}.{suffix}",
                f"Stage-164 {name} {suffix}",
            )
        metrics_doc = load(output / artifacts["metrics"]["path"], f"{name} metrics")
        metrics = metrics_doc["metrics"]
        require(
            metrics_doc.get("returncode") == 0
            and metrics_doc.get("timed_out") is False,
            f"development arm failed: {name}",
        )
        resource_rows.append(metrics)
        arm: dict[str, Any] = {
            "decision": decision(name),
            "metrics": metrics,
            "artifacts": artifacts,
        }
        stdout = (output / artifacts["stdout"]["path"]).read_text()
        if stdout.strip().startswith("{"):
            report = json.loads(stdout)
            if report.get("backend") == "native-f4":
                require(
                    report.get("status") == "sat"
                    and report.get("source_model_valid") is True
                    and report.get("source_witness_valid") is True
                    and report.get("solver_equations_blake3") == EQUATIONS,
                    f"development F4 terminal changed: {name}",
                )
                arm["cost"] = report["cost"]
                arm["solver_echelon_policy"] = report["solver_echelon_policy"]
        arms[name] = arm

    failed: dict[str, Any] = {}
    for index, source in enumerate(failed_builds, 1):
        destination = output / f"failed-build-{index:02d}"
        shutil.copytree(source / "evidence", destination)
        inventory = phase_b.all_regular_inventory(destination)
        rows = metric_rows(destination)
        resource_rows.extend(rows)
        vendor_stderr = (destination / "vendor.stderr").read_text()
        require(
            "redis v0.32.7" in vendor_stderr and "--offline was specified" in vendor_stderr,
            "failed build reason changed",
        )
        failed[destination.name] = {
            "status": "offline_vendor_cache_miss_retained",
            "scientific_terminal_promoted": False,
            "resources": sum_resources(rows),
            "inventory": inventory,
            "inventory_sha256": phase_b.canonical_sha256(inventory),
            "path": destination.name,
        }

    require(
        {"clean-streaming", "clean-batch-disabled", "clean-b7", "dependency-fetch"}
        <= set(arms),
        "required Stage-164 controls are absent",
    )
    manifest = {
        "schema": "koblitz_stage164_development.v1",
        "status": "optimization_ledger_and_controls_retained",
        "arms": arms,
        "failed_clean_builds": failed,
        "complete_development_cost": None,
        "reason_complete_is_null": (
            "incremental compiles, focused tests, one source-config vendor replay, "
            "composition and Git operations lack complete outer receipts"
        ),
    }
    phase_b.write_json_new(output / "manifest.json", manifest)
    return manifest, resource_rows


def write_markdown(result: dict[str, Any]) -> None:
    selected = result["selected_native_f4"]
    comparisons = result["comparisons"]
    metrics = selected["metrics"]
    extra = selected["report"]["cost"]["extra"]
    lines = [
        "# Stage 164: dense native-F4 pipeline optimization",
        "",
        result["claim_boundary"],
        "",
        "The repository native F4 now combines a measured block-8 M4RI kernel, cache-local and suffix-trimmed combination tables, consecutive-pivot pattern extraction, dense exact pair grouping and column indexing on the 18-variable systems, indexed symbolic reducers, and order-preserving batch basis installation. The generic hash/sort and streaming fallbacks remain available.",
        "",
        f"The clean blind F4 process takes {metrics['wall_seconds']:.6f} wall seconds, {metrics['total_core_seconds']:.6f} core-seconds, and {metrics['peak_rss_bytes']} bytes peak RSS. It performs {selected['report']['cost']['ops']:,} charged elimination-and-table word XORs, routes {extra['m4ri_matrices']} matrices through M4RI, and returns the same exact curve-verified witness.",
        "",
        f"The same clean binary takes {comparisons['same_binary_streaming_control']['wall_seconds']:.6f} seconds with M4RI disabled and {comparisons['same_binary_batch_disabled']['wall_seconds']:.6f} seconds with batch installation disabled. Relative to Stage 163, selected wall improves by {comparisons['selected_over_stage163']['wall_speedup']:.3f}x. Direct MITM still wins by {comparisons['selected_over_same_host_direct_mitm']['wall_ratio']:.2f}x wall and {comparisons['selected_over_same_host_direct_mitm']['core_ratio']:.2f}x CPU.",
        "",
        f"Dependency acquisition, a clean one-job build, and the custody run total {comparisons['selected_cold_dependency_build_plus_run']['wall_seconds']:.3f} charged wall seconds. The measured campaign lower bound is {result['campaign_accounting']['measured_lower_bound']['summed_wall_seconds']:.6f} sequential wall-seconds and {result['campaign_accounting']['measured_lower_bound']['total_core_seconds']:.6f} core-seconds; complete campaign cost remains null because some development commands were not outer-metered.",
        "",
        "Licensed Magma F4, the complete same-instance solver panel at all requested sizes, end-to-end index-calculus/rho crossover, and unaffiliated reproduction and novelty review remain open. This is not a Koblitz index-calculus SOTA result.",
        "",
    ]
    phase_b.write_new(STAGE / "RESULTS.md", ("\n".join(lines)).encode())


def compose(
    run_source: Path,
    score_source: Path,
    development_root: Path,
    failed_builds: list[Path],
) -> dict[str, Any]:
    require(not STAGE.exists(), f"refusing to overwrite {STAGE}")
    require(stage163_verify.verify().get("status") == "verified", "Stage-163 failed verification")
    baseline = load(STAGE163 / "result.json", "Stage-163 result")
    STAGE.mkdir()
    shutil.copytree(run_source, STAGE / "selected-run")
    shutil.copytree(score_source, STAGE / "selected-score")
    development, development_rows = copy_development(
        development_root, failed_builds, STAGE / "development"
    )
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
        and report["solver_echelon_policy"].startswith("shape_selected_m4ri_default_block8")
        and report["solver_pair_selector"].startswith("epoch_dense_lcm_groups")
        and report["solver_column_index"].startswith("dense_monomial_domain")
        and report["solver_equations_blake3"] == EQUATIONS
        and report["solver_terms_total"] == 3_864_601
        and extra["m4ri_matrices"] == 196
        and extra["m4ri_block_width_max"] == 8
        and extra["pair_dense_select_calls"] == 269_235
        and extra["dense_column_matrices"] == 391
        and report["cost"]["ops"] == 27_264_366_281,
        "selected dense-F4 result changed",
    )
    old = baseline["selected_native_f4"]
    require(
        report["solver_equations_blake3"] == old["report"]["solver_equations_blake3"]
        and report["solver_terms_total"] == old["report"]["solver_terms_total"]
        and report["witness_points"] == old["report"]["witness_points"],
        "Stage-164 changed equations or witness",
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
    table = copy.deepcopy(baseline["same_target_table"])
    f4 = next(row for row in table if row["backend"] == "native-f4-selected")
    f4.update(
        {
            "formulation": "fixed-X1 direct S4; dense exact pair/column indices; batch UPDATE; trimmed block-8 M4RI",
            "wall_seconds": metrics["wall_seconds"],
            "total_core_seconds": metrics["total_core_seconds"],
            "peak_rss_bytes": metrics["peak_rss_bytes"],
            "operations": {
                "word_xors_including_m4ri_tables": report["cost"]["ops"],
                "f4_calls": extra["f4_calls"],
                "m4ri_matrices": extra["m4ri_matrices"],
                "pair_dense_select_calls": extra["pair_dense_select_calls"],
            },
        }
    )
    mitm = next(row for row in table if row["backend"] == "direct-mitm")
    controls = development["arms"]
    streaming = controls["clean-streaming"]
    batch_disabled = controls["clean-batch-disabled"]
    block7 = controls["clean-b7"]
    fetch = controls["dependency-fetch"]
    require(
        streaming["cost"]["ops"] == 87_513_949_370
        and streaming["cost"]["extra"]["m4ri_matrices"] == 0
        and batch_disabled["cost"]["extra"]["batch_insert_groups"] == 0
        and block7["cost"]["extra"]["m4ri_block_width_max"] == 7,
        "same-binary controls changed",
    )
    comparisons = {
        "same_binary_streaming_control": {
            "wall_seconds": streaming["metrics"]["wall_seconds"],
            "total_core_seconds": streaming["metrics"]["total_core_seconds"],
            "peak_rss_bytes": streaming["metrics"]["peak_rss_bytes"],
            "selected_wall_ratio": metrics["wall_seconds"] / streaming["metrics"]["wall_seconds"],
            "selected_wall_speedup": streaming["metrics"]["wall_seconds"] / metrics["wall_seconds"],
            "selected_core_ratio": metrics["total_core_seconds"] / streaming["metrics"]["total_core_seconds"],
            "selected_word_xor_ratio": report["cost"]["ops"] / streaming["cost"]["ops"],
            "same_equation_fingerprint": True,
            "same_verified_witness": True,
        },
        "same_binary_batch_disabled": {
            "wall_seconds": batch_disabled["metrics"]["wall_seconds"],
            "total_core_seconds": batch_disabled["metrics"]["total_core_seconds"],
            "selected_wall_ratio": metrics["wall_seconds"] / batch_disabled["metrics"]["wall_seconds"],
            "selected_wall_speedup": batch_disabled["metrics"]["wall_seconds"] / metrics["wall_seconds"],
            "selected_core_ratio": metrics["total_core_seconds"] / batch_disabled["metrics"]["total_core_seconds"],
            "same_equation_fingerprint": True,
            "same_verified_witness": True,
        },
        "same_binary_block7_control": {
            "wall_seconds": block7["metrics"]["wall_seconds"],
            "total_core_seconds": block7["metrics"]["total_core_seconds"],
            "selected_wall_ratio": metrics["wall_seconds"] / block7["metrics"]["wall_seconds"],
            "selected_core_ratio": metrics["total_core_seconds"] / block7["metrics"]["total_core_seconds"],
            "selection_basis": "three-run development medians favored block 8; this clean control sample favored block 7",
        },
        "selected_over_stage163": {
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
        "selected_cold_dependency_build_plus_run": {
            "wall_seconds": fetch["metrics"]["wall_seconds"]
            + receipt["resources"]["summed_process_wall_seconds"]
            + run_summary["outer_resources"]["wall_seconds"],
            "total_core_seconds": fetch["metrics"]["total_core_seconds"]
            + receipt["resources"]["total_core_seconds"]
            + run_summary["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": max(
                fetch["metrics"]["peak_rss_bytes"],
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
        "stage163_measured_lower_bound": prior,
        "stage164_incremental_measured": incremental,
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
        "unmetered_development_commands_present": True,
        "unmetered_scope": [
            "incremental compiles and focused tests",
            "one source-config vendor replay",
            "composition, documentation and Git operations",
        ],
        "complete_campaign_cost": None,
        "reason_complete_is_null": "some development commands lack complete outer receipts",
    }
    result = {
        "schema": "koblitz_stage164_native_f4_dense.v1",
        "status": "complete_dense_f4_true_positive",
        "date": "2026-09-23",
        "claim_boundary": (
            "Strong internal engineering and one finite public toy-PDP native-F4 improvement. "
            "It is not a full solver panel, full index-calculus/rho crossover, independent "
            "reproduction, novelty review, or Koblitz index-calculus SOTA."
        ),
        "instance": baseline["instance"],
        "factor_base": baseline["factor_base"],
        "mechanisms": {
            "pair_selection": report["solver_pair_selector"],
            "pair_installation": report["solver_pair_installer"],
            "symbolic_reducers": report["solver_symbolic_reducer_selector"],
            "column_index": report["solver_column_index"],
            "echelon": report["solver_echelon_policy"],
            "same_binary_batch_disable": "PQ_F4_DISABLE_BATCH_INSERT=1",
            "same_binary_m4ri_disable": "PQ_F4_DISABLE_M4RI=1",
            "width_override": "PQ_F4_M4RI_BLOCK=2..10",
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
            "stage163_result_sha256": phase_b.sha256_file(
                STAGE163 / "result.json", "Stage-163 result"
            ),
            "stage163_seal_sha256": phase_b.sha256_file(
                STAGE163 / "result-seal.json", "Stage-163 seal"
            ),
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
        "schema": "koblitz_stage164_native_f4_dense_seal.v1",
        "status": "stage164_native_f4_dense_result_frozen",
        "result_sha256": phase_b.sha256_file(STAGE / "result.json", "Stage-164 result"),
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
    parser.add_argument("--failed-build", type=Path, action="append", default=[])
    args = parser.parse_args()
    try:
        result = compose(
            args.run.resolve(strict=True),
            args.score.resolve(strict=True),
            args.development_root.resolve(strict=True),
            [path.resolve(strict=True) for path in args.failed_build],
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, IndexError, Stage164Error, phase_b.PhaseBError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()

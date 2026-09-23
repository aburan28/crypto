#!/usr/bin/env python3
"""Compose the sealed Stage-167 deterministic parallel fixed-X1 result."""

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
import verify_koblitz_stage166_dense_mul as stage166_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE166 = GATES / "stage-166-dense-product-cancellation-20260923"
STAGE = GATES / "stage-167-parallel-x1-batches-20260923"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
SOURCE_ID = "fb0ea6bb30be3dea841c609a2668939128b4b799d4c0e051b0b52dc6da05f0fa"
COMMIT = "509e43c79bcf36425f4329a09c34fd20c75949ed"
FINGERPRINT = "bd95f10ffd0ab09a2384a2a283fbb85f886df79cd33776b3a3a5a852ccf3bd13"
BACKEND_SHA256 = "58767973f4bc1fcf30d5bbca231d1d5647e6a79f6282c55f1c3874a962820931"
WORD_XORS = 16_821_055_616


class Stage167Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage167Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be an object")
    return value


def only_task(root: Path, context: str) -> dict[str, Any]:
    paths = sorted((root / "tasks").glob("*/task-result.json"))
    require(len(paths) == 1, f"{context} must contain one task")
    task = load(paths[0], f"{context} task")
    require(
        task.get("blind_instance_id") == BLIND_ID
        and task.get("source_instance_id") == SOURCE_ID,
        f"{context} task identity changed",
    )
    return task


def verify_score(root: Path, context: str) -> None:
    seal = load(root / "score-seal.json", f"{context} score seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(root, {"score-seal.json"})
    require(
        seal.get("schema") == native_f4.SCORE_SEAL_SCHEMA
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and load(root / "score.json", context).get("classification_counts")
        == {"true_positive": 1},
        f"{context} score changed",
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
        "total_core_seconds": round(
            math.fsum(row["total_core_seconds"] for row in rows), 12
        ),
        "peak_rss_bytes": max(row["peak_rss_bytes"] for row in rows),
    }


def median(values: list[float]) -> float:
    require(len(values) == 3, "three-run median requires exactly three rows")
    return sorted(values)[1]


def reports_match_parallel_control(parallel: dict[str, Any], single: dict[str, Any]) -> bool:
    parallel = copy.deepcopy(parallel)
    single = copy.deepcopy(single)
    schedule_fields = {
        "solver_schedule",
        "solver_x1_batch_size",
        "solver_x1_batches_completed",
        "solver_x1_speculative_systems_completed",
        "solver_rayon_threads_requested",
        "single_thread_requested",
    }
    for report in (parallel, single):
        report.pop("timing_ns", None)
        for key in schedule_fields:
            report.pop(key, None)
        cost = report["cost"]
        cost.pop("wall_ns", None)
        for key in list(cost["extra"]):
            if key.endswith("_ns"):
                cost["extra"].pop(key)
    return parallel == single


def copy_development(root: Path, output: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    output.mkdir()
    sources = load(root / "sources.json", "Stage-167 source bindings")
    require(
        sources["binary_bindings"]["clean-selected"]["sha256"] == BACKEND_SHA256,
        "selected executable binding changed",
    )
    shutil.copyfile(root / "sources.json", output / "sources.json")
    names = sorted(path.name.removesuffix(".metrics") for path in root.glob("*.metrics"))
    expected = {
        *(f"pattern-tail-{arm}-r{index}" for arm in ("candidate", "control") for index in (1, 2, 3)),
        "pivot-tail-candidate-r1",
        "pivot-tail-control-r1",
        "table-copy-candidate-r1",
        "table-copy-control-r1",
        "x1-stale-batch2-r1",
        *(f"x1-confounded-batch{batch}-r1" for batch in (2, 4, 8, 10, 12, 13, 14)),
        *(f"x1-batch{batch}-r{index}" for batch in (1, 13) for index in (1, 2, 3)),
    }
    require(set(names) == expected, "Stage-167 development arm set changed")
    arms: dict[str, Any] = {}
    resource_rows: list[dict[str, Any]] = []
    for name in names:
        artifacts = {
            suffix: copy_file(
                root / f"{name}.{suffix}",
                output / f"{name}.{suffix}",
                f"Stage-167 {name} {suffix}",
            )
            for suffix in ("metrics", "stdout", "stderr")
        }
        process = load(output / artifacts["metrics"]["path"], f"{name} metrics")
        report = load(output / artifacts["stdout"]["path"], f"{name} stdout")
        metrics = process["metrics"]
        require(
            process.get("returncode") == 0
            and process.get("timed_out") is False
            and process.get("orphan_group_terminated") is False
            and report.get("status") == "sat"
            and report.get("source_witness_valid") is True,
            f"development arm changed: {name}",
        )
        if name.startswith("pattern-tail-"):
            decision = "rejected_pattern_tail_variable_loop_regression"
        elif name.startswith("pivot-tail-"):
            decision = "rejected_pivot_tail_negligible_work_and_wall_regression"
        elif name.startswith("table-copy-"):
            decision = "rejected_table_copy_negligible_work_and_wall_regression"
        elif name == "x1-stale-batch2-r1":
            decision = "rejected_stale_binary_batch_override_absent"
            require("solver_x1_batch_size" not in report, "stale binary unexpectedly accepted batch")
        elif name.startswith("x1-confounded-"):
            decision = "exploratory_width_sweep_rejected_table_copy_confound"
            require(
                report.get("solver_x1_batch_size") in {2, 4, 8, 10, 12, 13, 14}
                and report["cost"]["ops"] != WORD_XORS,
                f"confounded width arm changed: {name}",
            )
        elif name.startswith("x1-batch1-"):
            decision = "same_binary_single_thread_development_control"
            require(
                report.get("solver_x1_batch_size") == 1
                and report.get("solver_rayon_threads_requested") == 1
                and report["cost"]["ops"] == WORD_XORS,
                f"single development control changed: {name}",
            )
        elif name.startswith("x1-batch13-"):
            decision = "parallel_batch13_development_repeat"
            require(
                report.get("solver_x1_batch_size") == 13
                and report.get("solver_rayon_threads_requested") == 13
                and report.get("solver_x1_speculative_systems_completed") == 0
                and report["cost"]["ops"] == WORD_XORS,
                f"parallel development arm changed: {name}",
            )
        else:
            raise Stage167Error(f"unclassified development arm: {name}")
        resource_rows.append(metrics)
        arms[name] = {
            "decision": decision,
            "metrics": metrics,
            "cost": report["cost"],
            "x1_batch_size": report.get("solver_x1_batch_size"),
            "rayon_threads_requested": report.get("solver_rayon_threads_requested"),
            "systems_completed": report.get("fixed_x1_systems_completed"),
            "speculative_systems_completed": report.get(
                "solver_x1_speculative_systems_completed"
            ),
            "equation_fingerprint": report.get("solver_equations_blake3"),
            "witness_points": report.get("witness_points"),
            "artifacts": artifacts,
        }
    manifest = {
        "schema": "koblitz_stage167_development.v1",
        "status": "parallel_width_repeats_controls_and_rejections_retained",
        "source_bindings": sources,
        "arms": arms,
        "complete_development_cost": None,
        "reason_complete_is_null": (
            "incremental compilation, focused tests, profiling, composition, documentation "
            "and Git operations lack complete outer receipts"
        ),
    }
    phase_b.write_json_new(output / "manifest.json", manifest)
    return manifest, resource_rows


def write_markdown(result: dict[str, Any]) -> None:
    parallel = result["selected_parallel_native_f4"]
    single = result["single_thread_control"]
    comparisons = result["comparisons"]
    lines = [
        "# Stage 167: deterministic parallel fixed-X1 F4 batches",
        "",
        result["claim_boundary"],
        "",
        "The algebraic fixed-X1 systems are independent. The backend now gathers a deterministic batch of rational fixed-X1 systems, solves the batch through the bounded Rayon pool, waits for every launched system, charges every result, and then inspects results in the original schedule order. `PQ_F4_X1_BATCH=1` and one Rayon thread retain the prior single-thread path.",
        "",
        f"The clean 13-thread F4 process takes {parallel['metrics']['wall_seconds']:.6f} wall seconds, {parallel['metrics']['total_core_seconds']:.6f} core-seconds, and {parallel['metrics']['peak_rss_bytes']} bytes peak RSS. The same clean binary at batch one takes {single['metrics']['wall_seconds']:.6f} wall / {single['metrics']['total_core_seconds']:.6f} core-seconds / {single['metrics']['peak_rss_bytes']} bytes. Parallel wall improves {comparisons['parallel_over_single']['wall_speedup']:.3f}x while CPU rises {comparisons['parallel_over_single']['core_ratio']:.3f}x and RSS rises {comparisons['parallel_over_single']['rss_ratio']:.3f}x.",
        "",
        "Batch 13 was selected post hoc on this target: the known successful system is the 39th rational system, so three batches cover it without speculative systems. The execution mechanism does not use the witness, truth labels, subgroup enumeration, known scalar, or factor-base logarithms, but this batch-width result is target-specific and is not fresh-target evidence.",
        "",
        f"Direct MITM remains {comparisons['parallel_over_same_host_direct_mitm']['wall_ratio']:.2f}x faster by wall, {comparisons['parallel_over_same_host_direct_mitm']['core_ratio']:.2f}x cheaper by CPU, and {comparisons['parallel_over_same_host_direct_mitm']['rss_ratio']:.2f}x smaller by RSS. The parallel arm does not pass the full-cost gate.",
        "",
        "Licensed Magma F4, the complete same-instance solver panel, fresh-target validation, full index-calculus/rho crossover, and unaffiliated reproduction and novelty review remain open. This is not a Koblitz index-calculus SOTA result.",
        "",
    ]
    phase_b.write_new(STAGE / "RESULTS.md", ("\n".join(lines)).encode())


def compose(
    parallel_run_source: Path,
    parallel_score_source: Path,
    single_run_source: Path,
    single_score_source: Path,
    development_root: Path,
) -> dict[str, Any]:
    require(not STAGE.exists(), f"refusing to overwrite {STAGE}")
    require(stage166_verify.verify().get("status") == "verified", "Stage-166 failed verification")
    baseline = load(STAGE166 / "result.json", "Stage-166 result")
    STAGE.mkdir()
    shutil.copytree(parallel_run_source, STAGE / "selected-parallel-run")
    shutil.copytree(parallel_score_source, STAGE / "selected-parallel-score")
    shutil.copytree(single_run_source, STAGE / "single-thread-run")
    shutil.copytree(single_score_source, STAGE / "single-thread-score")
    development, development_rows = copy_development(
        development_root, STAGE / "development"
    )

    parallel_seal, parallel_plan = native_f4.validate_run_seal(
        STAGE / "selected-parallel-run"
    )
    single_seal, single_plan = native_f4.validate_run_seal(STAGE / "single-thread-run")
    parallel_result = only_task(STAGE / "selected-parallel-run", "parallel")["result"]
    single_result = only_task(STAGE / "single-thread-run", "single")["result"]
    verify_score(STAGE / "selected-parallel-score", "parallel")
    verify_score(STAGE / "single-thread-score", "single")
    parallel_report = parallel_result["backend_report"]
    single_report = single_result["backend_report"]
    parallel_metrics = copy.deepcopy(parallel_result["metrics"])
    parallel_metrics["single_core_seconds"] = None
    single_metrics = single_result["metrics"]
    common = lambda report: (
        report["fixed_x1_masks_visited"] == 85
        and report["fixed_x1_nonrational_skipped"] == 46
        and report["fixed_x1_systems_constructed"] == 39
        and report["fixed_x1_systems_completed"] == 39
        and report["solver_equations_total"] == 2_301
        and report["solver_terms_total"] == 2_295_166
        and report["solver_equations_blake3"] == FINGERPRINT
        and report["cost"]["ops"] == WORD_XORS
        and report["cost"]["extra"]["f4_calls"] == 39
        and report["solver_x1_speculative_systems_completed"] == 0
    )
    require(
        parallel_result["status"] == "sat"
        and parallel_result["source_model_valid"] is True
        and parallel_result["source_witness_valid"] is True
        and parallel_report["solver_x1_batch_size"] == 13
        and parallel_report["solver_rayon_threads_requested"] == 13
        and parallel_report["single_thread_requested"] is False
        and parallel_report["solver_x1_batches_completed"] == 3
        and common(parallel_report),
        "selected parallel result changed",
    )
    require(
        single_result["status"] == "sat"
        and single_result["source_model_valid"] is True
        and single_result["source_witness_valid"] is True
        and single_report["solver_x1_batch_size"] == 1
        and single_report["solver_rayon_threads_requested"] == 1
        and single_report["single_thread_requested"] is True
        and single_report["solver_x1_batches_completed"] == 39
        and common(single_report)
        and reports_match_parallel_control(parallel_report, single_report),
        "same-binary single-thread control changed",
    )
    require(
        parallel_report["witness_points"] == single_report["witness_points"]
        == baseline["selected_native_f4"]["report"]["witness_points"],
        "Stage-167 changed the exact witness",
    )
    require(
        parallel_plan["source_revision"]
        == single_plan["source_revision"]
        == {"commit": COMMIT, "dirty": False, "porcelain": []}
        and parallel_plan["x1_batch_size"] == 13
        and parallel_plan["rayon_threads_requested"] == 13
        and parallel_plan["single_thread_requested"] is False
        and single_plan["x1_batch_size"] == 1
        and single_plan["rayon_threads_requested"] == 1
        and single_plan["single_thread_requested"] is True,
        "Stage-167 execution plans changed",
    )

    receipt_path = STAGE / "selected-parallel-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path, "selected receipt")
    build_tool.validate_receipt(
        receipt_path, receipt["binaries"], receipt["rust_source_objects"]
    )
    require(
        receipt["source_commit"] == COMMIT
        and receipt["source_clean"] is True
        and receipt["binaries"]["backend"]["sha256"] == BACKEND_SHA256,
        "selected build changed",
    )
    parallel_summary = load(
        STAGE / "selected-parallel-run/run-summary.json", "parallel summary"
    )
    single_summary = load(STAGE / "single-thread-run/run-summary.json", "single summary")
    require(
        parallel_summary["f4_process_resources"]["single_core_seconds"] is None
        and math.isclose(
            single_summary["f4_process_resources"]["single_core_seconds"],
            single_metrics["total_core_seconds"],
            rel_tol=0,
            abs_tol=1e-9,
        ),
        "single-core accounting boundary changed",
    )

    candidate_development = [
        development["arms"][f"x1-batch13-r{index}"] for index in (1, 2, 3)
    ]
    control_development = [
        development["arms"][f"x1-batch1-r{index}"] for index in (1, 2, 3)
    ]
    candidate_wall = median([row["metrics"]["wall_seconds"] for row in candidate_development])
    control_wall = median([row["metrics"]["wall_seconds"] for row in control_development])
    candidate_core = median(
        [row["metrics"]["total_core_seconds"] for row in candidate_development]
    )
    control_core = median(
        [row["metrics"]["total_core_seconds"] for row in control_development]
    )
    candidate_rss = median(
        [row["metrics"]["peak_rss_bytes"] for row in candidate_development]
    )
    control_rss = median(
        [row["metrics"]["peak_rss_bytes"] for row in control_development]
    )

    table = copy.deepcopy(baseline["same_target_table"])
    single_row = next(row for row in table if row["backend"] == "native-f4-selected")
    single_row.update(
        {
            "formulation": "fixed-X1 direct S4; sparse schedule; dense-cancel F4; one thread",
            "wall_seconds": single_metrics["wall_seconds"],
            "total_core_seconds": single_metrics["total_core_seconds"],
            "peak_rss_bytes": single_metrics["peak_rss_bytes"],
        }
    )
    parallel_row = copy.deepcopy(single_row)
    parallel_row.update(
        {
            "backend": "native-f4-parallel-x1-batch13",
            "formulation": "same 39 fixed-X1 systems in three deterministic Rayon batches",
            "wall_seconds": parallel_metrics["wall_seconds"],
            "total_core_seconds": parallel_metrics["total_core_seconds"],
            "peak_rss_bytes": parallel_metrics["peak_rss_bytes"],
            "single_core_seconds": None,
            "requested_parallel_threads": 13,
        }
    )
    table.insert(table.index(single_row) + 1, parallel_row)
    mitm = next(row for row in table if row["backend"] == "direct-mitm")
    comparisons = {
        "parallel_over_single": {
            "wall_speedup": single_metrics["wall_seconds"] / parallel_metrics["wall_seconds"],
            "wall_ratio": parallel_metrics["wall_seconds"] / single_metrics["wall_seconds"],
            "core_ratio": parallel_metrics["total_core_seconds"] / single_metrics["total_core_seconds"],
            "rss_ratio": parallel_metrics["peak_rss_bytes"] / single_metrics["peak_rss_bytes"],
            "same_scientific_report_outside_schedule_and_timing_fields": True,
            "same_verified_witness": True,
        },
        "three_run_medians": {
            "parallel_wall_seconds": candidate_wall,
            "single_wall_seconds": control_wall,
            "wall_speedup": control_wall / candidate_wall,
            "parallel_core_seconds": candidate_core,
            "single_core_seconds": control_core,
            "core_ratio": candidate_core / control_core,
            "parallel_rss_bytes": candidate_rss,
            "single_rss_bytes": control_rss,
            "rss_ratio": candidate_rss / control_rss,
            "parallel_efficiency_vs_13_threads": (control_wall / candidate_wall) / 13,
        },
        "parallel_over_stage166_single": {
            "wall_speedup": baseline["selected_native_f4"]["metrics"]["wall_seconds"]
            / parallel_metrics["wall_seconds"],
            "core_ratio": parallel_metrics["total_core_seconds"]
            / baseline["selected_native_f4"]["metrics"]["total_core_seconds"],
            "rss_ratio": parallel_metrics["peak_rss_bytes"]
            / baseline["selected_native_f4"]["metrics"]["peak_rss_bytes"],
        },
        "parallel_over_same_host_direct_mitm": {
            "wall_ratio": parallel_metrics["wall_seconds"] / mitm["wall_seconds"],
            "core_ratio": parallel_metrics["total_core_seconds"] / mitm["total_core_seconds"],
            "rss_ratio": parallel_metrics["peak_rss_bytes"] / mitm["peak_rss_bytes"],
        },
        "single_over_same_host_direct_mitm": {
            "wall_ratio": single_metrics["wall_seconds"] / mitm["wall_seconds"],
            "core_ratio": single_metrics["total_core_seconds"] / mitm["total_core_seconds"],
            "rss_ratio": single_metrics["peak_rss_bytes"] / mitm["peak_rss_bytes"],
        },
        "cold_build_plus_parallel_run": {
            "wall_seconds": receipt["resources"]["summed_process_wall_seconds"]
            + parallel_summary["outer_resources"]["wall_seconds"],
            "total_core_seconds": receipt["resources"]["total_core_seconds"]
            + parallel_summary["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": max(
                receipt["resources"]["largest_child_peak_rss_bytes"],
                parallel_metrics["peak_rss_bytes"],
            ),
        },
    }
    require(
        comparisons["parallel_over_single"]["wall_speedup"] > 3.5
        and comparisons["three_run_medians"]["wall_speedup"] > 3.5
        and comparisons["parallel_over_same_host_direct_mitm"]["wall_ratio"] > 1.8,
        "parallel comparison boundary changed",
    )

    prior = baseline["campaign_accounting"]["measured_lower_bound"]
    incremental_rows = list(development_rows)
    incremental_rows.extend(row["metrics"] for row in receipt["build_processes"])
    incremental_rows.extend(
        [
            {
                "wall_seconds": parallel_summary["outer_resources"]["wall_seconds"],
                "total_core_seconds": parallel_summary["outer_resources"][
                    "total_core_seconds"
                ],
                "peak_rss_bytes": parallel_metrics["peak_rss_bytes"],
            },
            {
                "wall_seconds": single_summary["outer_resources"]["wall_seconds"],
                "total_core_seconds": single_summary["outer_resources"][
                    "total_core_seconds"
                ],
                "peak_rss_bytes": single_metrics["peak_rss_bytes"],
            },
        ]
    )
    incremental = sum_resources(incremental_rows)
    campaign = {
        "stage166_measured_lower_bound": prior,
        "stage167_incremental_measured": incremental,
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
            "incremental compilation, focused tests and local profiling",
            "score composition, documentation and Git operations",
        ],
        "complete_campaign_cost": None,
        "reason_complete_is_null": "some development commands lack complete outer receipts",
    }
    result = {
        "schema": "koblitz_stage167_parallel_x1_batches.v1",
        "status": "complete_parallel_and_single_true_positive",
        "date": "2026-09-23",
        "claim_boundary": (
            "Strong target-specific parallel engineering on one finite public toy PDP. Batch "
            "width 13 was selected post hoc because this target's first valid relation occurs "
            "in rational system 39. This is not fresh-target evidence, a full solver panel, a "
            "full index-calculus/rho crossover, independent reproduction, novelty review, or SOTA."
        ),
        "instance": baseline["instance"],
        "factor_base": baseline["factor_base"],
        "mechanism": {
            "name": "deterministic_parallel_fixed_x1_batches",
            "definition": (
                "Gather the next B rational fixed-X1 systems in the algebraic schedule, solve "
                "all B through a bounded Rayon pool, charge every launched result, then inspect "
                "the completed results in schedule order"
            ),
            "selected_batch_size": 13,
            "selected_rayon_threads": 13,
            "same_binary_control": "PQ_F4_X1_BATCH=1 and RAYON_NUM_THREADS=1",
            "uses_target_labels": False,
            "uses_discrete_log_labels": False,
            "uses_known_witness": False,
            "selection_caveat": (
                "batch width selected post hoc on this target; 39 rational systems divide into "
                "three batches of 13 without speculative systems"
            ),
        },
        "selected_parallel_native_f4": {
            "source_revision": parallel_plan["source_revision"],
            "run_inventory_sha256": parallel_seal["inventory_sha256"],
            "status": parallel_result["status"],
            "classification": "true_positive",
            "metrics": parallel_metrics,
            "report": parallel_report,
        },
        "single_thread_control": {
            "source_revision": single_plan["source_revision"],
            "run_inventory_sha256": single_seal["inventory_sha256"],
            "status": single_result["status"],
            "classification": "true_positive",
            "metrics": single_metrics,
            "report": single_report,
        },
        "same_target_table": table,
        "comparisons": comparisons,
        "development": development,
        "predecessor": {
            "stage166_result_sha256": phase_b.sha256_file(
                STAGE166 / "result.json", "Stage-166 result"
            ),
            "stage166_seal_sha256": phase_b.sha256_file(
                STAGE166 / "result-seal.json", "Stage-166 seal"
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
        "schema": "koblitz_stage167_parallel_x1_batches_seal.v1",
        "status": "stage167_parallel_x1_result_frozen",
        "result_sha256": phase_b.sha256_file(STAGE / "result.json", "Stage-167 result"),
        "inventory": inventory,
        "inventory_sha256": phase_b.canonical_sha256(inventory),
    }
    seal["seal_payload_sha256"] = phase_b.canonical_sha256(seal)
    phase_b.write_json_new(STAGE / "result-seal.json", seal)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parallel-run", type=Path, required=True)
    parser.add_argument("--parallel-score", type=Path, required=True)
    parser.add_argument("--single-run", type=Path, required=True)
    parser.add_argument("--single-score", type=Path, required=True)
    parser.add_argument("--development-root", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = compose(
            args.parallel_run.resolve(strict=True),
            args.parallel_score.resolve(strict=True),
            args.single_run.resolve(strict=True),
            args.single_score.resolve(strict=True),
            args.development_root.resolve(strict=True),
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (
        OSError,
        ValueError,
        KeyError,
        IndexError,
        Stage167Error,
        phase_b.PhaseBError,
    ) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()

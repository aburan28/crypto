#!/usr/bin/env python3
"""Compose the sealed Stage-168 one-batch fixed-X1 parallel result."""

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
import verify_koblitz_stage167_parallel_x1 as stage167_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE167 = GATES / "stage-167-parallel-x1-batches-20260923"
STAGE = GATES / "stage-168-single-batch-39-20260923"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
SOURCE_ID = "fb0ea6bb30be3dea841c609a2668939128b4b799d4c0e051b0b52dc6da05f0fa"
COMMIT = "509e43c79bcf36425f4329a09c34fd20c75949ed"
FINGERPRINT = "bd95f10ffd0ab09a2384a2a283fbb85f886df79cd33776b3a3a5a852ccf3bd13"
BACKEND_SHA256 = "58767973f4bc1fcf30d5bbca231d1d5647e6a79f6282c55f1c3874a962820931"
WORD_XORS = 16_821_055_616


class Stage168Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage168Error(message)


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
        "total_core_seconds": round(
            math.fsum(row["total_core_seconds"] for row in rows), 12
        ),
        "peak_rss_bytes": max(row["peak_rss_bytes"] for row in rows),
    }


def median(values: list[float]) -> float:
    require(len(values) == 3, "three-run median requires exactly three rows")
    return sorted(values)[1]


def reports_match(selected: dict[str, Any], reference: dict[str, Any]) -> bool:
    selected = copy.deepcopy(selected)
    reference = copy.deepcopy(reference)
    schedule_fields = {
        "solver_schedule",
        "solver_x1_batch_size",
        "solver_x1_batches_completed",
        "solver_x1_speculative_systems_completed",
        "solver_rayon_threads_requested",
        "single_thread_requested",
    }
    for report in (selected, reference):
        report.pop("timing_ns", None)
        for key in schedule_fields:
            report.pop(key, None)
        report["cost"].pop("wall_ns", None)
        for key in list(report["cost"]["extra"]):
            if key.endswith("_ns"):
                report["cost"]["extra"].pop(key)
    return selected == reference


def copy_development(root: Path, output: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    output.mkdir()
    sources = load(root / "sources.json", "development sources")
    require(
        sources["binary"]["sha256"] == BACKEND_SHA256
        and sources["binary"]["source_commit"] == COMMIT,
        "development source binding changed",
    )
    shutil.copyfile(root / "sources.json", output / "sources.json")
    names = sorted(path.name.removesuffix(".metrics") for path in root.glob("*.metrics"))
    expected = {
        "batch39-t10-r1",
        "batch39-t11-r1",
        "batch39-t12-r1",
        "batch39-t12-r2",
        "batch39-t12-r3",
        "batch39-t13-r1",
        "batch39-t14-r1",
    }
    require(set(names) == expected, "Stage-168 development arm set changed")
    arms: dict[str, Any] = {}
    resources: list[dict[str, Any]] = []
    for name in names:
        artifacts = {
            suffix: copy_file(
                root / f"{name}.{suffix}",
                output / f"{name}.{suffix}",
                f"Stage-168 {name} {suffix}",
            )
            for suffix in ("metrics", "stdout", "stderr")
        }
        process = load(output / artifacts["metrics"]["path"], f"{name} metrics")
        report = load(output / artifacts["stdout"]["path"], f"{name} stdout")
        metrics = process["metrics"]
        threads = int(name.split("-t", 1)[1].split("-", 1)[0])
        require(
            process.get("returncode") == 0
            and process.get("timed_out") is False
            and process.get("orphan_group_terminated") is False
            and report.get("status") == "sat"
            and report.get("source_witness_valid") is True
            and report.get("solver_x1_batch_size") == 39
            and report.get("solver_rayon_threads_requested") == threads
            and report.get("solver_x1_batches_completed") == 1
            and report.get("solver_x1_speculative_systems_completed") == 0
            and report.get("fixed_x1_systems_completed") == 39
            and report.get("solver_equations_blake3") == FINGERPRINT
            and report["cost"]["ops"] == WORD_XORS,
            f"development arm changed: {name}",
        )
        decision = (
            "selected_thread_count_repeat"
            if name.startswith("batch39-t12-")
            else "thread_count_sweep_control"
        )
        resources.append(metrics)
        arms[name] = {
            "decision": decision,
            "metrics": metrics,
            "valid_single_core_seconds": None,
            "threads": threads,
            "batch_size": 39,
            "cost": report["cost"],
            "equation_fingerprint": report["solver_equations_blake3"],
            "witness_points": report["witness_points"],
            "artifacts": artifacts,
        }
    manifest = {
        "schema": "koblitz_stage168_development.v1",
        "status": "one_batch_thread_sweep_retained",
        "source_binding": sources,
        "arms": arms,
        "selection_caveat": sources["selection_caveat"],
        "complete_development_cost": None,
        "reason_complete_is_null": (
            "worktree creation, score composition, documentation and Git operations lack "
            "complete outer receipts"
        ),
    }
    phase_b.write_json_new(output / "manifest.json", manifest)
    return manifest, resources


def write_markdown(result: dict[str, Any]) -> None:
    selected = result["selected_parallel_native_f4"]
    comparisons = result["comparisons"]
    lines = [
        "# Stage 168: one deterministic batch of 39 fixed-X1 systems",
        "",
        result["claim_boundary"],
        "",
        "Stage 168 removes the two barriers in Stage 167 by placing the same 39 rational fixed-X1 systems into one deterministic batch. Twelve Rayon workers dynamically execute that complete batch; every system is charged before results are inspected in algebraic schedule order.",
        "",
        f"The clean F4 process takes {selected['metrics']['wall_seconds']:.6f} wall seconds, {selected['metrics']['total_core_seconds']:.6f} total core-seconds, and {selected['metrics']['peak_rss_bytes']} bytes peak RSS. `single_core_seconds` is null. It completes the same 39 systems, charges 16,821,055,616 word XORs, and returns the same exact witness as the Stage-167 parallel and single-thread arms.",
        "",
        f"Relative to Stage 167 batch 13, clean wall improves {comparisons['selected_over_stage167_parallel']['wall_speedup']:.3f}x, CPU falls to {comparisons['selected_over_stage167_parallel']['core_ratio']:.3f}x, and RSS falls to {comparisons['selected_over_stage167_parallel']['rss_ratio']:.3f}x. Three 12-thread development runs have a {comparisons['three_run_median']['wall_seconds']:.6f}-second median.",
        "",
        "Batch size 39 is selected post hoc because it is exactly the successful prefix on this target. This is target-specific scheduling evidence and not a valid generic stopping rule or fresh-target result.",
        "",
        f"Direct MITM remains {comparisons['selected_over_same_host_direct_mitm']['wall_ratio']:.2f}x faster by wall, {comparisons['selected_over_same_host_direct_mitm']['core_ratio']:.2f}x cheaper by CPU, and {comparisons['selected_over_same_host_direct_mitm']['rss_ratio']:.2f}x smaller by RSS. The full-cost gate remains false.",
        "",
        "Licensed Magma F4, the complete same-instance solver panel, fresh-target validation, full index-calculus/rho crossover, and unaffiliated reproduction and novelty review remain open. This is not a Koblitz index-calculus SOTA result.",
        "",
    ]
    phase_b.write_new(STAGE / "RESULTS.md", ("\n".join(lines)).encode())


def compose(run_source: Path, score_source: Path, development_root: Path) -> dict[str, Any]:
    require(not STAGE.exists(), f"refusing to overwrite {STAGE}")
    require(stage167_verify.verify().get("status") == "verified", "Stage-167 failed verification")
    baseline = load(STAGE167 / "result.json", "Stage-167 result")
    STAGE.mkdir()
    shutil.copytree(run_source, STAGE / "selected-run")
    shutil.copytree(score_source, STAGE / "selected-score")
    development, development_rows = copy_development(
        development_root, STAGE / "development"
    )
    run_seal, plan = native_f4.validate_run_seal(STAGE / "selected-run")
    selected_result = only_task(STAGE / "selected-run")["result"]
    verify_score(STAGE / "selected-score")
    report = selected_result["backend_report"]
    raw_metrics = selected_result["metrics"]
    selected_metrics = copy.deepcopy(raw_metrics)
    selected_metrics["single_core_seconds"] = None
    require(
        selected_result["status"] == "sat"
        and selected_result["source_model_valid"] is True
        and selected_result["source_witness_valid"] is True
        and plan["source_revision"] == {"commit": COMMIT, "dirty": False, "porcelain": []}
        and plan["x1_batch_size"] == 39
        and plan["rayon_threads_requested"] == 12
        and plan["single_thread_requested"] is False
        and report["solver_x1_batch_size"] == 39
        and report["solver_rayon_threads_requested"] == 12
        and report["solver_x1_batches_completed"] == 1
        and report["solver_x1_speculative_systems_completed"] == 0
        and report["fixed_x1_masks_visited"] == 85
        and report["fixed_x1_systems_completed"] == 39
        and report["solver_equations_blake3"] == FINGERPRINT
        and report["cost"]["ops"] == WORD_XORS
        and report["cost"]["extra"]["f4_calls"] == 39,
        "selected one-batch result changed",
    )
    stage167_parallel = baseline["selected_parallel_native_f4"]
    stage167_single = baseline["single_thread_control"]
    require(
        report["witness_points"] == stage167_parallel["report"]["witness_points"]
        and reports_match(report, stage167_parallel["report"])
        and reports_match(report, stage167_single["report"]),
        "Stage-168 changed scientific output",
    )
    summary = load(STAGE / "selected-run/run-summary.json", "selected summary")
    require(
        summary["f4_process_resources"]["single_core_seconds"] is None,
        "parallel single-core field must be null",
    )
    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path, "build receipt")
    build_tool.validate_receipt(
        receipt_path, receipt["binaries"], receipt["rust_source_objects"]
    )
    require(
        receipt["source_commit"] == COMMIT
        and receipt["source_clean"] is True
        and receipt["binaries"]["backend"]["sha256"] == BACKEND_SHA256,
        "selected build changed",
    )
    repeats = [development["arms"][f"batch39-t12-r{index}"] for index in (1, 2, 3)]
    median_wall = median([row["metrics"]["wall_seconds"] for row in repeats])
    median_core = median([row["metrics"]["total_core_seconds"] for row in repeats])
    median_rss = median([row["metrics"]["peak_rss_bytes"] for row in repeats])
    table = copy.deepcopy(baseline["same_target_table"])
    previous_parallel = next(
        row for row in table if row["backend"] == "native-f4-parallel-x1-batch13"
    )
    selected_row = copy.deepcopy(previous_parallel)
    selected_row.update(
        {
            "backend": "native-f4-parallel-x1-batch39-t12",
            "formulation": "same 39 fixed-X1 systems in one deterministic batch on 12 threads",
            "wall_seconds": selected_metrics["wall_seconds"],
            "total_core_seconds": selected_metrics["total_core_seconds"],
            "peak_rss_bytes": selected_metrics["peak_rss_bytes"],
            "requested_parallel_threads": 12,
        }
    )
    table.insert(table.index(previous_parallel) + 1, selected_row)
    mitm = next(row for row in table if row["backend"] == "direct-mitm")
    comparisons = {
        "selected_over_stage167_parallel": {
            "wall_speedup": stage167_parallel["metrics"]["wall_seconds"]
            / selected_metrics["wall_seconds"],
            "core_ratio": selected_metrics["total_core_seconds"]
            / stage167_parallel["metrics"]["total_core_seconds"],
            "rss_ratio": selected_metrics["peak_rss_bytes"]
            / stage167_parallel["metrics"]["peak_rss_bytes"],
            "same_scientific_report": True,
        },
        "selected_over_single": {
            "wall_speedup": stage167_single["metrics"]["wall_seconds"]
            / selected_metrics["wall_seconds"],
            "core_ratio": selected_metrics["total_core_seconds"]
            / stage167_single["metrics"]["total_core_seconds"],
            "rss_ratio": selected_metrics["peak_rss_bytes"]
            / stage167_single["metrics"]["peak_rss_bytes"],
        },
        "three_run_median": {
            "wall_seconds": median_wall,
            "total_core_seconds": median_core,
            "peak_rss_bytes": median_rss,
        },
        "selected_over_same_host_direct_mitm": {
            "wall_ratio": selected_metrics["wall_seconds"] / mitm["wall_seconds"],
            "core_ratio": selected_metrics["total_core_seconds"] / mitm["total_core_seconds"],
            "rss_ratio": selected_metrics["peak_rss_bytes"] / mitm["peak_rss_bytes"],
        },
        "cold_build_plus_selected_run": {
            "wall_seconds": receipt["resources"]["summed_process_wall_seconds"]
            + summary["outer_resources"]["wall_seconds"],
            "total_core_seconds": receipt["resources"]["total_core_seconds"]
            + summary["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": max(
                receipt["resources"]["largest_child_peak_rss_bytes"],
                selected_metrics["peak_rss_bytes"],
            ),
            "build_already_charged_in_stage167_campaign": True,
        },
    }
    require(
        comparisons["selected_over_stage167_parallel"]["wall_speedup"] > 1.1
        and comparisons["selected_over_same_host_direct_mitm"]["wall_ratio"] > 1.6,
        "Stage-168 comparison boundary changed",
    )
    prior = baseline["campaign_accounting"]["measured_lower_bound"]
    incremental_rows = list(development_rows)
    incremental_rows.append(
        {
            "wall_seconds": summary["outer_resources"]["wall_seconds"],
            "total_core_seconds": summary["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": selected_metrics["peak_rss_bytes"],
        }
    )
    incremental = sum_resources(incremental_rows)
    campaign = {
        "stage167_measured_lower_bound": prior,
        "stage168_incremental_measured": incremental,
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
        "unmetered_development_commands_present": True,
        "unmetered_scope": [
            "detached worktree creation and failed preflight before process launch",
            "score composition, documentation and Git operations",
        ],
        "complete_campaign_cost": None,
        "reason_complete_is_null": "some development commands lack complete outer receipts",
    }
    result = {
        "schema": "koblitz_stage168_single_batch_39.v1",
        "status": "complete_one_batch_true_positive",
        "date": "2026-09-23",
        "claim_boundary": (
            "Strong target-specific parallel scheduling on one finite public toy PDP. Batch "
            "size 39 is exactly the known successful prefix and was selected post hoc. This "
            "is not a generic stopping rule, fresh-target evidence, a full solver panel, a "
            "full index-calculus/rho crossover, independent reproduction, novelty review, or SOTA."
        ),
        "instance": baseline["instance"],
        "factor_base": baseline["factor_base"],
        "mechanism": {
            "name": "one_complete_prefix_batch",
            "selected_batch_size": 39,
            "selected_rayon_threads": 12,
            "uses_target_labels": False,
            "uses_discrete_log_labels": False,
            "uses_known_witness": False,
            "selection_caveat": (
                "batch size equals the already known successful 39-system prefix; selected "
                "post hoc and invalid as a generic stopping rule"
            ),
        },
        "selected_parallel_native_f4": {
            "source_revision": plan["source_revision"],
            "run_inventory_sha256": run_seal["inventory_sha256"],
            "status": selected_result["status"],
            "classification": "true_positive",
            "metrics": selected_metrics,
            "report": report,
        },
        "same_target_table": table,
        "comparisons": comparisons,
        "development": development,
        "predecessor": {
            "stage167_result_sha256": phase_b.sha256_file(
                STAGE167 / "result.json", "Stage-167 result"
            ),
            "stage167_seal_sha256": phase_b.sha256_file(
                STAGE167 / "result-seal.json", "Stage-167 seal"
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
        "schema": "koblitz_stage168_single_batch_39_seal.v1",
        "status": "stage168_single_batch_result_frozen",
        "result_sha256": phase_b.sha256_file(STAGE / "result.json", "Stage-168 result"),
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
    except (
        OSError,
        ValueError,
        KeyError,
        IndexError,
        Stage168Error,
        phase_b.PhaseBError,
    ) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()

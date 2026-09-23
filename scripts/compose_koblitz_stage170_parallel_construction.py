#!/usr/bin/env python3
"""Compose the sealed Stage-170 parallel fixed-X1 construction result."""

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
import verify_koblitz_stage169_fresh_holdout as stage169_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE169 = GATES / "stage-169-fresh-target-holdout-20260923"
STAGE = GATES / "stage-170-parallel-fixed-x1-construction-20260923"
BLIND_ID = "b-421e22a9c1c3b9d56396c8bbd0e46185bbebc32de0306f2bee6e9585703c2be4"
SOURCE_SYSTEM_ID = "x-af85283c5728becd64890c6eaa654ccb164cb499fcc4cde7b5cfd1484b4e0efa"
SOURCE_INSTANCE_ID = "954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7"
COMMIT = "080d3c2bfe9fefaf86bd10c8279a7fb9ebe000a0"
BACKEND_SHA256 = "5515eb70f24c195172707dfaa0dce99bd80bb5e21958747af6b8221597ec8b80"
FINAL_DEVELOPMENT_SHA256 = "d691a28e858d7c8223fc6b393470ac0697c69fb7a7df7967284221fd17dfcab0"
INTERMEDIATE_PLAN_SHA256 = "b19ef110d5d4caf4992f6db06d1d09273babb078f11f4d6057d6eedef455d7d9"
STAGE169_BACKEND_SHA256 = "58767973f4bc1fcf30d5bbca231d1d5647e6a79f6282c55f1c3874a962820931"
FINGERPRINT = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
WORD_XORS = 99_199_976_264


class Stage170Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage170Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be an object")
    return value


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copy_file(source: Path, destination: Path, context: str) -> dict[str, Any]:
    require(source.is_file(), f"missing {context}: {source}")
    shutil.copyfile(source, destination)
    return {
        "path": destination.name,
        "bytes": destination.stat().st_size,
        "sha256": sha256(destination),
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


def only_task(root: Path) -> dict[str, Any]:
    paths = sorted((root / "tasks").glob("*/task-result.json"))
    require(len(paths) == 1, "selected run must contain one task")
    task = load(paths[0], "selected task")
    require(
        task.get("blind_instance_id") == BLIND_ID
        and task.get("source_system_id") == SOURCE_SYSTEM_ID
        and task.get("source_instance_id") == SOURCE_INSTANCE_ID,
        "selected task identity changed",
    )
    return task


def verify_score(root: Path) -> dict[str, Any]:
    seal = load(root / "score-seal.json", "selected score seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    inventory = phase_b.all_regular_inventory(root, {"score-seal.json"})
    score = load(root / "score.json", "selected score")
    require(
        seal.get("schema") == native_f4.SCORE_SEAL_SCHEMA
        and claimed == phase_b.canonical_sha256(payload)
        and inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256")
        and score.get("classification_counts") == {"true_negative": 1}
        and score["rows"][0]["blind_instance_id"] == BLIND_ID
        and score["rows"][0]["target_class"] == "nondecomposable"
        and score["rows"][0]["solver_status"] == "unsat",
        "selected score changed",
    )
    return score


def validate_report(
    report: dict[str, Any], context: str, expected_word_xors: int | None = WORD_XORS
) -> None:
    require(
        report.get("status") == "unsat"
        and report.get("exhaustive") is True
        and report.get("fixed_x1_masks_visited") == 512
        and report.get("fixed_x1_nonrational_skipped") == 270
        and report.get("fixed_x1_systems_constructed") == 242
        and report.get("fixed_x1_systems_completed") == 242
        and report.get("solver_equations_total") == 14_278
        and report.get("solver_terms_total") == 14_515_915
        and report.get("solver_equations_blake3") == FINGERPRINT
        and report.get("algebraic_roots") == 0
        and report.get("witness_points") is None
        and (
            report["cost"]["ops"] == expected_word_xors
            if expected_word_xors is not None
            else report["cost"]["ops"] > 0
        )
        and report["cost"]["extra"]["f4_calls"] == 242,
        f"{context} scientific result changed",
    )


def reports_match(left: dict[str, Any], right: dict[str, Any]) -> bool:
    left = copy.deepcopy(left)
    right = copy.deepcopy(right)
    schedule = {
        "solver_schedule",
        "solver_parallel_construction",
        "solver_parallel_construction_control",
        "solver_x1_batch_size",
        "solver_x1_batches_completed",
        "solver_rayon_threads_requested",
        "single_thread_requested",
    }
    for report in (left, right):
        report.pop("timing_ns", None)
        for key in schedule:
            report.pop(key, None)
        report["cost"].pop("wall_ns", None)
        for key in list(report["cost"]["extra"]):
            if key.endswith("_ns"):
                report["cost"]["extra"].pop(key)
    return left == right


def copy_arm(
    source: Path,
    output: Path,
    name: str,
    decision: str,
    binary_sha256: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    artifacts = {
        suffix: copy_file(
            source / f"native-f4.{suffix}",
            output / f"{name}.{suffix}",
            f"Stage-170 {name} {suffix}",
        )
        for suffix in ("metrics", "stdout", "stderr")
    }
    process = load(output / f"{name}.metrics", f"Stage-170 {name} metrics")
    report = load(output / f"{name}.stdout", f"Stage-170 {name} report")
    require(
        process.get("returncode") == 0
        and process.get("timed_out") is False
        and process.get("orphan_group_terminated") is False,
        f"Stage-170 {name} process changed",
    )
    validate_report(
        report,
        name,
        None if name in {"preparallel-m4ri-b7", "preparallel-m4ri-b9"} else WORD_XORS,
    )
    metrics = copy.deepcopy(process["metrics"])
    metrics["single_core_seconds"] = None
    return (
        {
            "decision": decision,
            "binary_sha256": binary_sha256,
            "metrics": metrics,
            "valid_single_core_seconds": None,
            "batch_size": report["solver_x1_batch_size"],
            "rayon_threads": report["solver_rayon_threads_requested"],
            "parallel_construction": report.get("solver_parallel_construction"),
            "construction_wall_seconds": report["timing_ns"][
                "fixed_x1_s4_construction"
            ]
            / 1e9,
            "cost": report["cost"],
            "equation_fingerprint": report["solver_equations_blake3"],
            "artifacts": artifacts,
        },
        process["metrics"],
    )


def copy_development(root: Path, output: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    output.mkdir()
    arms: dict[str, Any] = {}
    resources: list[dict[str, Any]] = []

    mappings = [
        ("preparallel-batch64", root / "koblitz-stage170-batch64-dev", "rejected_small_barrier_gain", STAGE169_BACKEND_SHA256),
        ("preparallel-m4ri-b7", root / "koblitz-stage170-m4ri-b7-dev", "rejected_m4ri_width_slower", STAGE169_BACKEND_SHA256),
        ("preparallel-m4ri-b9", root / "koblitz-stage170-m4ri-b9-dev", "rejected_m4ri_width_slower", STAGE169_BACKEND_SHA256),
        ("preparallel-t14", root / "koblitz-stage170-t14-dev", "retained_thread_width_probe", STAGE169_BACKEND_SHA256),
        ("preparallel-profile", root / "koblitz-stage170-profile-dev", "retained_metered_profile_subject", STAGE169_BACKEND_SHA256),
        ("plan-control", root / "koblitz-stage170-shared-plan-dev/control", "rejected_intermediate_plan_control", INTERMEDIATE_PLAN_SHA256),
        ("plan-serial", root / "koblitz-stage170-shared-plan-dev/plan-serial", "rejected_shared_plan_no_serial_gain", INTERMEDIATE_PLAN_SHA256),
        ("plan-parallel", root / "koblitz-stage170-shared-plan-dev/plan-parallel", "rejected_shared_plan_slower_than_reference_parallel", INTERMEDIATE_PLAN_SHA256),
        ("reference-parallel", root / "koblitz-stage170-shared-plan-dev/reference-parallel", "retained_mechanism_isolation", INTERMEDIATE_PLAN_SHA256),
    ]
    for prefix in ("control", "candidate"):
        for index in range(1, 4):
            mappings.append(
                (
                    f"{prefix}-r{index}",
                    root / f"koblitz-stage170-parallel-construction-dev/{prefix}-r{index}",
                    "same_binary_serial_repeat"
                    if prefix == "control"
                    else "same_binary_parallel_repeat",
                    FINAL_DEVELOPMENT_SHA256,
                )
            )
    mappings.extend(
        [
            (
                "posthoc-t14-b39",
                root / "koblitz-stage170-parallel-construction-dev/candidate-t14-b39",
                "posthoc_thread_width_probe",
                FINAL_DEVELOPMENT_SHA256,
            ),
            (
                "posthoc-t14-b64",
                root / "koblitz-stage170-parallel-construction-dev/candidate-t14-b64",
                "posthoc_target_tuned_ceiling",
                FINAL_DEVELOPMENT_SHA256,
            ),
        ]
    )
    for name, source, decision, binary in mappings:
        arm, metrics = copy_arm(source, output, name, decision, binary)
        arms[name] = arm
        resources.append(metrics)

    sample_root = root / "koblitz-stage170-sample-dev"
    sample_artifacts = {
        name: copy_file(
            sample_root / name,
            output / f"profile-{name}",
            f"Stage-170 profiler {name}",
        )
        for name in ("sample.txt", "native-f4.stdout", "native-f4.stderr", "exit-status")
    }
    require(
        (output / "profile-exit-status").read_text().strip() == "0"
        and "Sort by top of stack" in (output / "profile-sample.txt").read_text()
        and "pq_f4_f2::multiply_for_f4" in (output / "profile-sample.txt").read_text(),
        "Stage-170 profiler artifacts changed",
    )
    manifest = {
        "schema": "koblitz_stage170_development.v1",
        "status": "parallel_construction_controls_repeats_and_rejections_retained",
        "arms": arms,
        "profiler": {
            "decision": "read_only_stack_sample_used_to_select_next_mechanism",
            "metered_subject_arm": "preparallel-profile",
            "sample_process_resources": None,
            "artifacts": sample_artifacts,
        },
        "binary_bindings": {
            "stage169_clean": STAGE169_BACKEND_SHA256,
            "intermediate_shared_plan": INTERMEDIATE_PLAN_SHA256,
            "final_incremental": FINAL_DEVELOPMENT_SHA256,
        },
        "complete_development_cost": None,
        "reason_complete_is_null": (
            "incremental compilation, tests, stack sampling, composition, documentation and Git "
            "operations lack complete outer receipts"
        ),
    }
    phase_b.write_json_new(output / "manifest.json", manifest)
    return manifest, resources


def copy_failed_build(root: Path, output: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    evidence = root / "evidence"
    require(evidence.is_dir(), "failed build evidence is missing")
    shutil.copytree(evidence, output)
    roles = [
        "archive-source",
        "extract-source",
        "install-lock",
        "version-cargo",
        "version-rustc",
        "vendor",
    ]
    rows = []
    for role in roles:
        process = load(output / f"{role}.metrics.json", f"failed build {role}")
        rows.append(process["metrics"])
        require(process.get("timed_out") is False, f"failed build {role} timed out")
        if role == "vendor":
            require(process.get("returncode") == 101, "failed vendor terminal changed")
        else:
            require(process.get("returncode") == 0, f"failed build precursor changed: {role}")
    require(
        "failed to download `redis v0.32.7`" in (output / "vendor.stderr").read_text()
        and "--offline was specified" in (output / "vendor.stderr").read_text(),
        "failed build reason changed",
    )
    manifest = {
        "schema": "koblitz_stage170_rejected_build.v1",
        "status": "rejected_missing_dependency_in_default_offline_cache",
        "source_commit": COMMIT,
        "roles": roles,
        "resource_summary": sum_resources(rows),
        "reason": "default offline Cargo cache lacked redis 0.32.7; no network fallback occurred",
    }
    phase_b.write_json_new(output / "manifest.json", manifest)
    return manifest, rows


def median(values: list[float]) -> float:
    require(len(values) == 3, "three-run median requires exactly three values")
    return sorted(values)[1]


def write_markdown(result: dict[str, Any]) -> None:
    selected = result["selected_native_f4"]
    control = result["same_binary_serial_control"]
    comparison = result["comparisons"]
    lines = [
        "# Stage 170: parallel fixed-X1 S4 construction",
        "",
        result["claim_boundary"],
        "",
        "Each rational fixed-X1 S4 system in a deterministic batch is now constructed through the same bounded Rayon pool used by the repository-native F4 solver. Indexed parallel collection preserves algebraic mask order, and `PQ_F4_DISABLE_PARALLEL_CONSTRUCTION=1` restores serial construction in the same binary.",
        "",
        f"The clean selected batch-39/twelve-thread process completes the known Stage-169 true-negative target in {selected['metrics']['wall_seconds']:.6f} wall seconds / {selected['metrics']['total_core_seconds']:.6f} total core-seconds / {selected['metrics']['peak_rss_bytes']} bytes peak RSS. The same binary with serial construction takes {control['metrics']['wall_seconds']:.6f} seconds. This is a {comparison['clean_same_binary']['wall_speedup']:.3f}x wall speedup; construction falls from {control['construction_wall_seconds']:.6f} to {selected['construction_wall_seconds']:.6f} seconds.",
        "",
        f"Three alternating development pairs have a {comparison['three_run_medians']['wall_speedup']:.3f}x median wall speedup. Every arm preserves the 14,278 equations, 14,515,915 terms, equation fingerprint `{FINGERPRINT}`, 242 F4 calls, {WORD_XORS:,} word XORs, zero roots, and exhaustive UNSAT classification.",
        "",
        f"An explicitly post-hoc 14-thread/batch-64 clean arm reaches {result['posthoc_single_target_ceiling']['metrics']['wall_seconds']:.6f} seconds but raises peak RSS and can do extra speculative work on decomposable targets. It is retained as a target-tuned ceiling rather than the selected generic policy.",
        "",
        f"Direct MITM remains {comparison['selected_over_direct_mitm']['wall_ratio']:.2f}x faster by wall and {comparison['selected_over_direct_mitm']['core_ratio']:.2f}x cheaper by CPU on this target. The build plus selected run costs {comparison['cold_build_plus_selected']['wall_seconds']:.6f} wall seconds, so the full-cost gate remains false.",
        "",
        "This optimization was developed after the Stage-169 truth was opened. It strengthens implementation performance but is not a second fresh holdout, a full target distribution, a complete index-calculus/rho crossover, independent reproduction, novelty review, or SOTA.",
        "",
    ]
    phase_b.write_new(STAGE / "RESULTS.md", ("\n".join(lines)).encode())


def compose(
    selected_run: Path,
    selected_score: Path,
    clean_control_root: Path,
    clean_posthoc_root: Path,
    development_root: Path,
    failed_build_root: Path,
) -> dict[str, Any]:
    require(not STAGE.exists(), f"refusing to overwrite {STAGE}")
    require(stage169_verify.verify().get("status") == "verified", "Stage-169 failed verification")
    baseline = load(STAGE169 / "result.json", "Stage-169 result")
    STAGE.mkdir()
    shutil.copytree(selected_run, STAGE / "selected-run")
    shutil.copytree(selected_score, STAGE / "selected-score")
    clean = STAGE / "clean-controls"
    clean.mkdir()
    selected_control, control_row = copy_arm(
        clean_control_root,
        clean,
        "serial-construction",
        "same_clean_binary_serial_construction_control",
        BACKEND_SHA256,
    )
    posthoc, posthoc_row = copy_arm(
        clean_posthoc_root,
        clean,
        "posthoc-t14-b64",
        "posthoc_single_target_ceiling_not_selected_policy",
        BACKEND_SHA256,
    )
    development, development_rows = copy_development(
        development_root, STAGE / "development"
    )
    rejected_build, failed_build_rows = copy_failed_build(
        failed_build_root, STAGE / "rejected-build-default-cache"
    )

    run_seal, plan = native_f4.validate_run_seal(STAGE / "selected-run")
    selected_result = only_task(STAGE / "selected-run")["result"]
    verify_score(STAGE / "selected-score")
    report = selected_result["backend_report"]
    validate_report(report, "selected clean run")
    raw_metrics = selected_result["metrics"]
    selected_metrics = copy.deepcopy(raw_metrics)
    selected_metrics["single_core_seconds"] = None
    require(
        plan["source_revision"] == {"commit": COMMIT, "dirty": False, "porcelain": []}
        and plan["selected_blind_instance_ids"] == [BLIND_ID]
        and plan["x1_batch_size"] == 39
        and plan["rayon_threads_requested"] == 12
        and report["solver_parallel_construction"] is True
        and report["solver_x1_batch_size"] == 39
        and report["solver_rayon_threads_requested"] == 12
        and report["solver_x1_batches_completed"] == 7,
        "selected construction policy changed",
    )
    control_report = load(clean / "serial-construction.stdout", "clean serial report")
    posthoc_report = load(clean / "posthoc-t14-b64.stdout", "clean posthoc report")
    require(
        control_report["solver_parallel_construction"] is False
        and control_report["solver_x1_batch_size"] == 39
        and control_report["solver_rayon_threads_requested"] == 12
        and posthoc_report["solver_parallel_construction"] is True
        and posthoc_report["solver_x1_batch_size"] == 64
        and posthoc_report["solver_rayon_threads_requested"] == 14
        and reports_match(report, control_report)
        and reports_match(report, posthoc_report),
        "clean control changed scientific output",
    )

    summary = load(STAGE / "selected-run/run-summary.json", "selected summary")
    require(
        summary["f4_process_resources"]["single_core_seconds"] is None,
        "parallel selected run invented single-core time",
    )
    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path, "selected build receipt")
    build_tool.validate_receipt(
        receipt_path, receipt["binaries"], receipt["rust_source_objects"]
    )
    require(
        receipt["source_commit"] == COMMIT
        and receipt["source_clean"] is True
        and receipt["binaries"]["backend"]["sha256"] == BACKEND_SHA256,
        "selected build changed",
    )

    candidate_repeats = [
        development["arms"][f"candidate-r{index}"] for index in (1, 2, 3)
    ]
    control_repeats = [
        development["arms"][f"control-r{index}"] for index in (1, 2, 3)
    ]
    candidate_median = {
        key: median([row["metrics"][key] for row in candidate_repeats])
        for key in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
    }
    control_median = {
        key: median([row["metrics"][key] for row in control_repeats])
        for key in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
    }
    direct = next(row for row in baseline["same_target_table"] if row["backend"] == "direct-mitm")
    comparisons = {
        "clean_same_binary": {
            "wall_speedup": selected_control["metrics"]["wall_seconds"]
            / selected_metrics["wall_seconds"],
            "wall_reduction_fraction": 1
            - selected_metrics["wall_seconds"] / selected_control["metrics"]["wall_seconds"],
            "construction_speedup": selected_control["construction_wall_seconds"]
            / (report["timing_ns"]["fixed_x1_s4_construction"] / 1e9),
            "core_ratio": selected_metrics["total_core_seconds"]
            / selected_control["metrics"]["total_core_seconds"],
            "rss_ratio": selected_metrics["peak_rss_bytes"]
            / selected_control["metrics"]["peak_rss_bytes"],
            "same_scientific_report": True,
        },
        "three_run_medians": {
            "candidate": candidate_median,
            "control": control_median,
            "wall_speedup": control_median["wall_seconds"]
            / candidate_median["wall_seconds"],
        },
        "selected_over_stage169_preregistered_baseline": {
            "wall_speedup": baseline["fresh_holdout_native_f4"]["metrics"]["wall_seconds"]
            / selected_metrics["wall_seconds"],
            "core_ratio": selected_metrics["total_core_seconds"]
            / baseline["fresh_holdout_native_f4"]["metrics"]["total_core_seconds"],
            "rss_ratio": selected_metrics["peak_rss_bytes"]
            / baseline["fresh_holdout_native_f4"]["metrics"]["peak_rss_bytes"],
        },
        "selected_over_direct_mitm": {
            "wall_ratio": selected_metrics["wall_seconds"] / direct["wall_seconds"],
            "core_ratio": selected_metrics["total_core_seconds"]
            / direct["total_core_seconds"],
            "rss_ratio": selected_metrics["peak_rss_bytes"] / direct["peak_rss_bytes"],
        },
        "posthoc_over_direct_mitm": {
            "wall_ratio": posthoc["metrics"]["wall_seconds"] / direct["wall_seconds"],
            "core_ratio": posthoc["metrics"]["total_core_seconds"]
            / direct["total_core_seconds"],
            "rss_ratio": posthoc["metrics"]["peak_rss_bytes"] / direct["peak_rss_bytes"],
        },
        "cold_build_plus_selected": {
            "wall_seconds": receipt["resources"]["summed_process_wall_seconds"]
            + summary["outer_resources"]["wall_seconds"],
            "total_core_seconds": receipt["resources"]["total_core_seconds"]
            + summary["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": max(
                receipt["resources"]["largest_child_peak_rss_bytes"],
                selected_metrics["peak_rss_bytes"],
            ),
        },
    }
    require(
        comparisons["clean_same_binary"]["wall_speedup"] > 1.2
        and comparisons["three_run_medians"]["wall_speedup"] > 1.2
        and comparisons["selected_over_direct_mitm"]["wall_ratio"] > 5.7,
        "Stage-170 comparison boundary changed",
    )

    prior = baseline["campaign_accounting"]["measured_lower_bound"]
    incremental_rows = list(development_rows)
    incremental_rows.extend([control_row, posthoc_row])
    incremental_rows.extend(failed_build_rows)
    incremental_rows.extend(process["metrics"] for process in receipt["build_processes"])
    incremental_rows.append(
        {
            "wall_seconds": summary["outer_resources"]["wall_seconds"],
            "total_core_seconds": summary["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": selected_metrics["peak_rss_bytes"],
        }
    )
    incremental = sum_resources(incremental_rows)
    campaign = {
        "stage169_measured_lower_bound": prior,
        "stage170_incremental_measured": incremental,
        "measured_lower_bound": {
            "resource_components": prior["resource_components"]
            + incremental["resource_components"],
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
            "incremental compilation and focused tests",
            "read-only stack sampling process resources",
            "score composition, documentation and Git operations",
        ],
        "reason_complete_is_null": "some development and orchestration commands lack complete outer receipts",
    }

    result = {
        "schema": "koblitz_stage170_parallel_fixed_x1_construction.v1",
        "status": "complete_parallel_construction_true_negative",
        "date": "2026-09-23",
        "claim_boundary": (
            "Same-target engineering after the Stage-169 truth was opened. It establishes an "
            "exact deterministic construction-stage speedup for the repository-native F4 arm, "
            "not a fresh-target distribution, full index calculus, independent reproduction, "
            "novelty, or SOTA."
        ),
        "instance": baseline["instance"],
        "factor_base": baseline["factor_base"],
        "mechanism": {
            "name": "ordered_parallel_fixed_x1_construction",
            "target_independent": True,
            "uses_truth_or_witness": False,
            "preserves_batch_and_mask_order": True,
            "serial_control": "PQ_F4_DISABLE_PARALLEL_CONSTRUCTION=1",
            "selected_batch_size": 39,
            "selected_rayon_threads": 12,
        },
        "selected_native_f4": {
            "source_revision": plan["source_revision"],
            "run_inventory_sha256": run_seal["inventory_sha256"],
            "status": selected_result["status"],
            "classification": "true_negative",
            "metrics": selected_metrics,
            "construction_wall_seconds": report["timing_ns"][
                "fixed_x1_s4_construction"
            ]
            / 1e9,
            "report": report,
        },
        "same_binary_serial_control": selected_control,
        "posthoc_single_target_ceiling": posthoc,
        "comparisons": comparisons,
        "development": development,
        "rejected_build": rejected_build,
        "predecessor": {
            "stage169_result_sha256": sha256(STAGE169 / "result.json"),
            "stage169_seal_sha256": sha256(STAGE169 / "result-seal.json"),
            "stage169_preregistered_holdout_preserved": True,
        },
        "campaign_accounting": campaign,
        "gates": baseline["gates"],
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
        "schema": "koblitz_stage170_parallel_fixed_x1_construction_seal.v1",
        "status": "stage170_parallel_construction_result_frozen",
        "result_sha256": sha256(STAGE / "result.json"),
        "inventory": inventory,
        "inventory_sha256": phase_b.canonical_sha256(inventory),
    }
    seal["seal_payload_sha256"] = phase_b.canonical_sha256(seal)
    phase_b.write_json_new(STAGE / "result-seal.json", seal)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-run", type=Path, required=True)
    parser.add_argument("--selected-score", type=Path, required=True)
    parser.add_argument("--clean-control", type=Path, required=True)
    parser.add_argument("--clean-posthoc", type=Path, required=True)
    parser.add_argument("--development-root", type=Path, required=True)
    parser.add_argument("--failed-build", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = compose(
            args.selected_run.resolve(strict=True),
            args.selected_score.resolve(strict=True),
            args.clean_control.resolve(strict=True),
            args.clean_posthoc.resolve(strict=True),
            args.development_root.resolve(strict=True),
            args.failed_build.resolve(strict=True),
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (
        OSError,
        ValueError,
        KeyError,
        IndexError,
        Stage170Error,
        phase_b.PhaseBError,
    ) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()

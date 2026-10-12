#!/usr/bin/env python3
"""Compose the sealed Stage-166 dense F4 monomial-product result."""

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
import verify_koblitz_stage165_x1_order as stage165_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE165 = GATES / "stage-165-x1-weight-order-20260923"
STAGE = GATES / "stage-166-dense-product-cancellation-20260923"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
SOURCE_ID = "fb0ea6bb30be3dea841c609a2668939128b4b799d4c0e051b0b52dc6da05f0fa"
COMMIT = "ac3e385cec9c90b489017f6012ebe0b50ca06405"
FINGERPRINT = "bd95f10ffd0ab09a2384a2a283fbb85f886df79cd33776b3a3a5a852ccf3bd13"
BACKEND_SHA256 = "013d3525557a7e3c811f11b1b189d52aaafc782512c35e6c3a5f98947652ba00"
DEVELOPMENT_BINARY_SHA256 = "cc7afb5e856c48dc175249eaf446484256b4f78f09137a0524bbb29957755a92"
STAGE165_BACKEND_SHA256 = "9d2bbc2a0fc567dc771c91e9c1bf7d03473be55784fa1f3f300f03eaf5c46f8f"


class Stage166Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage166Error(message)


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


def equivalent_reports(selected: dict[str, Any], control: dict[str, Any]) -> bool:
    selected = copy.deepcopy(selected)
    control = copy.deepcopy(control)
    for report in (selected, control):
        report.pop("timing_ns", None)
        report.pop("solver_description", None)
        report.pop("solver_monomial_multiply", None)
        report.pop("solver_monomial_multiply_control", None)
        cost = report["cost"]
        cost.pop("wall_ns", None)
        extra = cost["extra"]
        for key in list(extra):
            if key.endswith("_ns") or key.startswith("dense_mul_"):
                extra.pop(key)
    return selected == control


def copy_development(root: Path, output: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    output.mkdir()
    source_bindings = load(root / "sources.json", "development source bindings")
    require(
        source_bindings.get("binary_bindings", {})
        .get("dense-final-pairs", {})
        .get("sha256")
        == DEVELOPMENT_BINARY_SHA256
        and source_bindings["binary_bindings"]["m4ri-block2-through10-subset"]["sha256"]
        == STAGE165_BACKEND_SHA256
        and source_bindings["binary_bindings"]["clean-control"]["sha256"]
        == BACKEND_SHA256,
        "development executable bindings changed",
    )
    shutil.copyfile(root / "sources.json", output / "sources.json")
    names = sorted(path.name.removesuffix(".metrics") for path in root.glob("*.metrics"))
    expected = {
        *(f"m4ri-block{block}" for block in (6, 7, 8, 9, 10)),
        "dense-v1-candidate-r1",
        "dense-v1-control-r1",
        *(f"dense-candidate-r{index}" for index in (1, 2, 3)),
        *(f"dense-control-r{index}" for index in (1, 2, 3)),
        "clean-control",
    }
    require(set(names) == expected, "Stage-166 development arm set changed")
    arms: dict[str, Any] = {}
    resource_rows: list[dict[str, Any]] = []
    for name in names:
        artifacts = {
            suffix: copy_file(
                root / f"{name}.{suffix}",
                output / f"{name}.{suffix}",
                f"Stage-166 {name} {suffix}",
            )
            for suffix in ("metrics", "stdout", "stderr")
        }
        process = load(output / artifacts["metrics"]["path"], f"{name} metrics")
        metrics = process["metrics"]
        report = load(output / artifacts["stdout"]["path"], f"{name} stdout")
        require(
            process.get("returncode") == 0
            and process.get("timed_out") is False
            and process.get("orphan_group_terminated") is False
            and report.get("status") == "sat"
            and report.get("source_witness_valid") is True
            and report.get("solver_equations_blake3") == FINGERPRINT,
            f"development arm changed: {name}",
        )
        if name.startswith("m4ri-block"):
            block = int(name.removeprefix("m4ri-block"))
            decision = "retained_block_width_sweep_no_new_default"
            require(
                report["cost"]["extra"]["m4ri_block_width_max"] == block,
                f"M4RI block width changed: {name}",
            )
        elif name == "dense-v1-candidate-r1":
            decision = "rejected_two_array_dense_parity_wall_regression"
            require(report["cost"]["extra"]["dense_mul_calls"] == 818_171, "bad v1 candidate")
        elif name == "dense-v1-control-r1":
            decision = "rejected_v1_same_binary_control"
            require(report["cost"]["extra"]["dense_mul_calls"] == 0, "bad v1 control")
        elif name.startswith("dense-candidate-"):
            decision = "encoded_key_dense_product_development_repeat"
            require(
                report["cost"]["extra"]["dense_mul_calls"] == 818_171
                and report["cost"]["extra"]["dense_mul_cancelled_terms"] == 101_593_414,
                f"bad dense candidate: {name}",
            )
        elif name.startswith("dense-control-") or name == "clean-control":
            decision = "same_binary_sorted_product_control"
            require(report["cost"]["extra"]["dense_mul_calls"] == 0, f"bad control: {name}")
        else:
            raise Stage166Error(f"unclassified development arm: {name}")
        resource_rows.append(metrics)
        arms[name] = {
            "decision": decision,
            "metrics": metrics,
            "cost": report["cost"],
            "solver_monomial_multiply": report.get("solver_monomial_multiply"),
            "equation_fingerprint": report["solver_equations_blake3"],
            "witness_points": report["witness_points"],
            "artifacts": artifacts,
        }
    manifest = {
        "schema": "koblitz_stage166_development.v1",
        "status": "dense_product_repeats_controls_and_rejections_retained",
        "source_bindings": source_bindings,
        "arms": arms,
        "complete_development_cost": None,
        "reason_complete_is_null": (
            "incremental compilation, tests, profiling, composition, documentation and Git "
            "operations lack complete outer receipts"
        ),
    }
    phase_b.write_json_new(output / "manifest.json", manifest)
    return manifest, resource_rows


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
        "schema": "koblitz_stage166_rejected_build.v1",
        "status": "rejected_missing_dependency_in_default_offline_cache",
        "source_commit": COMMIT,
        "roles": roles,
        "resource_summary": sum_resources(rows),
        "reason": "default offline Cargo cache lacked redis 0.32.7; no network fallback occurred",
    }
    phase_b.write_json_new(output / "manifest.json", manifest)
    return manifest, rows


def write_markdown(result: dict[str, Any]) -> None:
    metrics = result["selected_native_f4"]["metrics"]
    comparisons = result["comparisons"]
    extra = result["selected_native_f4"]["report"]["cost"]["extra"]
    lines = [
        "# Stage 166: dense cancellation before F4 product sorting",
        "",
        result["claim_boundary"],
        "",
        "For Boolean systems with at most 20 variables, monomial multiplication now toggles mapped output masks in one generation-tagged dense array. Each first-seen mask stores a 32-bit key encoding descending degree and ascending mask order. The solver sorts only odd-parity survivors; `PQ_F4_DISABLE_DENSE_MUL=1` restores the prior map-all, sort-all, then-cancel path.",
        "",
        f"The clean blind F4 process takes {metrics['wall_seconds']:.6f} wall seconds, {metrics['total_core_seconds']:.6f} core-seconds, and {metrics['peak_rss_bytes']} bytes peak RSS. Across {extra['dense_mul_calls']:,} products it maps {extra['dense_mul_input_terms']:,} input terms to {extra['dense_mul_output_terms']:,} survivors and cancels {extra['dense_mul_cancelled_terms']:,} terms before sorting.",
        "",
        f"The exact clean disabled control takes {comparisons['same_binary_sorted_product_control']['wall_seconds']:.6f} wall seconds. The selected clean speedup is {comparisons['same_binary_sorted_product_control']['selected_wall_speedup']:.3f}x; three interleaved development pairs give {comparisons['three_run_medians']['wall_speedup']:.3f}x. Equation fingerprint, all algebraic and matrix work counters outside the new product counters, roots, and the exact curve-verified witness agree.",
        "",
        f"Direct MITM remains {comparisons['selected_over_same_host_direct_mitm']['wall_ratio']:.2f}x faster by wall, {comparisons['selected_over_same_host_direct_mitm']['core_ratio']:.2f}x cheaper by CPU, and {comparisons['selected_over_same_host_direct_mitm']['rss_ratio']:.2f}x smaller by RSS.",
        "",
        "Licensed Magma F4, the complete same-instance solver panel, full index-calculus/rho crossover, and unaffiliated reproduction and novelty review remain open. This is not a Koblitz index-calculus SOTA result.",
        "",
    ]
    phase_b.write_new(STAGE / "RESULTS.md", ("\n".join(lines)).encode())


def compose(
    run_source: Path,
    score_source: Path,
    development_root: Path,
    failed_build_root: Path,
) -> dict[str, Any]:
    require(not STAGE.exists(), f"refusing to overwrite {STAGE}")
    require(stage165_verify.verify().get("status") == "verified", "Stage-165 failed verification")
    baseline = load(STAGE165 / "result.json", "Stage-165 result")
    STAGE.mkdir()
    shutil.copytree(run_source, STAGE / "selected-run")
    shutil.copytree(score_source, STAGE / "selected-score")
    development, development_rows = copy_development(
        development_root, STAGE / "development"
    )
    rejected_build, rejected_build_rows = copy_failed_build(
        failed_build_root, STAGE / "rejected-build-default-cache"
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
        and report["solver_x1_order"] == "hamming_weight_then_mask"
        and report["solver_monomial_multiply"]
        == "epoch_dense_parity_cancel_then_encoded_u32_order_sort"
        and report["solver_monomial_multiply_control"] == "PQ_F4_DISABLE_DENSE_MUL=1"
        and report["fixed_x1_masks_visited"] == 85
        and report["fixed_x1_nonrational_skipped"] == 46
        and report["fixed_x1_systems_constructed"] == 39
        and report["fixed_x1_systems_completed"] == 39
        and report["solver_equations_total"] == 2_301
        and report["solver_terms_total"] == 2_295_166
        and report["solver_equations_blake3"] == FINGERPRINT
        and report["cost"]["ops"] == 16_821_055_616
        and extra["f4_calls"] == 39
        and extra["dense_mul_calls"] == 818_171
        and extra["dense_mul_input_terms"] == 651_502_431
        and extra["dense_mul_output_terms"] == 549_909_017
        and extra["dense_mul_cancelled_terms"] == 101_593_414
        and extra["dense_mul_scratch_bytes_max"] == 1_067_584,
        "selected dense-product result changed",
    )
    old = baseline["selected_native_f4"]
    require(
        report["witness_points"] == old["report"]["witness_points"],
        "Stage-166 changed the exact witness",
    )

    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path, "selected receipt")
    build_tool.validate_receipt(
        receipt_path, receipt["binaries"], receipt["rust_source_objects"]
    )
    require(
        receipt["source_commit"] == COMMIT
        and receipt["source_clean"] is True
        and receipt["binaries"]["backend"]["sha256"] == BACKEND_SHA256,
        "selected source or binary changed",
    )
    run_summary = load(STAGE / "selected-run/run-summary.json", "run summary")
    control = development["arms"]["clean-control"]
    control_report = load(
        STAGE / "development/clean-control.stdout", "clean control report"
    )
    require(
        equivalent_reports(report, control_report)
        and control_report["solver_monomial_multiply"] == "sort_all_mapped_terms_then_cancel"
        and control["witness_points"] == report["witness_points"],
        "clean same-binary control changed scientific output",
    )

    candidates = [development["arms"][f"dense-candidate-r{i}"] for i in (1, 2, 3)]
    controls = [development["arms"][f"dense-control-r{i}"] for i in (1, 2, 3)]
    candidate_wall = median([row["metrics"]["wall_seconds"] for row in candidates])
    control_wall = median([row["metrics"]["wall_seconds"] for row in controls])
    candidate_core = median([row["metrics"]["total_core_seconds"] for row in candidates])
    control_core = median([row["metrics"]["total_core_seconds"] for row in controls])

    table = copy.deepcopy(baseline["same_target_table"])
    f4 = next(row for row in table if row["backend"] == "native-f4-selected")
    f4.update(
        {
            "formulation": "fixed-X1 direct S4; sparse schedule; dense-cancel product F4",
            "wall_seconds": metrics["wall_seconds"],
            "total_core_seconds": metrics["total_core_seconds"],
            "peak_rss_bytes": metrics["peak_rss_bytes"],
            "operations": {
                "word_xors_including_m4ri_tables": report["cost"]["ops"],
                "f4_calls": extra["f4_calls"],
                "dense_mul_calls": extra["dense_mul_calls"],
                "dense_mul_cancelled_terms": extra["dense_mul_cancelled_terms"],
            },
        }
    )
    mitm = next(row for row in table if row["backend"] == "direct-mitm")
    comparisons = {
        "same_binary_sorted_product_control": {
            "wall_seconds": control["metrics"]["wall_seconds"],
            "total_core_seconds": control["metrics"]["total_core_seconds"],
            "peak_rss_bytes": control["metrics"]["peak_rss_bytes"],
            "selected_wall_ratio": metrics["wall_seconds"] / control["metrics"]["wall_seconds"],
            "selected_wall_speedup": control["metrics"]["wall_seconds"] / metrics["wall_seconds"],
            "selected_core_ratio": metrics["total_core_seconds"] / control["metrics"]["total_core_seconds"],
            "selected_rss_ratio": metrics["peak_rss_bytes"] / control["metrics"]["peak_rss_bytes"],
            "same_scientific_report_outside_declared_product_and_timing_fields": True,
            "same_verified_witness": True,
        },
        "three_run_medians": {
            "candidate_wall_seconds": candidate_wall,
            "control_wall_seconds": control_wall,
            "wall_speedup": control_wall / candidate_wall,
            "candidate_core_seconds": candidate_core,
            "control_core_seconds": control_core,
            "core_speedup": control_core / candidate_core,
        },
        "selected_over_stage165": {
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
    require(
        comparisons["same_binary_sorted_product_control"]["selected_wall_speedup"] > 1.08
        and comparisons["three_run_medians"]["wall_speedup"] > 1.06,
        "dense product speedup threshold changed",
    )

    prior = baseline["campaign_accounting"]["measured_lower_bound"]
    incremental_rows = list(development_rows)
    incremental_rows.extend(rejected_build_rows)
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
        "stage165_measured_lower_bound": prior,
        "stage166_incremental_measured": incremental,
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
        "schema": "koblitz_stage166_dense_product_cancellation.v1",
        "status": "complete_dense_product_true_positive",
        "date": "2026-09-23",
        "claim_boundary": (
            "Strong single-target internal engineering on one finite public toy PDP. This is "
            "not a full solver panel, full index-calculus/rho crossover, independent "
            "reproduction, novelty review, or SOTA."
        ),
        "instance": baseline["instance"],
        "factor_base": baseline["factor_base"],
        "mechanism": {
            "name": "epoch_dense_parity_product_cancellation",
            "definition": (
                "For n_vars <= 20, toggle mapped monomial masks in one generation-tagged "
                "u32 domain array; store one encoded degree/mask key per first-seen output; "
                "sort only odd-parity survivors"
            ),
            "exact_reference": "F2BoolPoly::mul_mono map-all sort-all duplicate cancellation",
            "same_binary_control": "PQ_F4_DISABLE_DENSE_MUL=1",
            "uses_target_labels": False,
            "uses_discrete_log_labels": False,
            "uses_known_witness": False,
            "equivalence_test": "dense_monomial_multiply_matches_sorted_reference",
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
        "rejected_build": rejected_build,
        "predecessor": {
            "stage165_result_sha256": phase_b.sha256_file(
                STAGE165 / "result.json", "Stage-165 result"
            ),
            "stage165_seal_sha256": phase_b.sha256_file(
                STAGE165 / "result-seal.json", "Stage-165 seal"
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
        "schema": "koblitz_stage166_dense_product_cancellation_seal.v1",
        "status": "stage166_dense_product_result_frozen",
        "result_sha256": phase_b.sha256_file(STAGE / "result.json", "Stage-166 result"),
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
    parser.add_argument("--failed-build-root", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = compose(
            args.run.resolve(strict=True),
            args.score.resolve(strict=True),
            args.development_root.resolve(strict=True),
            args.failed_build_root.resolve(strict=True),
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (
        OSError,
        ValueError,
        KeyError,
        IndexError,
        Stage166Error,
        phase_b.PhaseBError,
    ) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()

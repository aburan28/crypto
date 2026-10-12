#!/usr/bin/env python3
"""Compose the sealed Stage-160 fixed-X1 S4 construction result."""

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
import verify_koblitz_stage159_native_f4 as stage159_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE159 = GATES / "stage-159-native-f4-single-target-20260922"
STAGE = GATES / "stage-160-fixed-x1-constant-specialisation-20260922"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
SOURCE_ID = "fb0ea6bb30be3dea841c609a2668939128b4b799d4c0e051b0b52dc6da05f0fa"
CELL = "n59-l9-m3-standard-a1-f0"
RESULT_SCHEMA = "koblitz_stage160_fixed_x1_constant_specialisation.v1"
SEAL_SCHEMA = "koblitz_stage160_fixed_x1_constant_specialisation_seal.v1"

DEVELOPMENT_ARMS = (
    (
        "ordered-matrix",
        "stage159-sorted-f4",
        "rejected",
        "Sort complete F4 matrices by leading monomial and sparse tie-break; lower RSS but slower wall and more XOR work.",
    ),
    (
        "trim-row-end",
        "stage159-trim-end-f4",
        "rejected",
        "Trim zero suffix words after row XOR; negligible work reduction and slower wall.",
    ),
    (
        "hashed-reducer-index",
        "stage159-subset-index-f4",
        "rejected",
        "Enumerate monomial divisors through a hash index; probes fell sharply but wall did not improve.",
    ),
    (
        "dense-reducer-index",
        "stage159-dense-index-f4",
        "rejected",
        "Replace the reducer hash with a bounded dense index; internal F4 improved slightly but whole-target wall regressed.",
    ),
    (
        "fixed-x1-specialised-provisional",
        "stage159-fixed-x1-specialized-f4",
        "promoted_after_clean_rebuild",
        "Apply fixed and target field elements as linear maps and form X2*X3 once; provisional dirty-tree timing only.",
    ),
)


class Stage160Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage160Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be an object")
    return value


def only_task(root: Path) -> dict[str, Any]:
    paths = sorted((root / "tasks").glob("*/task-result.json"))
    require(len(paths) == 1, "selected run must contain exactly one task")
    task = load(paths[0], "selected task")
    require(
        task.get("blind_instance_id") == BLIND_ID
        and task.get("source_instance_id") == SOURCE_ID
        and task.get("cell_id") == CELL,
        "selected task identity changed",
    )
    return task


def verify_score(root: Path) -> dict[str, Any]:
    seal = load(root / "score-seal.json", "selected score seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == native_f4.SCORE_SEAL_SCHEMA
        and claimed == phase_b.canonical_sha256(payload),
        "selected score seal changed",
    )
    inventory = phase_b.all_regular_inventory(root, {"score-seal.json"})
    require(
        inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256"),
        "selected score inventory changed",
    )
    score = load(root / "score.json", "selected score")
    require(
        score.get("classification_counts") == {"true_positive": 1}
        and score.get("rows", [{}])[0].get("blind_instance_id") == BLIND_ID,
        "selected score classification changed",
    )
    return score


def copy_development(root: Path, output: Path) -> dict[str, Any]:
    output.mkdir()
    arms = []
    for name, prefix, decision, description in DEVELOPMENT_ARMS:
        copied = {}
        for suffix in ("metrics", "stdout", "stderr"):
            source = root / f"{prefix}.{suffix}"
            require(source.is_file(), f"missing development receipt: {source}")
            destination = output / f"{name}.{suffix}"
            shutil.copyfile(source, destination)
            copied[suffix] = {
                "path": destination.name,
                "sha256": phase_b.sha256_file(destination, f"{name} {suffix}"),
                "bytes": destination.stat().st_size,
            }
        metrics = load(output / copied["metrics"]["path"], f"{name} metrics")["metrics"]
        report = load(output / copied["stdout"]["path"], f"{name} stdout")
        require(
            report.get("status") == "sat"
            and report.get("source_model_valid") is True
            and report.get("source_witness_valid") is True
            and report.get("solver_equations_blake3")
            == "d6969c5f89b5be6a7eacbdac45799fc734a273c9034949a482723b4422504ae7",
            f"{name} exact terminal changed",
        )
        arms.append(
            {
                "name": name,
                "decision": decision,
                "description": description,
                "source_patch_retained": False,
                "source_boundary": "dirty incremental development binary; useful only as a charged exploratory receipt",
                "metrics": metrics,
                "fixed_x1_s4_construction_ns": report["timing_ns"]["fixed_x1_s4_construction"],
                "f4_wall_ns": report["cost"]["wall_ns"],
                "word_xors": report["cost"]["ops"],
                "divisor_tests": report["cost"]["extra"]["divisor_tests"],
                "peak_matrix_bytes": report["cost"]["peak_bytes"],
                "artifacts": copied,
            }
        )
    manifest = {
        "schema": "koblitz_stage160_development_arms.v1",
        "status": "metered_development_arms_retained",
        "complete_development_cost": None,
        "reason_complete_is_null": "incremental compiles and focused test invocations lacked an outer meter",
        "arms": arms,
    }
    phase_b.write_json_new(output / "manifest.json", manifest)
    return manifest


def add_resources(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "resource_components": len(rows),
        "summed_wall_seconds": round(math.fsum(row["wall_seconds"] for row in rows), 12),
        "total_core_seconds": round(math.fsum(row["total_core_seconds"] for row in rows), 12),
        "peak_rss_bytes": max(row["peak_rss_bytes"] for row in rows),
    }


def write_markdown(result: dict[str, Any]) -> None:
    selected = result["selected_native_f4"]
    comparisons = result["comparisons"]
    rows = result["same_target_table"]
    lines = [
        "# Stage 160: constant-linear fixed-X1 S4 construction",
        "",
        result["claim_boundary"],
        "",
        "The selected constructor forms the shared symbolic `X2*X3` product once and applies the fixed `X1` and target powers as exact linear maps in the polynomial basis. A direct unit gate compares its complete Boolean polynomial vectors with the generic S4 expansion over three fields, five fixed values, and four targets per field.",
        "",
        "| Backend | Host | Terminal | Wall s | Core s | Peak RSS B | Conflicts / native work |",
        "|:--|:--|:--|--:|--:|--:|:--|",
    ]
    for row in rows:
        work = row["operations"] if row["operations"] is not None else row["conflicts"]
        lines.append(
            "| {backend} | {host} | {status} | {wall} | {core} | {rss} | {work} |".format(
                backend=row["backend"],
                host=row["host"] or "not run",
                status=row["status"],
                wall="null" if row["wall_seconds"] is None else f"{row['wall_seconds']:.6f}",
                core="null" if row["total_core_seconds"] is None else f"{row['total_core_seconds']:.6f}",
                rss="null" if row["peak_rss_bytes"] is None else row["peak_rss_bytes"],
                work="null" if work is None else work,
            )
        )
    lines += [
        "",
        f"The clean selected F4 run takes {selected['metrics']['wall_seconds']:.6f} wall seconds and {selected['metrics']['total_core_seconds']:.6f} core-seconds. Construction falls from 26.323585 s to {selected['report']['timing_ns']['fixed_x1_s4_construction'] / 1e9:.6f} s, a {comparisons['selected_over_stage159']['construction_speedup']:.2f}x improvement. Whole-target wall improves {comparisons['selected_over_stage159']['wall_speedup']:.2f}x with identical equations, F4 calls, word XORs, roots, and verified witness.",
        "",
        f"Direct MITM still wins decisively on the same host: selected F4 uses {comparisons['selected_over_same_host_direct_mitm']['wall_ratio']:.2f}x its wall time, {comparisons['selected_over_same_host_direct_mitm']['core_ratio']:.2f}x its CPU, and {comparisons['selected_over_same_host_direct_mitm']['rss_ratio']:.2f}x its RSS. The clean build-plus-run cost is {comparisons['selected_cold_build_plus_run']['wall_seconds']:.3f} wall seconds and {comparisons['selected_cold_build_plus_run']['total_core_seconds']:.3f} core-seconds.",
        "",
        "WDSat and CryptoMiniSat remain exact-target frozen Linux timeouts on a different host. Licensed Magma remains unexecuted, GGMP has no same n=59 cell, full IC/rho cost is unchanged, and independent reproduction is still outstanding. This does not establish a SOTA.",
        "",
    ]
    phase_b.write_new(STAGE / "RESULTS.md", ("\n".join(lines)).encode())


def compose(run_source: Path, score_source: Path, development_root: Path) -> dict[str, Any]:
    require(not STAGE.exists(), f"refusing to overwrite {STAGE}")
    stage159 = stage159_verify.verify()
    require(stage159.get("status") == "verified", "Stage-159 baseline failed verification")
    baseline = load(STAGE159 / "result.json", "Stage-159 result")
    STAGE.mkdir()
    shutil.copytree(run_source, STAGE / "selected-run")
    shutil.copytree(score_source, STAGE / "selected-score")
    development = copy_development(development_root, STAGE / "development")

    run_seal, plan = native_f4.validate_run_seal(STAGE / "selected-run")
    require(plan["selected_blind_instance_ids"] == [BLIND_ID], "selected run changed")
    task = only_task(STAGE / "selected-run")
    selected = task["result"]
    verify_score(STAGE / "selected-score")
    report = selected["backend_report"]
    metrics = selected["metrics"]
    require(
        selected["status"] == "sat"
        and selected["source_model_valid"] is True
        and selected["source_witness_valid"] is True
        and report["solver_constructor"] == "fixed_x1_constant_linear_specialisation"
        and report["fixed_x1_systems_completed"] == 65
        and report["cost"]["extra"]["f4_calls"] == 65
        and report["cost"]["ops"] == 87_513_949_370,
        "selected native-F4 result changed",
    )
    old = baseline["selected_native_f4"]
    require(
        report["solver_equations_blake3"] == old["report"]["solver_equations_blake3"]
        and report["solver_terms_total"] == old["report"]["solver_terms_total"]
        and report["witness_points"] == old["report"]["witness_points"],
        "specialised constructor changed equations or witness",
    )

    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path, "selected clean build")
    build_tool.validate_receipt(receipt_path, receipt["binaries"], receipt["rust_source_objects"])
    require(
        receipt["source_commit"] == "3cb24b82e8a059a293f6371178bd0d3647272a63"
        and receipt["source_clean"] is True,
        "selected build source changed",
    )
    run_summary = load(STAGE / "selected-run/run-summary.json", "selected run summary")

    table = copy.deepcopy(baseline["same_target_table"])
    f4_row = next(row for row in table if row["backend"] == "native-f4-selected")
    f4_row.update(
        {
            "formulation": "fixed-X1 direct S4; constant-linear constructor; non-rational X1 skipped",
            "status": "sat",
            "classification": "true_positive",
            "wall_seconds": metrics["wall_seconds"],
            "total_core_seconds": metrics["total_core_seconds"],
            "peak_rss_bytes": metrics["peak_rss_bytes"],
            "conflicts": None,
            "operations": {
                "word_xors_elimination_only": report["cost"]["ops"],
                "f4_calls": report["cost"]["extra"]["f4_calls"],
            },
        }
    )
    mitm = next(row for row in table if row["backend"] == "direct-mitm")

    old_construction = old["report"]["timing_ns"]["fixed_x1_s4_construction"]
    new_construction = report["timing_ns"]["fixed_x1_s4_construction"]
    comparisons = {
        "selected_over_stage159": {
            "wall_ratio": metrics["wall_seconds"] / old["metrics"]["wall_seconds"],
            "wall_speedup": old["metrics"]["wall_seconds"] / metrics["wall_seconds"],
            "core_ratio": metrics["total_core_seconds"] / old["metrics"]["total_core_seconds"],
            "rss_ratio": metrics["peak_rss_bytes"] / old["metrics"]["peak_rss_bytes"],
            "construction_ratio": new_construction / old_construction,
            "construction_speedup": old_construction / new_construction,
            "same_equation_fingerprint": True,
            "same_term_count": True,
            "same_f4_calls": True,
            "same_word_xors": True,
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
    development_rows = [arm["metrics"] for arm in development["arms"]]
    incremental = [
        {
            "wall_seconds": receipt["resources"]["summed_process_wall_seconds"],
            "total_core_seconds": receipt["resources"]["total_core_seconds"],
            "peak_rss_bytes": receipt["resources"]["largest_child_peak_rss_bytes"],
        },
        {
            "wall_seconds": run_summary["outer_resources"]["wall_seconds"],
            "total_core_seconds": run_summary["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": metrics["peak_rss_bytes"],
        },
        *development_rows,
    ]
    incremental_total = add_resources(incremental)
    campaign = {
        "stage159_measured_lower_bound": prior,
        "stage160_incremental_measured": incremental_total,
        "measured_lower_bound": {
            "resource_components": prior["resource_components"]
            + incremental_total["resource_components"],
            "summed_wall_seconds": round(
                prior["summed_wall_seconds"] + incremental_total["summed_wall_seconds"], 12
            ),
            "total_core_seconds": round(
                prior["total_core_seconds"] + incremental_total["total_core_seconds"], 12
            ),
            "peak_rss_bytes": max(prior["peak_rss_bytes"], incremental_total["peak_rss_bytes"]),
        },
        "unmetered_development_commands_present": True,
        "unmetered_scope": [
            "incremental release compiles for development binaries",
            "focused Rust tests and Python evidence composition",
            "score composition, Git operations, and issue/PR updates",
        ],
        "complete_campaign_cost": None,
        "reason_complete_is_null": "some development commands lack an outer meter; measured totals remain a lower bound",
    }

    result = {
        "schema": RESULT_SCHEMA,
        "status": "complete_fixed_x1_constant_specialisation_true_positive",
        "date": "2026-09-22",
        "claim_boundary": (
            "Strong internal engineering and one finite public toy-PDP native-F4 construction improvement. "
            "It is not licensed Magma F4, a full native-F4 panel, an end-to-end index-calculus speedup, "
            "a direct-MITM improvement, or a Koblitz index-calculus SOTA."
        ),
        "instance": baseline["instance"],
        "factor_base": baseline["factor_base"],
        "mechanism": {
            "name": "fixed_x1_constant_linear_specialisation",
            "shared_symbolic_product": "X2*X3 formed once",
            "known_field_factors": "applied as exact linear maps in the polynomial basis",
            "equivalence_test": "fixed_x1_specialisation_matches_the_generic_s4_system",
            "generic_equivalence_cases": 60,
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
            "stage159_result_sha256": phase_b.sha256_file(
                STAGE159 / "result.json", "Stage-159 result"
            ),
            "stage159_seal_sha256": phase_b.sha256_file(
                STAGE159 / "result-seal.json", "Stage-159 seal"
            ),
        },
        "campaign_accounting": campaign,
        "gates": {
            "all_seven_gates_passed": False,
            "1_all_stage_resource_charging": "partial_measured_lower_bound_complete_total_null",
            "2_same_instance_solver_matrix": "partial_one_native_f4_n59_target_licensed_magma_and_full_panel_missing",
            "3_single_core_core_memory_conflicts_wall": "partial_native_f4_complete_conflicts_null_magma_missing",
            "4_n31_n41_larger_pdp": "inherited_satisfied_finite_coverage",
            "5_unknown_scalar_without_known_base_logs": "inherited_satisfied_finite_controls",
            "6_full_cost_vs_automorphism_rho": "unchanged_not_tested_by_pdp_stage_and_current_full_cost_gate_false",
            "7_independent_external_reproduction": "requested_but_missing",
        },
        "licensed_magma_f4_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    phase_b.write_json_new(STAGE / "result.json", result)
    write_markdown(result)
    inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
    seal = {
        "schema": SEAL_SCHEMA,
        "status": "stage160_fixed_x1_result_frozen",
        "result_sha256": phase_b.sha256_file(STAGE / "result.json", "Stage-160 result"),
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
        Stage160Error,
        phase_b.PhaseBError,
    ) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()

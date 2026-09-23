#!/usr/bin/env python3
"""Compose the sealed Stage-161 trusted-mask hashing result."""

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
import verify_koblitz_stage160_fixed_x1 as stage160_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE160 = GATES / "stage-160-fixed-x1-constant-specialisation-20260922"
STAGE = GATES / "stage-161-fast-mask-hash-20260923"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
SOURCE_ID = "fb0ea6bb30be3dea841c609a2668939128b4b799d4c0e051b0b52dc6da05f0fa"
CELL = "n59-l9-m3-standard-a1-f0"
RESULT_SCHEMA = "koblitz_stage161_fast_mask_hash.v1"
SEAL_SCHEMA = "koblitz_stage161_fast_mask_hash_seal.v1"


class Stage161Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage161Error(message)


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


def verify_score(root: Path) -> None:
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
        "selected score changed",
    )


def copy_artifact(source: Path, destination: Path, context: str) -> dict[str, Any]:
    require(source.is_file(), f"missing {context}: {source}")
    shutil.copyfile(source, destination)
    return {
        "path": destination.name,
        "bytes": destination.stat().st_size,
        "sha256": phase_b.sha256_file(destination, context),
    }


def copy_development(root: Path, output: Path) -> dict[str, Any]:
    output.mkdir()
    profile_artifacts = {
        suffix: copy_artifact(
            root / f"stage161-profile2-f4.{suffix}",
            output / f"profile-baseline.{suffix}",
            f"profile {suffix}",
        )
        for suffix in ("metrics", "stdout", "stderr", "sample.txt")
    }
    candidate_artifacts = {
        suffix: copy_artifact(
            root / f"stage161-fast-hash-f4.{suffix}",
            output / f"fast-hash-provisional.{suffix}",
            f"candidate {suffix}",
        )
        for suffix in ("metrics", "stdout", "stderr")
    }
    profile_metrics = load(
        output / profile_artifacts["metrics"]["path"], "profile metrics"
    )["metrics"]
    profile_report = load(output / profile_artifacts["stdout"]["path"], "profile stdout")
    candidate_metrics = load(
        output / candidate_artifacts["metrics"]["path"], "candidate metrics"
    )["metrics"]
    candidate_report = load(
        output / candidate_artifacts["stdout"]["path"], "candidate stdout"
    )
    for name, report in (("profile", profile_report), ("candidate", candidate_report)):
        require(
            report.get("status") == "sat"
            and report.get("source_model_valid") is True
            and report.get("source_witness_valid") is True
            and report.get("solver_equations_blake3")
            == "d6969c5f89b5be6a7eacbdac45799fc734a273c9034949a482723b4422504ae7",
            f"{name} development terminal changed",
        )
    manifest = {
        "schema": "koblitz_stage161_development.v1",
        "status": "profile_and_provisional_candidate_retained",
        "profile": {
            "decision": "diagnostic_baseline",
            "metrics": profile_metrics,
            "sample_interval_ms": 5,
            "sample_duration_seconds": 20,
            "samples": 3334,
            "top_of_stack": {
                "echelon": 1143,
                "state_insert": 533,
                "columns_pack": 238,
                "siphash_write": 225,
                "vec_retain": 118,
            },
            "artifacts": profile_artifacts,
        },
        "candidate": {
            "decision": "promoted_after_clean_rebuild",
            "metrics": candidate_metrics,
            "f4_build_ns": candidate_report["cost"]["extra"]["build_ns"],
            "f4_eliminate_ns": candidate_report["cost"]["extra"]["eliminate_ns"],
            "word_xors": candidate_report["cost"]["ops"],
            "artifacts": candidate_artifacts,
        },
        "failed_profiler_attempt": {
            "status": "sandbox_process_enumeration_denied_then_interrupted",
            "resource_receipt": None,
            "charged": False,
        },
        "complete_development_cost": None,
        "reason_complete_is_null": "one interrupted sandboxed profiling attempt and incremental compiles lack complete outer receipts",
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
        "# Stage 161: fast hashing for trusted F4 monomial masks",
        "",
        result["claim_boundary"],
        "",
        "A charged 20-second sample of the Stage 160 binary put 1,143 of 3,334 top-of-stack samples in dense echelon elimination, 533 in critical-pair installation, 238 in column packing, and 225 directly in SipHash writes. Stage 161 retains exact-key hash maps but replaces SipHash only for solver-constructed `u64` monomial masks with deterministic SplitMix64 hashing.",
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
    old = result["predecessor_metrics"]
    lines += [
        "",
        f"The clean selected F4 process takes {selected['metrics']['wall_seconds']:.6f} wall seconds and {selected['metrics']['total_core_seconds']:.6f} core-seconds. Symbolic/matrix build falls from {old['f4_build_seconds']:.6f} to {selected['report']['cost']['extra']['build_ns'] / 1e9:.6f} seconds, a {comparisons['selected_over_stage160']['build_phase_speedup']:.2f}x improvement. Whole-target wall improves {comparisons['selected_over_stage160']['wall_speedup']:.2f}x with identical equations, calls, matrices, pair counters, word XORs, roots, and verified witness.",
        "",
        f"Peak RSS increases {100 * (comparisons['selected_over_stage160']['rss_ratio'] - 1):.2f}% to {selected['metrics']['peak_rss_bytes']} bytes. Direct MITM still wins: selected F4 uses {comparisons['selected_over_same_host_direct_mitm']['wall_ratio']:.2f}x its wall, {comparisons['selected_over_same_host_direct_mitm']['core_ratio']:.2f}x its CPU, and {comparisons['selected_over_same_host_direct_mitm']['rss_ratio']:.2f}x its RSS. Clean build plus run costs {comparisons['selected_cold_build_plus_run']['wall_seconds']:.3f} wall seconds and {comparisons['selected_cold_build_plus_run']['total_core_seconds']:.3f} core-seconds.",
        "",
        "Licensed Magma, a full native-F4 panel, end-to-end IC/rho cost, and independent external reproduction remain open. This one-target engineering result does not establish a SOTA.",
        "",
    ]
    phase_b.write_new(STAGE / "RESULTS.md", ("\n".join(lines)).encode())


def compose(run_source: Path, score_source: Path, development_root: Path) -> dict[str, Any]:
    require(not STAGE.exists(), f"refusing to overwrite {STAGE}")
    verification160 = stage160_verify.verify()
    require(verification160.get("status") == "verified", "Stage-160 failed verification")
    baseline = load(STAGE160 / "result.json", "Stage-160 result")
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
        and report["solver_internal_mask_hasher"]
        == "splitmix64_for_trusted_u64_masks_with_exact_key_equality"
        and report["fixed_x1_systems_completed"] == 65
        and report["cost"]["extra"]["f4_calls"] == 65
        and report["cost"]["ops"] == 87_513_949_370,
        "selected native-F4 result changed",
    )
    old = baseline["selected_native_f4"]
    structural = (
        "basis_len",
        "divisor_tests",
        "extraction_tests",
        "f4_calls",
        "field_pairs_reduced",
        "matrix_cols_max",
        "matrix_rows_max",
        "matrix_rows_sum",
        "new_elements",
        "oversize",
        "pairs_chain_skipped",
        "pairs_left",
        "pairs_product_skipped",
        "pairs_reduced",
        "reducer_rows",
        "steps",
        "symbolic_bytes_estimate_max",
        "symbolic_cap_hit",
    )
    require(
        report["solver_equations_blake3"] == old["report"]["solver_equations_blake3"]
        and report["solver_terms_total"] == old["report"]["solver_terms_total"]
        and report["witness_points"] == old["report"]["witness_points"]
        and all(
            report["cost"]["extra"][key] == old["report"]["cost"]["extra"][key]
            for key in structural
        ),
        "fast hash changed equations, F4 structure, or witness",
    )

    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path, "selected clean build")
    build_tool.validate_receipt(receipt_path, receipt["binaries"], receipt["rust_source_objects"])
    require(
        receipt["source_commit"] == "6906d0a26e4cb7e84a99e864c6f471fef36e5526"
        and receipt["source_clean"] is True,
        "selected build source changed",
    )
    run_summary = load(STAGE / "selected-run/run-summary.json", "selected run summary")

    table = copy.deepcopy(baseline["same_target_table"])
    f4_row = next(row for row in table if row["backend"] == "native-f4-selected")
    f4_row.update(
        {
            "formulation": "fixed-X1 direct S4; constant-linear constructor; fast trusted-mask hashing",
            "wall_seconds": metrics["wall_seconds"],
            "total_core_seconds": metrics["total_core_seconds"],
            "peak_rss_bytes": metrics["peak_rss_bytes"],
            "operations": {
                "word_xors_elimination_only": report["cost"]["ops"],
                "f4_calls": report["cost"]["extra"]["f4_calls"],
            },
        }
    )
    mitm = next(row for row in table if row["backend"] == "direct-mitm")

    old_build_ns = old["report"]["cost"]["extra"]["build_ns"]
    new_build_ns = report["cost"]["extra"]["build_ns"]
    comparisons = {
        "selected_over_stage160": {
            "wall_ratio": metrics["wall_seconds"] / old["metrics"]["wall_seconds"],
            "wall_speedup": old["metrics"]["wall_seconds"] / metrics["wall_seconds"],
            "core_ratio": metrics["total_core_seconds"] / old["metrics"]["total_core_seconds"],
            "rss_ratio": metrics["peak_rss_bytes"] / old["metrics"]["peak_rss_bytes"],
            "build_phase_ratio": new_build_ns / old_build_ns,
            "build_phase_speedup": old_build_ns / new_build_ns,
            "same_equation_fingerprint": True,
            "same_term_count": True,
            "same_f4_structure": True,
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
    incremental = [
        development["profile"]["metrics"],
        development["candidate"]["metrics"],
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
    ]
    incremental_total = add_resources(incremental)
    campaign = {
        "stage160_measured_lower_bound": prior,
        "stage161_incremental_measured": incremental_total,
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
            "one interrupted sandbox-denied profiler attempt",
            "incremental release compiles and focused tests",
            "score/evidence composition, Git operations, and issue/PR updates",
        ],
        "complete_campaign_cost": None,
        "reason_complete_is_null": "some development commands lack complete outer receipts; measured totals remain a lower bound",
    }

    result = {
        "schema": RESULT_SCHEMA,
        "status": "complete_fast_mask_hash_true_positive",
        "date": "2026-09-23",
        "claim_boundary": (
            "Strong internal engineering and one finite public toy-PDP native-F4 hashing improvement. "
            "It is not licensed Magma F4, a full native-F4 panel, an end-to-end index-calculus speedup, "
            "a direct-MITM improvement, or a Koblitz index-calculus SOTA."
        ),
        "instance": baseline["instance"],
        "factor_base": baseline["factor_base"],
        "mechanism": {
            "name": "splitmix64_trusted_f4_masks",
            "scope": "column indices and symbolic-preprocessing sets keyed only by solver-constructed u64 monomial masks",
            "collision_semantics": "hashbrown retains exact key equality",
            "untrusted_input_hashed_directly": False,
        },
        "predecessor_metrics": {
            "stage160_wall_seconds": old["metrics"]["wall_seconds"],
            "stage160_core_seconds": old["metrics"]["total_core_seconds"],
            "stage160_peak_rss_bytes": old["metrics"]["peak_rss_bytes"],
            "f4_build_seconds": old_build_ns / 1e9,
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
            "stage160_result_sha256": phase_b.sha256_file(
                STAGE160 / "result.json", "Stage-160 result"
            ),
            "stage160_seal_sha256": phase_b.sha256_file(
                STAGE160 / "result-seal.json", "Stage-160 seal"
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
        "status": "stage161_fast_hash_result_frozen",
        "result_sha256": phase_b.sha256_file(STAGE / "result.json", "Stage-161 result"),
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
        Stage161Error,
        phase_b.PhaseBError,
    ) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()

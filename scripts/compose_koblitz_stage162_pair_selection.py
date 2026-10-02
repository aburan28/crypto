#!/usr/bin/env python3
"""Compose the sealed Stage-162 critical-pair selection result."""

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
import verify_koblitz_stage161_fast_hash as stage161_verify


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE161 = GATES / "stage-161-fast-mask-hash-20260923"
STAGE = GATES / "stage-162-grouped-pair-selection-20260923"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
SOURCE_ID = "fb0ea6bb30be3dea841c609a2668939128b4b799d4c0e051b0b52dc6da05f0fa"
CELL = "n59-l9-m3-standard-a1-f0"


class Stage162Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage162Error(message)


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
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256"),
        "selected score seal changed",
    )
    score = load(root / "score.json", "score")
    require(score.get("classification_counts") == {"true_positive": 1}, "score changed")


def copy_file(source: Path, destination: Path, context: str) -> dict[str, Any]:
    require(source.is_file(), f"missing {context}: {source}")
    shutil.copyfile(source, destination)
    return {
        "path": destination.name,
        "bytes": destination.stat().st_size,
        "sha256": phase_b.sha256_file(destination, context),
    }


def copy_development(root: Path, output: Path) -> dict[str, Any]:
    output.mkdir()
    profile = {
        suffix: copy_file(
            root / f"stage162-profile-fast-f4.{suffix}",
            output / f"profile-fast-hash.{suffix}",
            f"profile {suffix}",
        )
        for suffix in ("metrics", "stdout", "stderr", "sample.txt")
    }
    candidate = {
        suffix: copy_file(
            root / f"stage162-grouped-pairs-f4.{suffix}",
            output / f"grouped-pairs-provisional.{suffix}",
            f"candidate {suffix}",
        )
        for suffix in ("metrics", "stdout", "stderr")
    }
    profile_metrics = load(output / profile["metrics"]["path"], "profile metrics")["metrics"]
    profile_report = load(output / profile["stdout"]["path"], "profile stdout")
    candidate_metrics = load(output / candidate["metrics"]["path"], "candidate metrics")["metrics"]
    candidate_report = load(output / candidate["stdout"]["path"], "candidate stdout")
    for name, report in (("profile", profile_report), ("candidate", candidate_report)):
        require(
            report.get("status") == "sat"
            and report.get("source_witness_valid") is True
            and report.get("solver_equations_blake3")
            == "d6969c5f89b5be6a7eacbdac45799fc734a273c9034949a482723b4422504ae7",
            f"{name} terminal changed",
        )
    manifest = {
        "schema": "koblitz_stage162_development.v1",
        "status": "post_hash_profile_and_candidate_retained",
        "profile": {
            "metrics": profile_metrics,
            "samples": 3327,
            "sample_interval_ms": 5,
            "sample_duration_seconds": 20,
            "top_of_stack": {
                "echelon": 1366,
                "state_insert": 695,
                "stable_quicksort": 251,
                "columns_pack": 153,
                "hashmap_insert": 121,
                "vec_retain": 109,
            },
            "artifacts": profile,
        },
        "candidate": {
            "decision": "promoted_after_clean_rebuild",
            "metrics": candidate_metrics,
            "word_xors": candidate_report["cost"]["ops"],
            "artifacts": candidate,
        },
        "complete_development_cost": None,
        "reason_complete_is_null": "incremental compile and focused tests lack outer receipts",
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
    s = result["selected_native_f4"]
    c = result["comparisons"]
    lines = [
        "# Stage 162: grouped exact F4 critical-pair selection",
        "",
        result["claim_boundary"],
        "",
        "A charged post-hash profile put 695 of 3,327 top-of-stack samples in `State::insert`. Disassembly located the dominant loop in the quadratic scan that asks whether any other new-pair LCM divides the current LCM. Stage 162 groups equal LCMs, finds proper divisor LCMs by exact submask lookup, and then emits the same lowest-index representative in the same order as the prior UPDATE algorithm.",
        "",
        f"The clean selected process takes {s['metrics']['wall_seconds']:.6f} wall seconds, {s['metrics']['total_core_seconds']:.6f} core-seconds, and {s['metrics']['peak_rss_bytes']} bytes peak RSS. Whole-target wall improves {c['selected_over_stage161']['wall_speedup']:.3f}x over Stage 161. Equations, term count, F4 calls, matrices, pair counters, 87,513,949,370 elimination word XORs, roots, and witness are identical.",
        "",
        f"Direct MITM still wins by {c['selected_over_same_host_direct_mitm']['wall_ratio']:.2f}x wall, {c['selected_over_same_host_direct_mitm']['core_ratio']:.2f}x CPU, and {c['selected_over_same_host_direct_mitm']['rss_ratio']:.2f}x RSS. Clean build plus run costs {c['selected_cold_build_plus_run']['wall_seconds']:.3f} wall seconds and {c['selected_cold_build_plus_run']['total_core_seconds']:.3f} core-seconds.",
        "",
        "Licensed Magma, a full native-F4 panel, end-to-end IC/rho cost, and independent reproduction remain open. This does not establish a SOTA.",
        "",
    ]
    phase_b.write_new(STAGE / "RESULTS.md", ("\n".join(lines)).encode())


def compose(run_source: Path, score_source: Path, development_root: Path) -> dict[str, Any]:
    require(not STAGE.exists(), f"refusing to overwrite {STAGE}")
    require(stage161_verify.verify().get("status") == "verified", "Stage-161 failed verification")
    baseline = load(STAGE161 / "result.json", "Stage-161 result")
    STAGE.mkdir()
    shutil.copytree(run_source, STAGE / "selected-run")
    shutil.copytree(score_source, STAGE / "selected-score")
    development = copy_development(development_root, STAGE / "development")
    run_seal, plan = native_f4.validate_run_seal(STAGE / "selected-run")
    require(plan["selected_blind_instance_ids"] == [BLIND_ID], "selection changed")
    selected = only_task(STAGE / "selected-run")["result"]
    verify_score(STAGE / "selected-score")
    report = selected["backend_report"]
    metrics = selected["metrics"]
    require(
        selected["status"] == "sat"
        and selected["source_model_valid"] is True
        and selected["source_witness_valid"] is True
        and report["solver_pair_selector"]
        == "grouped_lcm_exact_submask_equivalent_to_quadratic_update"
        and report["cost"]["ops"] == 87_513_949_370,
        "selected F4 result changed",
    )
    old = baseline["selected_native_f4"]
    structural = (
        "basis_len", "divisor_tests", "extraction_tests", "f4_calls",
        "field_pairs_reduced", "matrix_cols_max", "matrix_rows_max",
        "matrix_rows_sum", "new_elements", "oversize", "pairs_chain_skipped",
        "pairs_left", "pairs_product_skipped", "pairs_reduced", "reducer_rows",
        "steps", "symbolic_bytes_estimate_max", "symbolic_cap_hit",
    )
    require(
        report["solver_equations_blake3"] == old["report"]["solver_equations_blake3"]
        and report["solver_terms_total"] == old["report"]["solver_terms_total"]
        and report["witness_points"] == old["report"]["witness_points"]
        and all(report["cost"]["extra"][k] == old["report"]["cost"]["extra"][k] for k in structural),
        "pair selector changed F4 structure",
    )
    receipt_path = STAGE / "selected-run/tool-builds/rust/receipt.json"
    receipt = load(receipt_path, "build receipt")
    build_tool.validate_receipt(receipt_path, receipt["binaries"], receipt["rust_source_objects"])
    require(
        receipt["source_commit"] == "29524b34ef5bfcbb726a927c362e825c3203c390"
        and receipt["source_clean"] is True,
        "clean source changed",
    )
    run_summary = load(STAGE / "selected-run/run-summary.json", "run summary")
    table = copy.deepcopy(baseline["same_target_table"])
    f4 = next(row for row in table if row["backend"] == "native-f4-selected")
    f4.update(
        {
            "formulation": "fixed-X1 direct S4; fast trusted-mask hash; grouped exact pair selection",
            "wall_seconds": metrics["wall_seconds"],
            "total_core_seconds": metrics["total_core_seconds"],
            "peak_rss_bytes": metrics["peak_rss_bytes"],
        }
    )
    mitm = next(row for row in table if row["backend"] == "direct-mitm")
    comparisons = {
        "selected_over_stage161": {
            "wall_ratio": metrics["wall_seconds"] / old["metrics"]["wall_seconds"],
            "wall_speedup": old["metrics"]["wall_seconds"] / metrics["wall_seconds"],
            "core_ratio": metrics["total_core_seconds"] / old["metrics"]["total_core_seconds"],
            "rss_ratio": metrics["peak_rss_bytes"] / old["metrics"]["peak_rss_bytes"],
            "same_equation_fingerprint": True,
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
            "wall_seconds": receipt["resources"]["summed_process_wall_seconds"] + run_summary["outer_resources"]["wall_seconds"],
            "total_core_seconds": receipt["resources"]["total_core_seconds"] + run_summary["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": max(receipt["resources"]["largest_child_peak_rss_bytes"], metrics["peak_rss_bytes"]),
        },
    }
    prior = baseline["campaign_accounting"]["measured_lower_bound"]
    incremental = [
        development["profile"]["metrics"], development["candidate"]["metrics"],
        {"wall_seconds": receipt["resources"]["summed_process_wall_seconds"], "total_core_seconds": receipt["resources"]["total_core_seconds"], "peak_rss_bytes": receipt["resources"]["largest_child_peak_rss_bytes"]},
        {"wall_seconds": run_summary["outer_resources"]["wall_seconds"], "total_core_seconds": run_summary["outer_resources"]["total_core_seconds"], "peak_rss_bytes": metrics["peak_rss_bytes"]},
    ]
    inc = add_resources(incremental)
    campaign = {
        "stage161_measured_lower_bound": prior,
        "stage162_incremental_measured": inc,
        "measured_lower_bound": {
            "resource_components": prior["resource_components"] + inc["resource_components"],
            "summed_wall_seconds": round(prior["summed_wall_seconds"] + inc["summed_wall_seconds"], 12),
            "total_core_seconds": round(prior["total_core_seconds"] + inc["total_core_seconds"], 12),
            "peak_rss_bytes": max(prior["peak_rss_bytes"], inc["peak_rss_bytes"]),
        },
        "unmetered_development_commands_present": True,
        "unmetered_scope": ["incremental compile and focused tests", "score/evidence composition and Git operations"],
        "complete_campaign_cost": None,
        "reason_complete_is_null": "some development commands lack complete outer receipts",
    }
    result = {
        "schema": "koblitz_stage162_grouped_pair_selection.v1",
        "status": "complete_grouped_pair_selection_true_positive",
        "date": "2026-09-23",
        "claim_boundary": "Strong internal engineering and one finite public toy-PDP native-F4 pair-selection improvement. It is not a full solver panel, direct-MITM improvement, full IC/rho crossover, independent reproduction, or Koblitz index-calculus SOTA.",
        "instance": baseline["instance"],
        "factor_base": baseline["factor_base"],
        "mechanism": {
            "name": "grouped_lcm_exact_submask",
            "reference_equivalence_test": "grouped_pair_selection_matches_quadratic_update",
            "randomized_reference_cases": 3600,
            "preserves": ["lowest-index duplicate representative", "critical-pair order", "chain/product counters"],
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
            "stage161_result_sha256": phase_b.sha256_file(STAGE161 / "result.json", "Stage-161 result"),
            "stage161_seal_sha256": phase_b.sha256_file(STAGE161 / "result-seal.json", "Stage-161 seal"),
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
        "schema": "koblitz_stage162_grouped_pair_selection_seal.v1",
        "status": "stage162_grouped_pair_result_frozen",
        "result_sha256": phase_b.sha256_file(STAGE / "result.json", "Stage-162 result"),
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
        print(json.dumps(compose(args.run.resolve(strict=True), args.score.resolve(strict=True), args.development_root.resolve(strict=True)), indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, IndexError, Stage162Error, phase_b.PhaseBError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()

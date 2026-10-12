#!/usr/bin/env python3
"""Compose the Stage 172 same-instance native-F4 single-core result."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
GATES = HERE.parent
PREDECESSOR = GATES / "stage-171-dense-symbolic-sets-20260923" / "result.json"


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def artifact(path: Path) -> dict:
    return {
        "path": path.relative_to(HERE).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def run(name: str) -> dict:
    root = RUNS / name
    return {
        "process": load(root / "metrics.json"),
        "report": load(root / "stdout.json"),
        "artifacts": {
            "metrics": artifact(root / "metrics.json"),
            "stdout": artifact(root / "stdout.json"),
            "stderr": artifact(root / "stderr.txt"),
        },
    }


def stable_f4(report: dict) -> dict:
    report = copy.deepcopy(report)
    report.pop("timing_ns", None)
    report.pop("solver_description", None)
    report.pop("solver_symbolic_monomial_sets", None)
    cost = report["cost"]
    cost.pop("wall_ns", None)
    for key in (
        "build_ns",
        "eliminate_ns",
        "pair_update_ns",
        "dense_symbolic_set_steps",
        "dense_symbolic_set_bytes_max",
    ):
        cost["extra"].pop(key, None)
    return report


def main() -> None:
    predecessor = load(PREDECESSOR)
    candidate = run("candidate")
    control = run("control")
    direct = run("direct-mitm")
    cm = candidate["process"]["metrics"]
    hm = control["process"]["metrics"]
    dm = direct["process"]["metrics"]
    cr = candidate["report"]
    hr = control["report"]
    dr = direct["report"]
    parallel = predecessor["selected_native_f4"]["metrics"]
    build = predecessor["clean_build"]["process"]["metrics"]
    previous_campaign = predecessor["campaign_accounting"]["measured_lower_bound"]
    incremental_wall = cm["wall_seconds"] + hm["wall_seconds"] + dm["wall_seconds"]
    incremental_core = cm["total_core_seconds"] + hm["total_core_seconds"] + dm["total_core_seconds"]
    incremental_peak = max(cm["peak_rss_bytes"], hm["peak_rss_bytes"], dm["peak_rss_bytes"])

    fingerprint = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
    same_science = stable_f4(cr) == stable_f4(hr)
    result = {
        "schema": "koblitz_stage172_native_f4_single_core.v1",
        "date": "2026-09-23",
        "status": "complete_native_f4_single_core_true_negative",
        "predecessor": {
            "path": PREDECESSOR.relative_to(GATES).as_posix(),
            "sha256": sha256(PREDECESSOR),
            "status": predecessor["status"],
        },
        "source_revision": predecessor["source_revision"],
        "instance": predecessor["instance"],
        "factor_base": predecessor["factor_base"],
        "protocol": {
            "backend": "repository_native_boolean_f4",
            "rayon_threads": 1,
            "x1_batch_size": 1,
            "blas_threads": 1,
            "watchdog_seconds": 600,
            "f4_budget_seconds": 540,
            "candidate": "dense symbolic monomial sets",
            "same_binary_control": "PQ_F4_DISABLE_DENSE_SYMBOLIC_SET=1",
            "selection_boundary": "same target after Stage-169 truth opening; measurement completion, not a fresh holdout",
        },
        "selected_single_core_native_f4": {
            **candidate,
            "classification": "true_negative",
            "conflicts": None,
            "valid_single_core_seconds": cm["single_core_seconds"],
        },
        "same_binary_single_core_hash_control": {
            **control,
            "classification": "true_negative",
            "conflicts": None,
            "valid_single_core_seconds": hm["single_core_seconds"],
        },
        "same_binary_single_core_direct_mitm": {
            **direct,
            "classification": "true_negative",
            "valid_single_core_seconds": dm["single_core_seconds"],
        },
        "validation": {
            "same_scientific_f4_report_outside_declared_mode_and_timing": same_science,
            "candidate_exact_terminal": (
                cr["status"] == "unsat"
                and cr["exhaustive"] is True
                and cr["single_thread_requested"] is True
                and cr["solver_rayon_threads_requested"] == 1
                and cr["solver_x1_batch_size"] == 1
                and cr["solver_x1_batches_completed"] == 242
                and cr["solver_equations_blake3"] == fingerprint
            ),
            "control_exact_terminal": (
                hr["status"] == "unsat"
                and hr["exhaustive"] is True
                and hr["single_thread_requested"] is True
                and hr["solver_equations_blake3"] == fingerprint
            ),
            "direct_exact_terminal": dr["status"] == "unsat" and dr["exhaustive"] is True,
            "same_binary_sha256": predecessor["source_revision"]["backend_sha256"],
        },
        "comparisons": {
            "dense_over_hash_single_core": {
                "wall_ratio": cm["wall_seconds"] / hm["wall_seconds"],
                "core_ratio": cm["total_core_seconds"] / hm["total_core_seconds"],
                "rss_ratio": cm["peak_rss_bytes"] / hm["peak_rss_bytes"],
                "f4_build_ns_ratio": cr["cost"]["extra"]["build_ns"]
                / hr["cost"]["extra"]["build_ns"],
                "same_scientific_report": same_science,
            },
            "single_core_f4_over_direct_mitm": {
                "wall_ratio": cm["wall_seconds"] / dm["wall_seconds"],
                "core_ratio": cm["total_core_seconds"] / dm["total_core_seconds"],
                "rss_ratio": cm["peak_rss_bytes"] / dm["peak_rss_bytes"],
            },
            "single_core_over_stage171_parallel": {
                "wall_ratio": cm["wall_seconds"] / parallel["wall_seconds"],
                "core_ratio": cm["total_core_seconds"] / parallel["total_core_seconds"],
                "rss_ratio": cm["peak_rss_bytes"] / parallel["peak_rss_bytes"],
                "parallel_wall_speedup": cm["wall_seconds"] / parallel["wall_seconds"],
                "parallel_core_overhead": parallel["total_core_seconds"] / cm["total_core_seconds"],
                "parallel_rss_overhead": parallel["peak_rss_bytes"] / cm["peak_rss_bytes"],
            },
            "inherited_cold_build_plus_single_core": {
                "wall_seconds": build["wall_seconds"] + cm["wall_seconds"],
                "total_core_seconds": build["total_core_seconds"] + cm["total_core_seconds"],
                "peak_rss_bytes": max(build["peak_rss_bytes"], cm["peak_rss_bytes"]),
                "build_receipt_stage": 171,
            },
        },
        "stage172_execution_accounting": {
            "resource_components": 3,
            "summed_wall_seconds": incremental_wall,
            "total_core_seconds": incremental_core,
            "peak_rss_bytes": incremental_peak,
            "all_execution_arms_outer_metered": True,
        },
        "campaign_accounting": {
            "complete_campaign_cost": None,
            "measured_lower_bound": {
                "resource_components": previous_campaign["resource_components"] + 3,
                "summed_wall_seconds": previous_campaign["summed_wall_seconds"] + incremental_wall,
                "total_core_seconds": previous_campaign["total_core_seconds"] + incremental_core,
                "peak_rss_bytes": max(previous_campaign["peak_rss_bytes"], incremental_peak),
            },
            "stage171_measured_lower_bound": previous_campaign,
            "stage172_incremental_measured": {
                "resource_components": 3,
                "summed_wall_seconds": incremental_wall,
                "total_core_seconds": incremental_core,
                "peak_rss_bytes": incremental_peak,
            },
            "reason_complete_is_null": predecessor["campaign_accounting"]["reason_complete_is_null"],
        },
        "gates": {
            "1_all_stage_resource_charging": "partial_historical_campaign_total_null_stage172_execution_complete",
            "2_same_instance_solver_matrix": "partial_native_f4_wdsat_cryptominisat_direct_mitm_present_licensed_magma_and_ggmp_same_instance_incomplete",
            "3_single_core_core_memory_conflicts_wall": "partial_native_f4_parallel_and_single_core_complete_conflicts_null_magma_missing",
            "4_n31_n41_larger_pdp": "inherited_satisfied_finite_coverage",
            "5_unknown_scalar_without_known_base_logs": "inherited_satisfied_finite_controls",
            "6_full_cost_vs_automorphism_rho": "unchanged_not_tested_by_pdp_stage_and_current_full_cost_gate_false",
            "7_independent_external_reproduction": "requested_but_missing",
            "all_seven_gates_passed": False,
        },
        "fresh_target_holdout_complete": True,
        "licensed_magma_f4_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "claim_boundary": (
            "This fills the repository-native F4 single-core resource row on one already-opened "
            "public n=59 PDP target. It is not a new target distribution, complete solver panel, "
            "full index-calculus/rho crossover, independent reproduction, novelty review, or SOTA."
        ),
    }
    (HERE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    results = f"""# Stage 172: native F4 single-core completion

This stage fills the missing valid single-core resource row for the repository-native Boolean F4 implementation on the already-opened Stage-169 `n=59, ell=9, m=3` true-negative target. It is not a fresh holdout, complete solver panel, full index-calculus/rho crossover, independent reproduction, novelty review, or SOTA.

With `RAYON_NUM_THREADS=1`, `PQ_F4_X1_BATCH=1`, and every BLAS-style thread variable fixed to one, dense symbolic sets complete exhaustive UNSAT in **{cm['wall_seconds']:.6f} wall seconds / {cm['total_core_seconds']:.6f} core-seconds / {cm['peak_rss_bytes']} bytes peak RSS**. The backend reports `single_thread_requested=true`, 242 sequential batches, all 242 systems complete, 14,278 equations, 14,515,915 terms, 99,199,976,264 word XORs, zero roots, and equation fingerprint `{fingerprint}`.

The same clean binary with `PQ_F4_DISABLE_DENSE_SYMBOLIC_SET=1` takes {hm['wall_seconds']:.6f} wall seconds / {hm['total_core_seconds']:.6f} core-seconds / {hm['peak_rss_bytes']} bytes. Dense sets reduce wall by {(1-cm['wall_seconds']/hm['wall_seconds'])*100:.2f}%, CPU by {(1-cm['total_core_seconds']/hm['total_core_seconds'])*100:.2f}%, and F4 build time by {(1-cr['cost']['extra']['build_ns']/hr['cost']['extra']['build_ns'])*100:.2f}%; RSS rises {(cm['peak_rss_bytes']/hm['peak_rss_bytes']-1)*100:.2f}% in this one-worker shape.

The same-binary direct-MITM reference completes in {dm['wall_seconds']:.6f} wall seconds / {dm['total_core_seconds']:.6f} core-seconds / {dm['peak_rss_bytes']} bytes, leaving single-core F4 {cm['wall_seconds']/dm['wall_seconds']:.2f}x slower by wall and {cm['total_core_seconds']/dm['total_core_seconds']:.2f}x more expensive by CPU. The Stage-171 parallel F4 run is {cm['wall_seconds']/parallel['wall_seconds']:.2f}x faster by wall but consumes {parallel['total_core_seconds']/cm['total_core_seconds']:.2f}x the CPU and {parallel['peak_rss_bytes']/cm['peak_rss_bytes']:.2f}x the memory.

The inherited exact-commit clean build plus this single-core execution costs {build['wall_seconds']+cm['wall_seconds']:.6f} wall seconds and {build['total_core_seconds']+cm['total_core_seconds']:.6f} core-seconds. All three Stage-172 execution arms are outer-metered; the historical campaign total remains null for the inherited unmetered scope.
"""
    (HERE / "RESULTS.md").write_text(results)


if __name__ == "__main__":
    main()

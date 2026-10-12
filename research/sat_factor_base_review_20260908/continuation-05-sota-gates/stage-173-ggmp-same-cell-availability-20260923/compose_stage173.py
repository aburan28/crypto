#!/usr/bin/env python3
"""Compose the Stage 173 target-independent GGMP availability result."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEV = HERE / "development"
GATES = HERE.parent
PREDECESSOR = GATES / "stage-172-native-f4-single-core-20260923" / "result.json"
COMMIT = "c6205c5697cff517a6f46d0a9a96a0293d4af695"
TREE = "6bb1c1a4abdfe3d3aa2d9b1753b052c932a6e43c"
ARCHIVE_SHA256 = "b594668c9f04bfe8750f4d9feeadcd667b71ec1f53755a858159c1695ba51d4b"
BINARY_SHA256 = "7704956b1e739c559192108ff2bc106173d8abe73f82e2a8eb54eda29ca009eb"
LOCK_SHA256 = "10461fc25d189a7e1b533c0a59e1c66d7376ca04d3af5c5868a6f155efc8e5cf"


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


def process(root: Path) -> dict:
    return {
        "process": load(root / "metrics.json"),
        "report": load(root / "stdout.json"),
        "artifacts": {
            "metrics": artifact(root / "metrics.json"),
            "stdout": artifact(root / "stdout.json"),
            "stderr": artifact(root / "stderr.txt"),
        },
    }


def all_metered() -> dict:
    paths = sorted(set(DEV.rglob("metrics.json")) | set(DEV.rglob("*.metrics.json")))
    rows = []
    for path in paths:
        try:
            record = load(path)
        except (OSError, json.JSONDecodeError):
            continue
        metrics = record.get("metrics", record)
        if all(key in metrics for key in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")):
            rows.append(metrics)
    return {
        "resource_components": len(rows),
        "summed_wall_seconds": sum(row["wall_seconds"] for row in rows),
        "total_core_seconds": sum(row["total_core_seconds"] for row in rows),
        "peak_rss_bytes": max(row["peak_rss_bytes"] for row in rows),
        "complete_total": None,
        "reason_complete_total_is_null": (
            "the initial consumer environment failure and cached development compilation lack "
            "outer process receipts; every formal build, probe, and final consumer test is metered"
        ),
    }


def test_receipt(name: str) -> dict:
    root = DEV / "tests" / name
    return {
        "process": load(root / "metrics.json"),
        "artifacts": {
            "metrics": artifact(root / "metrics.json"),
            "stdout": artifact(root / "stdout"),
            "stderr": artifact(root / "stderr"),
        },
    }


def main() -> None:
    predecessor = load(PREDECESSOR)
    same_cell = process(DEV / "final-runs" / "n59-l9")
    cap_control = process(DEV / "final-runs" / "n59-l58")
    positive_control = process(DEV / "final-runs" / "n31-l5")
    build_root = DEV / "final-clean-build"
    build = load(build_root / "metrics.json")
    tests = {
        "unknown_scalar_consumers": test_receipt("unknown-scalar"),
        "production_consumers": test_receipt("production-consumers"),
    }
    formal_processes = [
        build["metrics"],
        same_cell["process"]["metrics"],
        cap_control["process"]["metrics"],
        positive_control["process"]["metrics"],
        tests["unknown_scalar_consumers"]["process"]["metrics"],
        tests["production_consumers"]["process"]["metrics"],
    ]
    formal = {
        "resource_components": len(formal_processes),
        "summed_wall_seconds": sum(row["wall_seconds"] for row in formal_processes),
        "total_core_seconds": sum(row["total_core_seconds"] for row in formal_processes),
        "peak_rss_bytes": max(row["peak_rss_bytes"] for row in formal_processes),
        "all_formal_components_outer_metered": True,
    }
    development = all_metered()
    previous_campaign = predecessor["campaign_accounting"]["measured_lower_bound"]
    same_report = same_cell["report"]
    cap_report = cap_control["report"]
    positive_report = positive_control["report"]

    result = {
        "schema": "koblitz_stage173_ggmp_same_cell_availability.v1",
        "date": "2026-09-23",
        "status": "complete_ggmp_same_cell_structurally_unavailable",
        "predecessor": {
            "path": PREDECESSOR.relative_to(GATES).as_posix(),
            "sha256": sha256(PREDECESSOR),
            "status": predecessor["status"],
        },
        "source_revision": {
            "commit": COMMIT,
            "tree": TREE,
            "dirty": False,
            "archive_sha256": ARCHIVE_SHA256,
            "frozen_cargo_lock_sha256": LOCK_SHA256,
            "discovery_binary_sha256": BINARY_SHA256,
        },
        "protocol": {
            "construction": "GGMP linearised-polynomial kernel from divisors of T^n-1",
            "same_cell": {"n": 59, "ell": 9, "m": 3, "curve_a": 1},
            "factor_enumeration_degree_cap": 24,
            "selection_rule": same_report["selection_rule"],
            "single_thread_requested": True,
            "forbidden_inputs": same_report["forbidden_inputs"],
        },
        "same_cell_availability": same_cell,
        "enumeration_cap_control": cap_control,
        "available_positive_control": positive_control,
        "validation": {
            "same_cell_complete_factor_degrees": same_report["complete_factor_degrees"],
            "same_cell_enumerated_factor_degrees": same_report["factor_degrees"],
            "same_cell_requested_divisor_exists": same_report["requested_divisor_exists"],
            "same_cell_terminal": same_report["status"],
            "cap_control_requested_divisor_exists": cap_report["requested_divisor_exists"],
            "cap_control_terminal": cap_report["status"],
            "positive_control_terminal": positive_report["status"],
            "positive_control_selected_divisor_indices": positive_report["selected"]["divisor_indices"],
            "consumer_tests": tests,
        },
        "clean_build": {
            "process": build,
            "artifacts": {
                "metrics": artifact(build_root / "metrics.json"),
                "stdout": artifact(build_root / "stdout"),
                "stderr": artifact(build_root / "stderr"),
                "script": artifact(DEV / "clean-build-c6205c56.sh"),
            },
            "cold_build_plus_same_cell_probe": {
                "wall_seconds": build["metrics"]["wall_seconds"]
                + same_cell["process"]["metrics"]["wall_seconds"],
                "total_core_seconds": build["metrics"]["total_core_seconds"]
                + same_cell["process"]["metrics"]["total_core_seconds"],
                "peak_rss_bytes": max(
                    build["metrics"]["peak_rss_bytes"],
                    same_cell["process"]["metrics"]["peak_rss_bytes"],
                ),
            },
        },
        "rejected_and_superseded": {
            "failed_missing_calibration_build": {
                "reason": "source archive omitted a compile-time calibration file; no binary produced",
                "metrics": artifact(DEV / "clean-build" / "failed-missing-calibration.metrics.json"),
                "stderr": artifact(DEV / "clean-build" / "failed-missing-calibration.stderr"),
            },
            "superseded_degree_reporting_build": {
                "reason": "enumerated degrees were not distinguished from the complete cyclotomic factor degrees",
                "metrics": artifact(DEV / "clean-build" / "metrics.json"),
                "n59_result": artifact(DEV / "clean-runs" / "n59-l9" / "stdout.json"),
            },
        },
        "formal_protocol_accounting": formal,
        "stage173_development_accounting": development,
        "campaign_accounting": {
            "complete_campaign_cost": None,
            "measured_lower_bound": {
                "resource_components": previous_campaign["resource_components"]
                + development["resource_components"],
                "summed_wall_seconds": previous_campaign["summed_wall_seconds"]
                + development["summed_wall_seconds"],
                "total_core_seconds": previous_campaign["total_core_seconds"]
                + development["total_core_seconds"],
                "peak_rss_bytes": max(previous_campaign["peak_rss_bytes"], development["peak_rss_bytes"]),
            },
            "stage172_measured_lower_bound": previous_campaign,
            "stage173_incremental_measured": {
                "resource_components": development["resource_components"],
                "summed_wall_seconds": development["summed_wall_seconds"],
                "total_core_seconds": development["total_core_seconds"],
                "peak_rss_bytes": development["peak_rss_bytes"],
            },
            "reason_complete_is_null": development["reason_complete_total_is_null"],
        },
        "gates": {
            "1_all_stage_resource_charging": "partial_historical_campaign_total_null_stage173_formal_protocol_complete",
            "2_same_instance_solver_matrix": "partial_native_f4_wdsat_cryptominisat_direct_mitm_present_ggmp_same_cell_proven_inapplicable_licensed_magma_missing",
            "3_single_core_core_memory_conflicts_wall": predecessor["gates"]["3_single_core_core_memory_conflicts_wall"],
            "4_n31_n41_larger_pdp": "inherited_satisfied_finite_coverage",
            "5_unknown_scalar_without_known_base_logs": "inherited_satisfied_finite_controls",
            "6_full_cost_vs_automorphism_rho": "unchanged_not_tested_by_pdp_stage_and_current_full_cost_gate_false",
            "7_independent_external_reproduction": "requested_but_missing",
            "all_seven_gates_passed": False,
        },
        "ggmp_same_cell_benchmark_applicable": False,
        "ggmp_same_cell_unavailability_reason": (
            "the complete irreducible-factor degrees of T^59-1 over F_2 are 1 and 58, "
            "so no divisor kernel has dimension ell=9"
        ),
        "licensed_magma_f4_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "claim_boundary": (
            "This proves that the GGMP divisor-kernel family has no factor base at the exact "
            "n=59, ell=9, m=3 cell and distinguishes that mathematical absence from the "
            "implementation's degree-24 factor-enumeration cap. It is not a GGMP speed, "
            "licensed-Magma result, full solver panel, crossover, novelty review, or SOTA."
        ),
    }
    (HERE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    metrics = same_cell["process"]["metrics"]
    build_metrics = build["metrics"]
    results = f"""# Stage 173: GGMP availability at the exact n59 cell

The GGMP divisor-kernel construction is **mathematically inapplicable** at the exact Stage-169/172 `n=59, ell=9, m=3, a=1` cell. This is a target-independent availability result, not a GGMP speed, full solver-panel completion, crossover, novelty review, or SOTA.

The complete 2-cyclotomic factor degrees of `T^59-1` over `F_2` are `[1,58]`. Their divisor dimensions are therefore `0,1,58,59`; dimension 9 does not exist. The clean probe reports `unavailable_no_divisor_of_requested_dimension`, `requested_divisor_exists=false`, no candidates, and `selected=null` in **{metrics['wall_seconds']:.6f} wall seconds / {metrics['total_core_seconds']:.6f} core-seconds / {metrics['peak_rss_bytes']} bytes peak RSS**.

The tool separately reports its operational enumeration cap. At `n=59, ell=58`, where the complete degree list proves a divisor exists, it returns `unavailable_factor_enumeration_degree_cap` because exhaustive factor enumeration is capped at degree 24. At `n=31, ell=5`, all complete factors are enumerable, six candidates are censused, and the public rule selects divisor index `[1]` as expected.

Every probe records that it constructed no target, enumerated no target subgroup, constructed no discrete-log labels, used no relation yield, and used no solver timing. The exact-commit cold build costs {build_metrics['wall_seconds']:.6f} wall seconds / {build_metrics['total_core_seconds']:.6f} core-seconds / {build_metrics['peak_rss_bytes']} bytes peak RSS; build plus the same-cell availability probe costs {build_metrics['wall_seconds'] + metrics['wall_seconds']:.6f} wall seconds.

This resolves the GGMP part of the same-cell comparison as structurally inapplicable rather than silently missing. Licensed Magma remains the unresolved comparator in gate 2.
"""
    (HERE / "RESULTS.md").write_text(results)


if __name__ == "__main__":
    main()

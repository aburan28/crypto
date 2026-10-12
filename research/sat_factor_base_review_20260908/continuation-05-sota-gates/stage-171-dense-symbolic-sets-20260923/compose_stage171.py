#!/usr/bin/env python3
"""Compose and seal the Stage 171 native-F4 single-target evidence."""

from __future__ import annotations

import copy
import hashlib
import json
import statistics
from pathlib import Path


HERE = Path(__file__).resolve().parent
GATES = HERE.parent
DEV = HERE / "development"
STAGE170 = GATES / "stage-170-parallel-fixed-x1-construction-20260923" / "result.json"
STAGE169 = GATES / "stage-169-fresh-target-holdout-20260923" / "result.json"
COMMIT = "6de26d8c7c1791e916f7a595440319f0a1a2c030"
BACKEND_SHA256 = "1c955e78d1885c400d731d2a69e669e096dc5fef4c06bcff9d0540273584dd15"
EXPORTER_SHA256 = "fed446b4c4483b79f34e2d79d1b42e57afe797dd270c94a2ecca0becfc4b4e8b"
SOURCE_ARCHIVE_SHA256 = "45a9c8e42311eebac690520983946f2a91a66b1b1371db72fce9747c9e9e9a68"
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


def run(name: str, group: str = "clean-pairs") -> dict:
    root = DEV / group / name
    process = load(root / "metrics.json")
    report = load(root / "stdout.json")
    metrics = copy.deepcopy(process["metrics"])
    # This process uses a Rayon pool.  process_meter's legacy field contains
    # total child CPU, not a valid one-thread measurement.
    metrics["single_core_seconds"] = None
    return {
        "name": name,
        "process": {**process, "metrics": metrics},
        "report": report,
        "artifacts": {
            "metrics": artifact(root / "metrics.json"),
            "stdout": artifact(root / "stdout.json"),
            "stderr": artifact(root / "stderr.txt"),
        },
    }


def ratio(numerator: float, denominator: float) -> float:
    return numerator / denominator


def stable_report(report: dict) -> dict:
    report = copy.deepcopy(report)
    report.pop("timing_ns", None)
    report.pop("solver_description", None)
    report.pop("solver_symbolic_monomial_sets", None)
    cost = report.get("cost", {})
    cost.pop("wall_ns", None)
    extra = cost.get("extra", {})
    for key in (
        "build_ns",
        "eliminate_ns",
        "pair_update_ns",
        "dense_symbolic_set_steps",
        "dense_symbolic_set_bytes_max",
        "m4ri_scratch_bytes_max",
    ):
        extra.pop(key, None)
    return report


def metered_development() -> dict:
    paths = sorted(set(DEV.rglob("metrics.json")) | set(DEV.rglob("*.metrics.json")))
    receipts = []
    for path in paths:
        if "development/clean-t14-b64/clean-t14-b64/" in path.relative_to(HERE).as_posix():
            continue
        try:
            record = load(path)
        except (OSError, json.JSONDecodeError):
            continue
        metrics = record.get("metrics", record)
        if not all(key in metrics for key in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")):
            continue
        receipts.append((path, metrics))
    return {
        "metered_components": len(receipts),
        "summed_wall_seconds": sum(metrics["wall_seconds"] for _, metrics in receipts),
        "total_core_seconds": sum(metrics["total_core_seconds"] for _, metrics in receipts),
        "maximum_peak_rss_bytes": max(metrics["peak_rss_bytes"] for _, metrics in receipts),
        "complete_total": None,
        "reason_complete_total_is_null": (
            "exploratory compiles and tests were not all outer-metered, and one early profiler "
            "receipt was overwritten before archival; the retained sum is a measured lower bound"
        ),
    }


def compact_arm(group: str, name: str) -> dict:
    item = run(name, group)
    report = item.pop("report")
    item["status"] = report["status"]
    item["equation_fingerprint"] = report["solver_equations_blake3"]
    item["cost"] = report["cost"]
    return item


def compact_flat(name: str) -> dict:
    root = DEV / name
    process = load(root / "metrics.json")
    report = load(root / "stdout.json")
    process["metrics"]["single_core_seconds"] = None
    return {
        "name": name,
        "process": process,
        "status": report["status"],
        "equation_fingerprint": report["solver_equations_blake3"],
        "cost": report["cost"],
        "artifacts": {
            "metrics": artifact(root / "metrics.json"),
            "stdout": artifact(root / "stdout.json"),
            "stderr": artifact(root / "stderr.txt"),
        },
    }


def main() -> None:
    predecessor = load(STAGE170)
    stage169 = load(STAGE169)
    selected = run("candidate-r1")
    control = run("control-r1")
    candidates = [run(f"candidate-r{i}") for i in (1, 2, 3)]
    controls = [run(f"control-r{i}") for i in (1, 2, 3)]

    paired_ratios = []
    for index, (candidate, reference) in enumerate(zip(candidates, controls), 1):
        cm = candidate["process"]["metrics"]
        rm = reference["process"]["metrics"]
        cx = candidate["report"]["cost"]["extra"]
        rx = reference["report"]["cost"]["extra"]
        paired_ratios.append(
            {
                "pair": index,
                "candidate_over_control_wall": ratio(cm["wall_seconds"], rm["wall_seconds"]),
                "candidate_over_control_core": ratio(cm["total_core_seconds"], rm["total_core_seconds"]),
                "candidate_over_control_rss": ratio(cm["peak_rss_bytes"], rm["peak_rss_bytes"]),
                "candidate_over_control_build_ns": ratio(cx["build_ns"], rx["build_ns"]),
                "candidate_over_control_eliminate_ns": ratio(
                    cx["eliminate_ns"], rx["eliminate_ns"]
                ),
            }
        )
    paired_medians = {
        key: statistics.median(pair[key] for pair in paired_ratios)
        for key in paired_ratios[0]
        if key != "pair"
    }

    selected_metrics = selected["process"]["metrics"]
    selected_report = selected["report"]
    control_metrics = control["process"]["metrics"]
    control_report = control["report"]
    old_metrics = predecessor["selected_native_f4"]["metrics"]
    posthoc_metrics = predecessor["posthoc_single_target_ceiling"]["metrics"]
    direct = stage169["controls"]["direct-mitm"]["process"]["metrics"]
    build_process = load(DEV / "clean-build" / "metrics.json")
    development = metered_development()
    predecessor_campaign = predecessor["campaign_accounting"]["measured_lower_bound"]

    same_science = all(
        stable_report(item["report"]) == stable_report(controls[0]["report"])
        for item in candidates + controls
    )
    expected_fingerprint = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
    exact_terminals = all(
        item["report"]["status"] == "unsat"
        and item["report"]["exhaustive"] is True
        and item["report"]["solver_equations_blake3"] == expected_fingerprint
        for item in candidates + controls
    )

    stage_campaign_wall = development["summed_wall_seconds"]
    stage_campaign_core = development["total_core_seconds"]
    result = {
        "schema": "koblitz_stage171_dense_symbolic_sets.v1",
        "date": "2026-09-23",
        "status": "complete_dense_symbolic_sets_true_negative",
        "predecessor": {
            "path": STAGE170.relative_to(GATES).as_posix(),
            "sha256": sha256(STAGE170),
            "status": predecessor["status"],
        },
        "source_revision": {
            "commit": COMMIT,
            "tree": "9e5b4a5e2ca11e1b055a69a71e1d8656736e37f4",
            "dirty": False,
            "archive_sha256": SOURCE_ARCHIVE_SHA256,
            "frozen_cargo_lock_sha256": LOCK_SHA256,
            "backend_sha256": BACKEND_SHA256,
            "exporter_sha256": EXPORTER_SHA256,
        },
        "instance": predecessor["instance"],
        "factor_base": predecessor["factor_base"],
        "mechanism": {
            "name": "reused_bit_packed_symbolic_monomial_sets",
            "definition": (
                "For at most twenty Boolean variables, retain one bit per complete-domain "
                "monomial plus an insertion-order vector for the F4 lcm, examined, and "
                "no-divisor sets; clear and reuse them across F4 steps. Larger systems retain "
                "the exact SplitMix64 hash-set path."
            ),
            "same_binary_control": "PQ_F4_DISABLE_DENSE_SYMBOLIC_SET=1",
            "selected_rayon_threads": 12,
            "selected_x1_batch_size": 512,
            "selected_rational_systems": 242,
            "selected_batches": 1,
            "batch_cap_changed_from": 64,
            "batch_cap_changed_to": "2^ell",
            "uses_truth_or_witness": False,
            "uses_discrete_log_labels": False,
        },
        "validation": {
            "boolean_f4_tests": {
                "passed": 13,
                "failed": 0,
                "artifacts": {
                    "metrics": artifact(DEV / "final-tests" / "f4.metrics.json"),
                    "stdout": artifact(DEV / "final-tests" / "f4.stdout"),
                    "stderr": artifact(DEV / "final-tests" / "f4.stderr"),
                },
            },
            "fixed_x1_backend_tests": {
                "passed": 3,
                "failed": 0,
                "artifacts": {
                    "metrics": artifact(DEV / "final-tests" / "backend.metrics.json"),
                    "stdout": artifact(DEV / "final-tests" / "backend.stdout"),
                    "stderr": artifact(DEV / "final-tests" / "backend.stderr"),
                },
            },
            "same_scientific_report_outside_declared_mode_allocation_and_timing": same_science,
            "all_clean_pair_terminals_exact": exact_terminals,
        },
        "clean_build": {
            "process": build_process,
            "artifacts": {
                "metrics": artifact(DEV / "clean-build" / "metrics.json"),
                "stdout": artifact(DEV / "clean-build" / "stdout"),
                "stderr": artifact(DEV / "clean-build" / "stderr"),
                "script": artifact(DEV / "clean-build-6de26d8c.sh"),
            },
            "failed_attempts": [
                {
                    "reason": "offline Cargo home lacked the registry index; no binary produced",
                    "metrics": artifact(DEV / "clean-build" / "failed-default-cache.metrics.json"),
                    "stderr": artifact(DEV / "clean-build" / "failed-default-cache.stderr"),
                },
                {
                    "reason": "Cargo ran outside the extracted source and did not discover vendor config",
                    "metrics": artifact(DEV / "clean-build" / "failed-vendor-wrong-cwd.metrics.json"),
                    "stderr": artifact(DEV / "clean-build" / "failed-vendor-wrong-cwd.stderr"),
                },
            ],
        },
        "selected_native_f4": {
            "decision": "first_clean_candidate_at_preselected_12_threads_one_complete_batch",
            "classification": "true_negative",
            "process": selected["process"],
            "metrics": selected_metrics,
            "report": selected_report,
            "artifacts": selected["artifacts"],
            "conflicts": None,
            "valid_single_core_seconds": None,
        },
        "same_binary_hash_control": {
            "decision": "first_clean_control_same_binary_and_schedule",
            "process": control["process"],
            "metrics": control_metrics,
            "report": control_report,
            "artifacts": control["artifacts"],
            "conflicts": None,
            "valid_single_core_seconds": None,
        },
        "clean_paired_repeats": {
            "candidates": candidates,
            "controls": controls,
            "pair_ratios": paired_ratios,
            "median_pair_ratios": paired_medians,
        },
        "comparisons": {
            "selected_over_same_binary_hash_control": {
                "wall_ratio": ratio(selected_metrics["wall_seconds"], control_metrics["wall_seconds"]),
                "core_ratio": ratio(
                    selected_metrics["total_core_seconds"], control_metrics["total_core_seconds"]
                ),
                "rss_ratio": ratio(
                    selected_metrics["peak_rss_bytes"], control_metrics["peak_rss_bytes"]
                ),
                "build_ns_ratio": ratio(
                    selected_report["cost"]["extra"]["build_ns"],
                    control_report["cost"]["extra"]["build_ns"],
                ),
                "same_scientific_report": same_science,
            },
            "paired_medians_candidate_over_control": paired_medians,
            "selected_over_stage170_selected": {
                "wall_speedup": ratio(old_metrics["wall_seconds"], selected_metrics["wall_seconds"]),
                "core_ratio": ratio(
                    selected_metrics["total_core_seconds"], old_metrics["total_core_seconds"]
                ),
                "rss_ratio": ratio(
                    selected_metrics["peak_rss_bytes"], old_metrics["peak_rss_bytes"]
                ),
            },
            "selected_over_stage170_posthoc_ceiling": {
                "wall_ratio": ratio(
                    selected_metrics["wall_seconds"], posthoc_metrics["wall_seconds"]
                ),
                "core_ratio": ratio(
                    selected_metrics["total_core_seconds"], posthoc_metrics["total_core_seconds"]
                ),
                "rss_ratio": ratio(
                    selected_metrics["peak_rss_bytes"], posthoc_metrics["peak_rss_bytes"]
                ),
                "interpretation": "the clean selected median-scale result does not beat the earlier one-off posthoc wall sample",
            },
            "selected_over_direct_mitm": {
                "wall_ratio": ratio(selected_metrics["wall_seconds"], direct["wall_seconds"]),
                "core_ratio": ratio(
                    selected_metrics["total_core_seconds"], direct["total_core_seconds"]
                ),
                "rss_ratio": ratio(selected_metrics["peak_rss_bytes"], direct["peak_rss_bytes"]),
            },
            "cold_build_plus_selected": {
                "wall_seconds": build_process["metrics"]["wall_seconds"]
                + selected_metrics["wall_seconds"],
                "total_core_seconds": build_process["metrics"]["total_core_seconds"]
                + selected_metrics["total_core_seconds"],
                "peak_rss_bytes": max(
                    build_process["metrics"]["peak_rss_bytes"],
                    selected_metrics["peak_rss_bytes"],
                ),
            },
        },
        "rejected_development": {
            "radix_product_sort": {
                "reason": "higher wall and core cost",
                "candidate": compact_flat("radix-candidate-r1"),
                "control": compact_flat("radix-control-r1"),
            },
            "dense_reducer_index": {
                "reason": "two MiB per concurrent solve and higher core cost",
                "candidate": compact_flat("dense-reducer-candidate-r1"),
                "control": compact_flat("dense-reducer-control-r1"),
            },
            "encoded_column_sort": {
                "reason": "saved no CPU after key-buffer construction",
                "candidate": compact_flat("encoded-columns-candidate-r1"),
                "control": compact_flat("encoded-columns-control-r1"),
            },
            "packed_seen_rows": {
                "reason": "sub-percent paired CPU effect and no clear whole-run gain; reverted for simplicity",
                "candidate_runs": [compact_flat(f"packed-seen-candidate-r{i}") for i in (1, 2, 3)],
                "control_runs": [compact_flat(f"packed-seen-control-r{i}") for i in (1, 2, 3)],
            },
            "fourteen_thread_batch64": {
                "reason": "repeatable wall stalls despite comparable process CPU",
                "candidate_runs": [compact_arm("clean-t14-b64", f"candidate-r{i}") for i in (1, 2)],
                "control": compact_arm("clean-t14-b64", "control-r1"),
            },
        },
        "development_accounting": development,
        "campaign_accounting": {
            "complete_campaign_cost": None,
            "measured_lower_bound": {
                "resource_components": predecessor_campaign["resource_components"]
                + development["metered_components"],
                "summed_wall_seconds": predecessor_campaign["summed_wall_seconds"]
                + stage_campaign_wall,
                "total_core_seconds": predecessor_campaign["total_core_seconds"]
                + stage_campaign_core,
                "peak_rss_bytes": max(
                    predecessor_campaign["peak_rss_bytes"],
                    development["maximum_peak_rss_bytes"],
                ),
            },
            "stage170_measured_lower_bound": predecessor_campaign,
            "stage171_incremental_measured": {
                "resource_components": development["metered_components"],
                "summed_wall_seconds": stage_campaign_wall,
                "total_core_seconds": stage_campaign_core,
                "peak_rss_bytes": development["maximum_peak_rss_bytes"],
            },
            "reason_complete_is_null": development["reason_complete_total_is_null"],
            "unmetered_scope": [
                "some exploratory compilation and focused tests",
                "one overwritten early profiler receipt",
                "score composition, documentation and Git operations",
            ],
        },
        "gates": {
            "1_all_stage_resource_charging": "partial_measured_lower_bound_complete_total_null",
            "2_same_instance_solver_matrix": "partial_two_native_f4_n59_targets_one_preregistered_holdout_licensed_magma_and_full_panel_missing",
            "3_single_core_core_memory_conflicts_wall": "partial_native_f4_complete_conflicts_null_magma_missing",
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
            "Same-target engineering after the Stage-169 truth was opened. It improves the "
            "repository-native F4 symbolic-preprocessing path and permits one complete fixed-X1 "
            "batch. It is not a fresh-target distribution, complete index calculus, rho crossover, "
            "independent reproduction, novelty review, or SOTA."
        ),
    }

    (HERE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    results_md = f"""# Stage 171: dense symbolic monomial sets in native F4

Same-target engineering after the Stage-169 truth was opened. This stage improves the repository-native Boolean F4 arm; it is not a fresh-target distribution, complete index calculus, rho crossover, independent reproduction, novelty review, or SOTA.

For systems with at most twenty Boolean variables, F4 now represents the symbolic-preprocessing LCM, examined, and no-divisor sets with one bit per complete-domain monomial plus an insertion-order vector. The allocations are reused across steps. `PQ_F4_DISABLE_DENSE_SYMBOLIC_SET=1` restores the same binary's hash-set path. The fixed-X1 batch cap now reaches the complete public `2^ell` domain, while the default remains one and `PQ_F4_X1_BATCH=1` is the serial scheduling control.

The first clean selected run uses twelve Rayon threads and one complete batch of 242 rational systems. It finishes the known Stage-169 true-negative target in {selected_metrics['wall_seconds']:.6f} wall seconds / {selected_metrics['total_core_seconds']:.6f} total core-seconds / {selected_metrics['peak_rss_bytes']} bytes peak RSS. It visits all 512 masks, preserves 14,278 equations and 14,515,915 terms, performs 242 F4 calls and 99,199,976,264 charged word XORs, finds zero roots, and returns exhaustive UNSAT with equation fingerprint `{expected_fingerprint}`.

Across three interleaved clean pairs, the median candidate/control ratios are {paired_medians['candidate_over_control_wall']:.6f} wall, {paired_medians['candidate_over_control_core']:.6f} CPU, {paired_medians['candidate_over_control_rss']:.6f} RSS, and {paired_medians['candidate_over_control_build_ns']:.6f} F4 matrix-build time. The scientific reports agree after removing the declared representation, allocation, and timing counters. The selected run is {old_metrics['wall_seconds']/selected_metrics['wall_seconds']:.3f}x faster in wall time than Stage 170's selected run, but it does not beat Stage 170's one-off 14.661977-second post-hoc wall sample.

The exact-commit clean archive/build costs {build_process['metrics']['wall_seconds']:.6f} wall seconds, {build_process['metrics']['total_core_seconds']:.6f} core-seconds, and {build_process['metrics']['peak_rss_bytes']} bytes peak RSS. Build plus selected execution costs {build_process['metrics']['wall_seconds'] + selected_metrics['wall_seconds']:.6f} wall seconds and {build_process['metrics']['total_core_seconds'] + selected_metrics['total_core_seconds']:.6f} core-seconds. Two failed offline setup attempts are preserved and charged.

Direct MITM remains {selected_metrics['wall_seconds']/direct['wall_seconds']:.2f}x faster by wall, {selected_metrics['total_core_seconds']/direct['total_core_seconds']:.2f}x cheaper by CPU, and {selected_metrics['peak_rss_bytes']/direct['peak_rss_bytes']:.2f}x smaller by RSS on this target. The broader campaign has a {development['metered_components']}-component Stage-171 measured lower bound, while complete development cost remains null because some exploratory compilation/test work was not outer-metered and one early profile receipt was overwritten.
"""
    (HERE / "RESULTS.md").write_text(results_md)


if __name__ == "__main__":
    main()

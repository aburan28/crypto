#!/usr/bin/env python3
"""Compose the Stage 174 size-gated contiguous-M4RI evidence."""

from __future__ import annotations

import copy
import hashlib
import json
import math
import statistics
from pathlib import Path


HERE = Path(__file__).resolve().parent
GATES = HERE.parent
DEV = HERE / "development"
FLAT = DEV / "flat-campaign"
STAGE171 = GATES / "stage-171-dense-symbolic-sets-20260923" / "result.json"
STAGE172 = GATES / "stage-172-native-f4-single-core-20260923" / "result.json"
STAGE173 = GATES / "stage-173-ggmp-same-cell-availability-20260923" / "result.json"
COMMIT = "e51efb219edd4511a08d545de21184928e5e771c"
TREE = "321c58f585735c3e2c1d430e4c39a2555b654838"
BACKEND_SHA256 = "74593e230909b899feb09c7f2e1fa3695c87da6673c80b7457f3e6724b74ca6e"
EXPORTER_SHA256 = "0405cc55bc9b78c9e86de5755016046d2d0db2402f3949ee84d8ff9ff3603c46"
SOURCE_ARCHIVE_SHA256 = "71a9ce55473225ff5cc27e08903240b651d5d98ddf06d439da0d6233fa8cbe80"
LOCK_SHA256 = "10461fc25d189a7e1b533c0a59e1c66d7376ca04d3af5c5868a6f155efc8e5cf"
EQUATION_FINGERPRINT = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"


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


def run(name: str, *, parallel: bool = True) -> dict:
    root = FLAT / name
    process = load(root / "metrics.json")
    report = load(root / "stdout.json")
    normalized = copy.deepcopy(process)
    if parallel:
        normalized["metrics"]["single_core_seconds"] = None
    return {
        "name": name,
        "process": normalized,
        "report": report,
        "artifacts": {
            "metrics": artifact(root / "metrics.json"),
            "stdout": artifact(root / "stdout.json"),
            "stderr": artifact(root / "stderr.txt"),
        },
    }


def metric(item: dict, key: str) -> float:
    return item["process"]["metrics"][key]


def ratio(numerator: float, denominator: float) -> float:
    return numerator / denominator


def stable_report(report: dict) -> dict:
    value = copy.deepcopy(report)
    for key in (
        "timing_ns",
        "solver_description",
        "solver_echelon_policy",
        "solver_flat_m4ri_control",
        "solver_flat_m4ri_min_words",
    ):
        value.pop(key, None)
    cost = value.get("cost", {})
    cost.pop("wall_ns", None)
    extra = cost.get("extra", {})
    for key in (
        "build_ns",
        "eliminate_ns",
        "pair_update_ns",
        "m4ri_scratch_bytes_max",
        "flat_m4ri_matrices",
        "flat_m4ri_order_bytes_max",
    ):
        extra.pop(key, None)
    return value


def paired(candidate: dict, control: dict, index: int) -> dict:
    cm = candidate["process"]["metrics"]
    rm = control["process"]["metrics"]
    cx = candidate["report"]["cost"]["extra"]
    rx = control["report"]["cost"]["extra"]
    return {
        "pair": index,
        "candidate_over_control_wall": ratio(cm["wall_seconds"], rm["wall_seconds"]),
        "candidate_over_control_core": ratio(
            cm["total_core_seconds"], rm["total_core_seconds"]
        ),
        "candidate_over_control_rss": ratio(
            cm["peak_rss_bytes"], rm["peak_rss_bytes"]
        ),
        "candidate_over_control_build_ns": ratio(cx["build_ns"], rx["build_ns"]),
        "candidate_over_control_eliminate_ns": ratio(
            cx["eliminate_ns"], rx["eliminate_ns"]
        ),
    }


def median_metrics(items: list[dict]) -> dict:
    return {
        key: statistics.median(metric(item, key) for item in items)
        for key in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
    }


def metered_development() -> dict:
    receipts = []
    for path in sorted(DEV.rglob("*")):
        if not path.is_file() or not (
            path.name == "metrics.json"
            or path.name.endswith(".metrics.json")
            or path.name.endswith(".metrics")
        ):
            continue
        try:
            record = load(path)
        except (OSError, json.JSONDecodeError):
            continue
        metrics = record.get("metrics", record)
        if all(
            key in metrics
            for key in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
        ):
            receipts.append((path, metrics))
    return {
        "resource_components": len(receipts),
        "summed_wall_seconds": math.fsum(m["wall_seconds"] for _, m in receipts),
        "total_core_seconds": math.fsum(m["total_core_seconds"] for _, m in receipts),
        "peak_rss_bytes": max(m["peak_rss_bytes"] for _, m in receipts),
        "complete_total": None,
        "reason_complete_total_is_null": (
            "some exploratory incremental compiles and tests were not outer-metered; all "
            "retained solver runs, clean-build attempts, PGO steps, and final tests are charged"
        ),
    }


def summarize_codegen() -> dict:
    root = DEV / "codegen"
    generic = [
        load(path)["metrics"]
        for path in sorted(root.glob("pair*-*-generic/metrics.json"))
    ]
    native = [
        load(path)["metrics"]
        for path in sorted(root.glob("pair*-*-native-lto/metrics.json"))
    ]
    pgo = [
        load(path)["metrics"] for path in sorted(root.glob("pgo-pair*-*-pgo/metrics.json"))
    ]
    pgo_controls = [
        load(path)["metrics"]
        for path in sorted(root.glob("pgo-pair*-*-native-lto/metrics.json"))
    ]
    return {
        "native_lto": {
            "generic_runs": len(generic),
            "native_runs": len(native),
            "native_over_generic_wall_median": ratio(
                statistics.median(m["wall_seconds"] for m in native),
                statistics.median(m["wall_seconds"] for m in generic),
            ),
            "native_over_generic_core_median": ratio(
                statistics.median(m["total_core_seconds"] for m in native),
                statistics.median(m["total_core_seconds"] for m in generic),
            ),
            "reason_not_selected": (
                "small CPU gain, larger cold-build cost, and host-specific target-cpu native"
            ),
        },
        "pgo": {
            "pgo_runs": len(pgo),
            "native_lto_controls": len(pgo_controls),
            "pgo_over_native_lto_wall_median": ratio(
                statistics.median(m["wall_seconds"] for m in pgo),
                statistics.median(m["wall_seconds"] for m in pgo_controls),
            ),
            "pgo_over_native_lto_core_median": ratio(
                statistics.median(m["total_core_seconds"] for m in pgo),
                statistics.median(m["total_core_seconds"] for m in pgo_controls),
            ),
            "llvm23_merge_failed": load(
                DEV / "pgo/merge/metrics.json"
            )["returncode"]
            != 0,
            "llvm21_merge_succeeded": load(
                DEV / "pgo/merge-llvm21/metrics.json"
            )["returncode"]
            == 0,
            "reason_not_selected": "profile-use runs were slower than native LTO controls",
        },
    }


def summarize_variable_order() -> dict:
    root = DEV / "variable-order"
    identities = [
        load(path)["metrics"]
        for path in sorted(root.glob("pair*-*-identity/metrics.json"))
    ]
    reverses = [
        load(path)["metrics"]
        for path in sorted(root.glob("pair*-*-reverse/metrics.json"))
    ]
    random_reports = [
        load(path) for path in sorted(root.glob("random*/stdout.json"))
    ]
    return {
        "identity_runs": len(identities),
        "reverse_runs": len(reverses),
        "reverse_over_identity_wall_median": ratio(
            statistics.median(m["wall_seconds"] for m in reverses),
            statistics.median(m["wall_seconds"] for m in identities),
        ),
        "reverse_over_identity_core_median": ratio(
            statistics.median(m["total_core_seconds"] for m in reverses),
            statistics.median(m["total_core_seconds"] for m in identities),
        ),
        "random_orders": len(random_reports),
        "minimum_random_word_xors": min(report["cost"]["ops"] for report in random_reports),
        "identity_word_xors": 99_199_976_264,
        "reason_not_selected": (
            "reverse had no stable CPU gain and every random mixed order increased F4 word XORs"
        ),
    }


def main() -> None:
    stage171 = load(STAGE171)
    stage172 = load(STAGE172)
    stage173 = load(STAGE173)

    candidates = [
        run("clean-pairs-01-flat-1048576"),
        run("clean-pairs-02-flat-1048576"),
        run("clean-pairs-05-flat-1048576"),
    ]
    controls = [
        run("clean-pairs-00-segmented"),
        run("clean-pairs-03-segmented"),
        run("clean-pairs-04-segmented"),
    ]
    pairs = [paired(candidate, control, i) for i, (candidate, control) in enumerate(zip(candidates, controls), 1)]
    pair_medians = {
        key: statistics.median(item[key] for item in pairs)
        for key in pairs[0]
        if key != "pair"
    }
    candidate_medians = median_metrics(candidates)
    control_medians = median_metrics(controls)
    selected = candidates[2]
    selected_control = controls[2]

    directs = [run(f"clean-direct-r{i}", parallel=False) for i in (1, 2, 3)]
    direct_medians = median_metrics(directs)
    single_core = {
        "segmented": run("unsat-single-core-segmented", parallel=False),
        "hybrid": run("unsat-single-core-hybrid", parallel=False),
        "all_flat": run("unsat-single-core-flat", parallel=False),
    }
    sat_witness = {
        "segmented": run("sat-witness-control", parallel=False),
        "all_flat": run("sat-witness", parallel=False),
    }

    all_reports = [item["report"] for item in candidates + controls]
    same_science = all(
        stable_report(report) == stable_report(all_reports[0]) for report in all_reports[1:]
    )
    exact = all(
        report["status"] == "unsat"
        and report["exhaustive"] is True
        and report["fixed_x1_values"] == 512
        and report["fixed_x1_masks_visited"] == 512
        and report["fixed_x1_systems_completed"] == 242
        and report["solver_equations_blake3"] == EQUATION_FINGERPRINT
        and report["cost"]["ops"] == 99_199_976_264
        for report in all_reports
    )
    layout_exact = all(
        item["report"]["cost"]["extra"]["flat_m4ri_matrices"] == 241
        and item["report"]["solver_flat_m4ri_min_words"] == 1_048_576
        for item in candidates
    ) and all(
        item["report"]["cost"]["extra"]["flat_m4ri_matrices"] == 0
        and item["report"]["solver_flat_m4ri_min_words"] is None
        for item in controls
    )

    incremental = metered_development()
    previous_campaign = stage173["campaign_accounting"]["measured_lower_bound"]
    campaign = {
        "resource_components": previous_campaign["resource_components"]
        + incremental["resource_components"],
        "summed_wall_seconds": previous_campaign["summed_wall_seconds"]
        + incremental["summed_wall_seconds"],
        "total_core_seconds": previous_campaign["total_core_seconds"]
        + incremental["total_core_seconds"],
        "peak_rss_bytes": max(previous_campaign["peak_rss_bytes"], incremental["peak_rss_bytes"]),
    }

    selected_metrics = selected["process"]["metrics"]
    direct_wall = direct_medians["wall_seconds"]
    direct_core = direct_medians["total_core_seconds"]
    direct_rss = direct_medians["peak_rss_bytes"]
    result = {
        "schema": "koblitz_stage174_size_gated_flat_m4ri.v1",
        "date": "2026-09-23",
        "status": "complete_same_target_parallel_f4_improvement",
        "claim_boundary": (
            "Same-target implementation engineering on the already-opened public n=59 true-negative. "
            "It is not a fresh holdout, a complete index-calculus run, a rho crossover, independent "
            "reproduction, novelty review, or Koblitz index-calculus SOTA."
        ),
        "predecessors": {
            "stage171": {"path": STAGE171.relative_to(GATES).as_posix(), "sha256": sha256(STAGE171)},
            "stage172": {"path": STAGE172.relative_to(GATES).as_posix(), "sha256": sha256(STAGE172)},
            "stage173": {"path": STAGE173.relative_to(GATES).as_posix(), "sha256": sha256(STAGE173)},
        },
        "source_revision": {
            "commit": COMMIT,
            "tree": TREE,
            "dirty": False,
            "source_archive_sha256": SOURCE_ARCHIVE_SHA256,
            "backend_sha256": BACKEND_SHA256,
            "exporter_sha256": EXPORTER_SHA256,
            "frozen_cargo_lock_sha256": LOCK_SHA256,
        },
        "instance": stage172["instance"],
        "factor_base": stage172["factor_base"],
        "mechanism": {
            "name": "size_gated_contiguous_m4ri_rows",
            "definition": (
                "Pack eligible Boolean F4 matrices into one contiguous u64 arena, retain a u32 "
                "logical-to-physical row order for O(1) pivot swaps, and keep the arena alive "
                "through pivot extraction. Matrices below the packed-word gate retain segmented rows."
            ),
            "selected_environment": {
                "PQ_F4_FLAT_M4RI": "1",
                "PQ_F4_FLAT_M4RI_MIN_WORDS": "1048576",
                "RAYON_NUM_THREADS": "12",
                "PQ_F4_X1_BATCH": "512",
            },
            "packed_matrix_min_bytes": 8_388_608,
            "selected_flat_matrices": 241,
            "total_m4ri_matrices": 723,
            "segmented_control": "omit PQ_F4_FLAT_M4RI",
            "default_behavior_changed": False,
            "uses_truth_or_witness": False,
            "uses_discrete_log_labels": False,
        },
        "validation": {
            "boolean_f4_tests": {
                "passed": 13,
                "failed": 0,
                "artifacts": {
                    suffix: artifact(FLAT / "final-tests-f4" / filename)
                    for suffix, filename in (
                        ("metrics", "metrics.json"),
                        ("stdout", "stdout.txt"),
                        ("stderr", "stderr.txt"),
                    )
                },
            },
            "backend_tests_selected_environment": {
                "passed": 3,
                "failed": 0,
                "artifacts": {
                    suffix: artifact(FLAT / "final-tests-backend" / filename)
                    for suffix, filename in (
                        ("metrics", "metrics.json"),
                        ("stdout", "stdout.txt"),
                        ("stderr", "stderr.txt"),
                    )
                },
            },
            "random_matrix_equivalence": (
                "streaming, segmented M4RI, and flat M4RI have identical pivot columns, charged "
                "XOR counts, and canonical row spaces on four deterministic shapes"
            ),
            "all_clean_terminals_exact": exact,
            "selected_layout_counts_exact": layout_exact,
            "same_scientific_report_outside_declared_layout_and_timing": same_science,
            "sat_witness_control": sat_witness,
        },
        "clean_build": {
            "process": load(FLAT / "clean-build/metrics.json"),
            "script": artifact(DEV / "clean-build-e51efb21.sh"),
            "artifacts": {
                suffix: artifact(FLAT / "clean-build" / filename)
                for suffix, filename in (
                    ("metrics", "metrics.json"),
                    ("stdout", "stdout.txt"),
                    ("stderr", "stderr.txt"),
                )
            },
            "failed_offline_vendor_attempt": {
                "reason": "default offline Cargo cache lacked redis 0.32.7",
                "metrics": artifact(DEV / "failed-clean-vendor/vendor.metrics.json"),
                "stderr": artifact(DEV / "failed-clean-vendor/vendor.stderr"),
            },
        },
        "selected_native_f4": {
            "decision": "median paired candidate by wall and CPU proximity",
            **selected,
            "conflicts": None,
            "valid_single_core_seconds": None,
        },
        "matched_segmented_control": {
            **selected_control,
            "conflicts": None,
            "valid_single_core_seconds": None,
        },
        "clean_paired_repeats": {
            "candidates": candidates,
            "controls": controls,
            "pair_ratios": pairs,
            "median_pair_ratios": pair_medians,
            "candidate_median_metrics": candidate_medians,
            "control_median_metrics": control_medians,
        },
        "single_core_ablation": {
            **single_core,
            "hybrid_over_segmented": {
                "wall_ratio": ratio(
                    metric(single_core["hybrid"], "wall_seconds"),
                    metric(single_core["segmented"], "wall_seconds"),
                ),
                "core_ratio": ratio(
                    metric(single_core["hybrid"], "total_core_seconds"),
                    metric(single_core["segmented"], "total_core_seconds"),
                ),
                "rss_ratio": ratio(
                    metric(single_core["hybrid"], "peak_rss_bytes"),
                    metric(single_core["segmented"], "peak_rss_bytes"),
                ),
                "interpretation": "neutral/slightly slower; keep segmented for single-core",
            },
        },
        "same_binary_direct_mitm": {
            "runs": directs,
            "median_metrics": direct_medians,
        },
        "comparisons": {
            "paired_median_candidate_over_control": pair_medians,
            "selected_over_direct_mitm_median": {
                "wall_ratio": selected_metrics["wall_seconds"] / direct_wall,
                "core_ratio": selected_metrics["total_core_seconds"] / direct_core,
                "rss_ratio": selected_metrics["peak_rss_bytes"] / direct_rss,
            },
            "selected_over_stage171_selected": {
                "wall_ratio": selected_metrics["wall_seconds"]
                / stage171["selected_native_f4"]["metrics"]["wall_seconds"],
                "core_ratio": selected_metrics["total_core_seconds"]
                / stage171["selected_native_f4"]["metrics"]["total_core_seconds"],
                "rss_ratio": selected_metrics["peak_rss_bytes"]
                / stage171["selected_native_f4"]["metrics"]["peak_rss_bytes"],
                "interpretation": (
                    "CPU and RSS improve, but the current host's heavily descheduled wall time is "
                    "not directly comparable to the earlier snapshot"
                ),
            },
        },
        "rejected_development": {
            "variable_order": summarize_variable_order(),
            "code_generation": summarize_codegen(),
            "all_flat_single_core": {
                "candidate_over_segmented_core": ratio(
                    metric(single_core["all_flat"], "total_core_seconds"),
                    metric(single_core["segmented"], "total_core_seconds"),
                ),
                "candidate_over_segmented_rss": ratio(
                    metric(single_core["all_flat"], "peak_rss_bytes"),
                    metric(single_core["segmented"], "peak_rss_bytes"),
                ),
                "reason_not_selected": "regressed the exhaustive single-core true-negative",
            },
        },
        "campaign_accounting": {
            "stage174_incremental_measured_lower_bound": incremental,
            "measured_lower_bound_through_stage174": campaign,
            "complete_campaign_cost": None,
            "reason_complete_is_null": incremental["reason_complete_total_is_null"],
        },
        "gates": {
            **stage173["gates"],
            "1_all_stage_resource_charging": (
                "partial_historical_campaign_total_null_stage174_retained_components_charged"
            ),
            "3_single_core_core_memory_conflicts_wall": (
                "partial_native_f4_parallel_and_single_core_complete_conflicts_null_magma_missing"
            ),
            "all_seven_gates_passed": False,
        },
        "fresh_target_holdout_complete": False,
        "licensed_magma_f4_complete": False,
        "full_cost_gate_passed": False,
        "independent_external_reproduction_satisfied": False,
        "koblitz_index_calculus_sota": False,
    }

    (HERE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    median = result["clean_paired_repeats"]["median_pair_ratios"]
    selected_direct = result["comparisons"]["selected_over_direct_mitm_median"]
    text = f"""# Stage 174: size-gated contiguous rows in native Boolean F4

This stage keeps the single-target work on the repository-native Boolean F4 backend. It adds an opt-in contiguous row arena for large M4RI matrices while retaining the segmented representation below an explicit packed-word threshold. The exact public target remains the already-opened `n=59, ell=9, m=3` true-negative, so this is same-target engineering rather than a fresh holdout.

The selected setting is `PQ_F4_FLAT_M4RI=1` with `PQ_F4_FLAT_M4RI_MIN_WORDS=1048576`, twelve Rayon workers, and one complete 512-value X1 batch. It routes 241 of 723 M4RI matrices through the contiguous arena. Across three interleaved exact-commit pairs, the median candidate/control ratios are `{median['candidate_over_control_wall']:.6f}` wall, `{median['candidate_over_control_core']:.6f}` CPU, and `{median['candidate_over_control_rss']:.6f}` peak RSS. The representative selected run takes `{selected_metrics['wall_seconds']:.6f}` wall seconds, `{selected_metrics['total_core_seconds']:.6f}` core-seconds, and `{selected_metrics['peak_rss_bytes']}` bytes peak RSS.

Every candidate and control visits all 512 masks, completes all 242 rational fixed-X1 systems, performs exactly 99,199,976,264 charged F4 word XORs, preserves equation fingerprint `{EQUATION_FINGERPRINT}`, finds no algebraic roots, and returns exhaustive UNSAT. The same flat output path also returns an exact polynomial-valid and curve-group-valid witness on the earlier SAT fixture. Thirteen Boolean-F4 tests and three backend tests pass.

The all-flat policy is rejected for single-core UNSAT. The 8 MiB hybrid is neutral/slightly slower there, so single-core retains segmented rows. Variable-order search, native/LTO code generation, and PGO are also rejected and charged in the development ledger.

The exact clean build costs `{result['clean_build']['process']['metrics']['wall_seconds']:.6f}` wall seconds and `{result['clean_build']['process']['metrics']['total_core_seconds']:.6f}` core-seconds. Direct MITM remains decisively faster: the selected F4 run is `{selected_direct['wall_ratio']:.2f}x` slower by wall, `{selected_direct['core_ratio']:.2f}x` more expensive by CPU, and `{selected_direct['rss_ratio']:.2f}x` larger by RSS than the same-binary direct median.

The seven-gate status is unchanged in substance: licensed Magma is missing, full cost does not beat automorphism-optimized rho, and unaffiliated reproduction and novelty review are missing. This is not Koblitz index-calculus SOTA.
"""
    (HERE / "RESULTS.md").write_text(text)


if __name__ == "__main__":
    main()

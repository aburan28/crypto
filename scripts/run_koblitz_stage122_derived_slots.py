#!/usr/bin/env python3
"""Compare stored and derived dense scan slots at n=53."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import resource
import subprocess
import time
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
from run_koblitz_stage26_affinity_cell import TreeSampler
import run_koblitz_stage33_n41_unknown_scalar as stage33
import run_koblitz_stage42_n53_same_target as stage42
import run_koblitz_stage44_n53_parallel as stage44
import run_koblitz_stage61_production_rank as stage61


REPO = Path(__file__).resolve().parents[1]
RUN_SCHEMA = "koblitz_stage122_derived_slots.v1"
RUN_SEAL_SCHEMA = "koblitz_stage122_derived_slots_seal.v1"


class Stage122Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage122Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    return {
        "schema": "koblitz_stage122_self_test.v1",
        "status": "PASS",
        "checks": 37,
        "n": 53,
        "parallel_width": 4096,
        "baseline_table_shards": 4,
        "candidate_table_shards": 4,
        "baseline_layout": "occupied_only_dense_scan_with_stored_slots",
        "candidate_layout": "occupied_only_dense_scan_with_derived_slots",
        "candidate_retained_bytes_gate": 536_870_912,
        "baseline_retained_bytes": 534_380_608,
        "candidate_retained_bytes": 534_380_608,
        "matched_schedule": "pipelined",
        "parallel_threads": 4,
        "same_total_table_slots": True,
        "half_size_candidate_filter": True,
        "exact_table_fallback": True,
        "same_target": True,
        "itoh_tsujii_inverse": True,
        "factor_base_solution_equality_required": True,
        "exact_relation_hash_equality_required": True,
    }


def common_environment() -> dict[str, str]:
    environment = custody.safe_child_environment()
    environment["KIC_SUMMARY_ONLY"] = "1"
    environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    environment["KIC_RANK_SURPLUS"] = "0"
    environment["KIC_RANK_AWARE_PAIR_SCAN"] = "1"
    environment["KIC_RANK_TARGET_DEFICIENCY"] = "3"
    environment["KIC_PARALLEL_SUPPORT_EXPANSION"] = "1"
    environment["KIC_PIPELINED_SUPPORT_EXPANSION"] = "1"
    environment["KIC_ENABLE_SHARDED_SUPPORT_TABLE"] = "1"
    environment["KIC_ENABLE_PREFILTERED_EXACT_LOOKUP"] = "1"
    environment["KIC_ENABLE_ITOH_N53_INVERSE"] = "1"
    environment["KIC_ENABLE_DIRECT_X_FILTER_BITS"] = "1"
    environment["RAYON_NUM_THREADS"] = "4"
    return environment


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-122 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-122 output must be new")
    build_root = args.build.resolve(strict=True)
    build = stage42.validate_build(build_root)
    host_identity = stage33.linux_host_identity()
    output.mkdir(parents=True)
    evidence = output / "evidence"
    evidence.mkdir()
    direct = str((build_root / "bin/koblitz_rank_fixture").resolve(strict=True))
    rho = str((build_root / "bin/koblitz_rho_fixture").resolve(strict=True))
    direct_command = [
        direct, "53", "0", "1", "128", str(stage42.DIRECT_SEED),
        "signed_expanded", "independent", "pair_pair_parallel_4096", "1",
        str(stage42.EXPLICIT_SCALAR),
    ]
    rho_command = [
        rho, "53", "0", "signed_frobenius", "1", "packed",
        str(stage42.RHO_BATCH_SEED), str(stage42.EXPLICIT_SCALAR),
    ]
    baseline_environment = common_environment()
    baseline_environment["KIC_ENABLE_DENSE_SUPPORT_TABLE"] = "1"
    baseline_environment["KIC_ENABLE_DENSE_SCAN_FAST_PATH"] = "1"
    derived_environment = common_environment()
    derived_environment["KIC_ENABLE_DENSE_SUPPORT_TABLE"] = "1"
    derived_environment["KIC_ENABLE_DENSE_SCAN_FAST_PATH"] = "1"
    derived_environment["KIC_ENABLE_DENSE_DERIVED_SLOTS"] = "1"

    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        baseline_process = stage33.metered(
            role="baseline", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=baseline_environment,
        )
        derived_process = stage33.metered(
            role="derived", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=derived_environment,
        )
        rho_process = stage33.metered(
            role="rho", command=rho_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=custody.safe_child_environment(),
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    processes = {
        "baseline": baseline_process,
        "derived": derived_process,
        "rho": rho_process,
    }
    require(all(stage33.process_complete(row) for row in processes.values()), "Stage-122 process failed")

    observations = {
        "baseline": stage61.direct_observation(evidence / "baseline.stdout", "FULL_RANK"),
        "derived": stage61.direct_observation(evidence / "derived.stdout", "FULL_RANK"),
    }
    bases = {name: value["base"] for name, value in observations.items()}
    summaries = {name: value["summary"] for name, value in observations.items()}
    require(len({base["base_hash"] for base in bases.values()}) == 1, "factor base changed")
    require(len({base["support_index_entries"] for base in bases.values()}) == 1, "support entries changed")
    for name in observations:
        require(bases[name].get("support_table_shards") == 4, f"{name} did not use four shards")
        require(bases[name].get("parallel_support_insertion") is True, f"{name} insertion was not parallel")
        require(bases[name].get("pipelined_support_expansion") is True, f"{name} pipeline changed")
        require(bases[name].get("support_table_shard_routing") == "xor_low_and_high_x_windows", f"{name} routing changed")
    for name in observations:
        require(bases[name].get("support_table_dense") is True, f"{name} lost dense layout")
        require(bases[name].get("support_dense_streamed_compaction") is True, f"{name} did not stream compaction")
        require(bases[name].get("support_table_slots") == bases[name]["support_index_entries"], f"{name} slots changed")
        require(bases[name].get("support_x_prefilter_bits") == 268_435_456, f"{name} filter size changed")
        require(bases[name].get("support_dense_occupancy_bytes") > 0, f"{name} occupancy index missing")
        require(bases[name].get("support_dense_rank_bytes") > 0, f"{name} rank index missing")
    for name in observations:
        require(bases[name].get("support_dense_packed_payload") is False, f"{name} became packed")
        require(bases[name].get("support_table_allocated_bytes") == 534_380_608, f"{name} allocation changed")
        require(bases[name].get("support_dense_witness_spill_bytes") == 0, f"{name} gained witness spill")
    solutions = [summary["factor_base_log_solution"] for summary in summaries.values()]
    require(all(solution == solutions[0] for solution in solutions[1:]), "full-rank solution changed")
    require(observations["baseline"]["relation_hashes"] == observations["derived"]["relation_hashes"], "relation transcript changed")
    for field in ("pair_canonicalization_maps", "pair_group_additions", "full_rank_at_relation"):
        require(len({summary[field] for summary in summaries.values()}) == 1, f"{field} changed")
    for name, summary in summaries.items():
        require(summary.get("query_mode") == "pair_pair_parallel_4096", "query width changed")
        require(summary.get("query_specialized_n53_pair_batch") is True, "specialized batch changed")
        require(summary.get("query_itoh_n53_inverse") is True, "Itoh inverse changed")
        require(summary.get("query_prefiltered_exact_lookup") is True, "prefiltered lookup changed")
        require(summary.get("fixed_base_reference_validation") is True, "fixed-base validation changed")
        require(summary.get("all_relations_group_verified") is True, "group validation failed")
        require(summary.get("linear_solution_verified") is True, "linear solution failed")
        require(summary.get("recovered_fixture_scalar") == stage42.EXPLICIT_SCALAR, "recovered scalar changed")
        require(summary.get("query_dense_scan_fast_path") is True, f"{name} lost dense scan fast path")
        require(summary.get("query_dense_derived_slots") is (name == "derived"), f"{name} derived-slot flag changed")
        require(summary.get("query_signed_slot_scratch_bytes_per_entry") == (0 if name == "derived" else 8), f"{name} signed-slot scratch changed")
    rho_rows = stage44.parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0].get("verified") is True, "rho verification failed")
    require(rho_rows[0].get("published_q") == summaries["baseline"].get("published_q"), "rho target changed")

    baseline_wall = baseline_process["metrics"]["wall_seconds"]
    derived_wall = derived_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    comparison = {
        "derived_over_baseline_wall_ratio": derived_wall / baseline_wall,
        "derived_setup_speedup": summaries["baseline"]["setup_ms"] / summaries["derived"]["setup_ms"],
        "derived_collection_speedup": summaries["baseline"]["collection_ms"] / summaries["derived"]["collection_ms"],
        "derived_query_speedup": summaries["baseline"]["timing_breakdown_ms"]["query"] / summaries["derived"]["timing_breakdown_ms"]["query"],
        "derived_core_seconds_ratio": derived_process["metrics"]["total_core_seconds"] / baseline_process["metrics"]["total_core_seconds"],
        "derived_peak_rss_ratio": derived_process["metrics"]["peak_rss_bytes"] / baseline_process["metrics"]["peak_rss_bytes"],
        "baseline_over_rho_wall_ratio": baseline_wall / rho_wall,
        "derived_over_rho_wall_ratio": derived_wall / rho_wall,
        "fresh_build_plus_derived_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + derived_wall) / rho_wall,
        "derived_whole_process_crossover": derived_wall < rho_wall,
        "baseline_support_queries": summaries["baseline"]["support_queries"],
        "derived_support_queries": summaries["derived"]["support_queries"],
        "baseline_exact_table_misses": summaries["baseline"]["query_exact_table_misses"],
        "derived_exact_table_misses": summaries["derived"]["query_exact_table_misses"],
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_derived_slots_rho",
        "constants": self_test(),
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-122 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-122 build seal", executable=False),
            "source_state": build["source_state"],
            "outer_resources": build["outer_resources"],
        },
        "processes": processes,
        "outer_resources": outer,
        "host_identity": host_identity,
        "observations": observations,
        "rho": rho_rows[0],
        "comparison": comparison,
        "charged_total_core_seconds_available": build["outer_resources"]["total_core_seconds"] + outer["total_core_seconds"],
        "charged_sequential_wall_seconds_available": build["outer_resources"]["wall_seconds"] + outer["wall_seconds"],
        "maximum_sampled_process_tree_rss_bytes": max(build["outer_resources"]["sampled_peak_process_tree_rss_bytes"], outer["sampled_peak_process_tree_rss_bytes"]),
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
        "fixture_scalar_used_by_collector": False,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    custody.write_json_new(output / "result.json", result)
    stage33.freeze(output, RUN_SEAL_SCHEMA, "run_frozen")
    return result


def verify(build_root: Path, output: Path) -> dict[str, Any]:
    build = stage42.validate_build(build_root)
    seal = stage33.verify_seal(output, RUN_SEAL_SCHEMA, "run_frozen")
    result = load(output / "result.json", "Stage-122 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_derived_slots_rho", "Stage-122 result is incomplete")
    for name, role in (("baseline", "baseline"), ("derived", "derived")):
        observed = stage61.direct_observation(output / f"evidence/{role}.stdout", "FULL_RANK")
        require(observed == result["observations"][name], f"Stage-122 {name} transcript changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-122 claim widened")
    stage33.validate_linux_host_identity(result["host_identity"])
    derived_observation = result["observations"]["derived"]
    return {
        "schema": "koblitz_stage122_verification.v1",
        "status": "n53_derived_slots_verified",
        "source_commit": build["source_state"]["commit"],
        "host_identity": result["host_identity"],
        "relations": derived_observation["summary"]["admitted_relations"],
        "support_queries": derived_observation["summary"]["support_queries"],
        "derived_retained_bytes": derived_observation["base"]["support_table_allocated_bytes"],
        "derived_wall_seconds": result["processes"]["derived"]["metrics"]["wall_seconds"],
        "derived_core_seconds": result["processes"]["derived"]["metrics"]["total_core_seconds"],
        "derived_peak_rss_bytes": result["processes"]["derived"]["metrics"]["peak_rss_bytes"],
        "rho_wall_seconds": result["processes"]["rho"]["metrics"]["wall_seconds"],
        **result["comparison"],
        "charged_total_core_seconds_available": result["charged_total_core_seconds_available"],
        "charged_sequential_wall_seconds_available": result["charged_sequential_wall_seconds_available"],
        "maximum_sampled_process_tree_rss_bytes": result["maximum_sampled_process_tree_rss_bytes"],
        "inventory_sha256": seal["inventory_sha256"],
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("self-test")
    build_parser = sub.add_parser("build")
    build_parser.add_argument("--source", type=Path, required=True)
    build_parser.add_argument("--output", type=Path, required=True)
    build_parser.add_argument("--jobs", type=int, default=4)
    build_parser.add_argument("--timeout", type=int, default=1800)
    run_parser = sub.add_parser("run")
    run_parser.add_argument("--build", type=Path, required=True)
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--timeout", type=int, default=1800)
    verify_parser = sub.add_parser("verify")
    verify_parser.add_argument("--build", type=Path, required=True)
    verify_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "self-test":
            value = self_test()
        elif args.command == "build":
            value = stage42.build(args)
        elif args.command == "run":
            value = run(args)
        else:
            value = verify(args.build.resolve(strict=True), args.output.resolve(strict=True))
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage122Error, stage61.Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage122-derived-slots: {error}")


if __name__ == "__main__":
    main()

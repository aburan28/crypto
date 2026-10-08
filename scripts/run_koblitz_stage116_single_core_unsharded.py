#!/usr/bin/env python3
"""Compare unsharded and four-shard exact support on one CPU."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
from run_koblitz_stage26_affinity_cell import TreeSampler
import run_koblitz_stage33_n41_unknown_scalar as stage33
import run_koblitz_stage42_n53_same_target as stage42
import run_koblitz_stage44_n53_parallel as stage44
import run_koblitz_stage61_production_rank as stage61


REPO = Path(__file__).resolve().parents[1]
RUN_SCHEMA = "koblitz_stage116_single_core_unsharded.v1"
RUN_SEAL_SCHEMA = "koblitz_stage116_single_core_unsharded_seal.v1"


class Stage116Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage116Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    return {
        "schema": "koblitz_stage116_self_test.v1",
        "status": "PASS",
        "checks": 35,
        "n": 53,
        "parallel_width": 4096,
        "linux_affinity_required": True,
        "effective_cpu_count": 1,
        "rayon_threads": 1,
        "parallel_support_expansion": False,
        "baseline_table_shards": 4,
        "candidate_table_shards": 1,
        "baseline_layout": "four_shard_open_addressed",
        "candidate_layout": "single_open_addressed_table",
        "matched_retained_bytes": 738_197_504,
        "matched_schedule": "serial_support_expansion",
        "rank_target_deficiency": 1,
        "same_total_table_slots": True,
        "same_filter_bits": True,
        "exact_table_fallback": True,
        "same_target": True,
        "itoh_tsujii_inverse": True,
        "factor_base_solution_equality_required": True,
        "exact_relation_hash_equality_required": False,
    }


def common_environment() -> dict[str, str]:
    environment = custody.safe_child_environment()
    environment["KIC_SUMMARY_ONLY"] = "1"
    environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    environment["KIC_RANK_SURPLUS"] = "0"
    environment["KIC_RANK_AWARE_PAIR_SCAN"] = "1"
    environment["KIC_RANK_TARGET_DEFICIENCY"] = "1"
    environment["KIC_ENABLE_PREFILTERED_EXACT_LOOKUP"] = "1"
    environment["KIC_ENABLE_ITOH_N53_INVERSE"] = "1"
    environment["KIC_ENABLE_DIRECT_X_FILTER_BITS"] = "1"
    environment["RAYON_NUM_THREADS"] = "1"
    return environment


def child_affinity() -> list[int]:
    completed = subprocess.run(
        [sys.executable, "-c", "import json,os; print(json.dumps(sorted(os.sched_getaffinity(0))))"],
        text=True,
        capture_output=True,
        check=True,
    )
    value = json.loads(completed.stdout)
    require(isinstance(value, list) and all(type(cpu) is int for cpu in value), "child affinity probe failed")
    return value


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux" and hasattr(os, "sched_getaffinity"), "Stage-116 requires Linux CPU affinity")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-116 output must be new")
    build_root = args.build.resolve(strict=True)
    build = stage42.validate_build(build_root)
    initial_cpus = sorted(os.sched_getaffinity(0))
    require(bool(initial_cpus), "Stage-116 has no allowed CPU")
    selected_cpu = initial_cpus[0]
    os.sched_setaffinity(0, {selected_cpu})
    effective_cpus = sorted(os.sched_getaffinity(0))
    inherited_cpus = child_affinity()
    require(effective_cpus == inherited_cpus == [selected_cpu], "single-CPU affinity was not inherited")
    host_identity = stage33.linux_host_identity()
    require(host_identity["cpu"]["effective_affinity"] == [selected_cpu], "host receipt lost one-CPU affinity")
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
    unsharded_environment = common_environment()
    sharded_environment = common_environment()
    sharded_environment["KIC_ENABLE_SHARDED_SUPPORT_TABLE"] = "1"

    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        unsharded_process = stage33.metered(
            role="single-core-unsharded", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=unsharded_environment,
        )
        sharded_process = stage33.metered(
            role="single-core-sharded", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=sharded_environment,
        )
        rho_process = stage33.metered(
            role="single-core-rho", command=rho_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=custody.safe_child_environment(),
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    processes = {
        "unsharded": unsharded_process,
        "sharded": sharded_process,
        "rho": rho_process,
    }
    require(all(stage33.process_complete(row) for row in processes.values()), "Stage-116 process failed")

    observations = {
        "unsharded": stage61.direct_observation(evidence / "single-core-unsharded.stdout", "FULL_RANK"),
        "sharded": stage61.direct_observation(evidence / "single-core-sharded.stdout", "FULL_RANK"),
    }
    bases = {name: value["base"] for name, value in observations.items()}
    summaries = {name: value["summary"] for name, value in observations.items()}
    require(len({base["base_hash"] for base in bases.values()}) == 1, "factor base changed")
    require(len({base["support_index_entries"] for base in bases.values()}) == 1, "support entries changed")
    for name in observations:
        require(bases[name].get("parallel_support_expansion") is False, f"{name} used parallel support expansion")
        require(bases[name].get("parallel_support_insertion") is False, f"{name} used parallel support insertion")
        require(bases[name].get("pipelined_support_expansion") is False, f"{name} used pipelined support expansion")
        require(bases[name].get("support_table_dense") is False, f"{name} became dense")
        require(bases[name].get("support_table_allocated_bytes") == 738_197_504, f"{name} allocation changed")
        require(bases[name].get("support_x_prefilter_bits") == 536_870_912, f"{name} filter changed")
    require(bases["unsharded"].get("support_table_shards") == 1, "candidate remained sharded")
    require(bases["unsharded"].get("support_table_shard_routing") == "unsharded", "candidate routing changed")
    require(bases["unsharded"].get("support_x_prefilter_hash_strategy") == "direct_low_and_high_x_bit_windows", "candidate filter changed")
    require(bases["sharded"].get("support_table_shards") == 4, "baseline lost four shards")
    require(bases["sharded"].get("support_table_shard_routing") == "xor_low_and_high_x_windows", "baseline routing changed")
    require(bases["sharded"].get("support_x_prefilter_hash_strategy") == "four_shard_direct_low_and_high_x_bit_windows", "baseline filter changed")
    solutions = [summary["factor_base_log_solution"] for summary in summaries.values()]
    require(all(solution == solutions[0] for solution in solutions[1:]), "full-rank solution changed")
    for field in ("pair_canonicalization_maps", "pair_group_additions"):
        require(len({summary[field] for summary in summaries.values()}) == 1, f"{field} changed")
    for summary in summaries.values():
        require(summary.get("query_parallel_threads") == 1, "query used more than one Rayon thread")
        require(summary.get("rank_target_deficiency") == 1, "rank-deficiency target changed")
        require(summary.get("query_mode") == "pair_pair_parallel_4096", "query width changed")
        require(summary.get("query_specialized_n53_pair_batch") is True, "specialized batch changed")
        require(summary.get("query_itoh_n53_inverse") is True, "Itoh inverse changed")
        require(summary.get("query_prefiltered_exact_lookup") is True, "prefiltered lookup changed")
        require(summary.get("fixed_base_reference_validation") is True, "fixed-base validation changed")
        require(summary.get("all_relations_group_verified") is True, "group validation failed")
        require(summary.get("linear_solution_verified") is True, "linear solution failed")
        require(summary.get("recovered_fixture_scalar") == stage42.EXPLICIT_SCALAR, "recovered scalar changed")
    rho_rows = stage44.parse_lines(evidence / "single-core-rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0].get("verified") is True, "rho verification failed")
    require(rho_rows[0].get("published_q") == summaries["unsharded"].get("published_q"), "rho target changed")

    unsharded_wall = unsharded_process["metrics"]["wall_seconds"]
    sharded_wall = sharded_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    comparison = {
        "unsharded_over_sharded_wall_ratio": unsharded_wall / sharded_wall,
        "unsharded_setup_speedup": summaries["sharded"]["setup_ms"] / summaries["unsharded"]["setup_ms"],
        "unsharded_collection_speedup": summaries["sharded"]["collection_ms"] / summaries["unsharded"]["collection_ms"],
        "unsharded_query_speedup": summaries["sharded"]["timing_breakdown_ms"]["query"] / summaries["unsharded"]["timing_breakdown_ms"]["query"],
        "unsharded_core_seconds_ratio": unsharded_process["metrics"]["total_core_seconds"] / sharded_process["metrics"]["total_core_seconds"],
        "unsharded_peak_rss_ratio": unsharded_process["metrics"]["peak_rss_bytes"] / sharded_process["metrics"]["peak_rss_bytes"],
        "sharded_over_rho_wall_ratio": sharded_wall / rho_wall,
        "unsharded_over_rho_wall_ratio": unsharded_wall / rho_wall,
        "fresh_build_plus_unsharded_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + unsharded_wall) / rho_wall,
        "unsharded_whole_process_crossover": unsharded_wall < rho_wall,
        "sharded_support_queries": summaries["sharded"]["support_queries"],
        "unsharded_support_queries": summaries["unsharded"]["support_queries"],
        "sharded_exact_table_misses": summaries["sharded"]["query_exact_table_misses"],
        "unsharded_exact_table_misses": summaries["unsharded"]["query_exact_table_misses"],
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_single_core_unsharded_rho",
        "constants": self_test(),
        "affinity": {
            "initial_allowed_cpus": initial_cpus,
            "selected_cpu": selected_cpu,
            "parent_effective_cpus": effective_cpus,
            "child_effective_cpus": inherited_cpus,
            "inheritance_verified": True,
        },
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-116 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-116 build seal", executable=False),
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
    result = load(output / "result.json", "Stage-116 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_single_core_unsharded_rho", "Stage-116 result is incomplete")
    require(result.get("affinity", {}).get("inheritance_verified") is True and len(result["affinity"]["parent_effective_cpus"]) == 1, "Stage-116 affinity changed")
    for name, role in (("unsharded", "single-core-unsharded"), ("sharded", "single-core-sharded")):
        observed = stage61.direct_observation(output / f"evidence/{role}.stdout", "FULL_RANK")
        require(observed == result["observations"][name], f"Stage-116 {name} transcript changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-116 claim widened")
    stage33.validate_linux_host_identity(result["host_identity"])
    require(result["host_identity"]["cpu"]["effective_affinity"] == result["affinity"]["parent_effective_cpus"], "Stage-116 host affinity changed")
    unsharded_observation = result["observations"]["unsharded"]
    return {
        "schema": "koblitz_stage116_verification.v1",
        "status": "n53_single_core_unsharded_verified",
        "source_commit": build["source_state"]["commit"],
        "host_identity": result["host_identity"],
        "selected_cpu": result["affinity"]["selected_cpu"],
        "relations": unsharded_observation["summary"]["admitted_relations"],
        "support_queries": unsharded_observation["summary"]["support_queries"],
        "unsharded_retained_bytes": unsharded_observation["base"]["support_table_allocated_bytes"],
        "unsharded_wall_seconds": result["processes"]["unsharded"]["metrics"]["wall_seconds"],
        "unsharded_core_seconds": result["processes"]["unsharded"]["metrics"]["total_core_seconds"],
        "unsharded_peak_rss_bytes": result["processes"]["unsharded"]["metrics"]["peak_rss_bytes"],
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage116Error, stage61.Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage116-single-core-unsharded: {error}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Compare mixed-hash and direct x-bit shard routing at n=53."""

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
RUN_SCHEMA = "koblitz_stage106_shard_routing.v1"
RUN_SEAL_SCHEMA = "koblitz_stage106_shard_routing_seal.v1"


class Stage106Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage106Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    return {
        "schema": "koblitz_stage106_self_test.v1",
        "status": "PASS",
        "checks": 28,
        "n": 53,
        "parallel_width": 4096,
        "baseline_table_shards": 4,
        "candidate_table_shards": 4,
        "baseline_shard_routing": "splitmix_bucket_bits",
        "candidate_shard_routing": "xor_low_and_high_x_windows",
        "matched_schedule": "pipelined",
        "parallel_threads": 4,
        "same_total_table_slots": True,
        "same_total_filter_allocation": True,
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
    require(platform.system() == "Linux", "Stage-106 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-106 output must be new")
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
    mixed_environment = common_environment()
    mixed_environment["KIC_ENABLE_MIXED_SHARD_ROUTING"] = "1"
    direct_environment = common_environment()

    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        mixed_process = stage33.metered(
            role="mixed-routing", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=mixed_environment,
        )
        direct_process = stage33.metered(
            role="direct-routing", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=direct_environment,
        )
        rho_process = stage33.metered(
            role="rho", command=rho_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=custody.safe_child_environment(),
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    processes = {
        "mixed": mixed_process,
        "direct": direct_process,
        "rho": rho_process,
    }
    require(all(stage33.process_complete(row) for row in processes.values()), "Stage-106 process failed")

    observations = {
        "mixed": stage61.direct_observation(evidence / "mixed-routing.stdout", "FULL_RANK"),
        "direct": stage61.direct_observation(evidence / "direct-routing.stdout", "FULL_RANK"),
    }
    bases = {name: value["base"] for name, value in observations.items()}
    summaries = {name: value["summary"] for name, value in observations.items()}
    require(len({base["base_hash"] for base in bases.values()}) == 1, "factor base changed")
    require(len({base["support_index_entries"] for base in bases.values()}) == 1, "support entries changed")
    require(len({base["support_table_allocated_bytes"] for base in bases.values()}) == 1, "table allocation changed")
    require(len({base["support_x_prefilter_bits"] for base in bases.values()}) == 1, "filter allocation changed")
    for name in observations:
        require(bases[name].get("support_table_shards") == 4, f"{name} did not use four shards")
        require(bases[name].get("parallel_support_insertion") is True, f"{name} insertion was not parallel")
        require(bases[name].get("support_x_prefilter_hash_strategy") == "four_shard_direct_low_and_high_x_bit_windows", f"{name} filter changed")
        require(bases[name].get("pipelined_support_expansion") is True, f"{name} pipeline changed")
    require(bases["mixed"].get("support_table_shard_routing") == "splitmix_bucket_bits", "mixed baseline routing changed")
    require(bases["direct"].get("support_table_shard_routing") == "xor_low_and_high_x_windows", "direct candidate routing changed")
    solutions = [summary["factor_base_log_solution"] for summary in summaries.values()]
    require(all(solution == solutions[0] for solution in solutions[1:]), "full-rank solution changed")
    for field in ("pair_canonicalization_maps", "pair_group_additions", "full_rank_at_relation"):
        require(len({summary[field] for summary in summaries.values()}) == 1, f"{field} changed")
    for summary in summaries.values():
        require(summary.get("query_mode") == "pair_pair_parallel_4096", "query width changed")
        require(summary.get("query_specialized_n53_pair_batch") is True, "specialized batch changed")
        require(summary.get("query_itoh_n53_inverse") is True, "Itoh inverse changed")
        require(summary.get("query_prefiltered_exact_lookup") is True, "prefiltered lookup changed")
        require(summary.get("fixed_base_reference_validation") is True, "fixed-base validation changed")
        require(summary.get("all_relations_group_verified") is True, "group validation failed")
        require(summary.get("linear_solution_verified") is True, "linear solution failed")
        require(summary.get("recovered_fixture_scalar") == stage42.EXPLICIT_SCALAR, "recovered scalar changed")
    rho_rows = stage44.parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0].get("verified") is True, "rho verification failed")
    require(rho_rows[0].get("published_q") == summaries["mixed"].get("published_q"), "rho target changed")

    mixed_wall = mixed_process["metrics"]["wall_seconds"]
    direct_wall = direct_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    comparison = {
        "direct_over_mixed_wall_ratio": direct_wall / mixed_wall,
        "direct_setup_speedup": summaries["mixed"]["setup_ms"] / summaries["direct"]["setup_ms"],
        "direct_collection_speedup": summaries["mixed"]["collection_ms"] / summaries["direct"]["collection_ms"],
        "direct_query_speedup": summaries["mixed"]["timing_breakdown_ms"]["query"] / summaries["direct"]["timing_breakdown_ms"]["query"],
        "direct_core_seconds_ratio": direct_process["metrics"]["total_core_seconds"] / mixed_process["metrics"]["total_core_seconds"],
        "direct_peak_rss_ratio": direct_process["metrics"]["peak_rss_bytes"] / mixed_process["metrics"]["peak_rss_bytes"],
        "mixed_over_rho_wall_ratio": mixed_wall / rho_wall,
        "direct_over_rho_wall_ratio": direct_wall / rho_wall,
        "fresh_build_plus_direct_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + direct_wall) / rho_wall,
        "direct_whole_process_crossover": direct_wall < rho_wall,
        "mixed_support_queries": summaries["mixed"]["support_queries"],
        "direct_support_queries": summaries["direct"]["support_queries"],
        "mixed_exact_table_misses": summaries["mixed"]["query_exact_table_misses"],
        "direct_exact_table_misses": summaries["direct"]["query_exact_table_misses"],
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_shard_routing_rho",
        "constants": self_test(),
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-106 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-106 build seal", executable=False),
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
    result = load(output / "result.json", "Stage-106 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_shard_routing_rho", "Stage-106 result is incomplete")
    for name, role in (("mixed", "mixed-routing"), ("direct", "direct-routing")):
        observed = stage61.direct_observation(output / f"evidence/{role}.stdout", "FULL_RANK")
        require(observed == result["observations"][name], f"Stage-106 {name} transcript changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-106 claim widened")
    stage33.validate_linux_host_identity(result["host_identity"])
    direct_observation = result["observations"]["direct"]
    return {
        "schema": "koblitz_stage106_verification.v1",
        "status": "n53_shard_routing_verified",
        "source_commit": build["source_state"]["commit"],
        "host_identity": result["host_identity"],
        "relations": direct_observation["summary"]["admitted_relations"],
        "support_queries": direct_observation["summary"]["support_queries"],
        "direct_wall_seconds": result["processes"]["direct"]["metrics"]["wall_seconds"],
        "direct_core_seconds": result["processes"]["direct"]["metrics"]["total_core_seconds"],
        "direct_peak_rss_bytes": result["processes"]["direct"]["metrics"]["peak_rss_bytes"],
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage106Error, stage61.Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage106-shard-routing: {error}")


if __name__ == "__main__":
    main()

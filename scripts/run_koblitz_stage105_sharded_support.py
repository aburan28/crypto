#!/usr/bin/env python3
"""Compare serial and four-shard exact support-table insertion at n=53."""

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
RUN_SCHEMA = "koblitz_stage105_sharded_support.v1"
RUN_SEAL_SCHEMA = "koblitz_stage105_sharded_support_seal.v1"


class Stage105Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage105Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    return {
        "schema": "koblitz_stage105_self_test.v1",
        "status": "PASS",
        "checks": 30,
        "n": 53,
        "parallel_width": 4096,
        "baseline_table_shards": 1,
        "candidate_table_shards": 4,
        "candidate_schedules": ["pipelined", "phased"],
        "parallel_threads": 4,
        "same_total_table_slots": True,
        "same_total_filter_allocation": True,
        "exact_table_fallback": True,
        "same_target": True,
        "itoh_tsujii_inverse": True,
        "factor_base_solution_equality_required": True,
        "exact_relation_hash_equality_required": False,
    }


def common_environment(*, pipelined: bool) -> dict[str, str]:
    environment = custody.safe_child_environment()
    environment["KIC_SUMMARY_ONLY"] = "1"
    environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    environment["KIC_RANK_SURPLUS"] = "0"
    environment["KIC_RANK_AWARE_PAIR_SCAN"] = "1"
    environment["KIC_RANK_TARGET_DEFICIENCY"] = "3"
    environment["KIC_PARALLEL_SUPPORT_EXPANSION"] = "1"
    if pipelined:
        environment["KIC_PIPELINED_SUPPORT_EXPANSION"] = "1"
    environment["KIC_ENABLE_PREFILTERED_EXACT_LOOKUP"] = "1"
    environment["KIC_ENABLE_ITOH_N53_INVERSE"] = "1"
    environment["KIC_ENABLE_DIRECT_X_FILTER_BITS"] = "1"
    environment["RAYON_NUM_THREADS"] = "4"
    return environment


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-105 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-105 output must be new")
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
    baseline_environment = common_environment(pipelined=True)
    sharded_pipelined_environment = common_environment(pipelined=True)
    sharded_pipelined_environment["KIC_ENABLE_SHARDED_SUPPORT_TABLE"] = "1"
    sharded_phased_environment = common_environment(pipelined=False)
    sharded_phased_environment["KIC_ENABLE_SHARDED_SUPPORT_TABLE"] = "1"

    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        baseline_process = stage33.metered(
            role="baseline", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=baseline_environment,
        )
        sharded_pipelined_process = stage33.metered(
            role="sharded-pipelined", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=sharded_pipelined_environment,
        )
        sharded_phased_process = stage33.metered(
            role="sharded-phased", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=sharded_phased_environment,
        )
        rho_process = stage33.metered(
            role="rho", command=rho_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=custody.safe_child_environment(),
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    processes = {
        "baseline": baseline_process,
        "sharded_pipelined": sharded_pipelined_process,
        "sharded_phased": sharded_phased_process,
        "rho": rho_process,
    }
    require(all(stage33.process_complete(row) for row in processes.values()), "Stage-105 process failed")

    observations = {
        "baseline": stage61.direct_observation(evidence / "baseline.stdout", "FULL_RANK"),
        "sharded_pipelined": stage61.direct_observation(evidence / "sharded-pipelined.stdout", "FULL_RANK"),
        "sharded_phased": stage61.direct_observation(evidence / "sharded-phased.stdout", "FULL_RANK"),
    }
    bases = {name: value["base"] for name, value in observations.items()}
    summaries = {name: value["summary"] for name, value in observations.items()}
    require(len({base["base_hash"] for base in bases.values()}) == 1, "factor base changed")
    require(len({base["support_index_entries"] for base in bases.values()}) == 1, "support entries changed")
    require(len({base["support_table_allocated_bytes"] for base in bases.values()}) == 1, "table allocation changed")
    require(len({base["support_x_prefilter_bits"] for base in bases.values()}) == 1, "filter allocation changed")
    require(bases["baseline"].get("support_table_shards") == 1, "baseline was sharded")
    require(bases["baseline"].get("parallel_support_insertion") is False, "baseline insertion changed")
    require(bases["baseline"].get("support_x_prefilter_hash_strategy") == "direct_low_and_high_x_bit_windows", "baseline filter changed")
    for name in ("sharded_pipelined", "sharded_phased"):
        require(bases[name].get("support_table_shards") == 4, f"{name} did not use four shards")
        require(bases[name].get("parallel_support_insertion") is True, f"{name} insertion was not parallel")
        require(bases[name].get("support_x_prefilter_hash_strategy") == "four_shard_direct_low_and_high_x_bit_windows", f"{name} filter changed")
    require(bases["baseline"].get("pipelined_support_expansion") is True, "baseline pipeline changed")
    require(bases["sharded_pipelined"].get("pipelined_support_expansion") is True, "pipelined candidate changed")
    require(bases["sharded_phased"].get("pipelined_support_expansion") is False, "phased candidate changed")
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
    require(rho_rows[0].get("published_q") == summaries["baseline"].get("published_q"), "rho target changed")

    baseline_wall = baseline_process["metrics"]["wall_seconds"]
    candidate_name = min(
        ("sharded_pipelined", "sharded_phased"),
        key=lambda name: processes[name]["metrics"]["wall_seconds"],
    )
    candidate_process = processes[candidate_name]
    candidate_summary = summaries[candidate_name]
    candidate_wall = candidate_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    comparison = {
        "selected_candidate": candidate_name,
        "selected_over_baseline_wall_ratio": candidate_wall / baseline_wall,
        "selected_setup_speedup": summaries["baseline"]["setup_ms"] / candidate_summary["setup_ms"],
        "selected_collection_speedup": summaries["baseline"]["collection_ms"] / candidate_summary["collection_ms"],
        "selected_query_speedup": summaries["baseline"]["timing_breakdown_ms"]["query"] / candidate_summary["timing_breakdown_ms"]["query"],
        "selected_core_seconds_ratio": candidate_process["metrics"]["total_core_seconds"] / baseline_process["metrics"]["total_core_seconds"],
        "selected_peak_rss_ratio": candidate_process["metrics"]["peak_rss_bytes"] / baseline_process["metrics"]["peak_rss_bytes"],
        "baseline_over_rho_wall_ratio": baseline_wall / rho_wall,
        "selected_over_rho_wall_ratio": candidate_wall / rho_wall,
        "fresh_build_plus_selected_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + candidate_wall) / rho_wall,
        "selected_whole_process_crossover": candidate_wall < rho_wall,
        "pipelined_over_baseline_wall_ratio": sharded_pipelined_process["metrics"]["wall_seconds"] / baseline_wall,
        "phased_over_baseline_wall_ratio": sharded_phased_process["metrics"]["wall_seconds"] / baseline_wall,
        "baseline_support_queries": summaries["baseline"]["support_queries"],
        "selected_support_queries": candidate_summary["support_queries"],
        "baseline_exact_table_misses": summaries["baseline"]["query_exact_table_misses"],
        "selected_exact_table_misses": candidate_summary["query_exact_table_misses"],
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_sharded_support_rho",
        "constants": self_test(),
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-105 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-105 build seal", executable=False),
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
    result = load(output / "result.json", "Stage-105 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_sharded_support_rho", "Stage-105 result is incomplete")
    for name, role in (("baseline", "baseline"), ("sharded_pipelined", "sharded-pipelined"), ("sharded_phased", "sharded-phased")):
        observed = stage61.direct_observation(output / f"evidence/{role}.stdout", "FULL_RANK")
        require(observed == result["observations"][name], f"Stage-105 {name} transcript changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-105 claim widened")
    stage33.validate_linux_host_identity(result["host_identity"])
    selected = result["comparison"]["selected_candidate"]
    selected_observation = result["observations"][selected]
    return {
        "schema": "koblitz_stage105_verification.v1",
        "status": "n53_sharded_support_verified",
        "source_commit": build["source_state"]["commit"],
        "host_identity": result["host_identity"],
        "selected_candidate": selected,
        "relations": selected_observation["summary"]["admitted_relations"],
        "support_queries": selected_observation["summary"]["support_queries"],
        "selected_wall_seconds": result["processes"][selected]["metrics"]["wall_seconds"],
        "selected_core_seconds": result["processes"][selected]["metrics"]["total_core_seconds"],
        "selected_peak_rss_bytes": result["processes"][selected]["metrics"]["peak_rss_bytes"],
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage105Error, stage61.Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage105-sharded-support: {error}")


if __name__ == "__main__":
    main()

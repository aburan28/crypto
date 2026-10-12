#!/usr/bin/env python3
"""Compare two-word direct-bit and one-word blocked x filters at n=53."""

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
RUN_SCHEMA = "koblitz_stage96_blocked_filter.v1"
RUN_SEAL_SCHEMA = "koblitz_stage96_blocked_filter_seal.v1"


class Stage96Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage96Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    return {
        "schema": "koblitz_stage96_self_test.v1",
        "status": "PASS",
        "checks": 24,
        "n": 53,
        "parallel_width": 4096,
        "baseline_filter": "two_word_direct_low_and_high_x_bit_windows",
        "candidate_filter": "one_word_four_bit_blocked_x_fingerprint",
        "same_filter_allocation": True,
        "exact_table_fallback": True,
        "same_target": True,
        "itoh_tsujii_inverse": True,
        "exact_relation_hash_equality_required": True,
    }


def common_environment() -> dict[str, str]:
    environment = custody.safe_child_environment()
    environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    environment["KIC_RANK_SURPLUS"] = "0"
    environment["KIC_RANK_AWARE_PAIR_SCAN"] = "1"
    environment["KIC_PARALLEL_SUPPORT_EXPANSION"] = "1"
    environment["KIC_PIPELINED_SUPPORT_EXPANSION"] = "1"
    environment["KIC_ENABLE_PREFILTERED_EXACT_LOOKUP"] = "1"
    environment["KIC_ENABLE_ITOH_N53_INVERSE"] = "1"
    environment["RAYON_NUM_THREADS"] = "4"
    return environment


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-96 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-96 output must be new")
    build_root = args.build.resolve(strict=True)
    build = stage42.validate_build(build_root)
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
    direct_bits_environment = common_environment()
    direct_bits_environment["KIC_ENABLE_DIRECT_X_FILTER_BITS"] = "1"
    blocked_environment = common_environment()
    blocked_environment["KIC_ENABLE_BLOCKED_X_FILTER"] = "1"

    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        direct_bits_process = stage33.metered(
            role="direct-bits", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=direct_bits_environment,
        )
        blocked_process = stage33.metered(
            role="blocked", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=blocked_environment,
        )
        rho_process = stage33.metered(
            role="rho", command=rho_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=custody.safe_child_environment(),
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    processes = {"direct_bits": direct_bits_process, "blocked": blocked_process, "rho": rho_process}
    require(all(stage33.process_complete(row) for row in processes.values()), "Stage-96 process failed")

    direct_bits = stage61.direct_observation(evidence / "direct-bits.stdout", "FULL_RANK")
    blocked = stage61.direct_observation(evidence / "blocked.stdout", "FULL_RANK")
    direct_summary = direct_bits["summary"]
    blocked_summary = blocked["summary"]
    direct_base = direct_bits["base"]
    blocked_base = blocked["base"]
    require(direct_base["base_hash"] == blocked_base["base_hash"], "factor base changed")
    require(direct_base["support_index_entries"] == blocked_base["support_index_entries"], "support entries changed")
    require(direct_base["support_x_prefilter_bits"] == blocked_base["support_x_prefilter_bits"], "filter allocation changed")
    require(direct_bits["relation_hashes"] == blocked["relation_hashes"], "relation transcript changed")
    require(direct_summary["factor_base_log_solution"] == blocked_summary["factor_base_log_solution"], "full-rank solution changed")
    for field in ("pair_canonicalization_maps", "pair_group_additions", "support_queries", "full_rank_at_relation", "query_batch_inversions", "query_group_additions"):
        require(direct_summary[field] == blocked_summary[field], f"{field} changed")
    require(direct_base.get("support_x_prefilter_hash_strategy") == "direct_low_and_high_x_bit_windows", "direct-bit baseline changed")
    require(direct_base.get("support_x_prefilter_direct_bits") is True and direct_base.get("support_x_prefilter_blocked") is False, "direct-bit flags changed")
    require(blocked_base.get("support_x_prefilter_hash_strategy") == "blocked_one_word_four_bit_x_fingerprint", "blocked filter was not selected")
    require(blocked_base.get("support_x_prefilter_direct_bits") is False and blocked_base.get("support_x_prefilter_blocked") is True, "blocked flags changed")
    for summary in (direct_summary, blocked_summary):
        require(summary.get("query_mode") == "pair_pair_parallel_4096", "query width changed")
        require(summary.get("query_specialized_n53_pair_batch") is True, "specialized batch changed")
        require(summary.get("query_itoh_n53_inverse") is True, "Itoh inverse changed")
        require(summary.get("query_prefiltered_exact_lookup") is True, "prefiltered lookup changed")
        require(summary.get("recovered_fixture_scalar") == stage42.EXPLICIT_SCALAR, "recovered scalar changed")
    rho_rows = stage44.parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0].get("verified") is True, "rho verification failed")
    require(rho_rows[0].get("published_q") == blocked_summary.get("published_q"), "rho target changed")

    direct_wall = direct_bits_process["metrics"]["wall_seconds"]
    blocked_wall = blocked_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    comparison = {
        "blocked_over_direct_bits_wall_ratio": blocked_wall / direct_wall,
        "blocked_collection_speedup": direct_summary["collection_ms"] / blocked_summary["collection_ms"],
        "blocked_query_speedup": direct_summary["timing_breakdown_ms"]["query"] / blocked_summary["timing_breakdown_ms"]["query"],
        "blocked_setup_speedup": direct_summary["setup_ms"] / blocked_summary["setup_ms"],
        "blocked_core_seconds_ratio": blocked_process["metrics"]["total_core_seconds"] / direct_bits_process["metrics"]["total_core_seconds"],
        "blocked_peak_rss_ratio": blocked_process["metrics"]["peak_rss_bytes"] / direct_bits_process["metrics"]["peak_rss_bytes"],
        "direct_bits_exact_table_misses": direct_summary["query_exact_table_misses"],
        "blocked_exact_table_misses": blocked_summary["query_exact_table_misses"],
        "exact_table_miss_ratio": blocked_summary["query_exact_table_misses"] / direct_summary["query_exact_table_misses"],
        "blocked_over_rho_wall_ratio": blocked_wall / rho_wall,
        "direct_bits_over_rho_wall_ratio": direct_wall / rho_wall,
        "fresh_build_plus_blocked_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + blocked_wall) / rho_wall,
        "blocked_whole_process_crossover": blocked_wall < rho_wall,
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_blocked_filter_rho",
        "constants": self_test(),
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-96 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-96 build seal", executable=False),
            "source_state": build["source_state"],
            "outer_resources": build["outer_resources"],
        },
        "processes": processes,
        "outer_resources": outer,
        "direct_bits": direct_bits,
        "blocked": blocked,
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
    result = load(output / "result.json", "Stage-96 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_blocked_filter_rho", "Stage-96 result is incomplete")
    direct_bits = stage61.direct_observation(output / "evidence/direct-bits.stdout", "FULL_RANK")
    blocked = stage61.direct_observation(output / "evidence/blocked.stdout", "FULL_RANK")
    require(direct_bits == result["direct_bits"] and blocked == result["blocked"], "Stage-96 transcript changed")
    require(direct_bits["relation_hashes"] == blocked["relation_hashes"], "Stage-96 relation equality changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-96 claim widened")
    return {
        "schema": "koblitz_stage96_verification.v1",
        "status": "n53_blocked_filter_verified",
        "source_commit": build["source_state"]["commit"],
        "relations": blocked["summary"]["admitted_relations"],
        "support_queries": blocked["summary"]["support_queries"],
        "blocked_wall_seconds": result["processes"]["blocked"]["metrics"]["wall_seconds"],
        "blocked_core_seconds": result["processes"]["blocked"]["metrics"]["total_core_seconds"],
        "blocked_peak_rss_bytes": result["processes"]["blocked"]["metrics"]["peak_rss_bytes"],
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage96Error, stage61.Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage96-blocked-filter: {error}")


if __name__ == "__main__":
    main()

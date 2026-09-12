#!/usr/bin/env python3
"""Compare ordinary and rank-aware n=53 full-rank collection."""

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
RUN_SCHEMA = "koblitz_stage62_rank_aware.v1"
RUN_SEAL_SCHEMA = "koblitz_stage62_rank_aware_seal.v1"


class Stage62Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage62Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    predecessor = stage61.self_test()
    require(predecessor.get("status") == "PASS", "Stage-61 control failed")
    return {
        "schema": "koblitz_stage62_self_test.v1",
        "status": "PASS",
        "checks": 14,
        "n": 53,
        "eta": [1, 128],
        "query_mode": "pair_pair_parallel_4096",
        "parallel_threads": 4,
        "ordinary_surplus_relations": 0,
        "rank_aware_surplus_relations": 0,
        "targeting_trigger": "rank_at_least_orbit_columns",
        "targeting_material": "public factor-base orbit column only",
        "discrete_log_labels_used": False,
        "same_target": True,
        "common_relation_prefix_required": True,
        "full_rank_solution_equality_required": True,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-62 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-62 output must be new")
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
    ordinary_environment = custody.safe_child_environment()
    ordinary_environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    ordinary_environment["KIC_RANK_SURPLUS"] = "0"
    ordinary_environment["RAYON_NUM_THREADS"] = "4"
    candidate_environment = dict(ordinary_environment)
    candidate_environment["KIC_RANK_AWARE_PAIR_SCAN"] = "1"
    rho_environment = custody.safe_child_environment()
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        ordinary_process = stage33.metered(
            role="ordinary-full-rank", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=ordinary_environment,
        )
        candidate_process = stage33.metered(
            role="rank-aware", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=candidate_environment,
        )
        rho_process = stage33.metered(
            role="rho", command=rho_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=rho_environment,
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    require(all(stage33.process_complete(row) for row in (ordinary_process, candidate_process, rho_process)), "Stage-62 process failed")
    ordinary = stage61.direct_observation(evidence / "ordinary-full-rank.stdout", "FULL_RANK")
    candidate = stage61.direct_observation(evidence / "rank-aware.stdout", "FULL_RANK")
    ordinary_summary = ordinary["summary"]
    candidate_summary = candidate["summary"]
    require(ordinary["base"]["base_hash"] == candidate["base"]["base_hash"], "factor base changed")
    require(ordinary_summary.get("rank_aware_pair_scan") is False, "ordinary mode changed")
    require(candidate_summary.get("rank_aware_pair_scan") is True, "rank-aware mode was not enabled")
    require(candidate_summary["admitted_relations"] == candidate_summary["full_rank_at_relation"], "candidate did not stop at full rank")
    prefix = min(ordinary_summary["orbit_columns"], len(candidate["relation_hashes"]))
    require(candidate["relation_hashes"][:prefix] == ordinary["relation_hashes"][:prefix], "pre-targeting relation prefix changed")
    require(candidate_summary["factor_base_log_solution"] == ordinary_summary["factor_base_log_solution"], "full-rank solution changed")
    require(candidate_summary["recovered_fixture_scalar"] == ordinary_summary["recovered_fixture_scalar"] == stage42.EXPLICIT_SCALAR, "recovered scalar changed")
    require(candidate_summary.get("rank_target_slot_index_builds") == 1, "candidate slot index build count changed")
    require(candidate_summary.get("rank_target_slot_index_allocated_bytes", 0) < 4 * 1024 * 1024, "candidate slot index exceeds 4 MiB")
    require(candidate_summary.get("field_mul_backend") == candidate_summary.get("field_square_backend") == "x86_64_pclmulqdq", "PCLMUL backend unavailable")
    rho_rows = stage44.parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0].get("verified") is True, "rho verification failed")
    require(rho_rows[0].get("published_q") == candidate_summary.get("published_q"), "rho target changed")
    ordinary_wall = ordinary_process["metrics"]["wall_seconds"]
    candidate_wall = candidate_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    comparison = {
        "rank_aware_direct_wall_speedup": ordinary_wall / candidate_wall,
        "rank_aware_collection_speedup": ordinary_summary["collection_ms"] / candidate_summary["collection_ms"],
        "rank_aware_core_seconds_ratio": candidate_process["metrics"]["total_core_seconds"] / ordinary_process["metrics"]["total_core_seconds"],
        "rank_aware_support_query_ratio": candidate_summary["support_queries"] / ordinary_summary["support_queries"],
        "rank_aware_over_rho_wall_ratio": candidate_wall / rho_wall,
        "ordinary_over_rho_wall_ratio": ordinary_wall / rho_wall,
        "fresh_build_plus_rank_aware_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + candidate_wall) / rho_wall,
        "rank_aware_whole_process_crossover": candidate_wall < rho_wall,
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_rank_aware_rho",
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-62 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-62 build seal", executable=False),
            "source_state": build["source_state"],
            "outer_resources": build["outer_resources"],
        },
        "constants": self_test(),
        "ordinary_process": ordinary_process,
        "candidate_process": candidate_process,
        "rho_process": rho_process,
        "outer_resources": outer,
        "ordinary": ordinary,
        "candidate": candidate,
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
    result = load(output / "result.json", "Stage-62 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_rank_aware_rho", "Stage-62 result is incomplete")
    ordinary = stage61.direct_observation(output / "evidence/ordinary-full-rank.stdout", "FULL_RANK")
    candidate = stage61.direct_observation(output / "evidence/rank-aware.stdout", "FULL_RANK")
    require(ordinary == result["ordinary"] and candidate == result["candidate"], "Stage-62 transcript changed")
    require(candidate["summary"]["factor_base_log_solution"] == ordinary["summary"]["factor_base_log_solution"], "Stage-62 solution changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-62 claim widened")
    return {
        "schema": "koblitz_stage62_verification.v1",
        "status": "n53_rank_aware_verified",
        "source_commit": build["source_state"]["commit"],
        "relations": candidate["summary"]["admitted_relations"],
        "target_trials": candidate["summary"]["target_trials"],
        "support_queries": candidate["summary"]["support_queries"],
        "rank_target_slot_index_allocated_bytes": candidate["summary"]["rank_target_slot_index_allocated_bytes"],
        "candidate_wall_seconds": result["candidate_process"]["metrics"]["wall_seconds"],
        "rho_wall_seconds": result["rho_process"]["metrics"]["wall_seconds"],
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage62Error, stage61.Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage62-rank-aware: {error}")


if __name__ == "__main__":
    main()

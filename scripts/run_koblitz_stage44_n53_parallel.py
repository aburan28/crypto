#!/usr/bin/env python3
"""Run the n=53 same-target sequential/parallel direct and rho control."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
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
import verify_koblitz_stage42_result as stage42_result


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
PREDECESSOR = STAGE / "stage-42-n53-same-target-result-20260912/verification.json"
RUN_SCHEMA = "koblitz_stage44_n53_parallel.v1"
RUN_SEAL_SCHEMA = "koblitz_stage44_n53_parallel_seal.v1"


class Stage44Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage44Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def predecessor_identity() -> dict[str, Any]:
    value = custody.executable_identity(PREDECESSOR, "Stage-42 hosted result", executable=False)
    value["path"] = str(PREDECESSOR.relative_to(REPO))
    return value


def self_test() -> dict[str, Any]:
    predecessor = stage42_result.verify(PREDECESSOR.parent)
    require(predecessor.get("status") == "hosted_n53_same_target_verified", "Stage-44 predecessor is not verified")
    require(predecessor["verification"]["published_fixture_scalar"] == stage42.EXPLICIT_SCALAR, "Stage-44 scalar changed")
    require(predecessor["verification"]["published_q"] == [2565091273463387, 5885236316843894], "Stage-44 target changed")
    return {
        "schema": "koblitz_stage44_self_test.v1",
        "status": "PASS",
        "checks": 11,
        "n": 53,
        "eta": [1, 128],
        "explicit_public_validation_scalar": stage42.EXPLICIT_SCALAR,
        "direct_seed": stage42.DIRECT_SEED,
        "rho_batch_seed": stage42.RHO_BATCH_SEED,
        "sequential_query_mode": "pair_pair_256",
        "parallel_query_mode": "pair_pair_parallel_512",
        "parallel_threads": 4,
        "same_target": True,
        "exact_relation_hash_equality_required": True,
    }


def parse_lines(path: Path) -> list[dict[str, Any]]:
    rows = []
    for line in path.read_text().splitlines():
        value = json.loads(line)
        require(isinstance(value, dict), f"non-object JSON line in {path}")
        rows.append(value)
    return rows


def direct_observation(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base = next(row for row in rows if row.get("kind") == "point_defined_factor_base")
    summary = next(row for row in rows if row.get("kind") == "relation_rank_summary")
    hashes = [row["relation_hash"] for row in rows if row.get("kind") == "relation_rank_receipt"]
    require(summary.get("status") == "RANK_PLUS_32" and summary.get("linear_solution_verified") is True, "Stage-44 direct endpoint failed")
    require(summary.get("all_relations_group_verified") is True and len(hashes) == summary.get("admitted_relations"), "Stage-44 relation transcript is incomplete")
    require(base.get("selection_uses_scalar_labels") is False, "Stage-44 factor base used scalar labels")
    return {"base": base, "summary": summary, "relation_hashes": hashes}


def semantic_direct(value: dict[str, Any]) -> dict[str, Any]:
    summary = value["summary"]
    return {
        "base_hash": value["base"]["base_hash"],
        "factor_base_points": summary["factor_base_points"],
        "orbit_columns": summary["orbit_columns"],
        "matrix_columns": summary["matrix_columns"],
        "published_q": summary["published_q"],
        "published_fixture_scalar": summary["published_fixture_scalar"],
        "recovered_fixture_scalar": summary["recovered_fixture_scalar"],
        "target_trials": summary["target_trials"],
        "admitted_relations": summary["admitted_relations"],
        "full_rank_at_relation": summary["full_rank_at_relation"],
        "factor_base_log_solution": summary["factor_base_log_solution"],
        "relation_hashes": value["relation_hashes"],
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-44 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-44 output must be new")
    build_root = args.build.resolve(strict=True)
    build = stage42.validate_build(build_root)
    output.mkdir(parents=True)
    evidence = output / "evidence"
    evidence.mkdir()
    direct = str((build_root / "bin/koblitz_rank_fixture").resolve(strict=True))
    rho = str((build_root / "bin/koblitz_rho_fixture").resolve(strict=True))
    common = ["53", "0", "1", "128", str(stage42.DIRECT_SEED), "signed_expanded", "independent"]
    sequential_command = [direct, *common, "pair_pair_256", "1", str(stage42.EXPLICIT_SCALAR)]
    parallel_command = [direct, *common, "pair_pair_parallel_512", "1", str(stage42.EXPLICIT_SCALAR)]
    rho_command = [rho, "53", "0", "signed_frobenius", "1", "packed", str(stage42.RHO_BATCH_SEED), str(stage42.EXPLICIT_SCALAR)]
    sequential_env = custody.safe_child_environment()
    sequential_env["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    parallel_env = dict(sequential_env)
    parallel_env["RAYON_NUM_THREADS"] = "4"
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        sequential_process = stage33.metered(role="sequential-direct", command=sequential_command, cwd=REPO, evidence=evidence, timeout=args.timeout, environment=sequential_env)
        parallel_process = stage33.metered(role="parallel-direct", command=parallel_command, cwd=REPO, evidence=evidence, timeout=args.timeout, environment=parallel_env)
        rho_process = stage33.metered(role="rho", command=rho_command, cwd=REPO, evidence=evidence, timeout=args.timeout, environment=sequential_env)
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    require(all(stage33.process_complete(row) for row in (sequential_process, parallel_process, rho_process)), "Stage-44 process failed")
    sequential = direct_observation(parse_lines(evidence / "sequential-direct.stdout"))
    parallel = direct_observation(parse_lines(evidence / "parallel-direct.stdout"))
    rho_rows = parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1, "Stage-44 rho row count changed")
    rho_row = rho_rows[0]
    sequential_semantic = semantic_direct(sequential)
    parallel_semantic = semantic_direct(parallel)
    require(sequential_semantic == parallel_semantic, "Stage-44 parallel relation transcript changed")
    require(sequential["summary"]["query_mode"] == "pair_pair_256", "Stage-44 sequential mode changed")
    require(parallel["summary"]["query_mode"] == "pair_pair_parallel_512" and parallel["summary"]["query_parallel_threads"] == 4, "Stage-44 parallel mode changed")
    require(sequential_semantic["published_q"] == rho_row.get("published_q"), "Stage-44 direct/rho target mismatch")
    require(rho_row.get("recovered_fixture_scalar") == sequential_semantic["recovered_fixture_scalar"] == stage42.EXPLICIT_SCALAR, "Stage-44 scalar recovery mismatch")
    require(rho_row.get("verified") is True and rho_row.get("reference_group_validation") is True, "Stage-44 rho verification failed")
    sequential_wall = sequential_process["metrics"]["wall_seconds"]
    parallel_wall = parallel_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    build_wall = build["outer_resources"]["wall_seconds"]
    comparison = {
        "parallel_direct_wall_speedup": sequential_wall / parallel_wall,
        "parallel_direct_core_seconds_ratio": parallel_process["metrics"]["total_core_seconds"] / sequential_process["metrics"]["total_core_seconds"],
        "parallel_collection_speedup": sequential["summary"]["collection_ms"] / parallel["summary"]["collection_ms"],
        "parallel_over_rho_wall_ratio": parallel_wall / rho_wall,
        "sequential_over_rho_wall_ratio": sequential_wall / rho_wall,
        "fresh_build_plus_parallel_over_rho_wall_ratio": (build_wall + parallel_wall) / rho_wall,
        "parallel_whole_process_crossover": parallel_wall < rho_wall,
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_sequential_parallel_rho",
        "build_binding": {"result": custody.executable_identity(build_root / "result.json", "Stage-44 build result", executable=False), "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-44 build seal", executable=False), "source_state": build["source_state"], "binaries": build["binaries"], "outer_resources": build["outer_resources"], "process": build["process"]},
        "constants": self_test(),
        "predecessor_stage42": {"identity": predecessor_identity(), "verification": load(PREDECESSOR, "Stage-42 result")["verification"]},
        "sequential_process": sequential_process,
        "parallel_process": parallel_process,
        "rho_process": rho_process,
        "outer_resources": outer,
        "semantic_direct": sequential_semantic,
        "sequential_summary": sequential["summary"],
        "parallel_summary": parallel["summary"],
        "rho": rho_row,
        "comparison": comparison,
        "charged_total_core_seconds_available": build["outer_resources"]["total_core_seconds"] + outer["total_core_seconds"],
        "charged_sequential_wall_seconds_available": build_wall + outer["wall_seconds"],
        "maximum_sampled_process_tree_rss_bytes": max(build["outer_resources"]["sampled_peak_process_tree_rss_bytes"], outer["sampled_peak_process_tree_rss_bytes"]),
        "same_public_target": True,
        "exact_relation_hash_equality": True,
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
        "fixture_scalar_used_by_collector": False,
        "sat_conflicts": None,
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
    result = load(output / "result.json", "Stage-44 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_sequential_parallel_rho", "Stage-44 result is incomplete")
    sequential = direct_observation(parse_lines(output / "evidence/sequential-direct.stdout"))
    parallel = direct_observation(parse_lines(output / "evidence/parallel-direct.stdout"))
    require(semantic_direct(sequential) == semantic_direct(parallel) == result["semantic_direct"], "Stage-44 relation equality changed")
    rho_rows = parse_lines(output / "evidence/rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0] == result["rho"], "Stage-44 rho evidence changed")
    sw = result["sequential_process"]["metrics"]["wall_seconds"]
    pw = result["parallel_process"]["metrics"]["wall_seconds"]
    rw = result["rho_process"]["metrics"]["wall_seconds"]
    expected = {
        "parallel_direct_wall_speedup": sw / pw,
        "parallel_direct_core_seconds_ratio": result["parallel_process"]["metrics"]["total_core_seconds"] / result["sequential_process"]["metrics"]["total_core_seconds"],
        "parallel_collection_speedup": sequential["summary"]["collection_ms"] / parallel["summary"]["collection_ms"],
        "parallel_over_rho_wall_ratio": pw / rw,
        "sequential_over_rho_wall_ratio": sw / rw,
        "fresh_build_plus_parallel_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + pw) / rw,
        "parallel_whole_process_crossover": pw < rw,
    }
    require(result["comparison"] == expected, "Stage-44 comparison changed")
    expected_core = build["outer_resources"]["total_core_seconds"] + result["outer_resources"]["total_core_seconds"]
    expected_wall = build["outer_resources"]["wall_seconds"] + result["outer_resources"]["wall_seconds"]
    expected_peak = max(build["outer_resources"]["sampled_peak_process_tree_rss_bytes"], result["outer_resources"]["sampled_peak_process_tree_rss_bytes"])
    require(result["charged_total_core_seconds_available"] == expected_core and result["charged_sequential_wall_seconds_available"] == expected_wall, "Stage-44 resource totals changed")
    require(result["maximum_sampled_process_tree_rss_bytes"] == expected_peak, "Stage-44 memory charge changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-44 claim widened")
    return {
        "schema": "koblitz_stage44_verification.v1",
        "status": "n53_parallel_same_target_verified",
        "source_commit": build["source_state"]["commit"],
        "published_q": result["semantic_direct"]["published_q"],
        "published_fixture_scalar": result["semantic_direct"]["published_fixture_scalar"],
        "relations": len(result["semantic_direct"]["relation_hashes"]),
        "sequential_direct_wall_seconds": sw,
        "parallel_direct_wall_seconds": pw,
        "rho_wall_seconds": rw,
        "parallel_direct_wall_speedup": expected["parallel_direct_wall_speedup"],
        "parallel_over_rho_wall_ratio": expected["parallel_over_rho_wall_ratio"],
        "parallel_direct_core_seconds_ratio": expected["parallel_direct_core_seconds_ratio"],
        "charged_total_core_seconds_available": result["charged_total_core_seconds_available"],
        "charged_sequential_wall_seconds_available": result["charged_sequential_wall_seconds_available"],
        "maximum_sampled_process_tree_rss_bytes": result["maximum_sampled_process_tree_rss_bytes"],
        "inventory_sha256": seal["inventory_sha256"],
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage44Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage44-n53-parallel: {error}")


if __name__ == "__main__":
    main()

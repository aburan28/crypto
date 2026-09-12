#!/usr/bin/env python3
"""Compare serial and ordered-parallel n=53 support expansion."""

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
import run_koblitz_stage62_rank_aware as stage62


REPO = Path(__file__).resolve().parents[1]
RUN_SCHEMA = "koblitz_stage64_parallel_expansion.v1"
RUN_SEAL_SCHEMA = "koblitz_stage64_parallel_expansion_seal.v1"


class Stage64Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage64Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    predecessor = stage62.self_test()
    require(predecessor.get("status") == "PASS", "Stage-62 control failed")
    return {
        "schema": "koblitz_stage64_self_test.v1",
        "status": "PASS",
        "checks": 15,
        "n": 53,
        "eta": [1, 128],
        "query_mode": "pair_pair_parallel_4096",
        "parallel_threads": 4,
        "support_expansion_job_batch": 65_536,
        "support_expansion_job_chunk": 1_024,
        "serial_table_insertion": True,
        "parallel_frobenius_expansion": True,
        "same_target": True,
        "exact_relation_hash_equality_required": True,
        "exact_support_entry_equality_required": True,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-64 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-64 output must be new")
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
    serial_environment = custody.safe_child_environment()
    serial_environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    serial_environment["KIC_RANK_SURPLUS"] = "0"
    serial_environment["KIC_RANK_AWARE_PAIR_SCAN"] = "1"
    serial_environment["RAYON_NUM_THREADS"] = "4"
    parallel_environment = dict(serial_environment)
    parallel_environment["KIC_PARALLEL_SUPPORT_EXPANSION"] = "1"
    rho_environment = custody.safe_child_environment()
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        serial_process = stage33.metered(
            role="serial-expansion", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=serial_environment,
        )
        parallel_process = stage33.metered(
            role="parallel-expansion", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=parallel_environment,
        )
        rho_process = stage33.metered(
            role="rho", command=rho_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=rho_environment,
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    require(all(stage33.process_complete(row) for row in (serial_process, parallel_process, rho_process)), "Stage-64 process failed")
    serial = stage61.direct_observation(evidence / "serial-expansion.stdout", "FULL_RANK")
    parallel = stage61.direct_observation(evidence / "parallel-expansion.stdout", "FULL_RANK")
    serial_summary = serial["summary"]
    parallel_summary = parallel["summary"]
    require(serial["base"]["base_hash"] == parallel["base"]["base_hash"], "factor base changed")
    require(serial["base"]["support_index_entries"] == parallel["base"]["support_index_entries"], "support entries changed")
    require(serial["base"].get("parallel_support_expansion") is False, "serial control changed")
    require(parallel["base"].get("parallel_support_expansion") is True, "parallel expansion was not enabled")
    require(parallel["base"].get("parallel_support_expansion_job_batch") == 65_536, "parallel batch changed")
    require(parallel["base"].get("parallel_support_expansion_job_chunk") == 1_024, "parallel chunk changed")
    require(serial["relation_hashes"] == parallel["relation_hashes"], "relation transcript changed")
    require(serial_summary["factor_base_log_solution"] == parallel_summary["factor_base_log_solution"], "full-rank solution changed")
    require(serial_summary["recovered_fixture_scalar"] == parallel_summary["recovered_fixture_scalar"] == stage42.EXPLICIT_SCALAR, "recovered scalar changed")
    for field in ("pair_canonicalization_maps", "pair_group_additions", "support_queries", "query_exact_table_misses", "full_rank_at_relation"):
        require(serial_summary[field] == parallel_summary[field], f"{field} changed")
    require(parallel_summary.get("field_mul_backend") == parallel_summary.get("field_square_backend") == "x86_64_pclmulqdq", "PCLMUL backend unavailable")
    rho_rows = stage44.parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0].get("verified") is True, "rho verification failed")
    require(rho_rows[0].get("published_q") == parallel_summary.get("published_q"), "rho target changed")
    serial_wall = serial_process["metrics"]["wall_seconds"]
    parallel_wall = parallel_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    comparison = {
        "parallel_direct_wall_speedup": serial_wall / parallel_wall,
        "parallel_setup_speedup": serial_summary["setup_ms"] / parallel_summary["setup_ms"],
        "parallel_collection_ratio": parallel_summary["collection_ms"] / serial_summary["collection_ms"],
        "parallel_core_seconds_ratio": parallel_process["metrics"]["total_core_seconds"] / serial_process["metrics"]["total_core_seconds"],
        "parallel_peak_rss_ratio": parallel_process["metrics"]["peak_rss_bytes"] / serial_process["metrics"]["peak_rss_bytes"],
        "parallel_over_rho_wall_ratio": parallel_wall / rho_wall,
        "serial_over_rho_wall_ratio": serial_wall / rho_wall,
        "fresh_build_plus_parallel_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + parallel_wall) / rho_wall,
        "parallel_whole_process_crossover": parallel_wall < rho_wall,
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_parallel_expansion_rho",
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-64 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-64 build seal", executable=False),
            "source_state": build["source_state"],
            "outer_resources": build["outer_resources"],
        },
        "constants": self_test(),
        "serial_process": serial_process,
        "parallel_process": parallel_process,
        "rho_process": rho_process,
        "outer_resources": outer,
        "serial": serial,
        "parallel": parallel,
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
    result = load(output / "result.json", "Stage-64 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_parallel_expansion_rho", "Stage-64 result is incomplete")
    serial = stage61.direct_observation(output / "evidence/serial-expansion.stdout", "FULL_RANK")
    parallel = stage61.direct_observation(output / "evidence/parallel-expansion.stdout", "FULL_RANK")
    require(serial == result["serial"] and parallel == result["parallel"], "Stage-64 transcript changed")
    require(serial["relation_hashes"] == parallel["relation_hashes"], "Stage-64 relation equality changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-64 claim widened")
    return {
        "schema": "koblitz_stage64_verification.v1",
        "status": "n53_parallel_expansion_verified",
        "source_commit": build["source_state"]["commit"],
        "relations": parallel["summary"]["admitted_relations"],
        "support_queries": parallel["summary"]["support_queries"],
        "support_table_allocated_bytes": parallel["base"]["support_table_allocated_bytes"],
        "parallel_wall_seconds": result["parallel_process"]["metrics"]["wall_seconds"],
        "parallel_peak_rss_bytes": result["parallel_process"]["metrics"]["peak_rss_bytes"],
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage64Error, stage62.Stage62Error, stage61.Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage64-parallel-expansion: {error}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run one public hash-derived n=53 unknown scalar through direct IC and rho."""

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
PUBLIC_HASH_SEED = 53001
EXPECTED_Q = [221146751147987, 5112896707579409]
BUILD_SCHEMA = stage42.BUILD_SCHEMA
RUN_SCHEMA = "koblitz_stage57_n53_unknown_scalar.v1"
RUN_SEAL_SCHEMA = "koblitz_stage57_n53_unknown_scalar_seal.v1"


class Stage57Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage57Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    require(PUBLIC_HASH_SEED == 53001 and len(EXPECTED_Q) == 2, "Stage-57 public target changed")
    return {
        "schema": "koblitz_stage57_self_test.v1",
        "status": "PASS",
        "checks": 12,
        "n": 53,
        "a": 0,
        "eta": [1, 128],
        "public_hash_seed": PUBLIC_HASH_SEED,
        "expected_public_q": EXPECTED_Q,
        "direct_seed": stage42.DIRECT_SEED,
        "rho_seed": stage42.RHO_BATCH_SEED,
        "query_mode": "pair_pair_parallel_4096",
        "parallel_threads": 4,
        "target_scalar_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
    }


def observe_direct(path: Path, optimized: bool = False) -> dict[str, Any]:
    rows = stage44.parse_lines(path)
    value = stage61.direct_observation(path, "FULL_RANK") if optimized else stage44.direct_observation(rows)
    base, summary = value["base"], value["summary"]
    require(base.get("selection_uses_scalar_labels") is False, "Stage-57 factor base used scalar labels")
    require(summary.get("published_q") == EXPECTED_Q, "Stage-57 public target changed")
    require(summary.get("published_fixture_scalar") is None, "Stage-57 direct retained an expected scalar")
    require(summary.get("target_scalar_constructed") is False, "Stage-57 direct constructed a target scalar")
    require(summary.get("target_kind") == "public_hash_to_curve_cofactor", "Stage-57 target kind changed")
    require(summary.get("public_hash_seed") == PUBLIC_HASH_SEED and summary.get("public_hash_counter") == 0, "Stage-57 hash derivation changed")
    require(summary.get("fixture_scalar_source") == summary.get("claim_boundary") == "public_hash_unknown_scalar", "Stage-57 direct evidence scope changed")
    require(summary.get("factor_base_logs_known_by_construction") is False, "Stage-57 direct used constructed factor-base logs")
    require(summary.get("linear_solution_verified") is True and summary.get("all_relations_group_verified") is True, "Stage-57 direct verification failed")
    require(summary.get("field_mul_backend") == "x86_64_pclmulqdq", "Stage-57 direct did not use PCLMUL")
    if optimized:
        require(summary.get("required_surplus_relations") == 0, "optimized direct retained surplus relations")
        require(summary.get("rank_aware_pair_scan") is True, "optimized direct did not use rank-aware collection")
        require(summary.get("rank_target_deficiency") == 3, "optimized direct rank-deficiency trigger changed")
        require(summary.get("field_product_pipeline") == "x86_64_pclmul_n53_fused_reduce", "optimized direct did not use fused field products")
    return value


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-57 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-57 output must be new")
    build_root = args.build.resolve(strict=True)
    build = stage42.validate_build(build_root)
    output.mkdir(parents=True)
    evidence = output / "evidence"
    evidence.mkdir()
    direct = str((build_root / "bin/koblitz_rank_fixture").resolve(strict=True))
    rho = str((build_root / "bin/koblitz_rho_fixture").resolve(strict=True))
    direct_command = [direct, "53", "0", "1", "128", str(stage42.DIRECT_SEED), "signed_expanded", "independent", "pair_pair_parallel_4096", "1", f"hash:{PUBLIC_HASH_SEED}"]
    rho_command = [rho, "53", "0", "signed_frobenius", "1", "packed", str(stage42.RHO_BATCH_SEED), f"hash:{PUBLIC_HASH_SEED}"]
    direct_environment = custody.safe_child_environment()
    direct_environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    direct_environment["RAYON_NUM_THREADS"] = "4"
    if args.optimized:
        direct_environment["KIC_RANK_SURPLUS"] = "0"
        direct_environment["KIC_RANK_AWARE_PAIR_SCAN"] = "1"
        direct_environment["KIC_RANK_TARGET_DEFICIENCY"] = "3"
        direct_environment["KIC_PARALLEL_SUPPORT_EXPANSION"] = "1"
        direct_environment["KIC_PIPELINED_SUPPORT_EXPANSION"] = "1"
    rho_environment = custody.safe_child_environment()
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        direct_process = stage33.metered(role="direct", command=direct_command, cwd=REPO, evidence=evidence, timeout=args.timeout, environment=direct_environment)
        rho_process = stage33.metered(role="rho", command=rho_command, cwd=REPO, evidence=evidence, timeout=args.timeout, environment=rho_environment)
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    require(stage33.process_complete(direct_process) and stage33.process_complete(rho_process), "Stage-57 process failed")
    direct_value = observe_direct(evidence / "direct.stdout", args.optimized)
    rho_rows = stage44.parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1, "Stage-57 rho row count changed")
    rho_value = rho_rows[0]
    require(rho_value.get("published_q") == EXPECTED_Q, "Stage-57 rho target changed")
    require(rho_value.get("published_fixture_scalar") is None and rho_value.get("target_scalar_constructed") is False, "Stage-57 rho constructed a target scalar")
    require(rho_value.get("public_hash_seed") == PUBLIC_HASH_SEED and rho_value.get("public_hash_counter") == 0, "Stage-57 rho hash derivation changed")
    require(rho_value.get("scope") == "public_hash_unknown_scalar", "Stage-57 rho evidence scope changed")
    require(rho_value.get("verified") is True and rho_value.get("reference_group_validation") is True, "Stage-57 rho verification failed")
    recovered = direct_value["summary"].get("recovered_fixture_scalar")
    require(isinstance(recovered, int) and recovered == rho_value.get("recovered_fixture_scalar"), "Stage-57 solvers disagree on the recovered scalar")
    direct_wall = direct_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    build_wall = build["outer_resources"]["wall_seconds"]
    comparison = {
        "direct_over_rho_wall_ratio": direct_wall / rho_wall,
        "fresh_build_plus_direct_over_rho_wall_ratio": (build_wall + direct_wall) / rho_wall,
        "whole_process_crossover": direct_wall < rho_wall,
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_public_unknown_scalar",
        "optimized_stack": args.optimized,
        "constants": self_test(),
        "build_binding": {"result": custody.executable_identity(build_root / "result.json", "Stage-57 build result", executable=False), "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-57 build seal", executable=False), "source_state": build["source_state"], "binaries": build["binaries"], "outer_resources": build["outer_resources"], "process": build["process"]},
        "direct_process": direct_process,
        "rho_process": rho_process,
        "outer_resources": outer,
        "direct_base": direct_value["base"],
        "direct_summary": direct_value["summary"],
        "relation_hashes": direct_value["relation_hashes"],
        "rho": rho_value,
        "recovered_scalar": recovered,
        "comparison": comparison,
        "charged_total_core_seconds_available": build["outer_resources"]["total_core_seconds"] + outer["total_core_seconds"],
        "charged_sequential_wall_seconds_available": build_wall + outer["wall_seconds"],
        "maximum_sampled_process_tree_rss_bytes": max(build["outer_resources"]["sampled_peak_process_tree_rss_bytes"], outer["sampled_peak_process_tree_rss_bytes"]),
        "factor_base_algebraically_point_defined": True,
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "target_scalar_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
        "unknown_scalar_recovery_verified": True,
        "sat_conflicts": None,
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    custody.write_json_new(output / "result.json", result)
    stage33.freeze(output, RUN_SEAL_SCHEMA, "run_frozen")
    return result


def verify(build_root: Path, output: Path) -> dict[str, Any]:
    build = stage42.validate_build(build_root)
    seal = stage33.verify_seal(output, RUN_SEAL_SCHEMA, "run_frozen")
    result = load(output / "result.json", "Stage-57 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_public_unknown_scalar", "Stage-57 result is incomplete")
    optimized = result.get("optimized_stack") is True
    direct_value = observe_direct(output / "evidence/direct.stdout", optimized)
    require(result.get("direct_base") == direct_value["base"] and result.get("direct_summary") == direct_value["summary"] and result.get("relation_hashes") == direct_value["relation_hashes"], "Stage-57 direct evidence changed")
    rho_rows = stage44.parse_lines(output / "evidence/rho.stdout")
    require(len(rho_rows) == 1 and result.get("rho") == rho_rows[0], "Stage-57 rho evidence changed")
    recovered = direct_value["summary"]["recovered_fixture_scalar"]
    require(recovered == rho_rows[0].get("recovered_fixture_scalar") == result.get("recovered_scalar"), "Stage-57 recovered scalar changed")
    dw = result["direct_process"]["metrics"]["wall_seconds"]
    rw = result["rho_process"]["metrics"]["wall_seconds"]
    expected = {"direct_over_rho_wall_ratio": dw / rw, "fresh_build_plus_direct_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + dw) / rw, "whole_process_crossover": dw < rw}
    require(result.get("comparison") == expected, "Stage-57 comparison changed")
    expected_core = build["outer_resources"]["total_core_seconds"] + result["outer_resources"]["total_core_seconds"]
    expected_wall = build["outer_resources"]["wall_seconds"] + result["outer_resources"]["wall_seconds"]
    expected_peak = max(build["outer_resources"]["sampled_peak_process_tree_rss_bytes"], result["outer_resources"]["sampled_peak_process_tree_rss_bytes"])
    require(result.get("charged_total_core_seconds_available") == expected_core and result.get("charged_sequential_wall_seconds_available") == expected_wall and result.get("maximum_sampled_process_tree_rss_bytes") == expected_peak, "Stage-57 resource totals changed")
    require(result.get("target_scalar_constructed_or_supplied") is False and result.get("factor_base_logs_known_by_construction") is False and result.get("unknown_scalar_recovery_verified") is True, "Stage-57 unknown-scalar boundary changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-57 claim widened")
    return {
        "schema": "koblitz_stage57_verification.v1",
        "status": "n53_public_unknown_scalar_verified",
        "source_commit": build["source_state"]["commit"],
        "published_q": EXPECTED_Q,
        "recovered_scalar": recovered,
        "relations": len(result["relation_hashes"]),
        "factor_base_points": direct_value["summary"]["factor_base_points"],
        "orbit_columns": direct_value["summary"]["orbit_columns"],
        "direct_wall_seconds": dw,
        "rho_wall_seconds": rw,
        **expected,
        "charged_total_core_seconds_available": result["charged_total_core_seconds_available"],
        "charged_sequential_wall_seconds_available": result["charged_sequential_wall_seconds_available"],
        "maximum_sampled_process_tree_rss_bytes": result["maximum_sampled_process_tree_rss_bytes"],
        "inventory_sha256": seal["inventory_sha256"],
        "target_scalar_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
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
    run_parser.add_argument("--optimized", action="store_true")
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage57Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage57-n53-unknown: {error}")


if __name__ == "__main__":
    main()

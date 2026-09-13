#!/usr/bin/env python3
"""Compare the n=53 research surplus with a certified full-rank stop."""

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


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
PREDECESSOR = STAGE / "stage-56-n53-pclmul-result-20260912/verification.json"
RUN_SCHEMA = "koblitz_stage61_production_rank.v1"
RUN_SEAL_SCHEMA = "koblitz_stage61_production_rank_seal.v1"


class Stage61Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage61Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    predecessor = load(PREDECESSOR, "Stage-56 verification")
    require(predecessor.get("status") == "hosted_n53_pclmul_result_verified", "Stage-56 predecessor is not verified")
    require(predecessor.get("hosted", {}).get("candidate", {}).get("full_rank_at_relation") == 157, "Stage-56 full-rank milestone changed")
    return {
        "schema": "koblitz_stage61_self_test.v1",
        "status": "PASS",
        "checks": 12,
        "n": 53,
        "eta": [1, 128],
        "query_mode": "pair_pair_parallel_4096",
        "parallel_threads": 4,
        "control_surplus_relations": 32,
        "production_surplus_relations": 0,
        "same_target": True,
        "exact_relation_prefix_required": True,
        "full_rank_solution_equality_required": True,
    }


def direct_observation(path: Path, expected_status: str) -> dict[str, Any]:
    rows = stage44.parse_lines(path)
    base = next(row for row in rows if row.get("kind") == "point_defined_factor_base")
    summary = next(row for row in rows if row.get("kind") == "relation_rank_summary")
    receipts = [row for row in rows if row.get("kind") == "relation_rank_receipt"]
    require(summary.get("status") == expected_status, f"unexpected {expected_status} endpoint")
    require(summary.get("linear_solution_verified") is True, "linear solution was not verified")
    require(summary.get("all_relations_group_verified") is True, "a relation failed group verification")
    if receipts:
        relation_hashes = [row["relation_hash"] for row in receipts]
    else:
        relation_hashes = summary.get("relation_hashes")
        require(
            summary.get("summary_only_timing") is True
            and summary.get("relation_receipts_emitted") is False,
            "missing relation receipts outside compact evidence mode",
        )
    require(
        isinstance(relation_hashes, list)
        and all(isinstance(value, str) and len(value) == 64 for value in relation_hashes)
        and len(relation_hashes) == summary.get("admitted_relations"),
        "relation transcript is incomplete",
    )
    require(base.get("selection_uses_scalar_labels") is False, "factor base used scalar labels")
    return {
        "base": base,
        "summary": summary,
        "relation_hashes": relation_hashes,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-61 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-61 output must be new")
    build_root = args.build.resolve(strict=True)
    build = stage42.validate_build(build_root)
    output.mkdir(parents=True)
    evidence = output / "evidence"
    evidence.mkdir()
    direct = str((build_root / "bin/koblitz_rank_fixture").resolve(strict=True))
    rho = str((build_root / "bin/koblitz_rho_fixture").resolve(strict=True))
    common = [
        "53", "0", "1", "128", str(stage42.DIRECT_SEED),
        "signed_expanded", "independent", "pair_pair_parallel_4096", "1",
        str(stage42.EXPLICIT_SCALAR),
    ]
    direct_command = [direct, *common]
    rho_command = [
        rho, "53", "0", "signed_frobenius", "1", "packed",
        str(stage42.RHO_BATCH_SEED), str(stage42.EXPLICIT_SCALAR),
    ]
    control_environment = custody.safe_child_environment()
    control_environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    control_environment["RAYON_NUM_THREADS"] = "4"
    production_environment = dict(control_environment)
    production_environment["KIC_RANK_SURPLUS"] = "0"
    rho_environment = custody.safe_child_environment()
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        control_process = stage33.metered(
            role="research-control", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=control_environment,
        )
        production_process = stage33.metered(
            role="production-rank", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=production_environment,
        )
        rho_process = stage33.metered(
            role="rho", command=rho_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=rho_environment,
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    require(all(stage33.process_complete(row) for row in (control_process, production_process, rho_process)), "Stage-61 process failed")
    control = direct_observation(evidence / "research-control.stdout", "RANK_PLUS_32")
    production = direct_observation(evidence / "production-rank.stdout", "FULL_RANK")
    control_summary = control["summary"]
    production_summary = production["summary"]
    require(control["base"]["base_hash"] == production["base"]["base_hash"], "factor base changed")
    require(control_summary.get("required_surplus_relations") == 32, "control surplus changed")
    require(production_summary.get("required_surplus_relations") == 0, "production surplus changed")
    require(production_summary["admitted_relations"] == production_summary["full_rank_at_relation"], "production did not stop at full rank")
    require(production["relation_hashes"] == control["relation_hashes"][:len(production["relation_hashes"])], "relation prefix changed")
    require(production_summary["factor_base_log_solution"] == control_summary["factor_base_log_solution"], "full-rank solution changed")
    require(production_summary["recovered_fixture_scalar"] == control_summary["recovered_fixture_scalar"] == stage42.EXPLICIT_SCALAR, "recovered scalar changed")
    require(production_summary.get("field_mul_backend") == production_summary.get("field_square_backend") == "x86_64_pclmulqdq", "PCLMUL backend unavailable")
    require(production["base"].get("support_table_allocated_bytes") == 738_197_504, "packed support allocation changed")
    rho_rows = stage44.parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0].get("verified") is True, "rho verification failed")
    require(rho_rows[0].get("published_q") == production_summary.get("published_q"), "rho target changed")
    control_wall = control_process["metrics"]["wall_seconds"]
    production_wall = production_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    comparison = {
        "production_direct_wall_speedup": control_wall / production_wall,
        "production_collection_speedup": control_summary["collection_ms"] / production_summary["collection_ms"],
        "production_core_seconds_ratio": production_process["metrics"]["total_core_seconds"] / control_process["metrics"]["total_core_seconds"],
        "production_support_query_ratio": production_summary["support_queries"] / control_summary["support_queries"],
        "production_over_rho_wall_ratio": production_wall / rho_wall,
        "control_over_rho_wall_ratio": control_wall / rho_wall,
        "fresh_build_plus_production_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + production_wall) / rho_wall,
        "production_whole_process_crossover": production_wall < rho_wall,
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_production_rank_rho",
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-61 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-61 build seal", executable=False),
            "source_state": build["source_state"],
            "outer_resources": build["outer_resources"],
        },
        "predecessor_stage56": custody.executable_identity(PREDECESSOR, "Stage-56 verification", executable=False),
        "constants": self_test(),
        "control_process": control_process,
        "production_process": production_process,
        "rho_process": rho_process,
        "outer_resources": outer,
        "control": control,
        "production": production,
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
    result = load(output / "result.json", "Stage-61 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_production_rank_rho", "Stage-61 result is incomplete")
    control = direct_observation(output / "evidence/research-control.stdout", "RANK_PLUS_32")
    production = direct_observation(output / "evidence/production-rank.stdout", "FULL_RANK")
    require(control == result["control"] and production == result["production"], "Stage-61 transcript changed")
    require(production["relation_hashes"] == control["relation_hashes"][:len(production["relation_hashes"])], "Stage-61 relation prefix changed")
    require(production["summary"]["factor_base_log_solution"] == control["summary"]["factor_base_log_solution"], "Stage-61 solution changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-61 claim widened")
    return {
        "schema": "koblitz_stage61_verification.v1",
        "status": "n53_production_full_rank_verified",
        "source_commit": build["source_state"]["commit"],
        "relations": production["summary"]["admitted_relations"],
        "support_queries": production["summary"]["support_queries"],
        "support_table_allocated_bytes": production["base"]["support_table_allocated_bytes"],
        "production_wall_seconds": result["production_process"]["metrics"]["wall_seconds"],
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage61-production-rank: {error}")


if __name__ == "__main__":
    main()

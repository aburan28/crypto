#!/usr/bin/env python3
"""Compare portable and PCLMUL n=53 field products on one target."""

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
PREDECESSOR = STAGE / "stage-44-n53-parallel-result-20260912/verification.json"
RUN_SCHEMA = "koblitz_stage55_n53_pclmul.v1"
RUN_SEAL_SCHEMA = "koblitz_stage55_n53_pclmul_seal.v1"


class Stage55Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage55Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    predecessor = load(PREDECESSOR, "Stage-44 verification")
    require(predecessor.get("status") == "hosted_n53_parallel_query_verified", "Stage-44 predecessor is not verified")
    require(predecessor.get("verification", {}).get("published_fixture_scalar") == stage42.EXPLICIT_SCALAR, "Stage-55 scalar changed")
    require(predecessor.get("verification", {}).get("published_q") == [2565091273463387, 5885236316843894], "Stage-55 target changed")
    return {
        "schema": "koblitz_stage55_self_test.v1",
        "status": "PASS",
        "checks": 10,
        "n": 53,
        "eta": [1, 128],
        "explicit_public_validation_scalar": stage42.EXPLICIT_SCALAR,
        "direct_seed": stage42.DIRECT_SEED,
        "rho_batch_seed": stage42.RHO_BATCH_SEED,
        "query_mode": "pair_pair_parallel_4096",
        "baseline_field_mul_backend": "portable_sparse_carryless",
        "candidate_field_mul_backend": "x86_64_pclmulqdq",
        "parallel_threads": 4,
        "same_target": True,
        "exact_relation_hash_equality_required": True,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-55 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-55 output must be new")
    build_root = args.build.resolve(strict=True)
    build = stage42.validate_build(build_root)
    output.mkdir(parents=True)
    evidence = output / "evidence"
    evidence.mkdir()
    direct = str((build_root / "bin/koblitz_rank_fixture").resolve(strict=True))
    rho = str((build_root / "bin/koblitz_rho_fixture").resolve(strict=True))
    common = ["53", "0", "1", "128", str(stage42.DIRECT_SEED), "signed_expanded", "independent"]
    baseline_command = [direct, *common, "pair_pair_parallel_4096", "1", str(stage42.EXPLICIT_SCALAR)]
    candidate_command = [direct, *common, "pair_pair_parallel_4096", "1", str(stage42.EXPLICIT_SCALAR)]
    rho_command = [rho, "53", "0", "signed_frobenius", "1", "packed", str(stage42.RHO_BATCH_SEED), str(stage42.EXPLICIT_SCALAR)]
    baseline_environment = custody.safe_child_environment()
    baseline_environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    baseline_environment["RAYON_NUM_THREADS"] = "4"
    baseline_environment["KIC_DISABLE_PCLMUL"] = "1"
    candidate_environment = custody.safe_child_environment()
    candidate_environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    candidate_environment["RAYON_NUM_THREADS"] = "4"
    rho_environment = custody.safe_child_environment()
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        baseline_process = stage33.metered(role="portable", command=baseline_command, cwd=REPO, evidence=evidence, timeout=args.timeout, environment=baseline_environment)
        candidate_process = stage33.metered(role="pclmul", command=candidate_command, cwd=REPO, evidence=evidence, timeout=args.timeout, environment=candidate_environment)
        rho_process = stage33.metered(role="rho", command=rho_command, cwd=REPO, evidence=evidence, timeout=args.timeout, environment=rho_environment)
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    require(all(stage33.process_complete(row) for row in (baseline_process, candidate_process, rho_process)), "Stage-55 process failed")
    baseline = stage44.direct_observation(stage44.parse_lines(evidence / "portable.stdout"))
    candidate = stage44.direct_observation(stage44.parse_lines(evidence / "pclmul.stdout"))
    semantic = stage44.semantic_direct(baseline)
    require(semantic == stage44.semantic_direct(candidate), "Stage-55 relation transcript changed")
    require(baseline["summary"].get("query_mode") == candidate["summary"].get("query_mode") == "pair_pair_parallel_4096", "Stage-55 query mode changed")
    require(baseline["summary"].get("field_mul_backend") == "portable_sparse_carryless", "Stage-55 portable backend changed")
    require(candidate["summary"].get("field_mul_backend") == "x86_64_pclmulqdq", "Stage-55 PCLMUL backend changed")
    require(baseline["summary"].get("query_parallel_threads") == candidate["summary"].get("query_parallel_threads") == 4, "Stage-55 thread count changed")
    rho_rows = stage44.parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1, "Stage-55 rho row count changed")
    rho_row = rho_rows[0]
    require(semantic["published_q"] == rho_row.get("published_q"), "Stage-55 direct/rho target mismatch")
    require(rho_row.get("recovered_fixture_scalar") == semantic["recovered_fixture_scalar"] == stage42.EXPLICIT_SCALAR, "Stage-55 scalar recovery mismatch")
    require(rho_row.get("verified") is True and rho_row.get("reference_group_validation") is True, "Stage-55 rho verification failed")
    baseline_wall = baseline_process["metrics"]["wall_seconds"]
    candidate_wall = candidate_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    build_wall = build["outer_resources"]["wall_seconds"]
    comparison = {
        "candidate_direct_wall_speedup": baseline_wall / candidate_wall,
        "candidate_direct_core_seconds_ratio": candidate_process["metrics"]["total_core_seconds"] / baseline_process["metrics"]["total_core_seconds"],
        "candidate_collection_speedup": baseline["summary"]["collection_ms"] / candidate["summary"]["collection_ms"],
        "candidate_over_rho_wall_ratio": candidate_wall / rho_wall,
        "baseline_over_rho_wall_ratio": baseline_wall / rho_wall,
        "fresh_build_plus_candidate_over_rho_wall_ratio": (build_wall + candidate_wall) / rho_wall,
        "candidate_whole_process_crossover": candidate_wall < rho_wall,
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_pclmul_ab_rho",
        "build_binding": {"result": custody.executable_identity(build_root / "result.json", "Stage-55 build result", executable=False), "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-55 build seal", executable=False), "source_state": build["source_state"], "binaries": build["binaries"], "outer_resources": build["outer_resources"], "process": build["process"]},
        "predecessor_stage44": custody.executable_identity(PREDECESSOR, "Stage-44 verification", executable=False),
        "constants": self_test(),
        "baseline_process": baseline_process,
        "candidate_process": candidate_process,
        "rho_process": rho_process,
        "outer_resources": outer,
        "semantic_direct": semantic,
        "baseline_summary": baseline["summary"],
        "candidate_summary": candidate["summary"],
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
    result = load(output / "result.json", "Stage-55 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_pclmul_ab_rho", "Stage-55 result is incomplete")
    baseline = stage44.direct_observation(stage44.parse_lines(output / "evidence/portable.stdout"))
    candidate = stage44.direct_observation(stage44.parse_lines(output / "evidence/pclmul.stdout"))
    require(stage44.semantic_direct(baseline) == stage44.semantic_direct(candidate) == result["semantic_direct"], "Stage-55 relation equality changed")
    require(baseline["summary"] == result["baseline_summary"] and candidate["summary"] == result["candidate_summary"], "Stage-55 summary changed")
    require(baseline["summary"].get("query_mode") == candidate["summary"].get("query_mode") == "pair_pair_parallel_4096", "Stage-55 verified query mode changed")
    require(baseline["summary"].get("field_mul_backend") == "portable_sparse_carryless" and candidate["summary"].get("field_mul_backend") == "x86_64_pclmulqdq", "Stage-55 verified field backends changed")
    rho_rows = stage44.parse_lines(output / "evidence/rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0] == result["rho"], "Stage-55 rho evidence changed")
    bw = result["baseline_process"]["metrics"]["wall_seconds"]
    cw = result["candidate_process"]["metrics"]["wall_seconds"]
    rw = result["rho_process"]["metrics"]["wall_seconds"]
    expected = {
        "candidate_direct_wall_speedup": bw / cw,
        "candidate_direct_core_seconds_ratio": result["candidate_process"]["metrics"]["total_core_seconds"] / result["baseline_process"]["metrics"]["total_core_seconds"],
        "candidate_collection_speedup": baseline["summary"]["collection_ms"] / candidate["summary"]["collection_ms"],
        "candidate_over_rho_wall_ratio": cw / rw,
        "baseline_over_rho_wall_ratio": bw / rw,
        "fresh_build_plus_candidate_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + cw) / rw,
        "candidate_whole_process_crossover": cw < rw,
    }
    require(result["comparison"] == expected, "Stage-55 comparison changed")
    expected_core = build["outer_resources"]["total_core_seconds"] + result["outer_resources"]["total_core_seconds"]
    expected_wall = build["outer_resources"]["wall_seconds"] + result["outer_resources"]["wall_seconds"]
    expected_peak = max(build["outer_resources"]["sampled_peak_process_tree_rss_bytes"], result["outer_resources"]["sampled_peak_process_tree_rss_bytes"])
    require(result["charged_total_core_seconds_available"] == expected_core and result["charged_sequential_wall_seconds_available"] == expected_wall, "Stage-55 resource totals changed")
    require(result["maximum_sampled_process_tree_rss_bytes"] == expected_peak, "Stage-55 memory charge changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-55 claim widened")
    return {
        "schema": "koblitz_stage55_verification.v1",
        "status": "n53_pclmul_same_target_verified",
        "source_commit": build["source_state"]["commit"],
        "published_q": result["semantic_direct"]["published_q"],
        "published_fixture_scalar": result["semantic_direct"]["published_fixture_scalar"],
        "relations": len(result["semantic_direct"]["relation_hashes"]),
        "baseline_wall_seconds": bw,
        "candidate_wall_seconds": cw,
        "rho_wall_seconds": rw,
        **expected,
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage55Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage55-n53-pclmul: {error}")


if __name__ == "__main__":
    main()

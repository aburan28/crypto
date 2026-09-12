#!/usr/bin/env python3
"""Compare generic and direct-PCLMUL n=53 initial support pair sums."""

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
import run_koblitz_stage64_parallel_expansion as stage64
import run_koblitz_stage65_n53_reduction as stage65


REPO = Path(__file__).resolve().parents[1]
RUN_SCHEMA = "koblitz_stage81_specialized_pair_sums.v1"
RUN_SEAL_SCHEMA = "koblitz_stage81_specialized_pair_sums_seal.v1"


class Stage81Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage81Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    predecessor = stage65.self_test()
    require(predecessor.get("status") == "PASS", "Stage-65 control failed")
    return {
        "schema": "koblitz_stage81_self_test.v1",
        "status": "PASS",
        "checks": 24,
        "n": 53,
        "field_modulus": "x^53 + x^6 + x^2 + x + 1",
        "baseline_pair_sums": "generic_batch_add",
        "candidate_pair_sums": "direct_pclmul_n53_batch_add",
        "parallel_threads": 4,
        "same_target": True,
        "exact_relation_hash_equality_required": True,
        "exact_support_counter_equality_required": True,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-81 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-81 output must be new")
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
    candidate_environment = custody.safe_child_environment()
    candidate_environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    candidate_environment["KIC_RANK_SURPLUS"] = "0"
    candidate_environment["KIC_RANK_AWARE_PAIR_SCAN"] = "1"
    candidate_environment["KIC_PARALLEL_SUPPORT_EXPANSION"] = "1"
    candidate_environment["KIC_PIPELINED_SUPPORT_EXPANSION"] = "1"
    candidate_environment["RAYON_NUM_THREADS"] = "4"
    baseline_environment = dict(candidate_environment)
    candidate_environment["KIC_ENABLE_SPECIALIZED_N53_PAIR_SUMS"] = "1"
    rho_environment = custody.safe_child_environment()
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        baseline_process = stage33.metered(
            role="generic-pair-sums", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=baseline_environment,
        )
        candidate_process = stage33.metered(
            role="specialized-pair-sums", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=candidate_environment,
        )
        rho_process = stage33.metered(
            role="rho", command=rho_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=rho_environment,
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    require(all(stage33.process_complete(row) for row in (baseline_process, candidate_process, rho_process)), "Stage-81 process failed")
    baseline = stage61.direct_observation(evidence / "generic-pair-sums.stdout", "FULL_RANK")
    candidate = stage61.direct_observation(evidence / "specialized-pair-sums.stdout", "FULL_RANK")
    baseline_summary = baseline["summary"]
    candidate_summary = candidate["summary"]
    require(baseline["base"]["base_hash"] == candidate["base"]["base_hash"], "factor base changed")
    require(baseline["base"]["support_index_entries"] == candidate["base"]["support_index_entries"], "support entries changed")
    require(baseline["base"].get("field_reduction_backend") == candidate["base"].get("field_reduction_backend") == "n53_fixed_two_fold", "fixed reducer changed")
    require(
        baseline["base"].get("field_product_pipeline")
        == candidate["base"].get("field_product_pipeline")
        == "x86_64_pclmul_n53_fused_reduce",
        "matched fused field pipeline was not selected",
    )
    require(baseline["relation_hashes"] == candidate["relation_hashes"], "relation transcript changed")
    require(baseline_summary["factor_base_log_solution"] == candidate_summary["factor_base_log_solution"], "full-rank solution changed")
    require(baseline_summary["recovered_fixture_scalar"] == candidate_summary["recovered_fixture_scalar"] == stage42.EXPLICIT_SCALAR, "recovered scalar changed")
    for field in ("pair_canonicalization_maps", "pair_group_additions", "support_queries", "full_rank_at_relation"):
        require(baseline_summary[field] == candidate_summary[field], f"{field} changed")
    require(
        baseline["base"].get("support_x_prefilter_hash_strategy")
        == candidate["base"].get("support_x_prefilter_hash_strategy")
        == "single_mix_split_29_bit_indices",
        "matched split-index filter changed",
    )
    require(
        baseline["base"].get("support_x_prefilter_insert_hash_reuse")
        == candidate["base"].get("support_x_prefilter_insert_hash_reuse")
        is True,
        "matched insert-hash reuse changed",
    )
    require(
        baseline_summary.get("query_dual_sign_denominator_sharing")
        == candidate_summary.get("query_dual_sign_denominator_sharing")
        is True,
        "matched dual-sign query kernel changed",
    )
    require(
        baseline_summary.get("query_shared_sign_denominator_inputs")
        == candidate_summary.get("query_shared_sign_denominator_inputs")
        > 0,
        "shared denominator count changed",
    )
    require(baseline_summary["query_batch_inversions"] == candidate_summary["query_batch_inversions"], "batch-inversion count changed")
    require(baseline_summary["query_exact_table_misses"] == candidate_summary["query_exact_table_misses"], "exact lookup count changed")
    require(
        baseline_summary.get("query_compact_pair_scratch")
        == candidate_summary.get("query_compact_pair_scratch")
        is True,
        "matched compact query scratch changed",
    )
    require(
        baseline_summary.get("query_raw_point_scratch_bytes")
        == candidate_summary.get("query_raw_point_scratch_bytes")
        == 24,
        "raw scratch element size changed",
    )
    require(
        baseline_summary.get("query_compact_point_scratch_bytes")
        == candidate_summary.get("query_compact_point_scratch_bytes")
        == 16,
        "compact scratch element size changed",
    )
    require(
        baseline["base"].get("field_product_dispatch")
        == candidate["base"].get("field_product_dispatch")
        == "inline_combined_n53_predicate",
        "matched combined field dispatch changed",
    )
    require(
        baseline_summary.get("query_specialized_n53_pair_batch")
        == candidate_summary.get("query_specialized_n53_pair_batch")
        is True,
        "matched specialized pair-query batch changed",
    )
    require(
        baseline["base"].get("specialized_n53_support_expansion")
        == candidate["base"].get("specialized_n53_support_expansion")
        is False,
        "matched generic support expansion changed",
    )
    require(
        baseline["base"].get("specialized_n53_pair_sums") is False,
        "generic pair-sum baseline was not selected",
    )
    require(
        candidate["base"].get("specialized_n53_pair_sums") is True,
        "specialized n53 pair sums were not selected",
    )
    require(candidate_summary.get("field_mul_backend") == candidate_summary.get("field_square_backend") == "x86_64_pclmulqdq", "PCLMUL backend unavailable")
    rho_rows = stage44.parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0].get("verified") is True, "rho verification failed")
    require(rho_rows[0].get("published_q") == candidate_summary.get("published_q"), "rho target changed")
    baseline_wall = baseline_process["metrics"]["wall_seconds"]
    candidate_wall = candidate_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    comparison = {
        "specialized_pair_sums_direct_wall_speedup": baseline_wall / candidate_wall,
        "specialized_pair_sums_setup_speedup": baseline_summary["setup_ms"] / candidate_summary["setup_ms"],
        "specialized_pair_sums_collection_speedup": baseline_summary["collection_ms"] / candidate_summary["collection_ms"],
        "specialized_pair_sums_core_seconds_ratio": candidate_process["metrics"]["total_core_seconds"] / baseline_process["metrics"]["total_core_seconds"],
        "specialized_pair_sums_peak_rss_ratio": candidate_process["metrics"]["peak_rss_bytes"] / baseline_process["metrics"]["peak_rss_bytes"],
        "specialized_pair_sums_over_rho_wall_ratio": candidate_wall / rho_wall,
        "generic_pair_sums_over_rho_wall_ratio": baseline_wall / rho_wall,
        "fresh_build_plus_specialized_pair_sums_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + candidate_wall) / rho_wall,
        "specialized_pair_sums_whole_process_crossover": candidate_wall < rho_wall,
        "baseline_exact_table_misses": baseline_summary["query_exact_table_misses"],
        "candidate_exact_table_misses": candidate_summary["query_exact_table_misses"],
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_specialized_pair_sums_rho",
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-81 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-81 build seal", executable=False),
            "source_state": build["source_state"],
            "outer_resources": build["outer_resources"],
        },
        "constants": self_test(),
        "baseline_process": baseline_process,
        "candidate_process": candidate_process,
        "rho_process": rho_process,
        "outer_resources": outer,
        "baseline": baseline,
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
    result = load(output / "result.json", "Stage-81 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_specialized_pair_sums_rho", "Stage-81 result is incomplete")
    baseline = stage61.direct_observation(output / "evidence/generic-pair-sums.stdout", "FULL_RANK")
    candidate = stage61.direct_observation(output / "evidence/specialized-pair-sums.stdout", "FULL_RANK")
    require(baseline == result["baseline"] and candidate == result["candidate"], "Stage-81 transcript changed")
    require(baseline["relation_hashes"] == candidate["relation_hashes"], "Stage-81 relation equality changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-81 claim widened")
    return {
        "schema": "koblitz_stage81_verification.v1",
        "status": "n53_specialized_pair_sums_verified",
        "source_commit": build["source_state"]["commit"],
        "relations": candidate["summary"]["admitted_relations"],
        "support_queries": candidate["summary"]["support_queries"],
        "candidate_wall_seconds": result["candidate_process"]["metrics"]["wall_seconds"],
        "candidate_peak_rss_bytes": result["candidate_process"]["metrics"]["peak_rss_bytes"],
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage81Error, stage65.Stage65Error, stage64.Stage64Error, stage61.Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage81-specialized-pair-sums: {error}")


if __name__ == "__main__":
    main()

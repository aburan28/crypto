#!/usr/bin/env python3
"""Compare n=53 exact lookup, x-prefilter, and batch-inversion strategies."""

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
RUN_SCHEMA = "koblitz_stage94_filter_inverse.v1"
RUN_SEAL_SCHEMA = "koblitz_stage94_filter_inverse_seal.v1"


class Stage94Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage94Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    return {
        "schema": "koblitz_stage94_self_test.v1",
        "status": "PASS",
        "checks": 26,
        "n": 53,
        "field_modulus": "x^53 + x^6 + x^2 + x + 1",
        "arms": [
            "mixed_hash_binary_inverse",
            "mixed_hash_prefiltered_lookup",
            "direct_x_bits_binary_inverse",
            "direct_x_bits_itoh_tsujii_inverse",
        ],
        "parallel_threads": 4,
        "same_target": True,
        "exact_relation_hash_equality_required": True,
        "exact_support_counter_equality_required": True,
        "itoh_addition_chain": [1, 2, 4, 8, 16, 32, 48, 52],
        "binary_inverse_field_products": 105,
        "itoh_inverse_field_products": 59,
    }


def selected_environment() -> dict[str, str]:
    environment = custody.safe_child_environment()
    environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    environment["KIC_RANK_SURPLUS"] = "0"
    environment["KIC_RANK_AWARE_PAIR_SCAN"] = "1"
    environment["KIC_PARALLEL_SUPPORT_EXPANSION"] = "1"
    environment["KIC_PIPELINED_SUPPORT_EXPANSION"] = "1"
    environment["RAYON_NUM_THREADS"] = "4"
    return environment


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-94 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-94 output must be new")
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
    mixed_environment = selected_environment()
    prefiltered_environment = dict(mixed_environment)
    prefiltered_environment["KIC_ENABLE_PREFILTERED_EXACT_LOOKUP"] = "1"
    direct_bits_environment = dict(prefiltered_environment)
    direct_bits_environment["KIC_ENABLE_DIRECT_X_FILTER_BITS"] = "1"
    itoh_environment = dict(direct_bits_environment)
    itoh_environment["KIC_ENABLE_ITOH_N53_INVERSE"] = "1"
    rho_environment = custody.safe_child_environment()

    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        mixed_process = stage33.metered(
            role="mixed-hash-binary", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=mixed_environment,
        )
        prefiltered_process = stage33.metered(
            role="mixed-hash-prefiltered", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=prefiltered_environment,
        )
        direct_bits_process = stage33.metered(
            role="direct-bits-binary", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=direct_bits_environment,
        )
        itoh_process = stage33.metered(
            role="direct-bits-itoh", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=itoh_environment,
        )
        rho_process = stage33.metered(
            role="rho", command=rho_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=rho_environment,
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    processes = {
        "mixed_hash_binary": mixed_process,
        "mixed_hash_prefiltered": prefiltered_process,
        "direct_bits_binary": direct_bits_process,
        "direct_bits_itoh": itoh_process,
        "rho": rho_process,
    }
    require(all(stage33.process_complete(row) for row in processes.values()), "Stage-94 process failed")

    arms = {
        "mixed_hash_binary": stage61.direct_observation(evidence / "mixed-hash-binary.stdout", "FULL_RANK"),
        "mixed_hash_prefiltered": stage61.direct_observation(evidence / "mixed-hash-prefiltered.stdout", "FULL_RANK"),
        "direct_bits_binary": stage61.direct_observation(evidence / "direct-bits-binary.stdout", "FULL_RANK"),
        "direct_bits_itoh": stage61.direct_observation(evidence / "direct-bits-itoh.stdout", "FULL_RANK"),
    }
    mixed = arms["mixed_hash_binary"]
    prefiltered = arms["mixed_hash_prefiltered"]
    direct_bits = arms["direct_bits_binary"]
    itoh = arms["direct_bits_itoh"]
    summaries = {name: value["summary"] for name, value in arms.items()}
    bases = {name: value["base"] for name, value in arms.items()}

    require(len({value["base_hash"] for value in bases.values()}) == 1, "factor base changed")
    require(len({value["support_index_entries"] for value in bases.values()}) == 1, "support entries changed")
    require(mixed["relation_hashes"] == prefiltered["relation_hashes"] == direct_bits["relation_hashes"] == itoh["relation_hashes"], "relation transcript changed")
    require(
        summaries["mixed_hash_binary"]["factor_base_log_solution"]
        == summaries["mixed_hash_prefiltered"]["factor_base_log_solution"]
        == summaries["direct_bits_binary"]["factor_base_log_solution"]
        == summaries["direct_bits_itoh"]["factor_base_log_solution"],
        "full-rank solution changed",
    )
    for field in (
        "pair_canonicalization_maps", "pair_group_additions", "support_queries",
        "full_rank_at_relation", "query_batch_inversions", "query_group_additions",
    ):
        require(len({value[field] for value in summaries.values()}) == 1, f"{field} changed")
    require(
        bases["mixed_hash_binary"].get("support_x_prefilter_hash_strategy")
        == "single_mix_split_29_bit_indices",
        "mixed-hash baseline changed",
    )
    require(
        bases["direct_bits_binary"].get("support_x_prefilter_hash_strategy")
        == bases["direct_bits_itoh"].get("support_x_prefilter_hash_strategy")
        == "direct_low_and_high_x_bit_windows",
        "direct-bit filter was not selected",
    )
    require(bases["mixed_hash_binary"].get("support_x_prefilter_direct_bits") is False, "baseline direct-bit flag changed")
    require(bases["mixed_hash_prefiltered"].get("support_x_prefilter_direct_bits") is False, "prefiltered direct-bit flag changed")
    require(bases["direct_bits_binary"].get("support_x_prefilter_direct_bits") is True, "direct-bit flag missing")
    require(summaries["mixed_hash_binary"].get("query_itoh_n53_inverse") is False, "baseline inverse changed")
    require(summaries["mixed_hash_prefiltered"].get("query_itoh_n53_inverse") is False, "prefiltered inverse changed")
    require(summaries["direct_bits_binary"].get("query_itoh_n53_inverse") is False, "direct-bit inverse changed")
    require(summaries["direct_bits_itoh"].get("query_itoh_n53_inverse") is True, "Itoh inverse was not selected")
    require(summaries["mixed_hash_binary"].get("query_prefiltered_exact_lookup") is False, "baseline lookup changed")
    require(summaries["mixed_hash_prefiltered"].get("query_prefiltered_exact_lookup") is True, "prefiltered lookup was not selected")
    require(summaries["direct_bits_binary"].get("query_prefiltered_exact_lookup") is True, "direct-bit lookup changed")
    require(summaries["direct_bits_itoh"].get("query_prefiltered_exact_lookup") is True, "Itoh lookup changed")
    require(
        all(value.get("query_specialized_n53_pair_batch") is True for value in summaries.values()),
        "specialized pair-query batch changed",
    )
    require(
        all(value.get("query_dual_sign_denominator_sharing") is True for value in summaries.values()),
        "dual-sign query kernel changed",
    )
    require(
        summaries["direct_bits_binary"]["query_exact_table_misses"]
        == summaries["direct_bits_itoh"]["query_exact_table_misses"],
        "matched direct-bit false-positive count changed",
    )
    require(
        all(value["recovered_fixture_scalar"] == stage42.EXPLICIT_SCALAR for value in summaries.values()),
        "recovered scalar changed",
    )
    rho_rows = stage44.parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0].get("verified") is True, "rho verification failed")
    require(rho_rows[0].get("published_q") == summaries["direct_bits_itoh"].get("published_q"), "rho target changed")

    mixed_wall = mixed_process["metrics"]["wall_seconds"]
    prefiltered_wall = prefiltered_process["metrics"]["wall_seconds"]
    direct_bits_wall = direct_bits_process["metrics"]["wall_seconds"]
    itoh_wall = itoh_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    comparison = {
        "prefiltered_over_baseline_wall_ratio": prefiltered_wall / mixed_wall,
        "prefiltered_collection_speedup": summaries["mixed_hash_binary"]["collection_ms"] / summaries["mixed_hash_prefiltered"]["collection_ms"],
        "direct_bits_over_prefiltered_wall_ratio": direct_bits_wall / prefiltered_wall,
        "direct_bits_collection_speedup": summaries["mixed_hash_prefiltered"]["collection_ms"] / summaries["direct_bits_binary"]["collection_ms"],
        "direct_bits_setup_speedup": summaries["mixed_hash_prefiltered"]["setup_ms"] / summaries["direct_bits_binary"]["setup_ms"],
        "itoh_over_direct_bits_wall_ratio": itoh_wall / direct_bits_wall,
        "itoh_collection_speedup": summaries["direct_bits_binary"]["collection_ms"] / summaries["direct_bits_itoh"]["collection_ms"],
        "combined_over_mixed_wall_ratio": itoh_wall / mixed_wall,
        "combined_collection_speedup": summaries["mixed_hash_binary"]["collection_ms"] / summaries["direct_bits_itoh"]["collection_ms"],
        "combined_setup_speedup": summaries["mixed_hash_binary"]["setup_ms"] / summaries["direct_bits_itoh"]["setup_ms"],
        "combined_core_seconds_ratio": itoh_process["metrics"]["total_core_seconds"] / mixed_process["metrics"]["total_core_seconds"],
        "combined_peak_rss_ratio": itoh_process["metrics"]["peak_rss_bytes"] / mixed_process["metrics"]["peak_rss_bytes"],
        "combined_over_rho_wall_ratio": itoh_wall / rho_wall,
        "mixed_over_rho_wall_ratio": mixed_wall / rho_wall,
        "fresh_build_plus_combined_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + itoh_wall) / rho_wall,
        "combined_whole_process_crossover": itoh_wall < rho_wall,
        "mixed_exact_table_misses": summaries["mixed_hash_binary"]["query_exact_table_misses"],
        "direct_bits_exact_table_misses": summaries["direct_bits_binary"]["query_exact_table_misses"],
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_filter_inverse_rho",
        "constants": self_test(),
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-94 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-94 build seal", executable=False),
            "source_state": build["source_state"],
            "outer_resources": build["outer_resources"],
        },
        "processes": processes,
        "outer_resources": outer,
        "arms": arms,
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
    result = load(output / "result.json", "Stage-94 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_filter_inverse_rho", "Stage-94 result is incomplete")
    files = {
        "mixed_hash_binary": "mixed-hash-binary.stdout",
        "mixed_hash_prefiltered": "mixed-hash-prefiltered.stdout",
        "direct_bits_binary": "direct-bits-binary.stdout",
        "direct_bits_itoh": "direct-bits-itoh.stdout",
    }
    replay = {name: stage61.direct_observation(output / "evidence" / path, "FULL_RANK") for name, path in files.items()}
    require(replay == result["arms"], "Stage-94 transcript changed")
    require(len({tuple(value["relation_hashes"]) for value in replay.values()}) == 1, "Stage-94 relation equality changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-94 claim widened")
    candidate = replay["direct_bits_itoh"]["summary"]
    return {
        "schema": "koblitz_stage94_verification.v1",
        "status": "n53_filter_inverse_verified",
        "source_commit": build["source_state"]["commit"],
        "relations": candidate["admitted_relations"],
        "support_queries": candidate["support_queries"],
        "candidate_wall_seconds": result["processes"]["direct_bits_itoh"]["metrics"]["wall_seconds"],
        "candidate_core_seconds": result["processes"]["direct_bits_itoh"]["metrics"]["total_core_seconds"],
        "candidate_peak_rss_bytes": result["processes"]["direct_bits_itoh"]["metrics"]["peak_rss_bytes"],
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage94Error, stage61.Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage94-filter-inverse: {error}")


if __name__ == "__main__":
    main()

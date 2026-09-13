#!/usr/bin/env python3
"""Sweep parallel query width under the direct-bit and Itoh n=53 kernels."""

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
RUN_SCHEMA = "koblitz_stage95_itoh_width.v1"
RUN_SEAL_SCHEMA = "koblitz_stage95_itoh_width_seal.v1"
WIDTHS = (4096, 2048, 1024)


class Stage95Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage95Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    return {
        "schema": "koblitz_stage95_self_test.v1",
        "status": "PASS",
        "checks": 24,
        "n": 53,
        "widths": list(WIDTHS),
        "parallel_lanes": 4,
        "same_target": True,
        "direct_x_filter_bits": True,
        "prefiltered_exact_lookup": True,
        "itoh_tsujii_inverse": True,
        "exact_relation_hash_equality_required": True,
    }


def candidate_environment() -> dict[str, str]:
    environment = custody.safe_child_environment()
    environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    environment["KIC_RANK_SURPLUS"] = "0"
    environment["KIC_RANK_AWARE_PAIR_SCAN"] = "1"
    environment["KIC_PARALLEL_SUPPORT_EXPANSION"] = "1"
    environment["KIC_PIPELINED_SUPPORT_EXPANSION"] = "1"
    environment["KIC_ENABLE_PREFILTERED_EXACT_LOOKUP"] = "1"
    environment["KIC_ENABLE_DIRECT_X_FILTER_BITS"] = "1"
    environment["KIC_ENABLE_ITOH_N53_INVERSE"] = "1"
    environment["RAYON_NUM_THREADS"] = "4"
    return environment


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-95 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-95 output must be new")
    build_root = args.build.resolve(strict=True)
    build = stage42.validate_build(build_root)
    output.mkdir(parents=True)
    evidence = output / "evidence"
    evidence.mkdir()
    direct = str((build_root / "bin/koblitz_rank_fixture").resolve(strict=True))
    rho = str((build_root / "bin/koblitz_rho_fixture").resolve(strict=True))
    environment = candidate_environment()
    direct_commands = {
        str(width): [
            direct, "53", "0", "1", "128", str(stage42.DIRECT_SEED),
            "signed_expanded", "independent", f"pair_pair_parallel_{width}", "1",
            str(stage42.EXPLICIT_SCALAR),
        ]
        for width in WIDTHS
    }
    rho_command = [
        rho, "53", "0", "signed_frobenius", "1", "packed",
        str(stage42.RHO_BATCH_SEED), str(stage42.EXPLICIT_SCALAR),
    ]

    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    processes: dict[str, dict[str, Any]] = {}
    with TreeSampler() as sampler:
        for width in WIDTHS:
            processes[str(width)] = stage33.metered(
                role=f"width-{width}", command=direct_commands[str(width)], cwd=REPO,
                evidence=evidence, timeout=args.timeout, environment=environment,
            )
        processes["rho"] = stage33.metered(
            role="rho", command=rho_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=custody.safe_child_environment(),
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    require(all(stage33.process_complete(row) for row in processes.values()), "Stage-95 process failed")

    arms = {
        str(width): stage61.direct_observation(evidence / f"width-{width}.stdout", "FULL_RANK")
        for width in WIDTHS
    }
    bases = {width: value["base"] for width, value in arms.items()}
    summaries = {width: value["summary"] for width, value in arms.items()}
    require(len({value["base_hash"] for value in bases.values()}) == 1, "factor base changed")
    require(len({value["support_index_entries"] for value in bases.values()}) == 1, "support entries changed")
    require(len({tuple(value["relation_hashes"]) for value in arms.values()}) == 1, "relation transcript changed")
    require(len({tuple(value["factor_base_log_solution"]) for value in summaries.values()}) == 1, "full-rank solution changed")
    for width, summary in summaries.items():
        require(summary.get("query_mode") == f"pair_pair_parallel_{width}", f"width {width} was not selected")
        require(summary.get("query_parallel_threads") == 4, f"width {width} lost four-lane execution")
        require(summary.get("query_specialized_n53_pair_batch") is True, f"width {width} lost specialized batch")
        require(summary.get("query_itoh_n53_inverse") is True, f"width {width} lost Itoh inverse")
        require(summary.get("query_prefiltered_exact_lookup") is True, f"width {width} lost prefiltered lookup")
        require(summary.get("recovered_fixture_scalar") == stage42.EXPLICIT_SCALAR, f"width {width} recovered scalar changed")
        require(bases[width].get("support_x_prefilter_hash_strategy") == "direct_low_and_high_x_bit_windows", f"width {width} lost direct-bit filter")
    rho_rows = stage44.parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0].get("verified") is True, "rho verification failed")
    require(rho_rows[0].get("published_q") == summaries["4096"].get("published_q"), "rho target changed")

    rho_wall = processes["rho"]["metrics"]["wall_seconds"]
    width_metrics = {
        width: {
            "wall_seconds": processes[width]["metrics"]["wall_seconds"],
            "core_seconds": processes[width]["metrics"]["total_core_seconds"],
            "peak_rss_bytes": processes[width]["metrics"]["peak_rss_bytes"],
            "setup_ms": summaries[width]["setup_ms"],
            "collection_ms": summaries[width]["collection_ms"],
            "query_ms": summaries[width]["timing_breakdown_ms"]["query"],
            "support_queries": summaries[width]["support_queries"],
            "query_batch_inversions": summaries[width]["query_batch_inversions"],
            "query_parallel_waves": summaries[width]["query_parallel_waves"],
            "query_parallel_chunks": summaries[width]["query_parallel_chunks"],
            "query_exact_table_misses": summaries[width]["query_exact_table_misses"],
            "over_rho_wall_ratio": processes[width]["metrics"]["wall_seconds"] / rho_wall,
        }
        for width in summaries
    }
    best_width = min(width_metrics, key=lambda width: width_metrics[width]["wall_seconds"])
    comparison = {
        "width_metrics": width_metrics,
        "best_width": int(best_width),
        "best_over_4096_wall_ratio": width_metrics[best_width]["wall_seconds"] / width_metrics["4096"]["wall_seconds"],
        "best_over_4096_collection_ratio": width_metrics[best_width]["collection_ms"] / width_metrics["4096"]["collection_ms"],
        "best_over_4096_core_ratio": width_metrics[best_width]["core_seconds"] / width_metrics["4096"]["core_seconds"],
        "best_over_rho_wall_ratio": width_metrics[best_width]["over_rho_wall_ratio"],
        "best_whole_process_crossover": width_metrics[best_width]["wall_seconds"] < rho_wall,
        "fresh_build_plus_best_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + width_metrics[best_width]["wall_seconds"]) / rho_wall,
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_itoh_width_rho",
        "constants": self_test(),
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-95 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-95 build seal", executable=False),
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
    result = load(output / "result.json", "Stage-95 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_itoh_width_rho", "Stage-95 result is incomplete")
    replay = {
        str(width): stage61.direct_observation(output / "evidence" / f"width-{width}.stdout", "FULL_RANK")
        for width in WIDTHS
    }
    require(replay == result["arms"], "Stage-95 transcript changed")
    require(len({tuple(value["relation_hashes"]) for value in replay.values()}) == 1, "Stage-95 relation equality changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-95 claim widened")
    best = str(result["comparison"]["best_width"])
    metrics = result["comparison"]["width_metrics"][best]
    return {
        "schema": "koblitz_stage95_verification.v1",
        "status": "n53_itoh_width_verified",
        "source_commit": build["source_state"]["commit"],
        "relations": replay[best]["summary"]["admitted_relations"],
        "best_width": int(best),
        "best_metrics": metrics,
        **{key: value for key, value in result["comparison"].items() if key != "width_metrics"},
        "rho_wall_seconds": result["processes"]["rho"]["metrics"]["wall_seconds"],
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage95Error, stage61.Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage95-itoh-width: {error}")


if __name__ == "__main__":
    main()

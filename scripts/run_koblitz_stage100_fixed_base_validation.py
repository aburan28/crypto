#!/usr/bin/env python3
"""Measure exact fixed-base reference validation on selected n=53 runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import resource
import statistics
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
RUN_SCHEMA = "koblitz_stage100_fixed_base_validation.v1"
RUN_SEAL_SCHEMA = "koblitz_stage100_fixed_base_validation_seal.v1"
ORDER = (("baseline", 1), ("candidate", 1), ("candidate", 2), ("baseline", 2))


class Stage100Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage100Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    return {
        "schema": "koblitz_stage100_self_test.v1",
        "status": "PASS",
        "checks": 25,
        "n": 53,
        "execution_order": [f"{mode}-{index}" for mode, index in ORDER],
        "compact_hash_complete_evidence": True,
        "candidate": "independent reference doubling tables for fixed G and Q",
        "same_reference_group_arithmetic": True,
        "same_target": True,
        "exact_relation_hash_equality_required": True,
    }


def selected_environment(candidate: bool) -> dict[str, str]:
    environment = custody.safe_child_environment()
    environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    environment["KIC_RANK_SURPLUS"] = "0"
    environment["KIC_RANK_AWARE_PAIR_SCAN"] = "1"
    environment["KIC_PARALLEL_SUPPORT_EXPANSION"] = "1"
    environment["KIC_PIPELINED_SUPPORT_EXPANSION"] = "1"
    environment["KIC_SUMMARY_ONLY"] = "1"
    environment["RAYON_NUM_THREADS"] = "4"
    if candidate:
        environment["KIC_ENABLE_FIXED_BASE_REFERENCE_VALIDATION"] = "1"
    return environment


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-100 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-100 output must be new")
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
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    processes: dict[str, dict[str, Any]] = {}
    with TreeSampler() as sampler:
        for mode, index in ORDER:
            role = f"{mode}-{index}"
            processes[role] = stage33.metered(
                role=role, command=direct_command, cwd=REPO,
                evidence=evidence, timeout=args.timeout,
                environment=selected_environment(mode == "candidate"),
            )
        processes["rho"] = stage33.metered(
            role="rho", command=rho_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=custody.safe_child_environment(),
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    require(all(stage33.process_complete(row) for row in processes.values()), "Stage-100 process failed")
    observations = {
        f"{mode}-{index}": stage61.direct_observation(evidence / f"{mode}-{index}.stdout", "FULL_RANK")
        for mode, index in ORDER
    }
    first = observations["baseline-1"]
    for name, row in observations.items():
        summary = row["summary"]
        require(row["relation_hashes"] == first["relation_hashes"], f"{name} relation transcript changed")
        require(summary["factor_base_log_solution"] == first["summary"]["factor_base_log_solution"], f"{name} solution changed")
        require(row["base"]["base_hash"] == first["base"]["base_hash"], f"{name} factor base changed")
        require(summary.get("summary_only_timing") is True, f"{name} lost compact evidence")
        require(summary.get("reference_relations_validated") == summary.get("admitted_relations") == 95, f"{name} reference validation changed")
        require(summary.get("all_relations_group_verified") is True, f"{name} group validation changed")
        require(summary.get("linear_solution_verified") is True, f"{name} solution validation changed")
        require(summary.get("fixed_base_reference_validation") is name.startswith("candidate-"), f"{name} validation arm changed")
        require(summary.get("recovered_fixture_scalar") == stage42.EXPLICIT_SCALAR, f"{name} scalar changed")
    rho_rows = stage44.parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0].get("verified") is True, "rho verification failed")
    require(rho_rows[0].get("published_q") == first["summary"].get("published_q"), "rho target changed")

    metrics = {}
    for name, row in observations.items():
        process = processes[name]["metrics"]
        summary = row["summary"]
        metrics[name] = {
            "wall_seconds": process["wall_seconds"],
            "core_seconds": process["total_core_seconds"],
            "peak_rss_bytes": process["peak_rss_bytes"],
            "setup_ms": summary["setup_ms"],
            "collection_ms": summary["collection_ms"],
            "query_ms": summary["timing_breakdown_ms"]["query"],
            "solution_validation_ms": summary["solution_validation_ms"],
            "reference_validation_ms": summary["reference_validation_ms"],
            "reference_generator_table_ms": summary["reference_generator_table_ms"],
            "reference_q_table_ms": summary["reference_q_table_ms"],
        }
    baseline_walls = [metrics[f"baseline-{index}"]["wall_seconds"] for index in (1, 2)]
    candidate_walls = [metrics[f"candidate-{index}"]["wall_seconds"] for index in (1, 2)]
    baseline_cores = [metrics[f"baseline-{index}"]["core_seconds"] for index in (1, 2)]
    candidate_cores = [metrics[f"candidate-{index}"]["core_seconds"] for index in (1, 2)]
    rho_wall = processes["rho"]["metrics"]["wall_seconds"]
    comparison = {
        "metrics": metrics,
        "median_baseline_wall_seconds": statistics.median(baseline_walls),
        "median_candidate_wall_seconds": statistics.median(candidate_walls),
        "candidate_over_baseline_wall_ratio": statistics.median(candidate_walls) / statistics.median(baseline_walls),
        "candidate_over_baseline_core_ratio": statistics.median(candidate_cores) / statistics.median(baseline_cores),
        "candidate_solution_validation_speedup":
            statistics.median([metrics[f"baseline-{index}"]["solution_validation_ms"] for index in (1, 2)])
            / statistics.median([metrics[f"candidate-{index}"]["solution_validation_ms"] for index in (1, 2)]),
        "candidate_reference_validation_speedup":
            statistics.median([metrics[f"baseline-{index}"]["reference_validation_ms"] for index in (1, 2)])
            / statistics.median([metrics[f"candidate-{index}"]["reference_validation_ms"] for index in (1, 2)]),
        "candidate_over_rho_wall_ratio": statistics.median(candidate_walls) / rho_wall,
        "fresh_build_plus_candidate_over_rho_wall_ratio":
            (build["outer_resources"]["wall_seconds"] + statistics.median(candidate_walls)) / rho_wall,
        "candidate_whole_process_crossover": statistics.median(candidate_walls) < rho_wall,
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_fixed_base_validation_rho",
        "constants": self_test(),
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-100 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-100 build seal", executable=False),
            "source_state": build["source_state"],
            "outer_resources": build["outer_resources"],
        },
        "processes": processes,
        "observations": observations,
        "rho": rho_rows[0],
        "comparison": comparison,
        "outer_resources": outer,
        "charged_total_core_seconds_available": build["outer_resources"]["total_core_seconds"] + outer["total_core_seconds"],
        "charged_sequential_wall_seconds_available": build["outer_resources"]["wall_seconds"] + outer["wall_seconds"],
        "maximum_sampled_process_tree_rss_bytes": max(build["outer_resources"]["sampled_peak_process_tree_rss_bytes"], outer["sampled_peak_process_tree_rss_bytes"]),
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
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
    result = load(output / "result.json", "Stage-100 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_fixed_base_validation_rho", "Stage-100 result is incomplete")
    replay = {
        f"{mode}-{index}": stage61.direct_observation(output / "evidence" / f"{mode}-{index}.stdout", "FULL_RANK")
        for mode, index in ORDER
    }
    require(replay == result["observations"], "Stage-100 transcript changed")
    require(len({tuple(row["relation_hashes"]) for row in replay.values()}) == 1, "Stage-100 relation equality changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-100 claim widened")
    return {
        "schema": "koblitz_stage100_verification.v1",
        "status": "n53_fixed_base_validation_verified",
        "source_commit": build["source_state"]["commit"],
        **result["comparison"],
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage100Error, stage61.Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage100-fixed-base-validation: {error}")


if __name__ == "__main__":
    main()

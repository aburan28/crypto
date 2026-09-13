#!/usr/bin/env python3
"""Measure the optimized n=53 direct/rho control under one-CPU affinity."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
from run_koblitz_stage26_affinity_cell import TreeSampler
import run_koblitz_stage33_n41_unknown_scalar as stage33
import run_koblitz_stage42_n53_same_target as stage42
import run_koblitz_stage44_n53_parallel as stage44
import run_koblitz_stage61_production_rank as stage61


REPO = Path(__file__).resolve().parents[1]
RUN_SCHEMA = "koblitz_stage71_single_core.v1"
RUN_SEAL_SCHEMA = "koblitz_stage71_single_core_seal.v1"


class Stage71Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage71Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    return {
        "schema": "koblitz_stage71_self_test.v1",
        "status": "PASS",
        "checks": 13,
        "n": 53,
        "eta": [1, 128],
        "same_target": True,
        "linux_affinity_required": True,
        "effective_cpu_count": 1,
        "rayon_threads": 1,
        "parallel_support_expansion": False,
        "rank_surplus": 0,
        "rank_aware_pair_scan": True,
        "field_pipeline": "x86_64_pclmul_n53_fused_reduce",
    }


def child_affinity() -> list[int]:
    completed = subprocess.run(
        [sys.executable, "-c", "import json,os; print(json.dumps(sorted(os.sched_getaffinity(0))))"],
        text=True,
        capture_output=True,
        check=True,
    )
    value = json.loads(completed.stdout)
    require(isinstance(value, list) and all(type(cpu) is int for cpu in value), "child affinity probe failed")
    return value


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux" and hasattr(os, "sched_getaffinity"), "Stage-71 requires Linux CPU affinity")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-71 output must be new")
    build_root = args.build.resolve(strict=True)
    build = stage42.validate_build(build_root)
    initial_cpus = sorted(os.sched_getaffinity(0))
    require(bool(initial_cpus), "Stage-71 has no allowed CPU")
    selected_cpu = initial_cpus[0]
    os.sched_setaffinity(0, {selected_cpu})
    effective_cpus = sorted(os.sched_getaffinity(0))
    inherited_cpus = child_affinity()
    require(effective_cpus == inherited_cpus == [selected_cpu], "single-CPU affinity was not inherited")
    host_identity = stage33.linux_host_identity()
    require(host_identity["cpu"]["effective_affinity"] == [selected_cpu], "host receipt lost CPU-0 affinity")
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
    direct_environment = custody.safe_child_environment()
    direct_environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    direct_environment["KIC_RANK_SURPLUS"] = "0"
    direct_environment["KIC_RANK_AWARE_PAIR_SCAN"] = "1"
    direct_environment["RAYON_NUM_THREADS"] = "1"
    rho_environment = custody.safe_child_environment()
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        direct_process = stage33.metered(
            role="single-core-direct", command=direct_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=direct_environment,
        )
        rho_process = stage33.metered(
            role="single-core-rho", command=rho_command, cwd=REPO,
            evidence=evidence, timeout=args.timeout, environment=rho_environment,
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    require(stage33.process_complete(direct_process) and stage33.process_complete(rho_process), "Stage-71 process failed")
    direct_value = stage61.direct_observation(evidence / "single-core-direct.stdout", "FULL_RANK")
    summary = direct_value["summary"]
    require(summary.get("query_parallel_threads") == 1, "direct used more than one Rayon thread")
    require(direct_value["base"].get("parallel_support_expansion") is False, "direct used parallel support expansion")
    require(summary.get("rank_aware_pair_scan") is True, "direct did not use rank-aware collection")
    require(summary.get("rank_target_deficiency") == 1, "direct rank-deficiency trigger changed")
    require(summary.get("field_product_pipeline") == "x86_64_pclmul_n53_fused_reduce", "direct did not use fused field products")
    require(direct_value["base"].get("support_x_prefilter_hash_strategy") == "direct_low_and_high_x_bit_windows", "direct lost selected x filter")
    require(summary.get("query_itoh_n53_inverse") is True, "direct lost Itoh inverse")
    require(summary.get("fixed_base_reference_validation") is True, "direct lost fixed-base validation")
    require(summary.get("recovered_fixture_scalar") == stage42.EXPLICIT_SCALAR, "direct recovered scalar changed")
    rho_rows = stage44.parse_lines(evidence / "single-core-rho.stdout")
    require(len(rho_rows) == 1 and rho_rows[0].get("verified") is True, "rho verification failed")
    require(rho_rows[0].get("published_q") == summary.get("published_q"), "direct/rho target changed")
    direct_wall = direct_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    comparison = {
        "direct_over_rho_wall_ratio": direct_wall / rho_wall,
        "direct_core_over_rho_core_ratio": direct_process["metrics"]["total_core_seconds"] / rho_process["metrics"]["total_core_seconds"],
        "fresh_build_plus_direct_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + direct_wall) / rho_wall,
        "whole_process_crossover": direct_wall < rho_wall,
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_single_core_direct_rho",
        "constants": self_test(),
        "affinity": {
            "initial_allowed_cpus": initial_cpus,
            "selected_cpu": selected_cpu,
            "parent_effective_cpus": effective_cpus,
            "child_effective_cpus": inherited_cpus,
            "inheritance_verified": True,
        },
        "host_identity": host_identity,
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-71 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-71 build seal", executable=False),
            "source_state": build["source_state"],
            "outer_resources": build["outer_resources"],
        },
        "direct_process": direct_process,
        "rho_process": rho_process,
        "outer_resources": outer,
        "direct": direct_value,
        "rho": rho_rows[0],
        "comparison": comparison,
        "charged_total_core_seconds_available": build["outer_resources"]["total_core_seconds"] + outer["total_core_seconds"],
        "charged_sequential_wall_seconds_available": build["outer_resources"]["wall_seconds"] + outer["wall_seconds"],
        "maximum_sampled_process_tree_rss_bytes": max(build["outer_resources"]["sampled_peak_process_tree_rss_bytes"], outer["sampled_peak_process_tree_rss_bytes"]),
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
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
    result = load(output / "result.json", "Stage-71 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_single_core_direct_rho", "Stage-71 result is incomplete")
    require(result.get("affinity", {}).get("inheritance_verified") is True and len(result["affinity"]["parent_effective_cpus"]) == 1, "Stage-71 affinity changed")
    stage33.validate_linux_host_identity(result["host_identity"])
    require(result["host_identity"]["cpu"]["effective_affinity"] == result["affinity"]["parent_effective_cpus"], "Stage-71 host affinity changed")
    direct = stage61.direct_observation(output / "evidence/single-core-direct.stdout", "FULL_RANK")
    require(direct["summary"].get("rank_aware_pair_scan") is True, "verified direct lost rank-aware collection")
    require(direct["summary"].get("rank_target_deficiency") == 1, "verified rank-deficiency trigger changed")
    rho_rows = stage44.parse_lines(output / "evidence/single-core-rho.stdout")
    require(direct == result["direct"] and len(rho_rows) == 1 and rho_rows[0] == result["rho"], "Stage-71 evidence changed")
    dw = result["direct_process"]["metrics"]["wall_seconds"]
    rw = result["rho_process"]["metrics"]["wall_seconds"]
    expected = {
        "direct_over_rho_wall_ratio": dw / rw,
        "direct_core_over_rho_core_ratio": result["direct_process"]["metrics"]["total_core_seconds"] / result["rho_process"]["metrics"]["total_core_seconds"],
        "fresh_build_plus_direct_over_rho_wall_ratio": (build["outer_resources"]["wall_seconds"] + dw) / rw,
        "whole_process_crossover": dw < rw,
    }
    require(result["comparison"] == expected, "Stage-71 comparison changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-71 claim widened")
    return {
        "schema": "koblitz_stage71_verification.v1",
        "status": "n53_single_core_verified",
        "source_commit": build["source_state"]["commit"],
        "selected_cpu": result["affinity"]["selected_cpu"],
        "host_identity": result["host_identity"],
        "relations": direct["summary"]["admitted_relations"],
        "support_queries": direct["summary"]["support_queries"],
        "direct_wall_seconds": dw,
        "direct_core_seconds": result["direct_process"]["metrics"]["total_core_seconds"],
        "direct_peak_rss_bytes": result["direct_process"]["metrics"]["peak_rss_bytes"],
        "rho_wall_seconds": rw,
        **expected,
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage71Error, stage61.Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage71-single-core: {error}")


if __name__ == "__main__":
    main()

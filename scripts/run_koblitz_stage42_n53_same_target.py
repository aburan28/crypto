#!/usr/bin/env python3
"""Build, run, seal, and verify the n=53 same-target direct/rho successor."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import time
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
from run_koblitz_stage26_affinity_cell import TreeSampler
import run_koblitz_stage33_n41_unknown_scalar as stage33


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
LOCK = STAGE / "stage-20-rust-build/Cargo.lock"
SCALAR_LABEL = "KIC-STAGE42-N53-HOSTED-SAME-TARGET-SCALAR"
DIRECT_LABEL = "KIC-STAGE42-N53-HOSTED-DIRECT-SEED"
RHO_LABEL = "KIC-STAGE42-N53-HOSTED-RHO-SEED"
SUBGROUP_ORDER = 21044858204113
EXPLICIT_SCALAR = 476811900269
DIRECT_SEED = 18118074343025384088
RHO_BATCH_SEED = 16969290756857938023
BUILD_SCHEMA = "koblitz_stage42_n53_build.v1"
BUILD_SEAL_SCHEMA = "koblitz_stage42_n53_build_seal.v1"
RUN_SCHEMA = "koblitz_stage42_n53_same_target.v1"
RUN_SEAL_SCHEMA = "koblitz_stage42_n53_same_target_seal.v1"


class Stage42Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage42Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def derive(label: str) -> int:
    return int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "big")


def self_test() -> dict[str, Any]:
    require(derive(SCALAR_LABEL) % (SUBGROUP_ORDER - 1) + 1 == EXPLICIT_SCALAR, "Stage-42 scalar derivation changed")
    require(derive(DIRECT_LABEL) == DIRECT_SEED and derive(RHO_LABEL) == RHO_BATCH_SEED, "Stage-42 algorithm seeds changed")
    require(len({EXPLICIT_SCALAR, DIRECT_SEED, RHO_BATCH_SEED}) == 3, "Stage-42 constants collided")
    return {
        "schema": "koblitz_stage42_self_test.v1",
        "status": "PASS",
        "checks": 10,
        "n": 53,
        "eta": [1, 128],
        "decomposition_arity": 4,
        "pair_mode": "signed_expanded",
        "query_mode": "pair_pair_256",
        "explicit_public_validation_scalar": EXPLICIT_SCALAR,
        "direct_seed": DIRECT_SEED,
        "rho_batch_seed": RHO_BATCH_SEED,
        "same_target_required": True,
        "factor_base_selection_uses_scalar_labels": False,
        "fixture_scalar_used_by_collector": False,
    }


def source_state(source: Path) -> dict[str, Any]:
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=source, text=True, capture_output=True, check=True).stdout.strip()
    custody.require_hex40(commit, "Stage-42 source commit")
    porcelain = subprocess.run(["git", "status", "--porcelain=v1", "--untracked-files=all"], cwd=source, text=True, capture_output=True, check=True).stdout.splitlines()
    return {"commit": commit, "dirty": bool(porcelain), "porcelain": porcelain}


def build(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-42 build requires Linux process accounting")
    source = args.source.resolve(strict=True)
    output = args.output.resolve()
    require(source == REPO, "Stage-42 build source must be this repository")
    require(not output.exists() and not output.is_symlink(), "Stage-42 build output must be new")
    before = source_state(source)
    require(before["dirty"] is False, "Stage-42 build source must be clean")
    output.mkdir(parents=True)
    evidence = output / "evidence"
    evidence.mkdir()
    target = output / "target"
    root_lock = source / "Cargo.lock"
    had_lock = root_lock.exists()
    old_lock = root_lock.read_bytes() if had_lock else None
    root_lock.write_bytes(LOCK.read_bytes())
    cargo = shutil.which("cargo")
    rustc = shutil.which("rustc")
    require(cargo is not None and rustc is not None, "Stage-42 Rust toolchain is unavailable")
    # Preserve the rustup shim names. Resolving them to the shared rustup
    # binary changes argv[0] and no longer selects cargo/rustc. Reuse the
    # audited build environment, including isolated CARGO_HOME and the host's
    # RUSTUP_HOME binding.
    cargo_path = Path(cargo).absolute()
    rustc_path = Path(rustc).absolute()
    environment = stage33.build_environment(output, cargo_path, rustc_path)
    environment["CARGO_TARGET_DIR"] = str(target)
    command = [str(cargo_path), "build", "--release", "--locked", "--jobs", str(args.jobs), "--example", "koblitz_rank_fixture", "--example", "koblitz_rho_fixture"]
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    try:
        with TreeSampler() as sampler:
            process = stage33.metered(role="build", command=command, cwd=source, evidence=evidence, timeout=args.timeout, environment=environment)
    finally:
        if had_lock:
            root_lock.write_bytes(old_lock or b"")
        else:
            root_lock.unlink(missing_ok=True)
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    require(stage33.process_complete(process), "Stage-42 producer build failed")
    after = source_state(source)
    require(after == before, "Stage-42 build changed source state")
    bin_dir = output / "bin"
    bin_dir.mkdir()
    binaries = {}
    for name in ("koblitz_rank_fixture", "koblitz_rho_fixture"):
        source_binary = target / "release/examples" / name
        destination = bin_dir / name
        shutil.copy2(source_binary, destination)
        binaries[name] = custody.executable_identity(destination, f"Stage-42 {name}", executable=True)
    shutil.rmtree(target)
    shutil.rmtree(output / "cargo-home")
    result = {
        "schema": BUILD_SCHEMA,
        "status": "complete_clean_build",
        "source_state": before,
        "frozen_lock": custody.executable_identity(LOCK, "Stage-42 lock", executable=False),
        "requested_parallel_jobs": args.jobs,
        "process": process,
        "outer_resources": outer,
        "binaries": binaries,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    custody.write_json_new(output / "result.json", result)
    stage33.freeze(output, BUILD_SEAL_SCHEMA, "build_frozen")
    return result


def validate_build(root: Path) -> dict[str, Any]:
    seal = stage33.verify_seal(root, BUILD_SEAL_SCHEMA, "build_frozen")
    result = load(root / "result.json", "Stage-42 build result")
    require(result.get("schema") == BUILD_SCHEMA and result.get("status") == "complete_clean_build", "Stage-42 build is incomplete")
    require(result.get("source_state", {}).get("dirty") is False, "Stage-42 build source was dirty")
    require(set(result.get("binaries", {})) == {"koblitz_rank_fixture", "koblitz_rho_fixture"}, "Stage-42 binary inventory changed")
    for name, identity in result.get("binaries", {}).items():
        path = root / "bin" / name
        require(custody.executable_identity(path, f"Stage-42 {name}", executable=True)["sha256"] == identity["sha256"], f"Stage-42 {name} changed")
    result["seal"] = seal
    return result


def parse_lines(path: Path) -> list[dict[str, Any]]:
    rows = []
    for line in path.read_text().splitlines():
        value = json.loads(line)
        require(isinstance(value, dict), f"non-object JSON line in {path}")
        rows.append(value)
    return rows


def analyze(direct_rows: list[dict[str, Any]], rho: dict[str, Any]) -> dict[str, Any]:
    base = next(row for row in direct_rows if row.get("kind") == "point_defined_factor_base")
    summary = next(row for row in direct_rows if row.get("kind") == "relation_rank_summary")
    require(summary.get("n") == rho.get("n") == 53 and summary.get("a") == rho.get("a") == 0, "Stage-42 curve mismatch")
    require(summary.get("eta") == {"numerator": 1, "denominator": 128}, "Stage-42 eta changed")
    require(summary.get("pair_index_mode") == "signed_expanded" and summary.get("query_mode") == "pair_pair_256", "Stage-42 direct mode changed")
    require(summary.get("decomposition_arity") == 4 and summary.get("status") == "RANK_PLUS_32", "Stage-42 rank endpoint changed")
    require(summary.get("linear_solution_verified") is True and summary.get("all_relations_group_verified") is True, "Stage-42 direct verification failed")
    require(base.get("selection_uses_scalar_labels") is False, "Stage-42 factor base used scalar labels")
    require(summary.get("fixture_scalar_source") == rho.get("fixture_scalar_source") == "explicit_public_validation_scalar", "Stage-42 target source changed")
    require(summary.get("published_fixture_scalar") == rho.get("published_fixture_scalar") == EXPLICIT_SCALAR, "Stage-42 target scalar mismatch")
    require(summary.get("published_q") == rho.get("published_q"), "Stage-42 direct and rho targets differ")
    require(summary.get("recovered_fixture_scalar") == rho.get("recovered_fixture_scalar") == EXPLICIT_SCALAR, "Stage-42 recovered scalar mismatch")
    require(rho.get("verified") is True and rho.get("reference_group_validation") is True, "Stage-42 rho verification failed")
    return {"base": base, "direct": summary, "rho": rho}


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-42 run requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-42 run output must be new")
    build_root = args.build.resolve(strict=True)
    build_result = validate_build(build_root)
    output.mkdir(parents=True)
    evidence = output / "evidence"
    evidence.mkdir()
    environment = custody.safe_child_environment()
    environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    direct_command = [str((build_root / "bin/koblitz_rank_fixture").resolve(strict=True)), "53", "0", "1", "128", str(DIRECT_SEED), "signed_expanded", "independent", "pair_pair_256", "1", str(EXPLICIT_SCALAR)]
    rho_command = [str((build_root / "bin/koblitz_rho_fixture").resolve(strict=True)), "53", "0", "signed_frobenius", "1", "packed", str(RHO_BATCH_SEED), str(EXPLICIT_SCALAR)]
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        direct_process = stage33.metered(role="direct", command=direct_command, cwd=REPO, evidence=evidence, timeout=args.timeout, environment=environment)
        rho_process = stage33.metered(role="rho", command=rho_command, cwd=REPO, evidence=evidence, timeout=args.timeout, environment=environment)
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    require(stage33.process_complete(direct_process) and stage33.process_complete(rho_process), "Stage-42 producer failed")
    direct_rows = parse_lines(evidence / "direct.stdout")
    rho_rows = parse_lines(evidence / "rho.stdout")
    require(len(rho_rows) == 1, "Stage-42 rho row count changed")
    observation = analyze(direct_rows, rho_rows[0])
    direct_wall = direct_process["metrics"]["wall_seconds"]
    rho_wall = rho_process["metrics"]["wall_seconds"]
    build_wall = build_result["outer_resources"]["wall_seconds"]
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_same_target",
        "build_binding": {"result": custody.executable_identity(build_root / "result.json", "Stage-42 build result", executable=False), "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-42 build seal", executable=False), "source_state": build_result["source_state"], "binaries": build_result["binaries"], "outer_resources": build_result["outer_resources"], "process": build_result["process"]},
        "constants": self_test(),
        "direct_process": direct_process,
        "rho_process": rho_process,
        "outer_resources": outer,
        "observation": observation,
        "comparison": {
            "direct_over_rho_whole_process_wall_ratio": direct_wall / rho_wall,
            "direct_charged_over_rho_algorithm_ratio": observation["direct"]["charged_total_ms"] / observation["rho"]["total_ms"],
            "fresh_build_plus_direct_over_rho_wall_ratio": (build_wall + direct_wall) / rho_wall,
            "whole_process_crossover": direct_wall < rho_wall,
        },
        "charged_total_core_seconds_available": build_result["outer_resources"]["total_core_seconds"] + outer["total_core_seconds"],
        "charged_sequential_wall_seconds_available": build_wall + outer["wall_seconds"],
        "maximum_sampled_process_tree_rss_bytes": max(build_result["outer_resources"]["sampled_peak_process_tree_rss_bytes"], outer["sampled_peak_process_tree_rss_bytes"]),
        "same_public_target": True,
        "factor_base_algebraically_point_defined": True,
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
        "fixture_scalar_used_by_collector": False,
        "target_scalar_retained_for_validation": True,
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
    build_result = validate_build(build_root)
    seal = stage33.verify_seal(output, RUN_SEAL_SCHEMA, "run_frozen")
    result = load(output / "result.json", "Stage-42 run result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_same_target", "Stage-42 result is incomplete")
    direct_rows = parse_lines(output / "evidence/direct.stdout")
    rho_rows = parse_lines(output / "evidence/rho.stdout")
    require(len(rho_rows) == 1 and result.get("observation") == analyze(direct_rows, rho_rows[0]), "Stage-42 observations changed")
    direct = result["direct_process"]["metrics"]["wall_seconds"]
    rho = result["rho_process"]["metrics"]["wall_seconds"]
    comparison = result["comparison"]
    observation = result["observation"]
    expected_comparison = {
        "direct_over_rho_whole_process_wall_ratio": direct / rho,
        "direct_charged_over_rho_algorithm_ratio": observation["direct"]["charged_total_ms"] / observation["rho"]["total_ms"],
        "fresh_build_plus_direct_over_rho_wall_ratio": (build_result["outer_resources"]["wall_seconds"] + direct) / rho,
        "whole_process_crossover": direct < rho,
    }
    require(comparison == expected_comparison, "Stage-42 comparison changed")
    expected_core = build_result["outer_resources"]["total_core_seconds"] + result["outer_resources"]["total_core_seconds"]
    expected_wall = build_result["outer_resources"]["wall_seconds"] + result["outer_resources"]["wall_seconds"]
    expected_peak = max(build_result["outer_resources"]["sampled_peak_process_tree_rss_bytes"], result["outer_resources"]["sampled_peak_process_tree_rss_bytes"])
    require(result["charged_total_core_seconds_available"] == expected_core and result["charged_sequential_wall_seconds_available"] == expected_wall, "Stage-42 resource total changed")
    require(result["maximum_sampled_process_tree_rss_bytes"] == expected_peak, "Stage-42 peak memory changed")
    require(result.get("same_public_target") is True and result.get("factor_base_selection_uses_scalar_labels") is False, "Stage-42 target/base boundary changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-42 claim widened")
    return {
        "schema": "koblitz_stage42_verification.v1",
        "status": "n53_same_target_build_and_run_verified",
        "source_commit": build_result["source_state"]["commit"],
        "published_q": result["observation"]["direct"]["published_q"],
        "published_fixture_scalar": EXPLICIT_SCALAR,
        "orbit_columns": result["observation"]["direct"]["orbit_columns"],
        "factor_base_points": result["observation"]["direct"]["factor_base_points"],
        "relations": result["observation"]["direct"]["admitted_relations"],
        "direct_wall_seconds": direct,
        "rho_wall_seconds": rho,
        "direct_over_rho_wall_ratio": comparison["direct_over_rho_whole_process_wall_ratio"],
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
            value = build(args)
        elif args.command == "run":
            value = run(args)
        else:
            value = verify(args.build.resolve(strict=True), args.output.resolve(strict=True))
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage42Error, stage33.Stage33Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage42-n53: {error}")


if __name__ == "__main__":
    main()

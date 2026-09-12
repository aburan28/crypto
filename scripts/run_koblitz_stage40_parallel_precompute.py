#!/usr/bin/env python3
"""Compare one- and four-core fixed-algebraic n=41 precomputation."""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
from run_koblitz_stage26_affinity_cell import Stage26Error, TreeSampler
import run_koblitz_stage33_n41_unknown_scalar as stage33
import run_koblitz_stage39_fixed_algebraic_holdout as stage39


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
BASE_PARAMS = STAGE / "stage-39-n41-fixed-algebraic-holdout-params.json"
PREDECESSOR = STAGE / "stage-39-fixed-algebraic-result-20260912/verification.json"
TARGET_SEEDS = [41401, 41402, 41403, 41404, 41405]
RELATION_SEED = 41431
CELLS = [("single-a", 1), ("four-a", 4), ("four-b", 4), ("single-b", 1)]
RUN_SCHEMA = "koblitz_stage40_parallel_precompute.v1"
RUN_SEAL_SCHEMA = "koblitz_stage40_parallel_precompute_seal.v1"


class Stage40Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage40Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def parameters() -> dict[str, Any]:
    value = copy.deepcopy(load(BASE_PARAMS, "Stage-39 parameters"))
    value["name"] = "k0n41-fixed-algebraic-four-core-control-stage40"
    value["seed"] = RELATION_SEED
    value["targets"] = [{"public_hash_seed": seed} for seed in TARGET_SEEDS]
    return value


def predecessor_identity() -> dict[str, Any]:
    value = custody.executable_identity(PREDECESSOR, "Stage-39 hosted result", executable=False)
    value["path"] = str(PREDECESSOR.relative_to(REPO))
    return value


def self_test() -> dict[str, Any]:
    base = load(BASE_PARAMS, "Stage-39 parameters")
    current = parameters()
    for key in set(base) | set(current):
        if key not in {"name", "seed", "targets"}:
            require(base.get(key) == current.get(key), f"Stage-40 changed {key}")
    require(base["seed"] == 41331 and current["seed"] == RELATION_SEED, "Stage-40 relation stream was reused")
    require(base["targets"] != current["targets"] and len(set(TARGET_SEEDS)) == 5, "Stage-40 target stream was reused")
    predecessor = load(PREDECESSOR, "Stage-39 hosted result")
    require(predecessor.get("status") == "hosted_fixed_algebraic_holdout_verified", "Stage-40 predecessor is not verified")
    require(predecessor.get("source_commit") == "6f7ff404d615e11df73c0a0e0607bf504b2e09af", "Stage-40 predecessor source changed")
    require(CELLS == [("single-a", 1), ("four-a", 4), ("four-b", 4), ("single-b", 1)], "Stage-40 ABBA order changed")
    return {
        "schema": "koblitz_stage40_self_test.v1",
        "status": "PASS",
        "checks": 13,
        "cells": CELLS,
        "relation_seed": RELATION_SEED,
        "target_seeds": TARGET_SEEDS,
        "same_instance_in_every_cell": True,
        "factor_base_discovery_target_samples": 0,
        "factor_base_discovery_scalar_labels": 0,
        "target_scalars_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
    }


def set_affinity(cpus: list[int], initial: list[int]) -> dict[str, Any]:
    require(platform.system() == "Linux" and hasattr(os, "sched_setaffinity"), "Stage-40 requires Linux CPU affinity")
    require(bool(cpus) and set(cpus).issubset(initial), "Stage-40 requested unavailable CPUs")
    os.sched_setaffinity(0, set(cpus))
    parent = sorted(os.sched_getaffinity(0))
    probe = subprocess.run(
        [sys.executable, "-c", "import json,os;print(json.dumps(sorted(os.sched_getaffinity(0))))"],
        text=True, capture_output=True, check=True,
    )
    child = json.loads(probe.stdout)
    require(parent == child == cpus, "Stage-40 CPU affinity was not inherited")
    return {"initial_allowed_cpus": initial, "requested": cpus, "parent_effective": parent, "child_effective": child}


def semantic_fingerprint(workflow: dict[str, Any], workflow_dir: Path) -> dict[str, Any]:
    baseline = next(row["vs_rho"] for row in workflow["stages"] if row["stage"] == "baseline")
    logs = next(row for row in workflow["stages"] if row["stage"] == "logs")
    log_doc = load(workflow_dir / "logs.json", "Stage-40 log database")
    solutions = workflow["solutions"]["items"]
    return {
        "factor_base": workflow["factor_base"],
        "log_columns": log_doc["columns"],
        "relations": logs["relations"],
        "relation_trials": logs["trials"],
        "summands_scanned": logs["summands_scanned"],
        "solutions": [{key: row[key] for key in ("index", "expected", "recovered", "descent_trials", "target", "verified")} for row in solutions],
        "rho": [{key: row[key] for key in ("index", "verified", "iterations", "restarts", "walk_group_additions", "recovered", "target")} for row in baseline["targets_detail"]],
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux" and hasattr(os, "sched_getaffinity"), "Stage-40 requires Linux affinity")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-40 output must be new")
    build_root = args.build.resolve(strict=True)
    build = stage33.validate_build(build_root)
    binary = (build_root / "bin/ic").resolve(strict=True)
    initial = sorted(os.sched_getaffinity(0))
    require(len(initial) >= 4, "Stage-40 needs at least four allowed CPUs")
    four = initial[:4]
    output.mkdir(parents=True)
    cells_root = output / "cells"
    cells_root.mkdir()
    rows: list[dict[str, Any]] = []
    fingerprints: list[dict[str, Any]] = []
    for label, threads in CELLS:
        cpus = [four[0]] if threads == 1 else four
        affinity = set_affinity(cpus, initial)
        root = cells_root / label
        root.mkdir()
        evidence = root / "evidence"
        evidence.mkdir()
        params_path = root / "params.json"
        custody.write_json_new(params_path, parameters())
        workflow_dir = root / "workflow"
        environment = custody.safe_child_environment()
        environment["RAYON_NUM_THREADS"] = str(threads)
        command = [str(binary), "--json", "workflow", "--params", str(params_path), "--dir", str(workflow_dir)]
        before_self = resource.getrusage(resource.RUSAGE_SELF)
        before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
        started = time.monotonic()
        with TreeSampler() as sampler:
            process = stage33.metered(
                role="workflow", command=command, cwd=REPO, evidence=evidence,
                timeout=args.timeout, environment=environment,
            )
        outer = stage33.outer_resources(started, before_self, before_children, sampler)
        require(stage33.process_complete(process), f"Stage-40 {label} process failed")
        workflow = load(evidence / "workflow.stdout", f"Stage-40 {label} workflow")
        baseline = stage39.validate_workflow(workflow, workflow_dir, TARGET_SEEDS)
        fingerprint = semantic_fingerprint(workflow, workflow_dir)
        fingerprints.append(fingerprint)
        row = {
            "label": label,
            "threads": threads,
            "affinity": affinity,
            "params": custody.executable_identity(params_path, f"Stage-40 {label} parameters", executable=False),
            "process": process,
            "outer_resources": outer,
            "workflow": workflow,
            "precompute_seconds": baseline["ic"]["precompute_seconds"],
            "descent_seconds_total": baseline["ic"]["descent_seconds_total"],
            "five_target_ic_seconds": baseline["ic"]["precompute_seconds"] + baseline["ic"]["descent_seconds_total"],
            "rho_seconds_total": baseline["rho"]["seconds_total"],
            "ic_over_rho_amortised_wall_ratio": 1.0 / baseline["ratio"]["amortised"],
            "rho_over_ic_online_wall_ratio": baseline["ratio"]["charged"],
        }
        custody.write_json_new(root / "cell-result.json", row)
        rows.append(row)
    require(all(value == fingerprints[0] for value in fingerprints[1:]), "Stage-40 thread counts changed mathematical output")
    singles = [row for row in rows if row["threads"] == 1]
    fours = [row for row in rows if row["threads"] == 4]
    med = lambda values: statistics.median(values)
    single_pre = med([row["precompute_seconds"] for row in singles])
    four_pre = med([row["precompute_seconds"] for row in fours])
    single_total = med([row["five_target_ic_seconds"] for row in singles])
    four_total = med([row["five_target_ic_seconds"] for row in fours])
    summary = {
        "single_core_precompute_seconds_median": single_pre,
        "four_core_precompute_seconds_median": four_pre,
        "precompute_wall_speedup": single_pre / four_pre,
        "single_core_five_target_ic_seconds_median": single_total,
        "four_core_five_target_ic_seconds_median": four_total,
        "five_target_ic_wall_speedup": single_total / four_total,
        "single_core_process_core_seconds_median": med([row["process"]["metrics"]["total_core_seconds"] for row in singles]),
        "four_core_process_core_seconds_median": med([row["process"]["metrics"]["total_core_seconds"] for row in fours]),
        "four_core_amortised_ic_over_rho_wall_ratio_median": med([row["ic_over_rho_amortised_wall_ratio"] for row in fours]),
    }
    build_outer = build["outer_resources"]
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_abba_parallel_precompute_control",
        "source_state": build["source_state"],
        "binary": build["binary"],
        "build": {"result": custody.executable_identity(build_root / "result.json", "Stage-40 build result", executable=False), "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-40 build seal", executable=False), "process_totals": build["process_totals"], "outer_resources": build_outer},
        "predecessor_stage39": {"identity": predecessor_identity(), "verification": load(PREDECESSOR, "Stage-39 hosted result")["verification"]},
        "cells": CELLS,
        "rows": rows,
        "semantic_fingerprint": fingerprints[0],
        "summary": summary,
        "charged_total_core_seconds_available": build_outer["total_core_seconds"] + math.fsum(row["outer_resources"]["total_core_seconds"] for row in rows),
        "charged_sequential_wall_seconds_available": build_outer["wall_seconds"] + math.fsum(row["outer_resources"]["wall_seconds"] for row in rows),
        "maximum_sampled_process_tree_rss_bytes": max(build_outer["sampled_peak_process_tree_rss_bytes"], *(row["outer_resources"]["sampled_peak_process_tree_rss_bytes"] for row in rows)),
        "same_instance_in_every_cell": True,
        "factor_base_algebraically_defined": True,
        "factor_base_discovery_target_samples": 0,
        "factor_base_discovery_scalar_labels": 0,
        "target_subgroup_enumerated_for_factor_base": False,
        "target_scalars_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
        "sat_conflicts": None,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    custody.write_json_new(output / "result.json", result)
    stage33.freeze(output, RUN_SEAL_SCHEMA, "parallel_control_frozen")
    return result


def verify(build_root: Path, output: Path) -> dict[str, Any]:
    build = stage33.validate_build(build_root)
    seal = stage33.verify_seal(output, RUN_SEAL_SCHEMA, "parallel_control_frozen")
    result = load(output / "result.json", "Stage-40 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_abba_parallel_precompute_control", "Stage-40 result is incomplete")
    require(result.get("source_state") == build.get("source_state") and result.get("binary") == build.get("binary"), "Stage-40 build binding changed")
    rows = result.get("rows")
    require(isinstance(rows, list) and [(row.get("label"), row.get("threads")) for row in rows] == CELLS, "Stage-40 cell order changed")
    fingerprints = []
    for row in rows:
        root = output / "cells" / row["label"]
        workflow = row["workflow"]
        baseline = stage39.validate_workflow(workflow, root / "workflow", TARGET_SEEDS)
        fingerprints.append(semantic_fingerprint(workflow, root / "workflow"))
        require(row["precompute_seconds"] == baseline["ic"]["precompute_seconds"], "Stage-40 precompute changed")
    require(all(value == fingerprints[0] for value in fingerprints) and result.get("semantic_fingerprint") == fingerprints[0], "Stage-40 mathematical outputs differ")
    singles = [row for row in rows if row["threads"] == 1]
    fours = [row for row in rows if row["threads"] == 4]
    summary = result["summary"]
    single_pre = statistics.median(row["precompute_seconds"] for row in singles)
    four_pre = statistics.median(row["precompute_seconds"] for row in fours)
    require(summary["single_core_precompute_seconds_median"] == single_pre and summary["four_core_precompute_seconds_median"] == four_pre, "Stage-40 precompute medians changed")
    require(summary["precompute_wall_speedup"] == single_pre / four_pre, "Stage-40 speedup changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-40 claim widened")
    return {
        "schema": "koblitz_stage40_verification.v1",
        "status": "abba_parallel_precompute_verified",
        "source_commit": build["source_state"]["commit"],
        "summary": summary,
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
    build_parser.add_argument("--timeout", type=int, default=3600)
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
            value = stage33.build(args)
        elif args.command == "run":
            value = run(args)
        else:
            value = verify(args.build.resolve(strict=True), args.output.resolve(strict=True))
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage26Error, Stage40Error, stage39.Stage39Error, stage33.Stage33Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage40-parallel: {error}")


if __name__ == "__main__":
    main()

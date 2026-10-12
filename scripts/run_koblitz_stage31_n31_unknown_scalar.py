#!/usr/bin/env python3
"""Run five public n=31 unknown-scalar targets from a validated factor base."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
from run_koblitz_stage26_affinity_cell import TreeSampler


SCHEMA = "koblitz_stage31_n31_unknown_scalar_result.v1"
PUBLIC_TARGET_SEEDS = [31001, 31002, 31003, 31004, 31005]


class Stage31Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage31Error(message)


def read_json(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be an object")
    return value


def self_test() -> dict[str, Any]:
    require(len(PUBLIC_TARGET_SEEDS) == 5, "target count changed")
    require(len(set(PUBLIC_TARGET_SEEDS)) == 5, "target seeds are not distinct")
    require(all(type(seed) is int and seed > 0 for seed in PUBLIC_TARGET_SEEDS), "target seed is invalid")
    saturated = {"kind": "two_torsion_saturated", "parent": {"kind": "divisor", "indices": [0, 1, 5]}}
    generated = params({"spec": saturated})
    require(generated["factor_base"] == {"mode": "spec", "spec": saturated}, "factor-base handoff changed")
    require(generated["targets"] == [{"public_hash_seed": seed} for seed in PUBLIC_TARGET_SEEDS], "public targets changed")
    require(all("known_log" not in target and "random_seed" not in target for target in generated["targets"]), "Stage-31 parameters construct target scalars")
    require(generated["baseline"]["rho"] is True and generated["linear_algebra"]["mode"] == "sparse", "rho or sparse linear algebra was disabled")
    return {"schema": "koblitz_stage31_self_test.v1", "status": "PASS", "checks": 7}


def validate_factor_base(path: Path) -> dict[str, Any]:
    value = read_json(path, "Stage-30 factor-base recipe")
    require(
        value.get("schema_version") == 1
        and value.get("degree") == 31
        and value.get("curve_a") == 0
        and value.get("subfield", 1) == 1
        and value.get("curve_b", 1) == 1,
        "Stage-30 factor-base recipe names the wrong curve",
    )
    spec = value.get("spec")
    require(isinstance(spec, dict), "Stage-30 factor-base specification is invalid")
    if spec.get("kind") == "divisor":
        divisor = spec
    else:
        require(spec.get("kind") == "two_torsion_saturated", "Stage-30 factor base is not a divisor or its two-torsion saturation")
        divisor = spec.get("parent")
        require(isinstance(divisor, dict) and divisor.get("kind") == "divisor", "Stage-30 saturated factor base has no divisor parent")
    indices = divisor.get("indices")
    require(isinstance(indices, list) and indices and all(type(index) is int and index >= 0 for index in indices), "Stage-30 divisor indices are invalid")
    return value


def params(recipe: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "name": "k0n31-public-unknown-scalar-stage31",
        "curve": {"degree": 31, "curve_a": 0},
        "summands": 3,
        "solver": "pair_table",
        "seed": 31031,
        "max_trials": 200000,
        "linear_algebra": {
            "mode": "sparse",
            "sparse": {
                "wiedemann": {"block_m": 4, "block_n": 4},
                "filter": {"target_excess": 32, "merge_max_weight": 8},
            },
        },
        "collection": {"unit_trials": 4096, "units": 4, "max_units": 64},
        "baseline": {
            "rho": True,
            "rho_seed": 5931241263826122831,
            "rho_max_iterations": 268435456,
        },
        "factor_base": {"mode": "spec", "spec": recipe["spec"]},
        "targets": [{"public_hash_seed": seed} for seed in PUBLIC_TARGET_SEEDS],
    }


def validate_workflow(value: dict[str, Any]) -> None:
    require(value.get("status") == "complete", "Stage-31 workflow did not complete")
    solutions = value.get("solutions", {})
    items = solutions.get("items") if isinstance(solutions, dict) else None
    require(solutions.get("count") == 5 and solutions.get("verified") == 5 and isinstance(items, list) and len(items) == 5, "Stage-31 did not verify five solutions")
    for index, item in enumerate(items):
        target = item.get("target", {})
        require(item.get("index") == index and item.get("verified") is True, "Stage-31 solution order or verification changed")
        require(item.get("expected") == "not_constructed" and isinstance(item.get("recovered"), str), "Stage-31 solution used an expected scalar")
        require(
            target.get("kind") == "public_hash_to_curve_cofactor"
            and target.get("target_scalar_constructed") is False
            and target.get("public_hash_seed") == PUBLIC_TARGET_SEEDS[index]
            and isinstance(target.get("public_hash_counter"), int)
            and isinstance(target.get("x"), str)
            and isinstance(target.get("y"), str),
            "Stage-31 public target record changed",
        )
    stages = value.get("stages")
    require(isinstance(stages, list), "Stage-31 stage report is missing")
    baseline = next((stage.get("vs_rho") for stage in stages if stage.get("stage") == "baseline"), None)
    require(isinstance(baseline, dict), "Stage-31 rho baseline is missing")
    require(
        baseline.get("claim_boundary") == "public_hash_unknown_scalar"
        and baseline.get("targets") == 5
        and baseline.get("ic", {}).get("verified") == 5
        and baseline.get("rho", {}).get("verified") == 5,
        "Stage-31 IC/rho comparison is incomplete",
    )


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(not args.output.exists() and not args.output.is_symlink(), "output must be new")
    ic = args.ic.resolve(strict=True)
    meter = args.meter.resolve(strict=True)
    recipe = validate_factor_base(args.factor_base.resolve(strict=True))
    args.output.mkdir(parents=True)
    workflow_dir = args.output / "workflow"
    parameter_path = args.output / "params.json"
    custody.write_json_new(parameter_path, params(recipe))
    command = [str(ic), "--json", "workflow", "--params", str(parameter_path), "--dir", str(workflow_dir)]
    metrics_path = args.output / "workflow.metrics.json"
    stdout_path = args.output / "workflow.stdout"
    stderr_path = args.output / "workflow.stderr"
    meter_command = [
        str(Path(sys.executable).resolve()), str(meter), "--cwd", str(Path.cwd().resolve()),
        "--timeout", str(args.timeout), "--stdout", str(stdout_path), "--stderr", str(stderr_path),
        "--metrics", str(metrics_path), "--", *command,
    ]
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        process = subprocess.run(meter_command, stdin=subprocess.DEVNULL, close_fds=True, check=False)
    elapsed = time.monotonic() - started
    after_self = resource.getrusage(resource.RUSAGE_SELF)
    after_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    outer_core = (
        after_self.ru_utime - before_self.ru_utime + after_self.ru_stime - before_self.ru_stime
        + after_children.ru_utime - before_children.ru_utime
        + after_children.ru_stime - before_children.ru_stime
    )
    metrics = read_json(metrics_path, "Stage-31 process receipt") if metrics_path.is_file() else None
    workflow = None
    try:
        workflow = json.loads(stdout_path.read_text()) if stdout_path.is_file() else None
        if isinstance(workflow, dict):
            validate_workflow(workflow)
    except (json.JSONDecodeError, Stage31Error):
        workflow = None
    clean = (
        process.returncode == 0 and isinstance(metrics, dict)
        and metrics.get("returncode") == 0 and metrics.get("timed_out") is False
        and metrics.get("orphan_group_terminated") is False and isinstance(workflow, dict)
    )
    result = {
        "schema": SCHEMA,
        "status": "complete_public_unknown_scalar_panel" if clean else "incomplete_public_unknown_scalar_panel",
        "factor_base_recipe": custody.executable_identity(args.factor_base.resolve(), "Stage-30 factor base", executable=False),
        "params": custody.executable_identity(parameter_path, "Stage-31 params", executable=False),
        "public_target_seeds": PUBLIC_TARGET_SEEDS,
        "target_scalars_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
        "factor_base_log_source": "relation collection and group-certified linear algebra inside the workflow",
        "workflow_result": workflow,
        "process_receipt": metrics,
        "outer_resources": {
            "wall_seconds": elapsed,
            "total_core_seconds": outer_core,
            "average_parallelism": outer_core / elapsed if elapsed else None,
            "sampled_peak_process_tree_rss_bytes": sampler.peak,
            "process_tree_rss_samples": sampler.samples,
            "sample_interval_seconds": 0.02,
        },
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    custody.write_json_new(args.output / "result.json", result)
    if not clean:
        raise Stage31Error("n=31 public unknown-scalar workflow did not complete; retain output as incomplete")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("self-test")
    run_parser = sub.add_parser("run")
    run_parser.add_argument("--ic", type=Path, required=True)
    run_parser.add_argument("--factor-base", type=Path, required=True)
    run_parser.add_argument("--meter", type=Path, default=Path(__file__).with_name("process_meter.py"))
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--timeout", type=int, default=3600)
    args = parser.parse_args()
    try:
        value = self_test() if args.command == "self-test" else run(args)
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, Stage31Error) as error:
        raise SystemExit(f"stage31-unknown-scalar: {error}")


if __name__ == "__main__":
    main()

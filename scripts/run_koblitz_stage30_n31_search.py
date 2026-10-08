#!/usr/bin/env python3
"""Run the frozen n=31 factor-base search with parallel resource sampling."""

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

import run_koblitz_blind_pdp_phase_b as custody
from run_koblitz_stage26_affinity_cell import TreeSampler


SCHEMA = "koblitz_stage30_n31_factor_base_search.v1"
COMMAND_ARGUMENTS = [
    "search", "--degree", "31", "--curve-a", "0", "--summands", "3",
    "--family", "divisor", "--min-dimension", "11", "--max-dimension", "11",
    "--max-abscissae", "2048", "--targets", "256", "--validate-top", "1",
    "--holdout", "2", "--max-trials", "10000", "--timeout-seconds", "300",
    "--solver", "pair-table", "--seed", "29031", "--json",
]


class Stage30Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage30Error(message)


def self_test() -> dict:
    require(COMMAND_ARGUMENTS.count("--targets") == 1, "target count is ambiguous")
    require(COMMAND_ARGUMENTS[COMMAND_ARGUMENTS.index("--targets") + 1] == "256", "target count changed")
    require(COMMAND_ARGUMENTS[COMMAND_ARGUMENTS.index("--validate-top") + 1] == "1", "validation count changed")
    require(COMMAND_ARGUMENTS[COMMAND_ARGUMENTS.index("--holdout") + 1] == "2", "holdout count changed")
    require(COMMAND_ARGUMENTS[COMMAND_ARGUMENTS.index("--min-dimension") + 1] == "11", "minimum dimension changed")
    require(COMMAND_ARGUMENTS[COMMAND_ARGUMENTS.index("--max-dimension") + 1] == "11", "maximum dimension changed")
    return {"schema": "koblitz_stage30_n31_search_self_test.v1", "status": "PASS", "checks": 6}


def run(args: argparse.Namespace) -> dict:
    require(not args.output.exists() and not args.output.is_symlink(), "output must be new")
    ic = args.ic.resolve(strict=True)
    meter = args.meter.resolve(strict=True)
    args.output.mkdir(parents=True)
    search_result = args.output / "search.json"
    factor_base = args.output / "factor-base.json"
    command = [str(ic), *COMMAND_ARGUMENTS, "--out", str(search_result), "--spec-out", str(factor_base)]
    metrics_path = args.output / "search.metrics.json"
    meter_command = [
        str(Path(sys.executable).resolve()), str(meter),
        "--cwd", str(Path.cwd().resolve()), "--timeout", str(args.timeout),
        "--stdout", str(args.output / "search.stdout"),
        "--stderr", str(args.output / "search.stderr"),
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
    metrics = json.loads(metrics_path.read_text()) if metrics_path.is_file() else None
    clean = (
        process.returncode == 0 and isinstance(metrics, dict)
        and metrics.get("returncode") == 0 and metrics.get("timed_out") is False
        and metrics.get("orphan_group_terminated") is False
        and search_result.is_file() and factor_base.is_file()
    )
    result = {
        "schema": SCHEMA,
        "status": "complete_validated_factor_base" if clean else "incomplete_search",
        "command": command,
        "fixed_parameters": {
            "degree": 31, "curve_a": 0, "summands": 3, "family": "divisor",
            "dimension": 11, "sampled_targets": 256, "validated_candidates": 1,
            "holdout_targets": 2, "search_seed": 29031,
        },
        "information_boundary": {
            "future_unknown_scalar_target_available": False,
            "target_subgroup_enumerated": False,
            "factor_base_discrete_log_labels_used": False,
            "relation_yield_and_holdouts": "public search samples and separate generated validation controls only",
        },
        "initial_allowed_cpus": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
        "process_meter_returncode": process.returncode,
        "search_process": metrics,
        "outer_resources": {
            "wall_seconds": elapsed,
            "total_core_seconds": outer_core,
            "average_parallelism": outer_core / elapsed if elapsed else None,
            "sampled_peak_process_tree_rss_bytes": sampler.peak,
            "process_tree_rss_samples": sampler.samples,
            "sample_interval_seconds": 0.02,
        },
        "artifacts": {
            "search": custody.executable_identity(search_result, "search result", executable=False) if search_result.is_file() else None,
            "factor_base": custody.executable_identity(factor_base, "factor-base recipe", executable=False) if factor_base.is_file() else None,
        },
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    custody.write_json_new(args.output / "wrapper-result.json", result)
    if not clean:
        raise Stage30Error("n=31 search did not complete; retain the output as incomplete evidence")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("self-test")
    run_parser = sub.add_parser("run")
    run_parser.add_argument("--ic", type=Path, required=True)
    run_parser.add_argument("--meter", type=Path, default=Path(__file__).with_name("process_meter.py"))
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--timeout", type=int, default=3600)
    args = parser.parse_args()
    try:
        result = self_test() if args.command == "self-test" else run(args)
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, ValueError, Stage30Error) as error:
        raise SystemExit(f"stage30-search: {error}")


if __name__ == "__main__":
    main()

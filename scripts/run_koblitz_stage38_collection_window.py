#!/usr/bin/env python3
"""Sweep relation-collection windows on the Stage-35 algebraic n=41 panel."""

from __future__ import annotations

import argparse
import copy
import json
import math
import platform
from pathlib import Path
import resource
import subprocess
import time
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
from run_koblitz_stage26_affinity_cell import Stage26Error, TreeSampler, singleton_affinity
import run_koblitz_stage33_n41_unknown_scalar as stage33
import run_koblitz_stage35_algebraic_walk as stage35


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
BASE_PARAMS = STAGE / "stage-35-n41-algebraic-walk-params.json"
WINDOWS: list[int | None] = [None, 595, 297, 149, 74, 37]
COLLECTION = {"unit_trials": 8192, "units": 1, "max_units": 64}
RUN_SCHEMA = "koblitz_stage38_algebraic_collection_window_sweep.v1"
RUN_SEAL_SCHEMA = "koblitz_stage38_algebraic_collection_window_sweep_seal.v1"


class Stage38Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage38Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def label(window: int | None) -> str:
    return "full" if window is None else f"window-{window}"


def parameters(window: int | None) -> dict[str, Any]:
    value = copy.deepcopy(load(BASE_PARAMS, "Stage-35 base parameters"))
    value["name"] = f"k0n41-algebraic-collection-{label(window)}-stage38"
    value["collection"] = COLLECTION
    if window is None:
        value.pop("collection_window", None)
    else:
        value["collection_window"] = window
    return value


def self_test() -> dict[str, Any]:
    base = load(BASE_PARAMS, "Stage-35 parameters")
    require(len(WINDOWS) == len(set(WINDOWS)) == 6 and WINDOWS[0] is None, "Stage-38 window set changed")
    require(all(window is None or 0 < window < 4759 for window in WINDOWS), "Stage-38 window is invalid")
    allowed = {"name", "collection", "collection_window"}
    for window in WINDOWS:
        value = parameters(window)
        for key in set(base) | set(value):
            if key not in allowed:
                require(base.get(key) == value.get(key), f"Stage-38 {label(window)} changed {key}")
        require(value["collection"] == COLLECTION, "Stage-38 work-unit boundary changed")
        require(value.get("collection_window") == window, "Stage-38 window label changed")
    return {
        "schema": "koblitz_stage38_self_test.v1",
        "status": "PASS",
        "checks": 10,
        "windows": WINDOWS,
        "collection": COLLECTION,
        "factor_base_algebraically_defined": True,
        "target_subgroup_enumerated": False,
        "target_scalars_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
    }


def validate_cell(value: dict[str, Any], root: Path, window: int | None) -> dict[str, Any]:
    baseline = stage35.validate_workflow(value, root / "workflow")
    require(value.get("collection_window") == window, f"Stage-38 {label(window)} report changed its window")
    require(value.get("summands") == 3 and value.get("descent_summands") == 2, "Stage-38 summand split changed")
    logs = next(row for row in value["stages"] if row["stage"] == "logs")
    collect = next(row for row in value["stages"] if row["stage"] == "collect")
    points = value["factor_base"]["points"]
    expected_per_trial = points if window is None else window
    require(logs.get("summands_scanned") == logs.get("trials") * expected_per_trial, f"Stage-38 {label(window)} lost its final lookup charge")
    require(collect.get("summands_scanned_total") == collect.get("trials_total") * expected_per_trial, f"Stage-38 {label(window)} lost its initial lookup charge")
    require(collect.get("summands_scanned_now") == collect.get("summands_scanned_total"), f"Stage-38 {label(window)} initial units were not fresh")
    require(logs["summands_scanned"] > 0 and logs["relations"] > 0, "Stage-38 collection charge is empty")
    return baseline


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-38 requires Linux singleton affinity")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-38 output must be new")
    build_root = args.build.resolve(strict=True)
    build = stage33.validate_build(build_root)
    binary = (build_root / "bin/ic").resolve(strict=True)
    affinity = singleton_affinity(args.cpu)
    output.mkdir(parents=True)
    cells_root = output / "cells"
    cells_root.mkdir()
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    rows: list[dict[str, Any]] = []
    with TreeSampler() as sampler:
        for window in WINDOWS:
            name = label(window)
            root = cells_root / name
            root.mkdir()
            evidence = root / "evidence"
            evidence.mkdir()
            params_path = root / "params.json"
            custody.write_json_new(params_path, parameters(window))
            workflow_root = root / "workflow"
            command = [str(binary), "--json", "workflow", "--params", str(params_path), "--dir", str(workflow_root)]
            process = stage33.metered(
                role="workflow", command=command, cwd=REPO, evidence=evidence,
                timeout=args.timeout, environment=custody.safe_child_environment(),
            )
            workflow: dict[str, Any] | None = None
            baseline: dict[str, Any] | None = None
            try:
                workflow = load(evidence / "workflow.stdout", f"Stage-38 {name} workflow")
                baseline = validate_cell(workflow, root, window)
            except (OSError, ValueError, stage35.Stage35Error, Stage38Error):
                workflow = None
                baseline = None
            require(stage33.process_complete(process) and workflow is not None and baseline is not None, f"Stage-38 {name} did not complete")
            logs = next(item for item in workflow["stages"] if item["stage"] == "logs")
            ic = baseline["ic"]
            collect_seconds = workflow["state"]["collect"]["elapsed_seconds"]
            logs_seconds = workflow["state"]["logs"]["elapsed_seconds"]
            row = {
                "label": name,
                "collection_window": window,
                "params": custody.executable_identity(params_path, f"Stage-38 {name} parameters", executable=False),
                "workflow": workflow,
                "process": process,
                "factor_base": workflow["factor_base"],
                "relations": logs["relations"],
                "trials": logs["trials"],
                "units": logs["units_used"],
                "summands_scanned": logs["summands_scanned"],
                "summands_scanned_per_relation": logs["summands_scanned"] / logs["relations"],
                "select_seconds": workflow["state"]["select"]["elapsed_seconds"],
                "collect_seconds": collect_seconds,
                "logs_seconds": logs_seconds,
                "collection_and_logs_seconds": collect_seconds + logs_seconds,
                "precompute_seconds": ic["precompute_seconds"],
                "descent_seconds_total": ic["descent_seconds_total"],
                "five_target_ic_seconds": ic["precompute_seconds"] + ic["descent_seconds_total"],
                "rho_seconds_total": baseline["rho"]["seconds_total"],
                "rho_over_ic_online_wall_ratio": baseline["ratio"]["charged"],
                "rho_over_ic_amortised_wall_ratio": baseline["ratio"]["amortised"],
                "targets_verified": workflow["solutions"]["verified"],
            }
            custody.write_json_new(root / "cell-result.json", row)
            rows.append(row)
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    # The window changes relation collection and the log solve fed by those
    # relations. Select on those stages only: repeating the identical factor-
    # base search and target descent in every cell is still charged, but their
    # timing noise must not choose the window.
    winner = min(rows, key=lambda row: row["collection_and_logs_seconds"])
    full = next(row for row in rows if row["collection_window"] is None)
    build_outer = build["outer_resources"]
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_six_window_algebraic_sweep",
        "source_state": build["source_state"],
        "binary": build["binary"],
        "build": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-38 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-38 build seal", executable=False),
            "process_totals": build["process_totals"],
            "outer_resources": build_outer,
        },
        "affinity": affinity,
        "windows": WINDOWS,
        "collection": COLLECTION,
        "rows": rows,
        "winner": {
            "label": winner["label"],
            "collection_window": winner["collection_window"],
            "selection_metric": "collection_and_logs_seconds",
            "collection_and_logs_seconds": winner["collection_and_logs_seconds"],
            "five_target_ic_seconds": winner["five_target_ic_seconds"],
            "precompute_seconds": winner["precompute_seconds"],
            "summands_scanned": winner["summands_scanned"],
            "summands_scanned_per_relation": winner["summands_scanned_per_relation"],
            "relative_collection_and_logs_speedup_over_full_scan": full["collection_and_logs_seconds"] / winner["collection_and_logs_seconds"],
            "relative_five_target_speedup_over_full_scan": full["five_target_ic_seconds"] / winner["five_target_ic_seconds"],
            "relative_lookup_reduction_from_full_scan": full["summands_scanned"] / winner["summands_scanned"],
        },
        "outer_resources": outer,
        "charged_total_core_seconds_available": build_outer["total_core_seconds"] + outer["total_core_seconds"],
        "charged_sequential_wall_seconds_available": build_outer["wall_seconds"] + outer["wall_seconds"],
        "maximum_sampled_process_tree_rss_bytes": max(build_outer["sampled_peak_process_tree_rss_bytes"], outer["sampled_peak_process_tree_rss_bytes"]),
        "factor_base_algebraically_defined": True,
        "target_subgroup_enumerated_for_factor_base": False,
        "target_scalars_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
        "same_public_targets_in_every_cell": True,
        "sat_conflicts": None,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    custody.write_json_new(output / "result.json", result)
    stage33.freeze(output, RUN_SEAL_SCHEMA, "sweep_frozen")
    return result


def verify(build_root: Path, output: Path) -> dict[str, Any]:
    build = stage33.validate_build(build_root)
    seal = stage33.verify_seal(output, RUN_SEAL_SCHEMA, "sweep_frozen")
    result = load(output / "result.json", "Stage-38 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_six_window_algebraic_sweep", "Stage-38 result is incomplete")
    require(result.get("source_state") == build.get("source_state") and result.get("binary") == build.get("binary"), "Stage-38 build binding changed")
    rows = result.get("rows")
    require(isinstance(rows, list) and len(rows) == len(WINDOWS), "Stage-38 row count changed")
    for expected, row in zip(WINDOWS, rows, strict=True):
        require(row.get("collection_window") == expected, "Stage-38 row order changed")
        baseline = validate_cell(row["workflow"], output / "cells" / label(expected), expected)
        require(row.get("rho_over_ic_online_wall_ratio") == baseline["ratio"]["charged"], "Stage-38 online ratio changed")
        require(row.get("rho_over_ic_amortised_wall_ratio") == baseline["ratio"]["amortised"], "Stage-38 amortised ratio changed")
    winner = min(rows, key=lambda row: row["collection_and_logs_seconds"])
    require(result.get("winner", {}).get("label") == winner["label"], "Stage-38 winner changed")
    require(result["winner"].get("selection_metric") == "collection_and_logs_seconds", "Stage-38 selection metric changed")
    require(result.get("factor_base_algebraically_defined") is True and result.get("target_subgroup_enumerated_for_factor_base") is False, "Stage-38 factor-base boundary changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-38 claim widened")
    return {
        "schema": "koblitz_stage38_verification.v1",
        "status": "six_algebraic_collection_windows_verified",
        "source_commit": build["source_state"]["commit"],
        "windows": WINDOWS,
        "winner": result["winner"],
        "rows": [{key: row[key] for key in (
            "label", "collection_window", "relations", "trials", "units",
            "summands_scanned", "summands_scanned_per_relation", "select_seconds",
            "collect_seconds", "logs_seconds", "collection_and_logs_seconds", "precompute_seconds",
            "descent_seconds_total", "five_target_ic_seconds", "rho_seconds_total",
            "rho_over_ic_online_wall_ratio", "rho_over_ic_amortised_wall_ratio",
        )} for row in rows],
        "sweep_outer_core_seconds": result["outer_resources"]["total_core_seconds"],
        "sweep_outer_wall_seconds": result["outer_resources"]["wall_seconds"],
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
    run_parser.add_argument("--cpu", type=int)
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
    except (
        OSError, ValueError, KeyError, subprocess.CalledProcessError,
        Stage26Error, Stage38Error, stage35.Stage35Error, stage33.Stage33Error,
        custody.PhaseBError,
    ) as error:
        raise SystemExit(f"stage38-window: {error}")


if __name__ == "__main__":
    main()

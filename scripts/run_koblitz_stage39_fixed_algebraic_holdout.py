#!/usr/bin/env python3
"""Run and verify the Stage-39 fixed algebraic n=41 holdout."""

from __future__ import annotations

import argparse
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


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
PARAMS = STAGE / "stage-39-n41-fixed-algebraic-holdout-params.json"
TARGET_SEEDS = [41201, 41202, 41203, 41204, 41205]
RELATION_SEED = 41231
WINDOW = 149
ALGEBRAIC_SPEC = {
    "kind": "two_torsion_saturated",
    "parent": {"kind": "frobenius_union", "seed_masks": [1, 2, 4, 8, 16, 32]},
}
RUN_SCHEMA = "koblitz_stage39_n41_fixed_algebraic_holdout.v1"
RUN_SEAL_SCHEMA = "koblitz_stage39_n41_fixed_algebraic_holdout_seal.v1"


class Stage39Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage39Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def validate_params(path: Path = PARAMS) -> dict[str, Any]:
    value = load(path, "Stage-39 parameters")
    require(value.get("schema_version") == 1, "Stage-39 schema changed")
    require(value.get("curve") == {"degree": 41, "curve_a": 0}, "Stage-39 curve changed")
    require(value.get("summands") == 3 and value.get("descent_summands") == 2, "Stage-39 summand split changed")
    require(value.get("solver") == "pair_table", "Stage-39 solver changed")
    require(value.get("collection_window") == WINDOW, "Stage-39 selected window changed")
    require(value.get("seed") == RELATION_SEED and value.get("max_trials") == 100_000_000, "Stage-39 seed or walked-probe cap changed")
    require(value.get("targets") == [{"public_hash_seed": seed} for seed in TARGET_SEEDS], "Stage-39 targets changed")
    require(all("known_log" not in row and "random_seed" not in row for row in value["targets"]), "Stage-39 constructs target logs")
    factor = value.get("factor_base", {})
    require(
        factor.get("mode") == "spec" and factor.get("spec") == ALGEBRAIC_SPEC,
        "Stage-39 factor base is not the fixed algebraic rule",
    )
    require(value.get("linear_algebra", {}).get("mode") == "sparse", "Stage-39 sparse linear algebra was disabled")
    require(value.get("baseline", {}).get("rho") is True, "Stage-39 rho baseline was disabled")
    return value


def self_test() -> dict[str, Any]:
    current = validate_params()
    require(current["collection"] == {"unit_trials": 8192, "units": 4, "max_units": 64}, "Stage-39 collection boundary changed")
    require(current["factor_base"].keys() == {"mode", "spec"}, "Stage-39 factor-base rule contains hidden search inputs")
    return {
        "schema": "koblitz_stage39_self_test.v1",
        "status": "PASS",
        "checks": 14,
        "factor_base_family": "fixed_two_torsion_saturated_frobenius_union",
        "factor_base_recipe": ALGEBRAIC_SPEC,
        "factor_base_discovery_target_samples": 0,
        "factor_base_discovery_scalar_labels": 0,
        "target_subgroup_enumerated": False,
        "target_scalars_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
        "collection_summands": 3,
        "descent_summands": 2,
        "walked_probe_cap": 100_000_000,
        "collection_window": WINDOW,
        "relation_seed_is_holdout": True,
        "public_target_seeds_are_holdout": True,
    }


def validate_workflow(value: dict[str, Any], workflow_dir: Path) -> dict[str, Any]:
    require(value.get("operation") == "workflow" and value.get("status") == "complete", "Stage-39 workflow is incomplete")
    require(value.get("evidence_scope") == "public_hash_unknown_scalar", "Stage-39 evidence label changed")
    require(value.get("degree") == 41 and value.get("summands") == 3 and value.get("descent_summands") == 2, "Stage-39 workflow summands changed")
    require(value.get("collection_window") == WINDOW, "Stage-39 workflow window changed")
    require(value.get("solver") == "pair_table", "Stage-39 workflow solver changed")
    factor = value.get("factor_base", {})
    require(
        factor.get("spec") == ALGEBRAIC_SPEC
        and factor.get("abscissae") == 2380
        and factor.get("points") == 4759
        and factor.get("signed_orbits") == 60
        and factor.get("columns") == 29,
        "Stage-39 fixed algebraic factor base changed",
    )
    stages = value.get("stages")
    require(isinstance(stages, list), "Stage-39 stage reports are missing")
    by_name = {row.get("stage"): row for row in stages}
    require(set(by_name) == {"select", "collect", "logs", "solve", "baseline"}, "Stage-39 stage inventory changed")
    require(
        by_name["select"].get("ran") is True
        and "search" not in by_name["select"]
        and value.get("state", {}).get("select", {}).get("elapsed_seconds", 0) > 0,
        "Stage-39 algebraic factor-base construction was not charged",
    )
    collect = by_name["collect"]
    logs = by_name["logs"]
    linear = logs.get("linear_algebra", {})
    require(
        logs.get("status") == "complete"
        and logs.get("columns") == factor["columns"]
        and logs.get("relations", 0) > 0
        and logs.get("rejected") == 0
        and linear.get("mode") == "sparse"
        and linear.get("sparse", {}).get("filter", {}).get("uncovered_columns") == 0,
        "Stage-39 relation/log solve is incomplete",
    )
    require(
        collect.get("summands_scanned_total") == collect.get("trials_total") * WINDOW
        and logs.get("summands_scanned") == logs.get("trials") * WINDOW,
        "Stage-39 third-summand lookup charge changed",
    )
    solutions = value.get("solutions", {})
    items = solutions.get("items") if isinstance(solutions, dict) else None
    require(solutions.get("count") == solutions.get("verified") == 5 and isinstance(items, list), "Stage-39 did not verify five solutions")
    trials = 0
    for index, item in enumerate(items):
        target = item.get("target", {})
        require(
            item.get("index") == index
            and item.get("expected") == "not_constructed"
            and item.get("verified") is True
            and isinstance(item.get("recovered"), str)
            and isinstance(item.get("descent_trials"), int)
            and 0 < item["descent_trials"] <= 100_000_000
            and target.get("kind") == "public_hash_to_curve_cofactor"
            and target.get("target_scalar_constructed") is False
            and target.get("public_hash_seed") == TARGET_SEEDS[index],
            f"Stage-39 target {index} changed or failed",
        )
        trials += item["descent_trials"]
    baseline = by_name["baseline"].get("vs_rho", {})
    require(
        baseline.get("claim_boundary") == "public_hash_unknown_scalar"
        and baseline.get("n") == 41
        and baseline.get("targets") == 5
        and baseline.get("ic", {}).get("verified") == 5
        and baseline.get("rho", {}).get("verified") == 5
        and baseline.get("ic", {}).get("descent_trials_total") == trials,
        "Stage-39 IC/rho comparison is incomplete",
    )
    for name in ("charged", "amortised"):
        require(isinstance(baseline.get("ratio", {}).get(name), (int, float)) and baseline["ratio"][name] > 0, f"Stage-39 {name} ratio is invalid")
    factor_doc = load(workflow_dir / "factor_base.json", "Stage-39 factor-base document")
    logs_doc = load(workflow_dir / "logs.json", "Stage-39 logs document")
    solutions_doc = load(workflow_dir / "solutions.json", "Stage-39 solutions document")
    baseline_doc = load(workflow_dir / "baseline.json", "Stage-39 baseline document")
    require(factor_doc.get("spec") == factor.get("spec"), "Stage-39 factor-base document differs")
    require(logs_doc.get("spec") == factor_doc.get("spec") and len(logs_doc.get("columns", [])) == factor["columns"], "Stage-39 log database differs")
    require(solutions_doc.get("solutions") == items and baseline_doc == baseline, "Stage-39 solution or baseline database differs")
    relation_files = sorted((workflow_dir / "relations").glob("unit-*.json"))
    require(len(relation_files) == logs.get("units_used"), "Stage-39 relation-unit inventory changed")
    for path in relation_files:
        unit = load(path, f"Stage-39 relation unit {path.name}")
        require(
            unit.get("summands_scanned") == unit.get("count") * WINDOW,
            f"Stage-39 relation unit {path.name} lost its lookup charge",
        )
    return baseline


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-39 requires Linux singleton affinity")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-39 output must be new")
    params = args.params.resolve(strict=True)
    validate_params(params)
    build_root = args.build.resolve(strict=True)
    build_result = stage33.validate_build(build_root)
    affinity = singleton_affinity(args.cpu)
    output.mkdir(parents=True)
    params_copy = output / "params.json"
    custody.write_new(params_copy, params.read_bytes())
    evidence = output / "evidence"
    evidence.mkdir()
    workflow_dir = output / "workflow"
    command = [str((build_root / "bin/ic").resolve(strict=True)), "--json", "workflow", "--params", str(params_copy), "--dir", str(workflow_dir)]
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        process = stage33.metered(
            role="workflow", command=command, cwd=REPO, evidence=evidence,
            timeout=args.timeout, environment=custody.safe_child_environment(),
        )
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    workflow: dict[str, Any] | None = None
    baseline: dict[str, Any] | None = None
    try:
        workflow = load(evidence / "workflow.stdout", "Stage-39 workflow output")
        baseline = validate_workflow(workflow, workflow_dir)
    except (OSError, ValueError, Stage39Error):
        workflow = None
        baseline = None
    complete = stage33.process_complete(process) and workflow is not None and baseline is not None
    build_outer = build_result["outer_resources"]
    rho_wall = baseline["rho"]["seconds_total"] if baseline else None
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n41_fixed_algebraic_holdout" if complete else "incomplete_n41_fixed_algebraic_holdout",
        "params": custody.executable_identity(params_copy, "Stage-39 parameters", executable=False),
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-39 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-39 build seal", executable=False),
            "binary": build_result["binary"],
            "source_state": build_result["source_state"],
            "requested_parallel_jobs": build_result["requested_parallel_jobs"],
            "process_totals": build_result["process_totals"],
            "outer_resources": build_outer,
        },
        "affinity": affinity,
        "workflow_process": process,
        "outer_resources": {
            **outer,
            "single_core_elapsed_seconds": outer["wall_seconds"],
            "one_cpu_utilization_percent": 100.0 * outer["total_core_seconds"] / outer["wall_seconds"],
        },
        "workflow_result": workflow,
        "factor_base_discovery": {
            "method": "fixed_algebraic_rule",
            "public_parameters": {"n": 41, "ell": 6, "collection_m": 3},
            "recipe": ALGEBRAIC_SPEC,
            "target_samples": 0,
            "target_scalar_labels": 0,
            "target_subgroup_enumerated": False,
            "predicate_materialization_charged": True,
            "construction_elapsed_seconds": workflow["state"]["select"]["elapsed_seconds"] if workflow else None,
        },
        "factor_base_algebraically_defined": True,
        "target_subgroup_enumerated_for_factor_base": False,
        "target_scalars_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
        "collection_summands": 3,
        "descent_summands": 2,
        "collection_window": WINDOW,
        "relation_seed": RELATION_SEED,
        "holdout_target_seeds": TARGET_SEEDS,
        "walked_probe_cap": 100_000_000,
        "sat_conflicts": None,
        "sat_conflict_semantics": "not applicable to the direct pair-table individual-log arm; Phase B retains the same-instance SAT conflict matrix",
        "charged_total_core_seconds_available": build_outer["total_core_seconds"] + outer["total_core_seconds"],
        "charged_sequential_wall_seconds_available": build_outer["wall_seconds"] + outer["wall_seconds"],
        "maximum_sampled_process_tree_rss_bytes": max(build_outer["sampled_peak_process_tree_rss_bytes"], outer["sampled_peak_process_tree_rss_bytes"]),
        "full_available_wall_over_same_targets_rho": (build_outer["wall_seconds"] + outer["wall_seconds"]) / rho_wall if rho_wall else None,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    custody.write_json_new(output / "result.json", result)
    stage33.freeze(output, RUN_SEAL_SCHEMA, "run_frozen")
    require(complete, "Stage-39 workflow did not complete; retain its sealed output")
    return result


def verify(build_root: Path, output: Path) -> dict[str, Any]:
    build_result = stage33.validate_build(build_root)
    seal = stage33.verify_seal(output, RUN_SEAL_SCHEMA, "run_frozen")
    result = load(output / "result.json", "Stage-39 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n41_fixed_algebraic_holdout", "Stage-39 result is incomplete")
    affinity = result.get("affinity", {})
    require(affinity.get("parent_effective") == affinity.get("child_effective") == [affinity.get("selected_cpu")], "Stage-39 affinity changed")
    workflow = result.get("workflow_result")
    require(isinstance(workflow, dict), "Stage-39 workflow result is absent")
    baseline = validate_workflow(workflow, output / "workflow")
    require(result.get("build_binding", {}).get("binary") == build_result.get("binary"), "Stage-39 build binding changed")
    require(result.get("factor_base_algebraically_defined") is True and result.get("target_subgroup_enumerated_for_factor_base") is False, "Stage-39 factor-base boundary changed")
    require(result.get("target_scalars_constructed_or_supplied") is False and result.get("factor_base_logs_known_by_construction") is False, "Stage-39 label boundary changed")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-39 claim widened")
    discovery = result.get("factor_base_discovery", {})
    require(
        discovery.get("method") == "fixed_algebraic_rule"
        and discovery.get("recipe") == ALGEBRAIC_SPEC
        and discovery.get("target_samples") == discovery.get("target_scalar_labels") == 0
        and discovery.get("target_subgroup_enumerated") is False
        and discovery.get("predicate_materialization_charged") is True
        and discovery.get("construction_elapsed_seconds") == workflow["state"]["select"]["elapsed_seconds"],
        "Stage-39 algebraic discovery accounting changed",
    )
    ratio = baseline["ratio"]
    return {
        "schema": "koblitz_stage39_verification.v1",
        "status": "n41_fixed_algebraic_holdout_build_and_run_verified",
        "source_commit": build_result["source_state"]["commit"],
        "factor_base": workflow["factor_base"],
        "factor_base_discovery": discovery,
        "relations": next(row["relations"] for row in workflow["stages"] if row["stage"] == "logs"),
        "relation_trials": next(row["trials"] for row in workflow["stages"] if row["stage"] == "logs"),
        "summands_scanned": next(row["summands_scanned"] for row in workflow["stages"] if row["stage"] == "logs"),
        "targets_verified": workflow["solutions"]["verified"],
        "descent_trials": baseline["ic"]["descent_trials_total"],
        "rho_iterations": baseline["rho"]["iterations_total"],
        "rho_over_ic_online_wall_ratio": ratio["charged"],
        "rho_over_ic_amortised_wall_ratio": ratio["amortised"],
        "ic_over_rho_online_wall_ratio": 1.0 / ratio["charged"],
        "ic_over_rho_amortised_wall_ratio": 1.0 / ratio["amortised"],
        "single_core_elapsed_seconds": result["outer_resources"]["single_core_elapsed_seconds"],
        "scientific_outer_core_seconds": result["outer_resources"]["total_core_seconds"],
        "charged_total_core_seconds_available": result["charged_total_core_seconds_available"],
        "charged_sequential_wall_seconds_available": result["charged_sequential_wall_seconds_available"],
        "maximum_sampled_process_tree_rss_bytes": result["maximum_sampled_process_tree_rss_bytes"],
        "full_available_wall_over_same_targets_rho": result["full_available_wall_over_same_targets_rho"],
        "sat_conflicts": None,
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
    run_parser.add_argument("--params", type=Path, default=PARAMS)
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--timeout", type=int, default=3600)
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
        OSError,
        ValueError,
        KeyError,
        subprocess.CalledProcessError,
        Stage26Error,
        Stage39Error,
        stage33.Stage33Error,
        custody.PhaseBError,
    ) as error:
        raise SystemExit(f"stage39-fixed-algebraic-holdout: {error}")


if __name__ == "__main__":
    main()

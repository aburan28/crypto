#!/usr/bin/env python3
"""Compose and verify optimized unknown-scalar and single-core n=53 results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage57_n53_unknown as stage57
import run_koblitz_stage71_single_core as stage71
import verify_koblitz_stage49_width_results as common


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
PREDECESSOR = STAGE / "stage-68-n53-fused-result-20260912/verification.json"
DEFAULT_OUTPUT = STAGE / "stage-72-current-results-20260912"
SCHEMA = "koblitz_stage72_current_results.v1"
SEAL_SCHEMA = "koblitz_stage72_current_results_seal.v1"
UNKNOWN_CONFIG = {
    "run_id": 34719465029,
    "commit": "7e693d9d778ca498ecddb273cfee1d2f13d05fc1",
    "artifact_id": 10306063039,
    "artifact_digest": "sha256:67983a67097d48fe5f69b06504b34f38068e441b089cf6808142fdea4ee5054c",
    "archive": "koblitz-stage70-n53-unknown-evidence-34719465029.tar.gz",
    "archive_root": "koblitz-stage70-n53-unknown-34719465029",
    "job": "n53-unknown-optimized",
}
SINGLE_CONFIG = {
    "run_id": 34719948884,
    "commit": "66b344f3a68a57b9fef6ec74e4d8c8d030753d4c",
    "artifact_id": 10305634413,
    "artifact_digest": "sha256:53ddc6434e4606c71bfdb5b642d2dddbea1c52bb08ab5880147befda7b3b1f78",
    "archive": "koblitz-stage71-single-core-evidence-34719948884.tar.gz",
    "archive_root": "koblitz-stage71-single-core-34719948884",
    "job": "n53-single-core",
}


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def unknown_result() -> dict[str, Any]:
    with tempfile.TemporaryDirectory() as directory:
        root = common.extract(70, UNKNOWN_CONFIG, Path(directory))
        verification = stage57.verify(root / "stage69-build", root / "stage69-run")
        embedded = load(root / "stage69-verification.json", "embedded Stage-70 verification")
        require(verification == embedded, "Stage-70 embedded verification changed")
        run = load(root / "stage69-run/result.json", "Stage-70 run")
        build_seal = custody.sha256_file(root / "stage69-build/result-seal.json", "Stage-70 build seal")
        run_seal = custody.sha256_file(root / "stage69-run/result-seal.json", "Stage-70 run seal")
    summary = run["direct_summary"]
    require(run.get("optimized_stack") is True, "Stage-70 optimized stack changed")
    require(verification.get("target_scalar_constructed_or_supplied") is False, "Stage-70 target scalar was supplied")
    require(verification.get("factor_base_logs_known_by_construction") is False, "Stage-70 factor-base logs were supplied")
    require(summary.get("rank_target_deficiency") == 3, "Stage-70 rank trigger changed")
    require(summary.get("admitted_relations") == summary.get("full_rank_at_relation") == 95, "Stage-70 rank endpoint changed")
    require(summary.get("target_trials") == 97 and summary.get("support_queries") == 50_656_466, "Stage-70 collection counters changed")
    require(run.get("recovered_scalar") == 7_892_094_459_170, "Stage-70 recovered scalar changed")
    require(run["direct_process"]["metrics"]["wall_seconds"] == 6.571509542999991, "Stage-70 direct wall changed")
    require(run["comparison"]["direct_over_rho_wall_ratio"] == 1.4676110945821472, "Stage-70 rho ratio changed")
    return {
        "source_commit": UNKNOWN_CONFIG["commit"],
        "archive": common.archive_identity(70, UNKNOWN_CONFIG),
        "workflow": common.workflow_accounting(70, UNKNOWN_CONFIG),
        "artifact": common.artifact_accounting(70, UNKNOWN_CONFIG),
        "build_seal_sha256": build_seal,
        "run_seal_sha256": run_seal,
        "verification": verification,
        "comparison": run["comparison"],
        "direct": summary,
        "rho": run["rho"],
        "processes": {"direct": run["direct_process"]["metrics"], "rho": run["rho_process"]["metrics"]},
        "resources": {
            "build_outer": run["build_binding"]["outer_resources"],
            "science_outer": run["outer_resources"],
            "charged_total_core_seconds_available": run["charged_total_core_seconds_available"],
            "charged_sequential_wall_seconds_available": run["charged_sequential_wall_seconds_available"],
            "maximum_sampled_process_tree_rss_bytes": run["maximum_sampled_process_tree_rss_bytes"],
        },
    }


def single_core_result() -> dict[str, Any]:
    with tempfile.TemporaryDirectory() as directory:
        root = common.extract(71, SINGLE_CONFIG, Path(directory))
        verification = stage71.verify(root / "stage71-build", root / "stage71-run")
        embedded = load(root / "stage71-verification.json", "embedded Stage-71 verification")
        require(verification == embedded, "Stage-71 embedded verification changed")
        run = load(root / "stage71-run/result.json", "Stage-71 run")
        build_seal = custody.sha256_file(root / "stage71-build/result-seal.json", "Stage-71 build seal")
        run_seal = custody.sha256_file(root / "stage71-run/result-seal.json", "Stage-71 run seal")
    affinity = run["affinity"]
    summary = run["direct"]["summary"]
    require(affinity.get("inheritance_verified") is True, "Stage-71 affinity inheritance changed")
    require(affinity.get("parent_effective_cpus") == affinity.get("child_effective_cpus") == [0], "Stage-71 effective CPU changed")
    require(summary.get("query_parallel_threads") == 1, "Stage-71 thread count changed")
    require(summary.get("rank_aware_pair_scan") is True and summary.get("rank_target_deficiency") == 1, "Stage-71 rank-aware controls changed")
    require(run["direct_process"]["metrics"]["wall_seconds"] == 10.945706588999997, "Stage-71 direct wall changed")
    require(run["direct_process"]["metrics"]["total_core_seconds"] == 8.576129, "Stage-71 direct CPU changed")
    require(run["comparison"]["direct_over_rho_wall_ratio"] == 1.993515234548992, "Stage-71 wall ratio changed")
    return {
        "source_commit": SINGLE_CONFIG["commit"],
        "archive": common.archive_identity(71, SINGLE_CONFIG),
        "workflow": common.workflow_accounting(71, SINGLE_CONFIG),
        "artifact": common.artifact_accounting(71, SINGLE_CONFIG),
        "build_seal_sha256": build_seal,
        "run_seal_sha256": run_seal,
        "verification": verification,
        "affinity": affinity,
        "comparison": run["comparison"],
        "direct": summary,
        "rho": run["rho"],
        "processes": {"direct": run["direct_process"]["metrics"], "rho": run["rho_process"]["metrics"]},
        "resources": {
            "build_outer": run["build_binding"]["outer_resources"],
            "science_outer": run["outer_resources"],
            "charged_total_core_seconds_available": run["charged_total_core_seconds_available"],
            "charged_sequential_wall_seconds_available": run["charged_sequential_wall_seconds_available"],
            "maximum_sampled_process_tree_rss_bytes": run["maximum_sampled_process_tree_rss_bytes"],
        },
    }


def result() -> dict[str, Any]:
    predecessor = load(PREDECESSOR, "Stage-68 predecessor")
    require(predecessor.get("status") == "hosted_n53_fused_result_verified", "Stage-68 predecessor changed")
    unknown = unknown_result()
    single = single_core_result()
    return {
        "schema": SCHEMA,
        "status": "current_unknown_and_single_core_results_verified",
        "predecessor_stage68": {"path": str(PREDECESSOR.relative_to(REPO)), "sha256": custody.sha256_file(PREDECESSOR, "Stage-68 result")},
        "unknown_scalar": unknown,
        "single_core": single,
        "factor_base_algebraically_point_defined": True,
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
        "unknown_scalar_recovery_verified": True,
        "single_core_affinity_verified": True,
        "sat_conflicts": None,
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "remaining_blockers": [
            "known-target parallel direct remains 1.182872 times rho wall",
            "optimized unknown-scalar direct remains 1.467611 times rho wall",
            "forced-single-core direct remains 1.993515 times rho wall",
            "all fresh-build comparisons remain far above rho",
            "licensed same-instance Magma and unaffiliated reproduction/novelty review remain absent",
        ],
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-72 output must be new")
    output.mkdir()
    custody.write_json_new(output / "verification.json", result())
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "hosted_results_frozen", "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "result-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "verification.json", "committed Stage-72 result")
    require(committed == result(), "committed Stage-72 result differs from source evidence")
    seal = load(output / "result-seal.json", "Stage-72 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-72 result seal changed")
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-72 result inventory changed")
    return committed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    create = sub.add_parser("compose")
    create.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    check = sub.add_parser("verify")
    check.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    try:
        value = freeze(args.output.resolve()) if args.command == "compose" else verify(args.output.resolve())
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, VerificationError, common.VerificationError, stage57.Stage57Error, stage71.Stage71Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage72-results: {error}")


if __name__ == "__main__":
    main()

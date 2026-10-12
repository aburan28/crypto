#!/usr/bin/env python3
"""Compose and verify the hosted Stage 68 n=53 fused-pipeline result."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage67_fused_pcl as stage67
import verify_koblitz_stage49_width_results as common


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
PREDECESSOR = STAGE / "stage-59-current-gate-audit-20260912/audit.json"
UNKNOWN = STAGE / "stage-58-n53-unknown-result-20260912/verification.json"
DEFAULT_OUTPUT = STAGE / "stage-68-n53-fused-result-20260912"
SCHEMA = "koblitz_stage68_fused_result.v1"
SEAL_SCHEMA = "koblitz_stage68_fused_result_seal.v1"
CONFIG = {
    "run_id": 34718443480,
    "commit": "14ce3e1a9bbea5f3682bbcb1553d67a1aaccee9a",
    "artifact_id": 10305982096,
    "artifact_digest": "sha256:6b00eff9bddb541f5de55dcca08669c2aca8de2e852df9eec03f09578689783d",
    "archive": "koblitz-stage68-n53-fused-evidence-34718443480.tar.gz",
    "archive_root": "koblitz-stage68-n53-fused-34718443480",
    "prefix": "stage67",
    "job": "n53-fused-pcl",
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


def compact(observation: dict[str, Any]) -> dict[str, Any]:
    base = observation["base"]
    summary = observation["summary"]
    return {
        "base_hash": base["base_hash"],
        "support_index_entries": base["support_index_entries"],
        "support_table_allocated_bytes": base["support_table_allocated_bytes"],
        "parallel_support_expansion": base["parallel_support_expansion"],
        "pipelined_support_expansion": base["pipelined_support_expansion"],
        "field_mul_backend": summary["field_mul_backend"],
        "field_square_backend": summary["field_square_backend"],
        "field_reduction_backend": summary["field_reduction_backend"],
        "field_product_pipeline": summary["field_product_pipeline"],
        "factor_base_points": summary["factor_base_points"],
        "orbit_columns": summary["orbit_columns"],
        "matrix_columns": summary["matrix_columns"],
        "admitted_relations": summary["admitted_relations"],
        "target_trials": summary["target_trials"],
        "full_rank_at_relation": summary["full_rank_at_relation"],
        "setup_ms": summary["setup_ms"],
        "collection_ms": summary["collection_ms"],
        "charged_total_ms": summary["charged_total_ms"],
        "support_queries": summary["support_queries"],
        "query_batch_inversions": summary["query_batch_inversions"],
        "query_exact_table_misses": summary["query_exact_table_misses"],
        "rank_target_slot_index_allocated_bytes": summary["rank_target_slot_index_allocated_bytes"],
        "published_q": summary["published_q"],
        "recovered_fixture_scalar": summary["recovered_fixture_scalar"],
        "linear_solution_verified": summary["linear_solution_verified"],
        "all_relations_group_verified": summary["all_relations_group_verified"],
    }


def result() -> dict[str, Any]:
    predecessor = load(PREDECESSOR, "Stage-59 predecessor")
    unknown = load(UNKNOWN, "Stage-58 unknown-scalar result")
    require(predecessor.get("schema") == "koblitz_stage59_current_gate_audit.v1", "Stage-59 predecessor changed")
    require(unknown.get("unknown_scalar_recovery_verified") is True, "Stage-58 unknown-scalar evidence changed")
    with tempfile.TemporaryDirectory() as directory:
        root = common.extract(68, CONFIG, Path(directory))
        verification = stage67.verify(root / "stage67-build", root / "stage67-run")
        embedded = load(root / "stage67-verification.json", "embedded Stage-68 verification")
        require(verification == embedded, "Stage-68 embedded verification changed")
        run = load(root / "stage67-run/result.json", "Stage-68 run")
        build_seal = custody.sha256_file(root / "stage67-build/result-seal.json", "Stage-68 build seal")
        run_seal = custody.sha256_file(root / "stage67-run/result-seal.json", "Stage-68 run seal")
    candidate = compact(run["candidate"])
    baseline = compact(run["baseline"])
    comparison = run["comparison"]
    require(verification.get("status") == "n53_fused_pcl_verified", "Stage-68 verification status changed")
    require(candidate["admitted_relations"] == candidate["full_rank_at_relation"] == 95, "Stage-68 rank endpoint changed")
    require(candidate["target_trials"] == 96 and candidate["support_queries"] == 48_531_878, "Stage-68 collection counters changed")
    require(candidate["support_index_entries"] == 24_805_379, "Stage-68 support entry count changed")
    require(candidate["support_table_allocated_bytes"] == 738_197_504, "Stage-68 support allocation changed")
    require(candidate["rank_target_slot_index_allocated_bytes"] == 2_111_100, "Stage-68 rank index allocation changed")
    require(candidate["field_product_pipeline"] == "x86_64_pclmul_n53_fused_reduce", "Stage-68 fused pipeline changed")
    require(baseline["field_product_pipeline"] == "separate_product_and_reduction", "Stage-68 baseline pipeline changed")
    require(comparison["fused_direct_wall_speedup"] == 1.0973798162569024, "Stage-68 wall speedup changed")
    require(comparison["fused_over_rho_wall_ratio"] == 1.1828715401273624, "Stage-68 rho ratio changed")
    require(comparison["fresh_build_plus_fused_over_rho_wall_ratio"] == 19.474746292593565, "Stage-68 full-cost ratio changed")
    require(comparison["fused_whole_process_crossover"] is False, "Stage-68 crossover boundary changed")
    require(run["candidate_process"]["metrics"]["wall_seconds"] == 6.039836534999949, "Stage-68 candidate wall changed")
    require(run["candidate_process"]["metrics"]["total_core_seconds"] == 11.26981, "Stage-68 candidate CPU changed")
    require(run["candidate_process"]["metrics"]["peak_rss_bytes"] == 1_051_987_968, "Stage-68 candidate RSS changed")
    return {
        "schema": SCHEMA,
        "status": "hosted_n53_fused_result_verified",
        "source_commit": CONFIG["commit"],
        "archive": common.archive_identity(68, CONFIG),
        "workflow": common.workflow_accounting(68, CONFIG),
        "artifact": common.artifact_accounting(68, CONFIG),
        "build_seal_sha256": build_seal,
        "run_seal_sha256": run_seal,
        "verification": verification,
        "comparison": comparison,
        "baseline": baseline,
        "candidate": candidate,
        "rho": run["rho"],
        "processes": {
            "baseline": run["baseline_process"]["metrics"],
            "candidate": run["candidate_process"]["metrics"],
            "rho": run["rho_process"]["metrics"],
        },
        "resources": {
            "build_outer": run["build_binding"]["outer_resources"],
            "science_outer": run["outer_resources"],
            "charged_total_core_seconds_available": run["charged_total_core_seconds_available"],
            "charged_sequential_wall_seconds_available": run["charged_sequential_wall_seconds_available"],
            "maximum_sampled_process_tree_rss_bytes": run["maximum_sampled_process_tree_rss_bytes"],
        },
        "predecessor_stage59": {"path": str(PREDECESSOR.relative_to(REPO)), "sha256": custody.sha256_file(PREDECESSOR, "Stage-59 audit")},
        "unknown_scalar_stage58": {"path": str(UNKNOWN.relative_to(REPO)), "sha256": custody.sha256_file(UNKNOWN, "Stage-58 result")},
        "same_public_target": True,
        "exact_relation_hash_equality": True,
        "factor_base_algebraically_point_defined": True,
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
        "fixture_scalar_used_by_collector": False,
        "unknown_scalar_n53_predecessor_verified": True,
        "sat_conflicts": None,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "remaining_blockers": [
            "preferred n=53 direct remains 1.182872 times same-target rho wall",
            "fresh build plus preferred direct remains 19.474746 times rho wall",
            "licensed same-instance Magma and unaffiliated reproduction/novelty review remain absent",
            "latest forced-single-core and optimized unknown-scalar reruns remain absent",
            "finite constant-factor improvements do not establish an asymptotic improvement",
        ],
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-68 output must be new")
    output.mkdir()
    custody.write_json_new(output / "verification.json", result())
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "hosted_result_frozen", "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "result-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "verification.json", "committed Stage-68 result")
    require(committed == result(), "committed Stage-68 result differs from source evidence")
    seal = load(output / "result-seal.json", "Stage-68 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-68 result seal changed")
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-68 result inventory changed")
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
    except (OSError, ValueError, KeyError, VerificationError, common.VerificationError, stage67.Stage67Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage68-result: {error}")


if __name__ == "__main__":
    main()

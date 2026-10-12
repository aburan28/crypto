#!/usr/bin/env python3
"""Compose Stage 35 into the current seven-gate Koblitz audit."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

import compose_koblitz_stage34_gate_audit as stage34
import run_koblitz_blind_pdp_phase_b as custody
import verify_koblitz_stage35_result as stage35


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE34 = STAGE / "stage-34-current-gate-audit-20260911/audit.json"
STAGE35 = STAGE / "stage-35-algebraic-walk-result-20260912/verification.json"
DEFAULT_OUTPUT = STAGE / "stage-36-current-gate-audit-20260912"
SCHEMA = "koblitz_stage36_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage36_current_gate_audit_seal.v1"


class AuditError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def compose() -> dict[str, Any]:
    predecessor = stage34.verify()
    walked = stage35.verify()
    require(predecessor == load(STAGE34, "Stage-34 audit"), "Stage-34 binding changed")
    require(walked == load(STAGE35, "Stage-35 result"), "Stage-35 binding changed")
    comparison = walked["comparison"]
    require(
        comparison["stage35_online_crossover"] is True
        and comparison["stage35_amortised_crossover"] is False
        and comparison["stage35_whole_process_crossover"] is False,
        "Stage-35 crossover boundary changed",
    )
    result = copy.deepcopy(predecessor)
    result.update({
        "schema": SCHEMA,
        "status": "current_evidence_audited_online_crossover_only_gates_incomplete",
        "predecessor_stage34": {
            "path": str(STAGE34.relative_to(REPO)),
            "sha256": custody.sha256_file(STAGE34, "Stage-34 audit"),
        },
        "stage35_binding": {
            "path": str(STAGE35.relative_to(REPO)),
            "sha256": custody.sha256_file(STAGE35, "Stage-35 result"),
            "verification": walked,
        },
        "stage35_algebraic_walk": {
            "factor_base_algebraically_defined": walked["factor_base"]["algebraically_defined"],
            "target_subgroup_enumerated": walked["factor_base"]["target_subgroup_enumerated"],
            "target_or_log_labels_used_for_selection": walked["factor_base"]["target_or_log_labels_used_for_selection"],
            "collection_summands": walked["collection_summands"],
            "descent_summands": walked["descent_summands"],
            "targets_verified": walked["targets_verified"],
            "descent_trials": walked["descent_trials"],
            "stage33_online_ic_over_rho_wall_ratio": comparison["stage33_online_ic_over_rho_wall_ratio"],
            "stage35_online_ic_over_rho_wall_ratio": comparison["stage35_online_ic_over_rho_wall_ratio"],
            "online_ratio_improvement_factor": comparison["online_ratio_improvement_factor"],
            "online_crossover": True,
            "amortised_ic_over_rho_wall_ratio": comparison["stage35_amortised_ic_over_rho_wall_ratio"],
            "amortised_crossover": False,
            "whole_process_crossover": False,
            "full_available_wall_over_same_targets_rho": comparison["stage35_full_available_wall_over_same_targets_rho"],
            "claim_boundary": comparison["claim_boundary"],
        },
        "all_seven_gates_passed": False,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "overall": "finite post-precomputation algebraic online crossover; full cost, Magma, external review, and SOTA gates remain false",
    })
    gate = result["gates"]["6_automorphism_rho"]
    gate["status"] = "finite_n23_n31_n41_complete_stage35_online_crossover_full_cost_no_crossover"
    gate["proved"].extend([
        f"Stage-35 algebraic walked descent is {1.0 / comparison['stage35_online_ic_over_rho_wall_ratio']:.12f} times faster than rho after precomputation",
        f"Stage-35 online IC/rho ratio improved {comparison['online_ratio_improvement_factor']:.12f} times over Stage 33",
        f"Stage-35 amortised IC remains {comparison['stage35_amortised_ic_over_rho_wall_ratio']:.12f} times slower than rho",
        f"Stage-35 fresh-build-plus-science wall remains {comparison['stage35_full_available_wall_over_same_targets_rho']:.12f} times the same-target rho wall",
    ])
    gate["limitations"] = [
        "online crossover assumes the factor-base log database already exists",
        "fresh build, discovery, relations, and logs do not cross for five targets",
        "finite 40-bit public synthetic subgroup only",
        "the optimization changes constants, not the asymptotic exponent",
    ]
    return result


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-36 output must be new")
    output.mkdir(parents=False)
    custody.write_json_new(output / "audit.json", compose())
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "audit_frozen", "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "audit-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "audit.json", "committed Stage-36 audit")
    require(committed == compose(), "committed Stage-36 audit differs from source evidence")
    seal = load(output / "audit-seal.json", "Stage-36 audit seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-36 seal is invalid")
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    require(inventory == seal.get("inventory") and custody.canonical_sha256(inventory) == seal.get("inventory_sha256"), "Stage-36 inventory changed")
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
    except (OSError, ValueError, KeyError, AuditError, custody.PhaseBError) as error:
        raise SystemExit(f"stage36-audit: {error}")


if __name__ == "__main__":
    main()

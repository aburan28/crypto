#!/usr/bin/env python3
"""Compose the current seven-gate Koblitz audit without widening its claims."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile
from typing import Any

import compose_koblitz_stage26_stage32 as phase_b
import run_koblitz_blind_pdp_phase_b as custody
import verify_koblitz_stage33_result as stage33


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
PHASE_B_SCORE = STAGE / "stage-26-affinity-matrix-stage32-corrected-result-20260911/score.json"
STAGE33_SCORE = STAGE / "stage-33-n41-public-result-20260911/verification.json"
DEFAULT_OUTPUT = STAGE / "stage-34-current-gate-audit-20260911"
SCHEMA = "koblitz_stage34_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage34_current_gate_audit_seal.v1"


class AuditError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def stage33_baseline() -> dict[str, Any]:
    with tempfile.TemporaryDirectory() as directory:
        root = stage33.extract(Path(directory))
        result = load(root / "stage33-run/result.json", "Stage-33 run result")
        workflow = result.get("workflow_result", {})
        stages = workflow.get("stages")
        require(isinstance(stages, list), "Stage-33 stages are missing")
        rows = [row.get("vs_rho") for row in stages if row.get("stage") == "baseline"]
        require(len(rows) == 1 and isinstance(rows[0], dict), "Stage-33 baseline is missing")
        return rows[0]


def compose() -> dict[str, Any]:
    phase_b_verification = phase_b.verify()
    stage33_verification = stage33.verify()
    phase_b_score = load(PHASE_B_SCORE, "corrected Phase-B score")
    stage33_score = load(STAGE33_SCORE, "Stage-33 hosted score")
    require(phase_b_score.get("koblitz_index_calculus_sota") is False, "Phase-B claim widened")
    require(stage33_score == stage33_verification, "Stage-33 score binding changed")
    baseline = stage33_baseline()
    rho_wall = baseline.get("rho", {}).get("seconds_total")
    require(isinstance(rho_wall, (int, float)) and rho_wall > 0, "Stage-33 rho wall is invalid")
    charged_wall = stage33_score["resources"]["charged_sequential_wall_seconds_available"]
    scientific_outer_wall = stage33_score["resources"]["scientific_single_core_elapsed_seconds"]
    full_path_ratio = charged_wall / rho_wall
    scientific_outer_ratio = scientific_outer_wall / rho_wall
    require(full_path_ratio > stage33_score["ic_over_rho_amortised_wall_ratio"] > 1.0, "Stage-33 full-cost ratio ordering changed")
    degree23 = phase_b_score["unknown_scalar_index_calculus_and_rho"]
    degree31 = phase_b_score["n31_unknown_scalar_index_calculus_and_rho"]
    gates = {
        "1_full_cost_accounting": {
            "status": "partial",
            "proved": [
                "Phase-B solver acquisition, builds, one-CPU cells, direct MITM, correction, conflicts, wall, and process-tree memory are charged",
                "Stage-33 fresh dependency acquisition, four-job build, algebraic discovery, predicate materialization, relations, sparse linear algebra, descent, rho, one-CPU time, and process-tree memory are charged",
            ],
            "missing": [
                "licensed Magma execution resources",
                "preinstalled operating-system and Rust-toolchain acquisition",
            ],
        },
        "2_same_instance_backend_matrix": {
            "status": "partial",
            "proved": [
                "native XOR SAT, WDSat, CryptoMiniSat, and direct MITM on the exact 160-input packet",
                "standard n=31/n=41 and GGMP n=31 cells",
            ],
            "missing": ["licensed Magma F4 on all 160 exact inputs"],
        },
        "3_required_resource_fields": {
            "status": "partial",
            "proved": [
                "single-core elapsed, total core-seconds, peak process and process-tree memory, conflicts or operation counters, and workflow wall are retained for executed arms"
            ],
            "missing": ["the same fields for the licensed Magma arm"],
        },
        "4_scaling": {
            "status": "finite_requested_regimes_complete",
            "proved": [
                "PDP matrix at n=31, n=41, and n=59",
                "end-to-end public unknown-scalar controls at n=31 and n=41",
            ],
            "limitations": ["n=59 remains PDP-only", "all subgroups are toy-sized"],
        },
        "5_unknown_scalar": {
            "status": "finite_n23_n31_n41_complete",
            "proved": [
                f"{degree23['completed_unknown_scalar_targets']} degree-23 public targets",
                f"{degree31['completed_unknown_scalar_targets']} degree-31 public hash targets",
                f"{stage33_score['targets_verified']} degree-41 public hash targets",
                "factor-base logs were not known by construction and every recovery was group-certified",
            ],
            "limitations": ["finite public synthetic controls only"],
        },
        "6_automorphism_rho": {
            "status": "finite_n23_n31_n41_complete_no_crossover",
            "proved": [
                "same-target signed-Frobenius or signed-Frobenius/negation rho at all three degrees",
                f"Stage-33 online descent is {stage33_score['ic_over_rho_charged_wall_ratio']:.12f} times slower than rho",
                f"Stage-33 discovery/relations/logs amortized path is {stage33_score['ic_over_rho_amortised_wall_ratio']:.12f} times slower than rho",
                f"Stage-33 complete available build-plus-science wall is {full_path_ratio:.12f} times the rho wall for the same five targets",
            ],
            "limitations": ["preinstalled system/toolchain acquisition is excluded", "no asymptotic inference"],
        },
        "7_external_reproduction_and_novelty": {
            "status": "missing",
            "proved": ["public issue 97 contains current source pins, artifacts, commands, and verdict format"],
            "missing": ["unaffiliated reproduction", "source-pinned external novelty verdict"],
        },
    }
    return {
        "schema": SCHEMA,
        "status": "current_evidence_audited_gates_incomplete",
        "phase_b_binding": {
            "score_sha256": custody.sha256_file(PHASE_B_SCORE, "Phase-B score"),
            "verification": phase_b_verification,
            "instances": phase_b_score["instances"],
            "backend_rows": phase_b_score["backend_rows"],
            "classification_counts": phase_b_score["classification_counts"],
            "charged_total_core_seconds_available": phase_b_score["charged_total_core_seconds_available"],
        },
        "stage33_binding": {
            "score_sha256": custody.sha256_file(STAGE33_SCORE, "Stage-33 score"),
            "verification": stage33_verification,
        },
        "stage33_full_cost_vs_rho": {
            "rho_same_five_targets_wall_seconds": rho_wall,
            "scientific_outer_wall_seconds": scientific_outer_wall,
            "scientific_outer_over_rho_wall_ratio": scientific_outer_ratio,
            "charged_build_plus_science_wall_seconds": charged_wall,
            "charged_build_plus_science_over_rho_wall_ratio": full_path_ratio,
            "online_descent_over_rho_wall_ratio": stage33_score["ic_over_rho_charged_wall_ratio"],
            "amortised_precompute_and_descent_over_rho_wall_ratio": stage33_score["ic_over_rho_amortised_wall_ratio"],
            "automorphism_discount": baseline["automorphism_discount"],
            "verdict": "rho_faster_under_every_stage33_cost_scope",
        },
        "gates": gates,
        "gates_passed": ["4_scaling", "5_unknown_scalar", "6_automorphism_rho"],
        "gates_partial": ["1_full_cost_accounting", "2_same_instance_backend_matrix", "3_required_resource_fields"],
        "gates_missing": ["7_external_reproduction_and_novelty"],
        "all_seven_gates_passed": False,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "overall": "strong internal engineering and finite public toy-research improvement; not a Koblitz index-calculus SOTA",
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-34 output must be new")
    output.mkdir(parents=False)
    custody.write_json_new(output / "audit.json", compose())
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    payload = {
        "schema": SEAL_SCHEMA,
        "status": "audit_frozen",
        "inventory": inventory,
        "inventory_sha256": custody.canonical_sha256(inventory),
    }
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "audit-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "audit.json", "committed Stage-34 audit")
    require(committed == compose(), "committed Stage-34 audit differs from current evidence")
    seal = load(output / "audit-seal.json", "Stage-34 audit seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-34 audit seal is invalid")
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    require(inventory == seal.get("inventory") and custody.canonical_sha256(inventory) == seal.get("inventory_sha256"), "Stage-34 audit inventory changed")
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
        raise SystemExit(f"stage34-audit: {error}")


if __name__ == "__main__":
    main()

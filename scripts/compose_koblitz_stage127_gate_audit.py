#!/usr/bin/env python3
"""Compose the current seven-gate audit over the frozen Phase-B matrix and
current scalar-blind n=53 full-cost evidence.

The two experiment families remain separate: the n31/n41/n59 packet is the
common Semaev solver matrix, while n53 is the optimized subgroup-orbit
MITM/rho workflow.  Stage 127 must never call them the same instance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any

REPO = Path(__file__).resolve().parents[1]
EVIDENCE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE126 = EVIDENCE / "stage-126-current-gate-audit-20260921"
PHASE_B = EVIDENCE / "stage-26-affinity-matrix-stage32-corrected-result-20260911"
N53 = REPO / "docs/ic/runs/koblitz-n53-one-unit-20260921.json"
N53_PARAMS = REPO / "docs/ic/params/k0n53-subgroup-one-unit-public-unknown.json"
GATE_STATUS = EVIDENCE / "GATE_STATUS.md"
SCHEMA = "koblitz_stage127_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage127_current_gate_audit_seal.v1"
GATE_STATUS_MARKER = "Current through Stage 127"

class Stage127Error(RuntimeError):
    pass

def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage127Error(message)

def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

def replay_stage126() -> dict[str, Any]:
    done = subprocess.run(
        [sys.executable, str(REPO / "scripts/compose_koblitz_stage126_gate_audit.py"),
         "verify", "--output", str(STAGE126)],
        cwd=REPO, text=True, capture_output=True, check=True,
    )
    value = json.loads(done.stdout)
    require(value.get("schema") == "koblitz_stage126_current_gate_audit.v1", "Stage-126 replay changed")
    return value

def phase_b_score() -> dict[str, Any]:
    seal = load(PHASE_B / "score-seal.json", "Phase-B score seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == "koblitz_stage26_stage32_composed_score_seal.v1", "Phase-B seal schema changed")
    require(claimed == canonical_sha256(payload), "Phase-B seal self-hash changed")
    inventory = seal.get("inventory")
    require(isinstance(inventory, list) and len(inventory) == 1, "Phase-B inventory changed")
    item = inventory[0]
    require(item.get("path") == "score.json", "Phase-B score path changed")
    require(item.get("sha256") == sha256(PHASE_B / "score.json"), "Phase-B score hash changed")
    require(item.get("bytes") == (PHASE_B / "score.json").stat().st_size, "Phase-B score size changed")
    require(seal.get("inventory_sha256") == canonical_sha256(inventory), "Phase-B inventory hash changed")
    score = load(PHASE_B / "score.json", "Phase-B score")
    require(score.get("status") == "complete_verified_stage32_corrected_four_cell_panel", "Phase-B status changed")
    require(score.get("factor_base_algebraic_without_target_subgroup_enumeration") is True, "Phase-B factor-base boundary changed")
    require(score.get("factor_base_logs_known_by_construction") is False, "Phase-B gained log labels")
    require(score.get("licensed_magma_f4_same_instance_panel_complete") is False, "licensed Magma unexpectedly marked complete")
    require(score.get("independent_external_reproduction_satisfied") is False, "external gate unexpectedly marked complete")
    audit = score["completion_gate_audit"]
    require(audit["2_same_instance_backend_matrix"]["status"] == "partial", "Phase-B matrix gate changed")
    require("native XOR SAT, WDSat, CryptoMiniSat, and direct MITM on the same 160 inputs" in audit["2_same_instance_backend_matrix"]["proved"], "Phase-B common matrix proof missing")
    require(score["matched_direct_mitm"]["instances"] == 160, "MITM instance count changed")
    return score

def compose() -> dict[str, Any]:
    predecessor = replay_stage126()
    matrix = phase_b_score()
    n53 = load(N53, "n53 evidence")
    params = load(N53_PARAMS, "n53 parameters")
    require(n53.get("schema_version") == 2, "n53 evidence schema changed")
    require(params.get("name") == "k0n53-subgroup-seed6-one-unit-17k-public-unknown", "n53 parameter identity changed")
    require(params["targets"] == [{"public_hash_seed": 53001}], "n53 public target changed")
    require(n53["instance"]["factor_base_logs_known_by_construction"] is False, "n53 gained factor-base labels")
    require(n53["instance"]["factor_base_points"] == 36464, "n53 factor-base width changed")
    require(n53["instance"]["signed_orbit_columns"] == 344, "n53 column count changed")
    fixed = n53["single_target_fixed_width_key_chunks"]
    one = fixed["matched_one_core_ab"]["medians"]
    default = fixed["matched_default_thread_ab"]["medians"]
    orbit = n53["single_target_orbit_derived_projection"]
    campaign = n53["science_process_accounting"]
    require(fixed["matched_one_core_ab"]["canonical_relations_sha256"] == "d8e18605b4a27f307101ecbb5ba9973b50fc01ba272d6f39d216c4fd6a615d37", "n53 relation hash changed")
    require(fixed["matched_one_core_ab"]["recovered_scalar_each"] == 7892094459170, "n53 recovered scalar changed")
    require(one["ic_wall"]["candidate"] > one["rho_wall"]["candidate"], "n53 one-core rho boundary changed")
    require(default["ic_wall"]["candidate"] < default["rho_wall"]["candidate"], "n53 online wall crossover changed")
    require(campaign["processes"] == 318, "n53 campaign process count changed")
    predicate = orbit["native_cost_current_head"]
    require(predicate["cofactor_multiplications"] == 344, "n53 predicate multiplication count changed")
    require(GATE_STATUS_MARKER in GATE_STATUS.read_text(), "gate-status marker changed")
    backends = matrix["backend_resources"]
    direct = matrix["matched_direct_mitm"]
    return {
        "schema": SCHEMA,
        "status": "current_seven_gate_audit_verified",
        "all_seven_gates_passed": False,
        "koblitz_index_calculus_sota": False,
        "claim_boundary": "strong internal engineering and finite public toy-research improvement; not a Koblitz index-calculus SOTA",
        "predecessor": {
            "stage126_audit_sha256": sha256(STAGE126 / "audit.json"),
            "stage126_seal_sha256": sha256(STAGE126 / "result-seal.json"),
            "stage126_status": predecessor["status"],
        },
        "evidence_pins": {
            "phase_b_score_sha256": sha256(PHASE_B / "score.json"),
            "phase_b_seal_sha256": sha256(PHASE_B / "score-seal.json"),
            "n53_evidence_sha256": sha256(N53),
            "n53_params_sha256": sha256(N53_PARAMS),
            "gate_status_sha256": sha256(GATE_STATUS),
        },
        "instance_separation": {
            "common_solver_matrix": "frozen n31-l5-m3 standard/GGMP, n41-l5-m3 standard, n59-l9-m3 standard Phase-B packet",
            "optimized_full_dlp": "n53 subgroup-orbit factor base with 344 projected columns and folded direct MITM",
            "same_instance_claim_across_these_families": False,
            "reason": "WDSat, CryptoMiniSat and Magma Semaev encodings require a named invariant subspace ell; the optimized n53 subgroup-orbit base has no equivalent ell encoding",
        },
        "phase_b_same_instance_matrix": {
            "instances": matrix["instances"],
            "cells": sorted(matrix["cells"]),
            "backends": {
                name: {
                    "total_core_seconds": row["total_core_seconds"],
                    "summed_process_wall_seconds": row["summed_process_wall_seconds"],
                    "peak_rss_bytes": row["peak_rss_bytes"],
                    "conflicts_sum": row["conflicts_sum"],
                    "conflicts_reported": row["conflicts_reported"],
                }
                for name, row in backends.items()
            },
            "direct_mitm": {
                "instances": direct["instances"],
                "true_positive": direct["classification_counts"]["true_positive"],
                "true_negative": direct["classification_counts"]["true_negative"],
                "group_additions": direct["operations"]["group_additions"],
                "pair_entries": direct["operations"]["pair_entries"],
                "summed_outer_total_core_seconds": direct["resources"]["summed_outer_total_core_seconds"],
                "summed_single_core_elapsed_seconds": direct["resources"]["summed_single_core_elapsed_seconds"],
                "maximum_sampled_process_tree_rss_bytes": direct["resources"]["maximum_sampled_process_tree_rss_bytes"],
                "workflow_wall_seconds": direct["resources"]["workflow_wall_seconds"],
                "conflicts": None,
            },
            "licensed_magma_complete": False,
        },
        "current_n53": {
            "factor_base": n53["instance"],
            "parameter_file": str(N53_PARAMS.relative_to(REPO)),
            "relation_hash": fixed["matched_one_core_ab"]["canonical_relations_sha256"],
            "recovered_scalar": fixed["matched_one_core_ab"]["recovered_scalar_each"],
            "one_core": {
                "ic_wall_seconds": one["ic_wall"]["candidate"],
                "rho_wall_seconds": one["rho_wall"]["candidate"],
                "ic_over_rho_wall_ratio": one["ic_wall"]["candidate"] / one["rho_wall"]["candidate"],
                "whole_process_cpu_seconds": one["whole_cpu"]["candidate"],
                "peak_rss_bytes": one["rss"]["candidate"],
            },
            "default_threads": {
                "ic_wall_seconds": default["ic_wall"]["candidate"],
                "rho_wall_seconds": default["rho_wall"]["candidate"],
                "rho_over_ic_wall_ratio": default["rho_over_ic"]["candidate"],
                "whole_process_cpu_seconds": default["whole_cpu"]["candidate"],
                "peak_rss_bytes": default["rss"]["candidate"],
            },
            "projected_predicate_native_cost": predicate,
            "campaign": campaign,
            "factor_base_logs_known_by_construction": False,
            "target_scalar_constructed_or_supplied": False,
        },
        "gates": {
            "1_all_stage_resource_charging": "partial_missing_licensed_magma_resources",
            "2_same_instance_solver_panel": "partial_phase_b_native_wdsat_cms_mitm_ggmp_complete_magma_missing",
            "3_single_core_core_memory_conflicts_wall": "partial_current_n53_refreshed_magma_missing",
            "4_n31_n41_larger_pdp_scaling": "satisfied_finite_execution_coverage",
            "5_unknown_scalar_without_constructed_base_logs": "satisfied_finite_n23_n31_n41_n53",
            "6_full_cost_vs_automorphism_rho": "partial_default_thread_wall_pass_single_core_and_core_cost_fail",
            "7_external_reproduction_novelty": "missing",
        },
        "next_targets": [
            "execute the frozen 160-input packet under licensed Magma F4 with one-core, total-core, wall, RSS and terminal receipts",
            "run one unified end-to-end algebraic-base IC/rho cost series at n41 and the larger PDP regime",
            "reduce the current n53 one-core IC/rho wall ratio below one without losing the default-thread wall crossover",
            "obtain unaffiliated reproduction and a source-pinned novelty/correctness review",
        ],
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
    }

def write_new(path: Path, value: dict[str, Any]) -> None:
    require(not path.exists(), f"refusing to overwrite {path}")
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")

def build(output: Path) -> dict[str, Any]:
    require(not output.exists(), f"refusing to overwrite {output}")
    output.mkdir(parents=True)
    audit = compose()
    write_new(output / "audit.json", audit)
    write_new(output / "result-seal.json", {
        "schema": SEAL_SCHEMA,
        "status": "audit_frozen",
        "audit_sha256": sha256(output / "audit.json"),
    })
    return audit

def verify(output: Path) -> dict[str, Any]:
    seal = load(output / "result-seal.json", "Stage-127 seal")
    require(seal.get("schema") == SEAL_SCHEMA, "seal schema changed")
    require(sha256(output / "audit.json") == seal.get("audit_sha256"), "audit seal changed")
    current = compose()
    require(current == load(output / "audit.json", "Stage-127 audit"), "current audit changed")
    return current

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build"); b.add_argument("--output", type=Path, required=True)
    v = sub.add_parser("verify"); v.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        value = build(args.output.resolve()) if args.command == "build" else verify(args.output.resolve(strict=True))
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage127Error) as error:
        raise SystemExit(f"stage127-gate-audit: {error}")

if __name__ == "__main__":
    main()

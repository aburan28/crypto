#!/usr/bin/env python3
"""Compose the current seven-gate audit through the Stage 108 selection,
re-sealed over the boundary ledger as updated on 2026-09-21.

Stage 109 sealed this audit on 2026-09-13 by pinning the mutable gate
documents (the boundary ledger and its Markdown twin, the gate status, the
autolab protocol) by hash.  The index-calculus boundary ledger of 2026-09-21
(research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md) then moved
the ledger and its twin -- operation-counted rows for the prime and binary
regimes, oracle pricing and the whole-process count on the Koblitz rows --
without moving any Koblitz gate fact.  Stage 124 recomposes the same seven
gates over the current documents and chains to the Stage 109 seal; Stage 109
itself is verified as an immutable historical snapshot by its own script.

The ledger moved again later on 2026-09-21 (the boundary ledger's Round 2),
so Stage 125 re-seals the audit and this script's verify is historical:
seal and frozen audit only.  `build` is kept for the record of how the
Stage 124 audit was composed and will not reproduce it on the moved
documents.
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
STAGE99 = EVIDENCE / "stage-99-optimization-chain-20260913"
STAGE108 = EVIDENCE / "stage-108-routing-selection-archive-20260913"
STAGE109 = EVIDENCE / "stage-109-current-gate-audit-20260913"
LEDGER = REPO / "docs/ic/boundary_targets.json"
SCOREBOARD = REPO / "docs/ic/BOUNDARY_TARGETS.md"
GATE_STATUS = EVIDENCE / "GATE_STATUS.md"
AUTOLAB_PROTOCOL = REPO / "research/sat_factor_base_review_20260908/autolab/protocol.json"

SCHEMA = "koblitz_stage124_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage124_current_gate_audit_seal.v1"
LEDGER_UPDATED = "2026-09-21"
GATE_STATUS_MARKER = "Current through Stage 124"


class Stage124Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage124Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay(script: str, output: Path) -> dict[str, Any]:
    completed = subprocess.run(
        [sys.executable, str(REPO / "scripts" / script), "verify", "--output", str(output)],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=True,
    )
    value = json.loads(completed.stdout)
    require(isinstance(value, dict), f"{script} replay was not an object")
    return value


def compose() -> dict[str, Any]:
    historical = replay("verify_koblitz_stage99_optimization_chain.py", STAGE99)
    selected = replay("verify_koblitz_stage108_routing_archive.py", STAGE108)
    # Stage 109 verifies as a historical snapshot: its own script checks the
    # seal and the frozen audit, never the documents that have moved since.
    stage109 = replay("compose_koblitz_stage109_gate_audit.py", STAGE109)
    ledger = load(LEDGER, "boundary ledger")
    require(
        ledger.get("schema_version") == 2 and ledger.get("updated") == LEDGER_UPDATED,
        "ledger version/date changed",
    )
    koblitz = ledger["regimes"]["koblitz"]["records"]
    factor = koblitz["factor_base"]["current"]["metrics"]
    vs_rho = koblitz["vs_rho"]["current"]["metrics"]
    unknown = koblitz["end_to_end_dlp"]["current"]["metrics"]
    require(factor["selection_uses_scalar_labels"] is False, "factor-base labels widened")
    require(factor["target_subgroup_enumerated"] is False, "target subgroup enumeration appeared")
    require(factor["support_table_shards"] == 4, "selected shard count changed")
    require(factor["support_table_shard_routing"] == "xor_low_and_high_x_windows", "selected routing changed")
    require(factor["retained_bytes"] == 738_197_504, "selected retained bytes changed")
    require(vs_rho["n53_selected_wins"] == 5, "selected panel wins changed")
    require(vs_rho["n53_selected_median_direct_over_rho_wall_ratio"] < 0.95, "online median gate changed")
    require(vs_rho["n53_current_direct_over_rho_core_ratio"] > 1.0, "core-cost boundary changed")
    require(vs_rho["n53_fresh_build_plus_current_over_rho_wall_ratio"] > 1.0, "build boundary changed")
    require(unknown["target_scalar_constructed_or_supplied"] is False, "unknown target gained scalar input")
    require(unknown["factor_base_logs_known_by_construction"] is False, "unknown run gained base logs")
    require(unknown["direct_over_rho_wall_ratio"] < 1.0, "unknown online crossover changed")
    require(selected["stage106"]["direct_route_wins"] == 5, "routing panel changed")
    require(selected["selected_panel"]["direct_wins"] == 5, "selected archive panel changed")
    # The gate facts Stage 109 sealed must be the ones this audit re-seals.
    sealed = stage109["selected_n53"]
    require(sealed["base_hash"] == factor["base_hash"], "selected base hash moved since Stage 109")
    require(sealed["support_entries"] == factor["support_index_entries"], "support entries moved since Stage 109")
    require(sealed["paired_wall_ratios"] == vs_rho["n53_selected_ratios"], "selected ratios moved since Stage 109")
    require(
        sealed["median_core_ratio"] == vs_rho["n53_current_direct_over_rho_core_ratio"]
        and sealed["fresh_build_ratio"] == vs_rho["n53_fresh_build_plus_current_over_rho_wall_ratio"]
        and sealed["unknown_scalar_ratio"] == unknown["direct_over_rho_wall_ratio"],
        "Koblitz gate ratios moved since Stage 109",
    )
    gate_text = GATE_STATUS.read_text()
    score_text = SCOREBOARD.read_text()
    require(GATE_STATUS_MARKER in gate_text, "gate-status version changed")
    require("17.636692 times rho" in gate_text, "gate-status full-cost boundary missing")
    require("0.8392x" in score_text and "below 512 MiB" in score_text, "scoreboard boundary missing")
    protocol = load(AUTOLAB_PROTOCOL, "autolab protocol")
    n53 = protocol["beats"]["koblitz.vs_rho.n53_probe"]
    env = n53.get("env", {})
    require(env.get("KIC_RANK_TARGET_DEFICIENCY") == "1", "autolab rank-target deficiency changed")
    require(env.get("KIC_RANK_AWARE_PAIR_SCAN") == "1", "autolab rank-aware scan changed")
    require(env.get("KIC_ENABLE_SHARDED_SUPPORT_TABLE") == "1", "autolab sharded support changed")
    require(n53.get("fixed_process_repetitions_required") == 5, "autolab repetition gate changed")
    selected_default = n53.get("selected_default", {})
    require(selected_default.get("support_table_shards") == 4, "autolab shard count changed")
    require(
        selected_default.get("support_table_shard_routing") == "xor_low_and_high_x_windows",
        "autolab shard routing changed",
    )
    require(
        selected_default.get("specialized_n53_pair_batch") is True
        and selected_default.get("specialized_n53_support_expansion") is False
        and selected_default.get("specialized_n53_pair_sums") is False,
        "autolab selected stack changed",
    )
    panel = n53.get("current_selected_panel", {})
    require(panel.get("direct_wins") == 5, "autolab selected panel wins changed")
    require(
        panel.get("median_ratio") == vs_rho["n53_selected_median_direct_over_rho_wall_ratio"],
        "autolab selected median changed",
    )
    require(
        panel.get("ratios") == vs_rho["n53_selected_ratios"],
        "autolab selected ratios changed",
    )
    require(
        n53.get("acceptance_gates")
        == [
            "algebraic factor base uses no scalar labels or target-subgroup enumeration",
            "same public target in every direct/rho pair",
            "at least five separately metered fixed-protocol process pairs",
            "all build/setup/relation/LA/wall/core/RSS resources charged",
            "measurement_schema.vs_rho required fields present",
            "median whole-process direct/rho ratio < 0.95",
            "not a ledger promotion until independent validation",
        ],
        "autolab acceptance gates changed",
    )
    return {
        "schema": SCHEMA,
        "status": "current_seven_gate_audit_verified",
        "predecessors": {
            "stage99_verification_sha256": sha256(STAGE99 / "verification.json"),
            "stage108_verification_sha256": sha256(STAGE108 / "verification.json"),
            "stage109_audit_sha256": sha256(STAGE109 / "audit.json"),
            "stage99_status": historical["status"],
            "stage108_status": selected["status"],
            "stage109_status": stage109["status"],
        },
        "since_stage109": {
            "what_moved": "docs/ic/boundary_targets.json, docs/ic/BOUNDARY_TARGETS.md and GATE_STATUS.md: operation-counted relation_yield / rank / end_to_end_dlp / vs_rho rows for the prime and binary regimes, oracle pricing and the whole-process operation count on the Koblitz decomposition and vs_rho rows, agent priority 8, and the research-notes move of 2026-09-19",
            "what_did_not": "every Koblitz gate fact Stage 109 sealed: the selected n=53 base hash, support entries and bytes, the five paired wall ratios, the core and fresh-build ratios, the unknown-scalar ratio, and the seven gate statuses",
            "reference": "research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md",
        },
        "gates": {
            "1_all_stage_resource_charging": "partial_missing_licensed_magma",
            "2_same_instance_solver_panel": "partial_missing_licensed_magma_full_panel",
            "3_single_core_core_memory_conflicts_wall": "partial_missing_magma_and_selected_single_core_refresh",
            "4_n31_n41_larger_pdp_scaling": "satisfied_finite_execution_coverage",
            "5_unknown_scalar_without_constructed_base_logs": "satisfied_finite_n23_n31_n41_n53",
            "6_full_cost_vs_automorphism_rho": "partial_online_wall_pass_full_cost_core_memory_fail",
            "7_external_reproduction_novelty": "missing",
        },
        "selected_n53": {
            "base_hash": factor["base_hash"],
            "support_entries": factor["support_index_entries"],
            "support_bytes": factor["retained_bytes"],
            "support_table_shards": factor["support_table_shards"],
            "routing": factor["support_table_shard_routing"],
            "paired_wall_ratios": vs_rho["n53_selected_ratios"],
            "median_wall_ratio": vs_rho["n53_selected_median_direct_over_rho_wall_ratio"],
            "median_core_ratio": vs_rho["n53_current_direct_over_rho_core_ratio"],
            "fresh_build_ratio": vs_rho["n53_fresh_build_plus_current_over_rho_wall_ratio"],
            "unknown_scalar_ratio": unknown["direct_over_rho_wall_ratio"],
        },
        "next_targets": [
            "reduce exact selected n=53 support below 512 MiB without changing the base or rank-95 solve",
            "reduce selected direct/rho core ratio to at most one and retain the online wall crossover",
            "refresh forced-single-core evidence for the selected four-shard route",
            "execute the frozen 160-input panel under licensed Magma F4",
            "obtain unaffiliated reproduction and novelty review",
        ],
        "ledger_sha256": sha256(LEDGER),
        "scoreboard_sha256": sha256(SCOREBOARD),
        "gate_status_sha256": sha256(GATE_STATUS),
        "autolab_protocol_sha256": sha256(AUTOLAB_PROTOCOL),
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


def write_new(path: Path, value: dict[str, Any]) -> None:
    require(not path.exists(), f"refusing to overwrite {path}")
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def build(output: Path) -> dict[str, Any]:
    require(not output.exists(), f"refusing to overwrite {output}")
    output.mkdir(parents=True)
    audit = compose()
    audit_path = output / "audit.json"
    write_new(audit_path, audit)
    write_new(
        output / "result-seal.json",
        {
            "schema": SEAL_SCHEMA,
            "status": "audit_frozen",
            "audit_sha256": sha256(audit_path),
        },
    )
    return audit


def verify(output: Path) -> dict[str, Any]:
    # Stage 124 is an immutable historical snapshot.  The mutable documents
    # it pinned by hash (the boundary ledger and its Markdown twin, the gate
    # status) moved again on 2026-09-21 with the boundary ledger's Round 2
    # (folded pair tables, walk targets, the exact counting ceiling, the
    # balanced Koblitz base), so recomposing here would compare two points
    # in time.  Stage 125 recomposes the current evidence and chains to this
    # seal (compose_koblitz_stage125_gate_audit.py); this verify checks the
    # seal and the frozen audit only, as Stage 109's does.
    seal = load(output / "result-seal.json", "Stage-124 seal")
    require(seal.get("schema") == SEAL_SCHEMA, "seal schema changed")
    require(sha256(output / "audit.json") == seal.get("audit_sha256"), "audit seal changed")
    committed = load(output / "audit.json", "Stage-124 audit")
    require(committed.get("schema") == SCHEMA, "committed audit schema changed")
    require(committed.get("status") == "current_seven_gate_audit_verified", "committed audit status changed")
    return committed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    build_parser = sub.add_parser("build")
    build_parser.add_argument("--output", type=Path, required=True)
    verify_parser = sub.add_parser("verify")
    verify_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        value = build(args.output.resolve()) if args.command == "build" else verify(args.output.resolve(strict=True))
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage124Error) as error:
        raise SystemExit(f"stage124-gate-audit: {error}")


if __name__ == "__main__":
    main()

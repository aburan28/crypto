"""Compose the current seven-gate audit through the Stage 108 selection,
re-sealed over the boundary ledger as updated by Round 5 on 2026-09-22.

Stage 109 sealed this audit on 2026-09-13 by pinning the mutable gate
documents (the boundary ledger and its Markdown twin, the gate status, the
autolab protocol) by hash; Stages 124, 125, 126 and 127 re-sealed it as the
index-calculus boundary ledger's first four rounds moved them.  Round 5
(research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md section 13)
moved them once more: the walk's sixteen restart offsets are now drawn on
first use instead of at setup, so a row is charged for the offsets it took
rather than for a pool it may never reach.  The saving is confined to the
small rungs and no Koblitz gate fact moved.  Stage 128 recomposes the same
seven gates over the current documents and chains to the Stage 127, 126,
125, 124 and 109 seals, all of which verify as immutable historical
snapshots by their own scripts.
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
STAGE124 = EVIDENCE / "stage-124-current-gate-audit-20260921"
STAGE125 = EVIDENCE / "stage-125-current-gate-audit-20260921"
STAGE126 = EVIDENCE / "stage-126-current-gate-audit-20260921"
STAGE127 = EVIDENCE / "stage-127-current-gate-audit-20260921"
LEDGER = REPO / "docs/ic/boundary_targets.json"
SCOREBOARD = REPO / "docs/ic/BOUNDARY_TARGETS.md"
GATE_STATUS = EVIDENCE / "GATE_STATUS.md"
AUTOLAB_PROTOCOL = REPO / "research/sat_factor_base_review_20260908/autolab/protocol.json"

SCHEMA = "koblitz_stage128_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage128_current_gate_audit_seal.v1"
LEDGER_UPDATED = "2026-09-22"
GATE_STATUS_MARKER = "Current through Stage 128"


class Stage128Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage128Error(message)


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
    # Stages 109, 124, 125, 126 and 127 verify as historical snapshots: their
    # own scripts check the seal and the frozen audit, never the documents
    # that have moved since.
    stage109 = replay("compose_koblitz_stage109_gate_audit.py", STAGE109)
    stage124 = replay("compose_koblitz_stage124_gate_audit.py", STAGE124)
    stage125 = replay("compose_koblitz_stage125_gate_audit.py", STAGE125)
    stage126 = replay("compose_koblitz_stage126_gate_audit.py", STAGE126)
    stage127 = replay("compose_koblitz_stage127_gate_audit.py", STAGE127)
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
    # The gate facts Stage 109 sealed, and every stage since re-sealed, must
    # be the ones this audit re-seals.
    require(stage124["selected_n53"] == stage109["selected_n53"], "Stage 124 and Stage 109 sealed different facts")
    require(stage125["selected_n53"] == stage109["selected_n53"], "Stage 125 and Stage 109 sealed different facts")
    require(stage126["selected_n53"] == stage109["selected_n53"], "Stage 126 and Stage 109 sealed different facts")
    require(stage127["selected_n53"] == stage109["selected_n53"], "Stage 127 and Stage 109 sealed different facts")
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
            "stage124_audit_sha256": sha256(STAGE124 / "audit.json"),
            "stage125_audit_sha256": sha256(STAGE125 / "audit.json"),
            "stage126_audit_sha256": sha256(STAGE126 / "audit.json"),
            "stage127_audit_sha256": sha256(STAGE127 / "audit.json"),
            "stage99_status": historical["status"],
            "stage108_status": selected["status"],
            "stage109_status": stage109["status"],
            "stage124_status": stage124["status"],
            "stage125_status": stage125["status"],
            "stage126_status": stage126["status"],
            "stage127_status": stage127["status"],
        },
        "since_stage127": {
            "what_moved": "docs/ic/boundary_targets.json, docs/ic/BOUNDARY_TARGETS.md and GATE_STATUS.md: the boundary ledger's Round 5, which draws the walk's sixteen restart offsets on first use instead of at setup. On the Round-4 ladder 110 of 204 walk rows had bought all sixteen and taken none. The offsets now come off their own random stream and are taken in a fixed cycle, so the eager and the lazy arm walk bit-identical trajectories: 537 of 537 rows identical in every counter but the pool's, zero trajectories moved, and zero rows whose saving differs from the offsets they stopped drawing, with the same on a holdout seed and an eight-repeat headline. The saving is confined to the small rungs, matching (48 log2 r - 48)/sqrt r to within 1 to 9 per cent. The same round measures that re-seeding the walk alone, at no change of cost, moves 204 of 537 rows with a median of 1.0000 and a range of 0.54 to 1.83, so the headline rows read lower than Round 4 for a reason classed accounting rather than progress", "what_did_not": "every Koblitz gate fact Stage 109 sealed and Stages 124, 125, 126 and 127 re-sealed: the selected n=53 base hash, support entries and bytes, the five paired wall ratios, the core and fresh-build ratios, the unknown-scalar ratio, and the seven gate statuses",
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
    # Stage 128 is the current audit: recompose from the live documents and
    # require the result to equal the sealed one.  When the documents it pins
    # move again, the next stage re-seals and this verify becomes historical,
    # exactly as Stage 109's, 124's, 125's, 126's and 127's did.
    seal = load(output / "result-seal.json", "Stage-128 seal")
    require(seal.get("schema") == SEAL_SCHEMA, "seal schema changed")
    require(sha256(output / "audit.json") == seal.get("audit_sha256"), "audit seal changed")
    current = compose()
    require(current == load(output / "audit.json", "Stage-128 audit"), "current audit changed")
    return current


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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage128Error) as error:
        raise SystemExit(f"stage128-gate-audit: {error}")


if __name__ == "__main__":
    main()

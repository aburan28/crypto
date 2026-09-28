"""Compose the current seven-gate audit through the Stage 108 selection,
re-sealed after the gate-6 rho control was re-measured on 2026-09-22.

Stage 109 sealed this audit on 2026-09-13 by pinning the mutable gate
documents by hash; Stages 124 to 128 re-sealed it as the index-calculus
boundary ledger moved them, and no Koblitz gate fact moved.  Stage 129
moves one: the selected n=53 route passed the online wall gate against
koblitz_rho_fixture, a correct but unbatched control (one inversion per
step, squaring-chain canonicalisation, every point stored, one thread).
Against a batched signed-Frobenius rho on the same target and host
(aburan28/cryptanalysis suite/examples/koblitz_batched_rho.rs) the selected
direct takes a median 17.3 s against 0.49 s (1 thread) and 0.20 s (4
threads), 46.3M support queries against about 0.56M rho steps.  Gate 6 is
therefore false outright.  The sealed Stage 108 measurements themselves are
unchanged, and this audit chains to the Stage 128, 127, 126, 125, 124 and
109 seals, all of which verify as immutable historical snapshots.
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
STAGE128 = EVIDENCE / "stage-128-current-gate-audit-20260922"
LEDGER = REPO / "docs/ic/boundary_targets.json"
SCOREBOARD = REPO / "docs/ic/BOUNDARY_TARGETS.md"
GATE_STATUS = EVIDENCE / "GATE_STATUS.md"
AUTOLAB_PROTOCOL = REPO / "research/sat_factor_base_review_20260908/autolab/protocol.json"

SCHEMA = "koblitz_stage129_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage129_current_gate_audit_seal.v1"
LEDGER_UPDATED = "2026-09-22"
GATE_STATUS_MARKER = "Current through Stage 129"
BATCHED_RHO_MARKER = "## The n=53 crossover against a batched rho"


class Stage129Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage129Error(message)


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
    # Stages 109 and 124 to 128 verify as historical snapshots: their
    # own scripts check the seal and the frozen audit, never the documents
    # that have moved since.
    stage109 = replay("compose_koblitz_stage109_gate_audit.py", STAGE109)
    stage124 = replay("compose_koblitz_stage124_gate_audit.py", STAGE124)
    stage125 = replay("compose_koblitz_stage125_gate_audit.py", STAGE125)
    stage126 = replay("compose_koblitz_stage126_gate_audit.py", STAGE126)
    stage127 = replay("compose_koblitz_stage127_gate_audit.py", STAGE127)
    stage128 = replay("compose_koblitz_stage128_gate_audit.py", STAGE128)
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
    require(stage128["selected_n53"] == stage109["selected_n53"], "Stage 128 and Stage 109 sealed different facts")
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
    require(BATCHED_RHO_MARKER in gate_text, "gate-status batched-rho control missing")
    require(
        "**False. The online wall gate passed only against the fixture rho control" in gate_text,
        "gate-status gate-6 verdict changed",
    )
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
            "stage128_audit_sha256": sha256(STAGE128 / "audit.json"),
            "stage99_status": historical["status"],
            "stage108_status": selected["status"],
            "stage109_status": stage109["status"],
            "stage124_status": stage124["status"],
            "stage125_status": stage125["status"],
            "stage126_status": stage126["status"],
            "stage127_status": stage127["status"],
            "stage128_status": stage128["status"],
        },
        "since_stage128": {
            "what_moved": "GATE_STATUS.md: the gate-6 rho control. The Stage 108 online wall crossing was measured against koblitz_rho_fixture, which inverts once per step, canonicalises by squaring chains, stores every point and runs on one thread (about 9.5 us a step). A batched signed-Frobenius rho (one inversion shared across 64 walks, normal-basis orbit key, distinguished points only; aburan28/cryptanalysis suite/examples/koblitz_batched_rho.rs, PR #54) on the same target and host solves in a median 0.49 s on 1 thread and 0.20 s on 4 threads against 17.3 s for the selected direct, every run verified by [d]G=Q. That is 46.3M support queries against about 0.56M rho steps; the query count r/(n|F|) grows as r^(2/3) against rho's r^(1/2). The host is an M4 Pro without the x86 PCLMUL paths, so the EPYC ratios need their own rerun",
            "what_did_not": "the Stage 108 measurements Stage 109 sealed and Stages 124 to 128 re-sealed: the selected n=53 base hash, support entries and bytes, the five paired wall ratios against the fixture control, the core and fresh-build ratios, and the unknown-scalar ratio; gates 1 to 5 and 7",
            "reference": "suite/docs/ic/runs/koblitz-n53-rho-control-20260922.json in aburan28/cryptanalysis",
        },
        "gates": {
            "1_all_stage_resource_charging": "partial_missing_licensed_magma",
            "2_same_instance_solver_panel": "partial_missing_licensed_magma_full_panel",
            "3_single_core_core_memory_conflicts_wall": "partial_missing_magma_and_selected_single_core_refresh",
            "4_n31_n41_larger_pdp_scaling": "satisfied_finite_execution_coverage",
            "5_unknown_scalar_without_constructed_base_logs": "satisfied_finite_n23_n31_n41_n53",
            "6_full_cost_vs_automorphism_rho": "false_online_wall_pass_only_vs_fixture_rho_loses_to_batched_rho",
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
            "beat a batched signed-Frobenius rho, not the fixture control, on the EPYC gate host",
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
    # Stage 129 is the current audit: recompose from the live documents and
    # require the result to equal the sealed one.  When the documents it pins
    # move again, the next stage re-seals and this verify becomes historical,
    # exactly as Stage 109's and 124's to 128's did.
    seal = load(output / "result-seal.json", "Stage-129 seal")
    require(seal.get("schema") == SEAL_SCHEMA, "seal schema changed")
    require(sha256(output / "audit.json") == seal.get("audit_sha256"), "audit seal changed")
    current = compose()
    require(current == load(output / "audit.json", "Stage-129 audit"), "current audit changed")
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage129Error) as error:
        raise SystemExit(f"stage129-gate-audit: {error}")


if __name__ == "__main__":
    main()

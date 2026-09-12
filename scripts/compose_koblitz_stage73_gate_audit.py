#!/usr/bin/env python3
"""Compose the current audit from retained fused, unknown, and single-core results."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import verify_koblitz_stage68_fused_result as stage68
import verify_koblitz_stage72_current_results as stage72


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE59 = STAGE / "stage-59-current-gate-audit-20260912/audit.json"
STAGE68 = STAGE / "stage-68-n53-fused-result-20260912/verification.json"
STAGE72 = STAGE / "stage-72-current-results-20260912/verification.json"
LEDGER = REPO / "docs/ic/boundary_targets.json"
SCOREBOARD = REPO / "docs/ic/BOUNDARY_TARGETS.md"
STATUS = STAGE / "GATE_STATUS.md"
AUTOLAB_PROTOCOL = REPO / "research/sat_factor_base_review_20260908/autolab/protocol.json"
DEFAULT_OUTPUT = STAGE / "stage-73-current-gate-audit-20260912"
SCHEMA = "koblitz_stage73_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage73_current_gate_audit_seal.v1"


class AuditError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def validate_documents() -> dict[str, Any]:
    ledger = load(LEDGER, "IC boundary ledger")
    records = ledger["regimes"]["koblitz"]["records"]
    current = records["vs_rho"]["current"]
    metrics = current["metrics"]
    require(current.get("verdict") == "N41_ONLINE_CHARGED_CROSSOVER_N53_FUSED_NEAR_WALL_LOSS", "Stage-73 ledger verdict changed")
    require(metrics.get("n53_parallel_direct_over_rho_wall_ratio") == 1.1828715401273624, "Stage-73 known ratio changed")
    require(metrics.get("n53_unknown_scalar_direct_over_rho_wall_ratio") == 1.4676110945821472, "Stage-73 unknown ratio changed")
    require(metrics.get("n53_single_core_direct_over_rho_wall_ratio") == 1.993515234548992, "Stage-73 single-core ratio changed")
    require(metrics.get("n53_support_table_allocated_bytes") == 738_197_504, "Stage-73 support allocation changed")
    require(records["rank"]["current"]["metrics"].get("relations_collected") == 95, "Stage-73 rank row changed")
    require(records["end_to_end_dlp"]["current"]["metrics"].get("target_scalar_constructed_or_supplied") is False, "Stage-73 unknown target boundary changed")
    protocol = load(AUTOLAB_PROTOCOL, "IC boundary autolab protocol")
    n53 = protocol["beats"]["koblitz.vs_rho.n53_probe"]
    require(n53.get("eta") == [1, 128] and n53.get("direct", {}).get("query_mode") == "pair_pair_parallel_4096", "Stage-73 autolab n=53 beat changed")
    require(n53.get("env", {}).get("KIC_PIPELINED_SUPPORT_EXPANSION") == "1", "Stage-73 autolab optimized stack changed")
    documents = SCOREBOARD.read_text() + STATUS.read_text()
    for marker in ("1.183x", "1.468x", "1.994x", "Satisfied for finite degrees 23, 31, 41, and 53"):
        require(marker in documents, f"Stage-73 documents are missing {marker!r}")
    return {
        "ledger": {"path": str(LEDGER.relative_to(REPO)), "sha256": custody.sha256_file(LEDGER, "IC boundary ledger")},
        "scoreboard": {"path": str(SCOREBOARD.relative_to(REPO)), "sha256": custody.sha256_file(SCOREBOARD, "IC boundary scoreboard")},
        "gate_status": {"path": str(STATUS.relative_to(REPO)), "sha256": custody.sha256_file(STATUS, "gate status")},
        "autolab_protocol": {"path": str(AUTOLAB_PROTOCOL.relative_to(REPO)), "sha256": custody.sha256_file(AUTOLAB_PROTOCOL, "autolab protocol")},
    }


def compose() -> dict[str, Any]:
    fused = stage68.verify(STAGE68.parent)
    current = stage72.verify(STAGE72.parent)
    require(fused == load(STAGE68, "Stage-68 result"), "Stage-68 binding changed")
    require(current == load(STAGE72, "Stage-72 result"), "Stage-72 binding changed")
    predecessor = load(STAGE59, "Stage-59 predecessor")
    require(predecessor.get("schema") == "koblitz_stage59_current_gate_audit.v1", "Stage-59 predecessor changed")
    result = copy.deepcopy(predecessor)
    result.update({
        "schema": SCHEMA,
        "status": "current_evidence_audited_n53_fused_unknown_single_core_losses_gates_incomplete",
        "predecessor_stage59": {"path": str(STAGE59.relative_to(REPO)), "sha256": custody.sha256_file(STAGE59, "Stage-59 audit")},
        "stage68_binding": {"path": str(STAGE68.relative_to(REPO)), "sha256": custody.sha256_file(STAGE68, "Stage-68 result"), "verification": fused},
        "stage72_binding": {"path": str(STAGE72.relative_to(REPO)), "sha256": custody.sha256_file(STAGE72, "Stage-72 result"), "verification": current},
        "boundary_ledger_binding": validate_documents(),
        "all_seven_gates_passed": False,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "overall": "n=53 fused known-target, optimized unknown-scalar, and kernel-enforced single-core results are retained; licensed Magma, whole-process crossover, and unaffiliated review remain false",
    })
    measurement = {
        "n": 53,
        "factor_base_points": 9964,
        "orbit_columns": 94,
        "published_q": fused["candidate"]["published_q"],
        "preferred_query_mode": "pair_pair_parallel_4096",
        "preferred_field_mul_backend": "x86_64_pclmulqdq",
        "preferred_field_reduction_backend": "n53_fixed_two_fold",
        "preferred_field_product_pipeline": "x86_64_pclmul_n53_fused_reduce",
        "parallel_direct_wall_seconds": fused["processes"]["candidate"]["wall_seconds"],
        "parallel_direct_core_seconds": fused["processes"]["candidate"]["total_core_seconds"],
        "parallel_direct_peak_rss_bytes": fused["processes"]["candidate"]["peak_rss_bytes"],
        "parallel_direct_over_rho_wall_ratio": fused["comparison"]["fused_over_rho_wall_ratio"],
        "rho_wall_seconds": fused["processes"]["rho"]["wall_seconds"],
        "fused_same_binary_wall_speedup": fused["comparison"]["fused_direct_wall_speedup"],
        "fused_same_binary_core_seconds_ratio": fused["comparison"]["fused_core_seconds_ratio"],
        "fresh_build_plus_parallel_over_rho": fused["comparison"]["fresh_build_plus_fused_over_rho_wall_ratio"],
        "support_table_allocated_bytes": fused["candidate"]["support_table_allocated_bytes"],
        "support_queries": fused["candidate"]["support_queries"],
        "relations": fused["candidate"]["admitted_relations"],
        "unknown_scalar": {
            "published_q": current["unknown_scalar"]["verification"]["published_q"],
            "recovered_scalar": current["unknown_scalar"]["verification"]["recovered_scalar"],
            "relations": current["unknown_scalar"]["verification"]["relations"],
            "support_queries": current["unknown_scalar"]["direct"]["support_queries"],
            "direct_wall_seconds": current["unknown_scalar"]["processes"]["direct"]["wall_seconds"],
            "rho_wall_seconds": current["unknown_scalar"]["processes"]["rho"]["wall_seconds"],
            "direct_over_rho_wall_ratio": current["unknown_scalar"]["comparison"]["direct_over_rho_wall_ratio"],
            "fresh_build_plus_direct_over_rho_wall_ratio": current["unknown_scalar"]["comparison"]["fresh_build_plus_direct_over_rho_wall_ratio"],
            "target_scalar_constructed_or_supplied": False,
            "factor_base_logs_known_by_construction": False,
        },
        "single_core": {
            "affinity": current["single_core"]["affinity"],
            "direct_wall_seconds": current["single_core"]["processes"]["direct"]["wall_seconds"],
            "direct_core_seconds": current["single_core"]["processes"]["direct"]["total_core_seconds"],
            "direct_peak_rss_bytes": current["single_core"]["processes"]["direct"]["peak_rss_bytes"],
            "rho_wall_seconds": current["single_core"]["processes"]["rho"]["wall_seconds"],
            "direct_over_rho_wall_ratio": current["single_core"]["comparison"]["direct_over_rho_wall_ratio"],
            "direct_over_rho_core_ratio": current["single_core"]["comparison"]["direct_core_over_rho_core_ratio"],
        },
    }
    result["current_measurements"]["n53_same_target"] = measurement
    gate1 = result["gates"]["1_full_cost_accounting"]
    gate1["proved"].extend([
        "Stage-68 charges clean build, separate/fused direct, rho, wall, core-seconds, process and tree RSS",
        "Stage-72 charges optimized public unknown-scalar and kernel-enforced single-core direct/rho runs with full build and workflow envelopes",
    ])
    gate3 = result["gates"]["3_required_resource_fields"]
    gate3["proved"].append("Stage-72 retains parent/child CPU-0 affinity plus direct/rho wall, core-seconds, RSS, and operation counts")
    gate5 = result["gates"]["5_unknown_scalar"]
    gate5["status"] = "finite_n23_n31_n41_n53_unknown_complete"
    gate5["proved"] = [
        "5 degree-23 public targets",
        "5 degree-31 public hash targets",
        "5 degree-41 public hash targets",
        "Stage-72 n=53 target derives from public seed 53001 without constructing or supplying its scalar",
        "Stage-72 derives all 94 factor-base logs from 95 verified relations",
        "Stage-72 direct IC and signed-Frobenius rho independently recover the same scalar and verify [d]G=Q",
    ]
    gate5["limitations"] = ["finite public synthetic controls only", "one retained optimized n=53 public hash target; independent seed panel and unaffiliated replay remain absent"]
    gate6 = result["gates"]["6_automorphism_rho"]
    gate6["status"] = "n41_online_crossover_full_cost_loss_n53_fused_unknown_single_core_wall_losses"
    gate6["proved"] = [
        "Stage-39 n=41 online IC/rho ratio is 0.285545560549 while amortized IC remains 5.025825 times slower",
        "Stage-68 fused n=53 direct remains 1.182872 times same-target signed-Frobenius rho",
        "Stage-68 fresh build plus fused direct remains 19.474746 times rho wall",
        "Stage-72 optimized public unknown-scalar direct remains 1.467611 times rho and 25.359925 times rho with fresh build",
        "Stage-72 CPU-0 direct remains 1.993515 times rho wall and 1.978676 times rho core-seconds",
    ]
    gate6["limitations"] = [
        "n=41 crossover is online only and assumes a factor-base log database",
        "no n=53 whole-process or fresh-build cell crosses",
        "hosted absolute timings vary between runners; optimization attribution uses same-run A/B controls",
        "all improvements are finite constants and do not change the asymptotic exponent",
    ]
    return result


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-73 output must be new")
    output.mkdir()
    custody.write_json_new(output / "audit.json", compose())
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "audit_frozen", "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "audit-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "audit.json", "committed Stage-73 audit")
    require(committed == compose(), "committed Stage-73 audit differs from source evidence")
    seal = load(output / "audit-seal.json", "Stage-73 audit seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-73 audit seal changed")
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-73 audit inventory changed")
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
        raise SystemExit(f"stage73-audit: {error}")


if __name__ == "__main__":
    main()

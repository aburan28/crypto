#!/usr/bin/env python3
"""Compose the Stage 56 PCLMUL and Stage 58 unknown-scalar results."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import compose_koblitz_stage53_gate_audit as stage53
import verify_koblitz_stage56_pclmul_result as stage56
import verify_koblitz_stage58_unknown_result as stage58


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE53 = STAGE / "stage-53-current-gate-audit-20260912/audit.json"
STAGE56 = STAGE / "stage-56-n53-pclmul-result-20260912/verification.json"
STAGE58 = STAGE / "stage-58-n53-unknown-result-20260912/verification.json"
LEDGER = REPO / "docs/ic/boundary_targets.json"
SCOREBOARD = REPO / "docs/ic/BOUNDARY_TARGETS.md"
STATUS = STAGE / "GATE_STATUS.md"
DEFAULT_OUTPUT = STAGE / "stage-59-current-gate-audit-20260912"
SCHEMA = "koblitz_stage59_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage59_current_gate_audit_seal.v1"


class AuditError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def validate_current_documents() -> dict[str, Any]:
    ledger = load(LEDGER, "IC boundary ledger")
    current = ledger["regimes"]["koblitz"]["records"]["vs_rho"]["current"]
    require(current.get("verdict") == "N41_ONLINE_CHARGED_CROSSOVER_N53_PARALLEL_WALL_LOSS", "Stage-59 ledger verdict changed")
    metrics = current["metrics"]
    require(metrics.get("n53_preferred_query_mode") == "pair_pair_parallel_4096", "Stage-59 preferred mode changed")
    require(metrics.get("n53_parallel_direct_over_rho_wall_ratio") == 3.688473082311116, "Stage-59 ratio changed")
    require(metrics.get("n53_unknown_scalar_completed") is True, "Stage-59 unknown-scalar status changed")
    require(metrics.get("n53_unknown_scalar_direct_over_rho_wall_ratio") == 4.428199126138434, "Stage-59 unknown-scalar ratio changed")
    require(metrics.get("n53_parallel_4096_support_table_allocated_bytes") == 1_275_068_416, "Stage-59 support allocation changed")
    documents = SCOREBOARD.read_text() + STATUS.read_text()
    for marker in ("3.688x slower", "1.419x faster", "0.614x its CPU", "Satisfied for finite degrees 23, 31, 41, and 53"):
        require(marker in documents, f"Stage-59 documents are missing {marker!r}")
    return {
        "ledger": {"path": str(LEDGER.relative_to(REPO)), "sha256": custody.sha256_file(LEDGER, "IC boundary ledger")},
        "scoreboard": {"path": str(SCOREBOARD.relative_to(REPO)), "sha256": custody.sha256_file(SCOREBOARD, "IC boundary scoreboard")},
        "gate_status": {"path": str(STATUS.relative_to(REPO)), "sha256": custody.sha256_file(STATUS, "gate status")},
    }


def compose() -> dict[str, Any]:
    predecessor = stage53.verify(STAGE53.parent)
    pclmul = stage56.verify(STAGE56.parent)
    unknown = stage58.verify(STAGE58.parent)
    require(predecessor == load(STAGE53, "Stage-53 predecessor"), "Stage-53 predecessor binding changed")
    require(pclmul == load(STAGE56, "Stage-56 result"), "Stage-56 result binding changed")
    require(unknown == load(STAGE58, "Stage-58 result"), "Stage-58 result binding changed")
    pclmul_verification = pclmul["hosted"]["verification"]
    unknown_verification = unknown["verification"]
    result = copy.deepcopy(predecessor)
    result.update({
        "schema": SCHEMA,
        "status": "current_evidence_audited_n53_pclmul_unknown_wall_loss_gates_incomplete",
        "predecessor_stage53": {"path": str(STAGE53.relative_to(REPO)), "sha256": custody.sha256_file(STAGE53, "Stage-53 audit")},
        "stage56_binding": {"path": str(STAGE56.relative_to(REPO)), "sha256": custody.sha256_file(STAGE56, "Stage-56 result"), "verification": pclmul},
        "stage58_binding": {"path": str(STAGE58.relative_to(REPO)), "sha256": custody.sha256_file(STAGE58, "Stage-58 result"), "verification": unknown},
        "boundary_ledger_binding": validate_current_documents(),
        "all_seven_gates_passed": False,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "overall": "n=53 PCLMUL and public unknown-scalar completion are verified; licensed Magma, whole-process crossover, and unaffiliated review remain false",
    })
    measurement = result["current_measurements"]["n53_same_target"]
    measurement.update({
        "preferred_query_mode": "pair_pair_parallel_4096",
        "preferred_field_mul_backend": "x86_64_pclmulqdq",
        "portable_wall_seconds_stage56": pclmul_verification["baseline_wall_seconds"],
        "pclmul_wall_seconds_stage56": pclmul_verification["candidate_wall_seconds"],
        "rho_wall_seconds_stage56": pclmul_verification["rho_wall_seconds"],
        "parallel_direct_wall_seconds": pclmul_verification["candidate_wall_seconds"],
        "parallel_direct_over_rho_wall_ratio": pclmul_verification["candidate_over_rho_wall_ratio"],
        "parallel_direct_wall_speedup": pclmul_verification["candidate_direct_wall_speedup"],
        "parallel_direct_wall_speedup_baseline": "portable_sparse_carryless same-binary Stage 56",
        "parallel_direct_core_seconds_ratio": pclmul_verification["candidate_direct_core_seconds_ratio"],
        "parallel_direct_core_seconds_ratio_baseline": "portable_sparse_carryless same-binary Stage 56",
        "pclmul_wall_speedup": pclmul_verification["candidate_direct_wall_speedup"],
        "pclmul_core_seconds_ratio": pclmul_verification["candidate_direct_core_seconds_ratio"],
        "pclmul_collection_speedup": pclmul_verification["candidate_collection_speedup"],
        "fresh_build_plus_parallel_over_rho": pclmul_verification["fresh_build_plus_candidate_over_rho_wall_ratio"],
        "unknown_scalar": {
            "published_q": unknown_verification["published_q"],
            "recovered_scalar": unknown_verification["recovered_scalar"],
            "relations": unknown_verification["relations"],
            "direct_wall_seconds": unknown_verification["direct_wall_seconds"],
            "rho_wall_seconds": unknown_verification["rho_wall_seconds"],
            "direct_over_rho_wall_ratio": unknown_verification["direct_over_rho_wall_ratio"],
            "fresh_build_plus_direct_over_rho_wall_ratio": unknown_verification["fresh_build_plus_direct_over_rho_wall_ratio"],
            "target_scalar_constructed_or_supplied": False,
            "factor_base_logs_known_by_construction": False,
        },
    })
    gate1 = result["gates"]["1_full_cost_accounting"]
    gate1["proved"].append("Stage-56 charges clean build, portable/PCLMUL direct, rho, wall, core-seconds, and tree RSS")
    gate1["proved"].append("Stage-58 charges clean build, public unknown-scalar direct/rho, target derivation, relation collection, LA, wall, core-seconds, and tree RSS")
    gate3 = result["gates"]["3_required_resource_fields"]
    gate3["proved"].append("Stages 56 and 58 retain process wall, user/system and total core-seconds, RSS, support operations, target derivation, and workflow resources")
    gate5 = result["gates"]["5_unknown_scalar"]
    gate5["status"] = "finite_n23_n31_n41_n53_unknown_complete"
    gate5["proved"] = [
        "5 degree-23 public targets",
        "5 degree-31 public hash targets",
        "5 degree-41 public hash targets",
        "Stage-58 n=53 target derives from public seed 53001 without constructing or supplying its scalar",
        "Stage-58 derives all 94 factor-base logs from 165 verified relations",
        "Stage-58 direct IC and signed-Frobenius rho independently recover the same scalar and verify [d]G=Q",
    ]
    gate5["limitations"] = [
        "finite public synthetic controls only",
        "one n=53 public hash target; independent additional seeds and unaffiliated replay would strengthen the result",
    ]
    gate6 = result["gates"]["6_automorphism_rho"]
    gate6["status"] = "n41_online_crossover_n41_amortised_loss_n53_pclmul_and_unknown_wall_loss"
    gate6["proved"] = [
        "Stage-39 n=41 online IC/rho ratio is 0.285545560549 while amortized IC remains 5.025825 times slower",
        "Stage-40 four-core n=41 amortized IC remains 3.605523 times slower",
        "Stage-56 PCLMUL is 1.418695 times faster in wall and uses 0.614352 times portable CPU",
        "Stage-56 PCLMUL direct remains 3.688473 times slower than same-target signed-Frobenius rho",
        "Stage-56 fresh build plus PCLMUL direct remains 21.182086 times rho wall",
        "Stage-58 public unknown-scalar direct remains 4.428199 times same-target rho and 25.722944 times rho with fresh build",
    ]
    gate6["limitations"] = [
        "n=41 crossover is online only and assumes a factor-base log database",
        "n=41 amortized and full available costs do not cross",
        "n=53 known-answer and unknown-scalar whole-process direct do not cross",
        "hosted absolute timings vary between runners, so PCLMUL attribution uses the same-binary A/B",
        "all improvements are finite constants and do not change the asymptotic exponent",
    ]
    return result


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-59 output must be new")
    output.mkdir()
    custody.write_json_new(output / "audit.json", compose())
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "audit_frozen", "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "audit-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "audit.json", "committed Stage-59 audit")
    require(committed == compose(), "committed Stage-59 audit differs from source evidence")
    seal = load(output / "audit-seal.json", "Stage-59 audit seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-59 audit seal changed")
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-59 audit inventory changed")
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
        raise SystemExit(f"stage59-audit: {error}")


if __name__ == "__main__":
    main()

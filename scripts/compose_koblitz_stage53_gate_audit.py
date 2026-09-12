#!/usr/bin/env python3
"""Compose the Stage 52 lazy-label result into the current gate audit."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import compose_koblitz_stage50_gate_audit as stage50
import verify_koblitz_stage52_lazy_result as stage52


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE50 = STAGE / "stage-50-current-gate-audit-20260912/audit.json"
STAGE52 = STAGE / "stage-52-n53-lazy-result-20260912/verification.json"
LEDGER = REPO / "docs/ic/boundary_targets.json"
SCOREBOARD = REPO / "docs/ic/BOUNDARY_TARGETS.md"
STATUS = STAGE / "GATE_STATUS.md"
DEFAULT_OUTPUT = STAGE / "stage-53-current-gate-audit-20260912"
SCHEMA = "koblitz_stage53_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage53_current_gate_audit_seal.v1"


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
    require(current.get("verdict") == "N41_ONLINE_CHARGED_CROSSOVER_N53_PARALLEL_WALL_LOSS", "Stage-53 ledger verdict changed")
    metrics = current["metrics"]
    require(metrics.get("n53_preferred_query_mode") == "pair_pair_parallel_4096", "Stage-53 preferred mode changed")
    require(metrics.get("n53_parallel_direct_over_rho_wall_ratio") == 5.16632272222771, "Stage-53 ratio changed")
    require(metrics.get("n53_parallel_4096_support_table_allocated_bytes") == 1_275_068_416, "Stage-53 support allocation changed")
    documents = SCOREBOARD.read_text() + STATUS.read_text()
    for marker in ("5.166x slower", "1.074x faster", "0.956x the CPU", "134,217,728"):
        require(marker in documents, f"Stage-53 documents are missing {marker!r}")
    return {
        "ledger": {"path": str(LEDGER.relative_to(REPO)), "sha256": custody.sha256_file(LEDGER, "IC boundary ledger")},
        "scoreboard": {"path": str(SCOREBOARD.relative_to(REPO)), "sha256": custody.sha256_file(SCOREBOARD, "IC boundary scoreboard")},
        "gate_status": {"path": str(STATUS.relative_to(REPO)), "sha256": custody.sha256_file(STATUS, "gate status")},
    }


def compose() -> dict[str, Any]:
    predecessor = stage50.verify(STAGE50.parent)
    latest = stage52.verify(STAGE52.parent)
    require(predecessor == load(STAGE50, "Stage-50 predecessor"), "Stage-50 predecessor binding changed")
    require(latest == load(STAGE52, "Stage-52 result"), "Stage-52 result binding changed")
    verification = latest["hosted"]["verification"]
    result = copy.deepcopy(predecessor)
    result.update({
        "schema": SCHEMA,
        "status": "current_evidence_audited_n53_lazy_4096_wall_loss_gates_incomplete",
        "predecessor_stage50": {"path": str(STAGE50.relative_to(REPO)), "sha256": custody.sha256_file(STAGE50, "Stage-50 audit")},
        "stage52_binding": {"path": str(STAGE52.relative_to(REPO)), "sha256": custody.sha256_file(STAGE52, "Stage-52 result"), "verification": latest},
        "boundary_ledger_binding": validate_current_documents(),
        "all_seven_gates_passed": False,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "overall": "n=41 fixed-algebraic online crossover and n=53 scratch/width/lazy-label improvements are verified; n=53 whole-process, amortized/full cost, Magma, and external gates remain false",
    })
    measurement = result["current_measurements"]["n53_same_target"]
    measurement.update({
        "preferred_query_mode": "pair_pair_parallel_4096",
        "parallel_1024_wall_seconds_stage52": verification["baseline_wall_seconds"],
        "parallel_4096_wall_seconds_stage52": verification["candidate_wall_seconds"],
        "rho_wall_seconds_stage52": verification["rho_wall_seconds"],
        "parallel_direct_wall_seconds": verification["candidate_wall_seconds"],
        "parallel_direct_over_rho_wall_ratio": verification["candidate_over_rho_wall_ratio"],
        "parallel_direct_wall_speedup": verification["candidate_direct_wall_speedup"],
        "parallel_direct_core_seconds_ratio": verification["candidate_direct_core_seconds_ratio"],
        "fresh_build_plus_parallel_over_rho": verification["fresh_build_plus_candidate_over_rho_wall_ratio"],
        "support_table_allocated_bytes": latest["support_table_allocated_bytes"],
        "support_table_bytes_saved": latest["support_table_bytes_saved"],
    })
    gate1 = result["gates"]["1_full_cost_accounting"]
    gate1["proved"].append("Stage-52 charges the clean build, lazy-label 1,024/4,096 direct arms, rho, wall, core-seconds, and process-tree RSS")
    gate3 = result["gates"]["3_required_resource_fields"]
    gate3["proved"].append("Stage-52 reports a 1,392,087,040-byte preferred-direct peak and 1,275,068,416-byte retained support allocation")
    gate6 = result["gates"]["6_automorphism_rho"]
    gate6["status"] = "n41_online_crossover_n41_amortised_loss_n53_lazy_4096_same_target_wall_loss"
    gate6["proved"] = [
        "Stage-39 n=41 online IC/rho ratio is 0.285545560549 while amortized IC remains 5.025825 times slower",
        "Stage-40 four-core n=41 amortized IC remains 3.605523 times slower and spends 1.286948 times one-core process CPU",
        "Stages 46–52 preserve the same 189 n=53 relation hashes, public target, rank endpoint, and recovered scalar",
        "Stage-52 4,096 cursors are 1.073641 times faster in wall and use 0.956319 times the CPU of 1,024 in the same hosted run",
        "Stage-52 lazy-label preferred direct remains 5.166323 times slower than same-target signed-Frobenius rho",
        "Stage-52 fresh build plus preferred direct remains 22.668683 times rho wall",
    ]
    gate6["limitations"] = [
        "n=41 crossover is online only and assumes a factor-base log database",
        "n=41 amortized and full available costs do not cross",
        "n=53 preferred same-target whole-process direct does not cross",
        "hosted absolute timings vary between runners, so width choices use same-run comparisons",
        "all improvements are finite constants and do not change the asymptotic exponent",
    ]
    return result


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-53 output must be new")
    output.mkdir()
    custody.write_json_new(output / "audit.json", compose())
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "audit_frozen", "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "audit-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "audit.json", "committed Stage-53 audit")
    require(committed == compose(), "committed Stage-53 audit differs from source evidence")
    seal = load(output / "audit-seal.json", "Stage-53 audit seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-53 audit seal changed")
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-53 audit inventory changed")
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
        raise SystemExit(f"stage53-audit: {error}")


if __name__ == "__main__":
    main()

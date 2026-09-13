#!/usr/bin/env python3
"""Compose the current seven-gate audit with the selected CPU-0 result."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import verify_koblitz_stage92_current_single_core as stage92


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE90 = STAGE / "stage-90-current-gate-audit-20260912/audit.json"
STAGE92 = STAGE / "stage-92-current-single-core-20260913/verification.json"
LEDGER = REPO / "docs/ic/boundary_targets.json"
SCOREBOARD = REPO / "docs/ic/BOUNDARY_TARGETS.md"
STATUS = STAGE / "GATE_STATUS.md"
AUTOLAB_PROTOCOL = REPO / "research/sat_factor_base_review_20260908/autolab/protocol.json"
DEFAULT_OUTPUT = STAGE / "stage-93-current-gate-audit-20260913"
SCHEMA = "koblitz_stage93_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage93_current_gate_audit_seal.v1"


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
    current = ledger["regimes"]["koblitz"]["records"]["vs_rho"]["current"]
    metrics = current["metrics"]
    require(
        metrics.get("n53_single_core_evidence_current_stack") is True,
        "Stage-93 current single-core marker changed",
    )
    require(
        metrics.get("n53_single_core_direct_over_rho_wall_ratio")
        == 1.3543032952267562,
        "Stage-93 single-core wall ratio changed",
    )
    require(
        metrics.get("n53_single_core_direct_over_rho_core_ratio")
        == 1.3401796163922106,
        "Stage-93 single-core CPU ratio changed",
    )
    protocol = load(AUTOLAB_PROTOCOL, "IC boundary autolab protocol")
    cpu0 = protocol["beats"]["koblitz.vs_rho.n53_probe"]["current_single_core"]
    require(
        cpu0.get("affinity_cpu") == 0
        and cpu0.get("direct_over_rho_wall_ratio") == 1.3543032952267562,
        "Stage-93 autolab CPU-0 binding changed",
    )
    documents = SCOREBOARD.read_text() + STATUS.read_text()
    for marker in ("1.354303", "1.340180", "863,047,680", "1.0188x"):
        require(marker in documents, f"Stage-93 documents are missing {marker!r}")
    return {
        "ledger": {
            "path": str(LEDGER.relative_to(REPO)),
            "sha256": custody.sha256_file(LEDGER, "IC boundary ledger"),
        },
        "scoreboard": {
            "path": str(SCOREBOARD.relative_to(REPO)),
            "sha256": custody.sha256_file(SCOREBOARD, "IC boundary scoreboard"),
        },
        "gate_status": {
            "path": str(STATUS.relative_to(REPO)),
            "sha256": custody.sha256_file(STATUS, "gate status"),
        },
        "autolab_protocol": {
            "path": str(AUTOLAB_PROTOCOL.relative_to(REPO)),
            "sha256": custody.sha256_file(AUTOLAB_PROTOCOL, "autolab protocol"),
        },
    }


def compose() -> dict[str, Any]:
    predecessor = load(STAGE90, "Stage-90 predecessor")
    require(
        predecessor.get("schema") == "koblitz_stage90_current_gate_audit.v1",
        "Stage-90 predecessor changed",
    )
    current = stage92.verify(STAGE92.parent)
    require(current == load(STAGE92, "Stage-92 result"), "Stage-92 binding changed")
    result = copy.deepcopy(predecessor)
    result.update(
        {
            "schema": SCHEMA,
            "status": "current_selected_panel_unknown_and_cpu0_archived_gates_incomplete",
            "predecessor_stage90": {
                "path": str(STAGE90.relative_to(REPO)),
                "sha256": custody.sha256_file(STAGE90, "Stage-90 audit"),
            },
            "stage92_binding": {
                "path": str(STAGE92.relative_to(REPO)),
                "sha256": custody.sha256_file(STAGE92, "Stage-92 result"),
                "verification": current,
            },
            "boundary_ledger_binding": validate_documents(),
            "all_seven_gates_passed": False,
            "independent_external_reproduction_satisfied": False,
            "licensed_magma_complete": False,
            "full_cost_gate_passed": False,
            "koblitz_index_calculus_sota": False,
            "overall": "selected n=53 five-run, optimized unknown-scalar, and current CPU-0 evidence are archived; median/fresh-build crossover, licensed Magma, complete backend scaling, and unaffiliated review remain incomplete",
        }
    )
    single = {
        "affinity": current["affinity"],
        "source_commit": current["source_commit"],
        "query_specialized_n53_pair_batch": current["summary"][
            "query_specialized_n53_pair_batch"
        ],
        "support_queries": current["summary"]["support_queries"],
        "relations": current["summary"]["admitted_relations"],
        "direct_wall_seconds": current["direct_process"]["wall_seconds"],
        "direct_core_seconds": current["direct_process"]["total_core_seconds"],
        "direct_peak_rss_bytes": current["direct_process"]["peak_rss_bytes"],
        "rho_wall_seconds": current["rho_process"]["wall_seconds"],
        "direct_over_rho_wall_ratio": current["comparison"][
            "direct_over_rho_wall_ratio"
        ],
        "direct_over_rho_core_ratio": current["comparison"][
            "direct_core_over_rho_core_ratio"
        ],
        "fresh_build_plus_direct_over_rho_wall_ratio": current["comparison"][
            "fresh_build_plus_direct_over_rho_wall_ratio"
        ],
    }
    result["current_measurements"]["n53_same_target"]["single_core"] = single
    gate1 = result["gates"]["1_full_cost_accounting"]
    gate1["status"] = "partial_current_selected_unknown_cpu0_archived_magma_missing"
    gate1["proved"].append(
        "Stage-92 archives current CPU-0 clean build, direct/rho children, wall, core-seconds, process/tree RSS, setup, collection, LA, validation, and affinity receipts"
    )
    gate3 = result["gates"]["3_required_resource_fields"]
    gate3["status"] = "partial_current_parallel_and_single_core_archived_magma_missing"
    gate3["proved"].append(
        "Stage-92 pins parent and child to CPU 0 and retains selected direct 9.405716 wall / 6.868316 core-seconds / 863,047,680 B RSS versus rho 6.945058 wall"
    )
    gate3["limitations"] = [
        "SAT conflicts are inapplicable to the exact pair-support selected arm and remain null",
        "licensed Magma wall, core-seconds, RSS, and conflict/split fields remain absent",
    ]
    gate6 = result["gates"]["6_automorphism_rho"]
    gate6["status"] = "n41_online_crossover_n53_two_of_five_wins_cpu0_and_median_losses"
    gate6["proved"].append(
        "Stage-92 selected CPU-0 direct remains 1.354303 times rho wall and 1.340180 times rho core-seconds; fresh build plus direct remains 14.691007 times rho"
    )
    gate6["limitations"] = [
        "n=41 crossover is online only and assumes a factor-base log database",
        "n=53 four-thread wins are not repeatable under every selected run and the five-run median remains above one",
        "selected CPU-0 direct remains slower than rho in wall and core-seconds",
        "all fresh-build comparisons remain far above one rho target",
        "all improvements are finite constants and do not change the asymptotic exponent",
    ]
    return result


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-93 output must be new")
    output.mkdir()
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
    committed = load(output / "audit.json", "committed Stage-93 audit")
    require(committed == compose(), "committed Stage-93 audit differs from source evidence")
    seal = load(output / "audit-seal.json", "Stage-93 audit seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == SEAL_SCHEMA
        and claimed == custody.canonical_sha256(payload),
        "Stage-93 audit seal changed",
    )
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    require(
        inventory == seal["inventory"]
        and custody.canonical_sha256(inventory) == seal["inventory_sha256"],
        "Stage-93 audit inventory changed",
    )
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
        value = (
            freeze(args.output.resolve())
            if args.command == "compose"
            else verify(args.output.resolve())
        )
        print(json.dumps(value, indent=2, sort_keys=True))
    except (
        OSError,
        ValueError,
        KeyError,
        AuditError,
        stage92.VerificationError,
        custody.PhaseBError,
    ) as error:
        raise SystemExit(f"stage93-audit: {error}")


if __name__ == "__main__":
    main()

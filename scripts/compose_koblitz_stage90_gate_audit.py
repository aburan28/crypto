#!/usr/bin/env python3
"""Compose the current seven-gate audit from the durable Stage-89 panel."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import verify_koblitz_stage89_current_panel as stage89


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE73 = STAGE / "stage-73-current-gate-audit-20260912/audit.json"
STAGE89 = STAGE / "stage-89-current-selected-panel-20260912/verification.json"
LEDGER = REPO / "docs/ic/boundary_targets.json"
SCOREBOARD = REPO / "docs/ic/BOUNDARY_TARGETS.md"
STATUS = STAGE / "GATE_STATUS.md"
AUTOLAB_PROTOCOL = REPO / "research/sat_factor_base_review_20260908/autolab/protocol.json"
DEFAULT_OUTPUT = STAGE / "stage-90-current-gate-audit-20260912"
SCHEMA = "koblitz_stage90_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage90_current_gate_audit_seal.v1"


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
    require(
        current.get("verdict")
        == "N41_ONLINE_CROSSOVER_N53_TWO_OF_FIVE_PROCESS_WINS_MEDIAN_NEAR_LOSS",
        "Stage-90 ledger verdict changed",
    )
    require(
        metrics.get("n53_selected_median_direct_over_rho_wall_ratio")
        == 1.0187756063132052,
        "Stage-90 selected median changed",
    )
    require(
        metrics.get("n53_selected_wins") == 2
        and metrics.get("n53_selected_losses") == 3,
        "Stage-90 selected win/loss panel changed",
    )
    require(
        metrics.get("n53_unknown_scalar_direct_over_rho_wall_ratio")
        == 1.279112162549439,
        "Stage-90 unknown-scalar ratio changed",
    )
    require(
        records["factor_base"]["current"]["metrics"].get("selection_uses_scalar_labels")
        is False
        and records["factor_base"]["current"]["metrics"].get(
            "target_subgroup_enumerated"
        )
        is False,
        "Stage-90 factor-base boundary changed",
    )
    protocol = load(AUTOLAB_PROTOCOL, "IC boundary autolab protocol")
    n53 = protocol["beats"]["koblitz.vs_rho.n53_probe"]
    require(
        n53.get("fixed_process_repetitions_required") == 5,
        "Stage-90 autolab repetition gate changed",
    )
    require(
        n53.get("selected_default", {}).get("specialized_n53_pair_batch") is True
        and n53.get("selected_default", {}).get(
            "specialized_n53_support_expansion"
        )
        is False
        and n53.get("selected_default", {}).get("specialized_n53_pair_sums")
        is False,
        "Stage-90 autolab selected stack changed",
    )
    documents = SCOREBOARD.read_text() + STATUS.read_text()
    for marker in (
        "1.0188x",
        "two wins",
        "1.279112",
        "Satisfied for finite degrees 23, 31, 41, and 53",
    ):
        require(marker in documents, f"Stage-90 documents are missing {marker!r}")
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
    current = stage89.verify(STAGE89.parent)
    require(current == load(STAGE89, "Stage-89 result"), "Stage-89 binding changed")
    predecessor = load(STAGE73, "Stage-73 predecessor")
    require(
        predecessor.get("schema") == "koblitz_stage73_current_gate_audit.v1",
        "Stage-73 predecessor changed",
    )
    result = copy.deepcopy(predecessor)
    panel = current["selected_panel"]
    selected = current["current_selected"]
    unknown = current["unknown_scalar"]
    summary = selected["summary"]
    process = selected["process"]
    rho = selected["rho_process"]
    result.update(
        {
            "schema": SCHEMA,
            "status": "current_selected_panel_unknown_archived_gates_incomplete",
            "predecessor_stage73": {
                "path": str(STAGE73.relative_to(REPO)),
                "sha256": custody.sha256_file(STAGE73, "Stage-73 audit"),
            },
            "stage89_binding": {
                "path": str(STAGE89.relative_to(REPO)),
                "sha256": custody.sha256_file(STAGE89, "Stage-89 result"),
                "verification": current,
            },
            "boundary_ledger_binding": validate_documents(),
            "all_seven_gates_passed": False,
            "independent_external_reproduction_satisfied": False,
            "licensed_magma_complete": False,
            "full_cost_gate_passed": False,
            "koblitz_index_calculus_sota": False,
            "overall": "five selected n=53 runs and the optimized unknown-scalar result are archived; selected median, fresh-build, refreshed single-core, licensed Magma, and unaffiliated review gates remain incomplete",
        }
    )
    result["current_measurements"]["n53_same_target"] = {
        "n": 53,
        "factor_base_points": 9964,
        "orbit_columns": 94,
        "published_q": summary["published_q"],
        "preferred_query_mode": summary["query_mode"],
        "selected_source_commit": selected["source_commit"],
        "selected_run_count": panel["run_count"],
        "selected_ratios": panel["ratios"],
        "selected_wins": panel["wins"],
        "selected_losses": panel["losses"],
        "selected_median_direct_over_rho_wall_ratio": panel[
            "median_direct_over_rho_wall_ratio"
        ],
        "repeatable_whole_process_crossover": False,
        "current_direct_wall_seconds": process["wall_seconds"],
        "current_direct_core_seconds": process["total_core_seconds"],
        "current_direct_peak_rss_bytes": process["peak_rss_bytes"],
        "current_rho_wall_seconds": rho["wall_seconds"],
        "current_direct_over_rho_wall_ratio": selected[
            "direct_over_rho_wall_ratio"
        ],
        "current_fresh_build_plus_direct_over_rho": selected["verification"][
            "fresh_build_plus_specialized_pair_batch_over_rho_wall_ratio"
        ],
        "support_table_allocated_bytes": selected["base"][
            "support_table_allocated_bytes"
        ],
        "support_queries": summary["support_queries"],
        "relations": summary["admitted_relations"],
        "query_specialized_n53_pair_batch": summary[
            "query_specialized_n53_pair_batch"
        ],
        "unknown_scalar": {
            "published_q": unknown["verification"]["published_q"],
            "recovered_scalar": unknown["verification"]["recovered_scalar"],
            "relations": unknown["verification"]["relations"],
            "support_queries": unknown["summary"]["support_queries"],
            "direct_wall_seconds": unknown["direct_process"]["wall_seconds"],
            "rho_wall_seconds": unknown["rho_process"]["wall_seconds"],
            "direct_over_rho_wall_ratio": unknown["comparison"][
                "direct_over_rho_wall_ratio"
            ],
            "target_scalar_constructed_or_supplied": False,
            "factor_base_logs_known_by_construction": False,
        },
        "single_core": {
            "current_selected_stack": False,
            "predecessor_stage72_retained": True,
            "direct_over_rho_wall_ratio": 1.993515234548992,
            "direct_over_rho_core_ratio": 1.978675797601307,
        },
    }
    gate1 = result["gates"]["1_full_cost_accounting"]
    gate1["status"] = "partial_current_selected_and_unknown_archived_magma_missing"
    gate1["proved"].extend(
        [
            "Stage-89 archives five clean-build selected-stack direct/rho process pairs with workflow, artifact, executable-seal, wall, core-seconds, process RSS, and tree RSS bindings",
            "Stage-89 archives the current optimized public unknown-scalar clean-build/direct/rho process series",
        ]
    )
    gate3 = result["gates"]["3_required_resource_fields"]
    gate3["status"] = "partial_current_parallel_archived_current_single_core_and_magma_missing"
    gate3["proved"].append(
        "Stage-89 retains current selected direct/rho wall, total core-seconds, RSS, operation counts, setup, collection, LA, and validation timings"
    )
    gate3["limitations"] = [
        "SAT conflicts are inapplicable to the exact pair-support selected arm and remain null",
        "the forced-single-core receipt predates the current specialized query batch",
        "licensed Magma resource fields remain absent",
    ]
    gate5 = result["gates"]["5_unknown_scalar"]
    gate5["status"] = "finite_n23_n31_n41_n53_unknown_complete_current_n53_archived"
    gate5["proved"].extend(
        [
            "Stage-89 archives the optimized n=53 public hash target with no supplied scalar or factor-base logs",
            "Stage-89 direct IC and signed-Frobenius rho recover 7892094459170 and verify [d]G=Q",
        ]
    )
    gate6 = result["gates"]["6_automorphism_rho"]
    gate6["status"] = "n41_online_crossover_n53_two_of_five_wins_median_and_full_build_losses"
    gate6["proved"] = [
        "Stage-39 n=41 online IC/rho ratio is 0.285545560549 while amortized IC remains 5.025825 times slower",
        "Stage-89 archives five selected n=53 whole-process ratios with two direct wins and three losses",
        "Stage-89 selected median and current-head direct/rho ratios are both 1.018775606313",
        "Stage-89 current fresh-build plus direct remains 20.262763701788 times rho",
        "Stage-89 optimized public unknown-scalar direct remains 1.279112162549 times rho",
    ]
    gate6["limitations"] = [
        "n=41 crossover is online only and assumes a factor-base log database",
        "n=53 whole-process wins are not repeatable under every fixed selected run and the five-run median remains above one",
        "the current selected-stack forced-single-core comparison remains absent",
        "all improvements are finite constants and do not change the asymptotic exponent",
    ]
    gate7 = result["gates"]["7_external_reproduction_and_novelty"]
    gate7["status"] = "missing_requests_updated_with_stage89_panel"
    gate7["proved"] = [
        "project and upstream reproduction requests contain the selected five-run panel, current source pin, and unknown-scalar evidence",
        "all project-authored archives and CI remain explicitly non-independent",
    ]
    gate7["limitations"] = [
        "no unaffiliated sealed reproduction received",
        "no unaffiliated source-pinned novelty/correctness assessment received",
    ]
    return result


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-90 output must be new")
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
    committed = load(output / "audit.json", "committed Stage-90 audit")
    require(committed == compose(), "committed Stage-90 audit differs from source evidence")
    seal = load(output / "audit-seal.json", "Stage-90 audit seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == SEAL_SCHEMA
        and claimed == custody.canonical_sha256(payload),
        "Stage-90 audit seal changed",
    )
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    require(
        inventory == seal["inventory"]
        and custody.canonical_sha256(inventory) == seal["inventory_sha256"],
        "Stage-90 audit inventory changed",
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
        stage89.VerificationError,
        custody.PhaseBError,
    ) as error:
        raise SystemExit(f"stage90-audit: {error}")


if __name__ == "__main__":
    main()

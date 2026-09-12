#!/usr/bin/env python3
"""Compose Stages 39, 40, and 42 into the current seven-gate audit."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

import compose_koblitz_stage36_gate_audit as stage36
import run_koblitz_blind_pdp_phase_b as custody
import verify_koblitz_stage39_result as stage39
import verify_koblitz_stage40_result as stage40
import verify_koblitz_stage42_result as stage42


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE36 = STAGE / "stage-36-current-gate-audit-20260912/audit.json"
STAGE39 = STAGE / "stage-39-fixed-algebraic-result-20260912/verification.json"
STAGE40 = STAGE / "stage-40-parallel-precompute-result-20260912/verification.json"
STAGE42 = STAGE / "stage-42-n53-same-target-result-20260912/verification.json"
LEDGER = REPO / "docs/ic/boundary_targets.json"
SCOREBOARD = REPO / "docs/ic/BOUNDARY_TARGETS.md"
DEFAULT_OUTPUT = STAGE / "stage-43-current-gate-audit-20260912"
SCHEMA = "koblitz_stage43_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage43_current_gate_audit_seal.v1"


class AuditError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def binding(path: Path, value: dict[str, Any]) -> dict[str, Any]:
    require(value == load(path, path.name), f"{path.name} binding changed")
    return {"path": str(path.relative_to(REPO)), "sha256": custody.sha256_file(path, path.name), "verification": value}


def validate_boundary_ledger() -> dict[str, Any]:
    ledger = load(LEDGER, "IC boundary ledger")
    require(ledger.get("schema_version") == 2, "IC boundary ledger schema changed")
    records = ledger["regimes"]["koblitz"]["records"]
    require(records["factor_base"]["current"]["metrics"]["n"] == 53, "n=53 factor-base row is stale")
    require(records["rank"]["current"]["metrics"]["terminal_rank"] == 95, "n=53 rank row is stale")
    require(records["end_to_end_dlp"]["current"]["metrics"]["recovered_d"] == 476811900269, "n=53 DLP row is stale")
    vs_rho = records["vs_rho"]["current"]
    require(vs_rho.get("verdict") == "N41_ONLINE_CHARGED_CROSSOVER_N53_WALL_LOSS", "vs-rho ledger verdict changed")
    require(vs_rho["metrics"]["n53_direct_over_rho_wall_ratio"] == 12.131346003064046, "n=53 ledger ratio changed")
    scoreboard = SCOREBOARD.read_text()
    for marker in ("9,964 points", "rank 95", "12.131x loss", "same-target whole-process"):
        require(marker in scoreboard, f"human boundary scoreboard is missing {marker!r}")
    return {
        "ledger": {"path": str(LEDGER.relative_to(REPO)), "sha256": custody.sha256_file(LEDGER, "IC boundary ledger")},
        "scoreboard": {"path": str(SCOREBOARD.relative_to(REPO)), "sha256": custody.sha256_file(SCOREBOARD, "IC boundary scoreboard")},
    }


def compose() -> dict[str, Any]:
    predecessor = stage36.verify()
    fixed = stage39.verify_result()
    parallel = stage40.verify_result()
    n53 = stage42.verify(STAGE42.parent)
    ledger_binding = validate_boundary_ledger()
    require(predecessor == load(STAGE36, "Stage-36 audit"), "Stage-36 audit binding changed")
    result = copy.deepcopy(predecessor)
    result.update({
        "schema": SCHEMA,
        "status": "current_evidence_audited_n53_same_target_loss_gates_incomplete",
        "predecessor_stage36": {"path": str(STAGE36.relative_to(REPO)), "sha256": custody.sha256_file(STAGE36, "Stage-36 audit")},
        "stage39_binding": binding(STAGE39, fixed),
        "stage40_binding": binding(STAGE40, parallel),
        "stage42_binding": binding(STAGE42, n53),
        "boundary_ledger_binding": ledger_binding,
        "current_measurements": {
            "n41_fixed_algebraic": {
                "factor_base_discovery_target_samples": 0,
                "factor_base_discovery_scalar_labels": 0,
                "precompute_seconds": fixed["verification"]["precompute_seconds"],
                "ic_over_rho_online_wall_ratio": fixed["verification"]["ic_over_rho_online_wall_ratio"],
                "ic_over_rho_amortised_wall_ratio": fixed["verification"]["ic_over_rho_amortised_wall_ratio"],
                "fresh_build_plus_science_over_rho": fixed["verification"]["full_available_wall_over_same_targets_rho"],
            },
            "n41_four_core": parallel["comparison"],
            "n53_same_target": {
                "published_q": n53["verification"]["published_q"],
                "factor_base_points": n53["verification"]["factor_base_points"],
                "orbit_columns": n53["verification"]["orbit_columns"],
                "relations": n53["verification"]["relations"],
                "direct_wall_seconds": n53["verification"]["direct_wall_seconds"],
                "rho_wall_seconds": n53["verification"]["rho_wall_seconds"],
                "direct_over_rho_wall_ratio": n53["verification"]["direct_over_rho_wall_ratio"],
                "fresh_build_plus_direct_over_rho": n53["hosted"]["comparison"]["fresh_build_plus_direct_over_rho_wall_ratio"],
            },
        },
        "all_seven_gates_passed": False,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "overall": "n=41 fixed-algebraic online crossover and n=53 same-target construction/rank are verified; amortized/full cost, Magma, and external gates remain false",
    })

    gate1 = result["gates"]["1_full_cost_accounting"]
    gate1["proved"].extend([
        "Stage-39 fixed algebraic construction, relations, LA, unknown-scalar descent, rho, build, core-seconds, wall, and tree RSS are charged",
        "Stage-40 charges the complete one/four/four/one parallel control and reports wall speedup plus added CPU cost",
        "Stage-42 charges clean build and same-target n=53 direct/rho arms with watchdog, CPU, wall, and process-tree RSS",
    ])
    gate1["missing"] = ["licensed Magma execution resources", "acquisition of the preinstalled operating system and Rust toolchain"]

    gate3 = result["gates"]["3_required_resource_fields"]
    gate3["proved"].append("Stage-42 n=53 retains single-process wall, total core-seconds, direct/rho peak RSS, support/query operations, build resources, and workflow wall")

    gate4 = result["gates"]["4_scaling"]
    gate4["status"] = "finite_n31_n41_n53_end_to_end_and_n59_pdp_complete"
    gate4["proved"].extend([
        "Stage-42 constructs, ranks, solves, and verifies a same-target known-answer n=53 direct instance",
        "Stage-42 measures a 9,964-point/94-column n=53 base and exact rho comparison",
    ])
    gate4["limitations"] = ["n=59 remains PDP-only", "all completed DLP controls are public synthetic toy subgroups", "finite measurements do not establish asymptotic scaling"]

    gate5 = result["gates"]["5_unknown_scalar"]
    gate5["status"] = "finite_n23_n31_n41_unknown_complete_n53_known_answer_only"
    gate5["proved"].append("Stage-39 n=41 uses zero factor-base discovery target samples or scalar labels and derives every log from relations")
    gate5["limitations"] = ["n=53 is explicit public known-answer validation, not an unknown-scalar run", "finite public synthetic controls only"]

    gate6 = result["gates"]["6_automorphism_rho"]
    gate6["status"] = "n41_online_crossover_n41_amortised_loss_n53_same_target_wall_loss"
    gate6["proved"] = [
        "Stage-39 n=41 online IC/rho ratio is 0.285545560549 (3.502068 times faster after precomputation)",
        "Stage-39 n=41 amortized IC remains 5.025825 times slower and fresh build plus science remains 110.791382 times rho",
        "Stage-40 four-core n=41 precompute wall is 2.655626 times faster while process CPU costs 1.286948 times the one-core median",
        "Stage-40 four-core amortized IC remains 3.605523 times slower than same-instance rho",
        "Stage-42 n=53 direct and packed signed-Frobenius rho use exactly the same public target and recovered scalar",
        "Stage-42 n=53 direct remains 12.131346 times slower by whole-process wall and fresh build plus direct remains 37.698325 times rho",
    ]
    gate6["limitations"] = [
        "n=41 crossover is online only and assumes the factor-base log database exists",
        "n=41 amortized and full available costs do not cross",
        "n=53 same-target whole-process direct does not cross",
        "parallel n=41 latency improvement spends more total CPU",
        "all improvements are finite constants and do not change the asymptotic exponent",
    ]

    gate7 = result["gates"]["7_external_reproduction_and_novelty"]
    gate7["proved"] = [
        "project issue 97 contains Stage-39, Stage-40, and Stage-42 source pins, artifacts, commands, and gate boundaries",
        "upstream EC-Index-Calculus-Benchmarks issue 1 contains the same reproduction request and current n=53 result",
    ]
    gate7["missing"] = ["unaffiliated sealed reproduction", "source-pinned external novelty and correctness verdict"]
    gate7["status"] = "missing"
    return result


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-43 output must be new")
    output.mkdir()
    custody.write_json_new(output / "audit.json", compose())
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "audit_frozen", "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "audit-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "audit.json", "committed Stage-43 audit")
    require(committed == compose(), "committed Stage-43 audit differs from source evidence")
    seal = load(output / "audit-seal.json", "Stage-43 audit seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-43 seal is invalid")
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-43 inventory changed")
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
        raise SystemExit(f"stage43-audit: {error}")


if __name__ == "__main__":
    main()

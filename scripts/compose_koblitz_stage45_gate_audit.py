#!/usr/bin/env python3
"""Compose the Stage 44 parallel n=53 result into the current gate audit."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import verify_koblitz_stage44_result as stage44


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE43 = STAGE / "stage-43-current-gate-audit-20260912/audit.json"
STAGE44 = STAGE / "stage-44-n53-parallel-result-20260912/verification.json"
LEDGER = REPO / "docs/ic/boundary_targets.json"
SCOREBOARD = REPO / "docs/ic/BOUNDARY_TARGETS.md"
DEFAULT_OUTPUT = STAGE / "stage-45-current-gate-audit-20260912"
SCHEMA = "koblitz_stage45_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage45_current_gate_audit_seal.v1"


class AuditError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def validate_ledger() -> dict[str, Any]:
    ledger = load(LEDGER, "IC boundary ledger")
    current = ledger["regimes"]["koblitz"]["records"]["vs_rho"]["current"]
    require(current.get("verdict") == "N41_ONLINE_CHARGED_CROSSOVER_N53_PARALLEL_WALL_LOSS", "Stage-45 ledger verdict changed")
    metrics = current["metrics"]
    require(metrics["n53_parallel_direct_over_rho_wall_ratio"] == 5.9746498402870625, "Stage-45 n=53 ratio changed")
    require(metrics["n53_parallel_direct_core_seconds_ratio"] == 1.4829302318645805, "Stage-45 n=53 CPU ratio changed")
    scoreboard = SCOREBOARD.read_text()
    for marker in ("5.975x slower", "1.579x faster", "1.483x costlier", "Stage 39, 40, 42, and 44"):
        require(marker in scoreboard, f"Stage-45 scoreboard is missing {marker!r}")
    return {
        "ledger": {"path": str(LEDGER.relative_to(REPO)), "sha256": custody.sha256_file(LEDGER, "IC boundary ledger")},
        "scoreboard": {"path": str(SCOREBOARD.relative_to(REPO)), "sha256": custody.sha256_file(SCOREBOARD, "IC boundary scoreboard")},
    }


def verify_stage43_predecessor() -> dict[str, Any]:
    audit = load(STAGE43, "Stage-43 predecessor audit")
    seal_path = STAGE43.parent / "audit-seal.json"
    seal = load(seal_path, "Stage-43 predecessor seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == "koblitz_stage43_current_gate_audit_seal.v1" and claimed == custody.canonical_sha256(payload), "Stage-43 predecessor seal changed")
    inventory = custody.all_regular_inventory(STAGE43.parent, {"audit-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-43 predecessor inventory changed")
    return audit


def compose() -> dict[str, Any]:
    predecessor = verify_stage43_predecessor()
    parallel = stage44.verify()
    require(parallel == load(STAGE44, "Stage-44 result"), "Stage-44 result binding changed")
    result = copy.deepcopy(predecessor)
    result.update({
        "schema": SCHEMA,
        "status": "current_evidence_audited_n53_parallel_wall_loss_gates_incomplete",
        "predecessor_stage43": {"path": str(STAGE43.relative_to(REPO)), "sha256": custody.sha256_file(STAGE43, "Stage-43 audit")},
        "stage44_binding": {"path": str(STAGE44.relative_to(REPO)), "sha256": custody.sha256_file(STAGE44, "Stage-44 result"), "verification": parallel},
        "boundary_ledger_binding": validate_ledger(),
        "all_seven_gates_passed": False,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "overall": "n=41 fixed-algebraic online crossover and n=53 deterministic parallel latency improvement are verified; amortized/full cost, Magma, and external gates remain false",
    })
    measurement = result["current_measurements"]["n53_same_target"]
    measurement.update({
        "sequential_direct_wall_seconds_stage44": parallel["verification"]["sequential_direct_wall_seconds"],
        "parallel_direct_wall_seconds": parallel["verification"]["parallel_direct_wall_seconds"],
        "rho_wall_seconds_stage44": parallel["verification"]["rho_wall_seconds"],
        "parallel_direct_over_rho_wall_ratio": parallel["verification"]["parallel_over_rho_wall_ratio"],
        "parallel_direct_wall_speedup": parallel["verification"]["parallel_direct_wall_speedup"],
        "parallel_direct_core_seconds_ratio": parallel["verification"]["parallel_direct_core_seconds_ratio"],
        "fresh_build_plus_parallel_over_rho": parallel["comparison"]["fresh_build_plus_parallel_over_rho_wall_ratio"],
    })
    gate1 = result["gates"]["1_full_cost_accounting"]
    gate1["proved"].append("Stage-44 charges clean build, sequential direct, deterministic four-thread direct, rho, parallel overshoot, wall, core-seconds, and process-tree RSS")
    gate3 = result["gates"]["3_required_resource_fields"]
    gate3["proved"].append("Stage-44 reports sequential/parallel process CPU and RSS, query operations, wave/chunk counts, workflow wall, and the exact same rho target")
    gate6 = result["gates"]["6_automorphism_rho"]
    gate6["status"] = "n41_online_crossover_n41_amortised_loss_n53_parallel_same_target_wall_loss"
    gate6["proved"] = [
        "Stage-39 n=41 online IC/rho ratio is 0.285545560549 while amortized IC remains 5.025825 times slower",
        "Stage-40 four-core n=41 amortized IC remains 3.605523 times slower and spends 1.286948 times one-core process CPU",
        "Stage-44 sequential and deterministic four-thread n=53 direct arms emit the same 189 relation hashes and recovered scalar",
        "Stage-44 four-thread n=53 direct is 1.578535 times faster in wall and costs 1.482930 times sequential process CPU",
        "Stage-44 four-thread n=53 direct remains 5.974650 times slower than same-target signed-Frobenius rho",
        "Stage-44 fresh build plus parallel direct remains 25.281709 times rho wall",
    ]
    gate6["limitations"] = [
        "n=41 crossover is online only and assumes a factor-base log database",
        "n=41 amortized and full available costs do not cross",
        "n=53 parallel same-target whole-process direct does not cross",
        "parallel n=53 latency improvement spends 48.3 percent more process CPU",
        "all improvements are finite constants and do not change the asymptotic exponent",
    ]
    return result


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-45 output must be new")
    output.mkdir()
    custody.write_json_new(output / "audit.json", compose())
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "audit_frozen", "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "audit-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "audit.json", "committed Stage-45 audit")
    # A successor audit may advance the mutable boundary ledger. Historical
    # verification therefore authenticates the frozen Stage-45 bytes and seal;
    # the current successor recomputes its own bindings against today's ledger.
    seal = load(output / "audit-seal.json", "Stage-45 audit seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-45 audit seal changed")
    inventory = custody.all_regular_inventory(output, {"audit-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-45 audit inventory changed")
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
        raise SystemExit(f"stage45-audit: {error}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Compose the current Koblitz gates with the completed n59
cofactor-projected public unknown-scalar workflow."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
EVIDENCE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE130 = EVIDENCE / "stage-130-current-gate-audit-20260921"
N59 = REPO / "docs/ic/runs/koblitz-n59-cofactor-projected-complete-20260921.json"
N59_PARAMS = REPO / "docs/ic/params/k1n59-cofactor-projected-l14-full-public-59001.json"
GATE_STATUS = EVIDENCE / "GATE_STATUS.md"
SCHEMA = "koblitz_stage131_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage131_current_gate_audit_seal.v1"
MARKER = "Current through Stage 131"


class Stage131Error(RuntimeError):
    pass


def require(value: bool, message: str) -> None:
    if not value:
        raise Stage131Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be an object")
    return value


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay_stage130() -> dict[str, Any]:
    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts/compose_koblitz_stage130_gate_audit.py"),
            "verify",
            "--output",
            str(STAGE130),
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=True,
    )
    value = json.loads(result.stdout)
    require(
        value.get("schema") == "koblitz_stage130_current_gate_audit.v1",
        "Stage-130 replay changed",
    )
    return value


def compose() -> dict[str, Any]:
    predecessor = replay_stage130()
    run = load(N59, "n59 completion evidence")
    params = load(N59_PARAMS, "n59 completion params")
    require(run.get("schema_version") == 1, "n59 evidence schema changed")
    require(
        run.get("operation") == "koblitz_n59_cofactor_projected_unknown_scalar_completion",
        "n59 evidence identity changed",
    )
    require(
        params.get("name") == "k1n59-cofactor-projected-standard-l14-full-public-59001-stage131",
        "n59 parameter identity changed",
    )
    require(params.get("targets") == [{"public_hash_seed": 59001}], "n59 target changed")
    require(run["params"]["sha256"] == sha256(N59_PARAMS), "n59 parameter pin changed")

    instance = run["instance"]
    factor_base = run["factor_base"]
    relations = run["relation_stream"]
    require(
        instance["target_scalar_constructed_or_supplied"] is False,
        "n59 target gained a constructed scalar",
    )
    require(
        factor_base["kind"] == "cofactor_projected"
        and factor_base["parent"] == {"kind": "standard_subspace", "ell": 14},
        "n59 factor-base recipe changed",
    )
    require(
        factor_base["factor_base_logs_known_by_construction"] is False
        and factor_base["target_subgroup_enumerated"] is False
        and factor_base["scalar_preimages_retained"] is False,
        "n59 factor-base knowledge boundary changed",
    )
    require(
        factor_base["points"] == 16308 and factor_base["projected_columns"] == 8094,
        "n59 factor-base dimensions changed",
    )
    require(
        relations["relations"] == 26796
        and relations["trials"] == 5300000
        and relations["canonical_relations_sha256"]
        == "f6b3234bfa93dab2cd05a3609dabf06ce4e5735f55dbcd1efb051983cfd0fc02",
        "n59 relation stream changed",
    )

    default = run["default_thread"]
    one_worker = run["one_worker"]
    for label, result in (("default", default), ("one-worker", one_worker)):
        require(result["verified"] is True, f"{label} n59 solution is unverified")
        require(
            result["recovered_scalar"] == 17861472351607,
            f"{label} n59 scalar changed",
        )
        require(result["relations"] == 26796, f"{label} relation count changed")
        require(
            result["ic_over_rho_wall_ratio"] > 1,
            f"{label} full-cost comparison boundary changed",
        )
        require(
            result["online_rho_over_descent_ratio"] > 1,
            f"{label} online comparison boundary changed",
        )
    require(
        default["peak_rss_bytes"] == 895205376
        and one_worker["peak_rss_bytes"] == 850935808,
        "n59 memory receipts changed",
    )
    accounting = run["process_accounting"]
    require(
        accounting["processes"] == 4
        and accounting["total_core_seconds"] > 2671
        and accounting["peak_rss_bytes"] == 895205376,
        "n59 process accounting changed",
    )
    require(MARKER in GATE_STATUS.read_text(), "gate-status marker changed")

    gates = dict(predecessor["gates"])
    gates["4_n31_n41_larger_pdp_scaling"] = (
        "satisfied_finite_n59_complete_single_public_target"
    )
    gates["3_single_core_core_memory_conflicts_wall"] = (
        "partial_current_n59_complete_magma_missing"
    )
    gates["5_unknown_scalar_without_constructed_base_logs"] = (
        "satisfied_finite_degrees_23_31_41_53_59"
    )
    gates["6_full_cost_vs_automorphism_rho"] = (
        "failed_current_n59_default_and_one_worker_online_descent_only"
    )
    return {
        "schema": SCHEMA,
        "status": "current_seven_gate_audit_verified",
        "all_seven_gates_passed": False,
        "koblitz_index_calculus_sota": False,
        "claim_boundary": (
            "finite n59 public unknown-scalar completion with a target-independent "
            "cofactor-projected algebraic base; full charged cost loses to rho, so this "
            "is not a Koblitz index-calculus SOTA"
        ),
        "predecessor": {
            "stage130_audit_sha256": sha256(STAGE130 / "audit.json"),
            "stage130_seal_sha256": sha256(STAGE130 / "result-seal.json"),
            "stage130_status": predecessor["status"],
        },
        "evidence_pins": {
            "n59_evidence_sha256": sha256(N59),
            "n59_params_sha256": sha256(N59_PARAMS),
            "gate_status_sha256": sha256(GATE_STATUS),
        },
        "inherited_phase_b_same_instance_matrix": predecessor[
            "inherited_phase_b_same_instance_matrix"
        ],
        "inherited_current_n53": predecessor["inherited_current_n53"],
        "inherited_current_n41": predecessor["inherited_current_n41"],
        "inherited_current_n59_standard_cap": predecessor["inherited_current_n59"],
        "inherited_current_n59_standard_frontier": predecessor[
            "current_n59_standard_frontier"
        ],
        "current_n59_cofactor_projected": run,
        "gates": gates,
        "next_targets": [
            "execute the frozen 160-input packet under licensed Magma F4 with complete resource and terminal receipts",
            "optimize the same n59 public target's factor-base width, relation collector and sparse solve under full-cost accounting",
            "obtain unaffiliated reproduction and a source-pinned novelty/correctness review",
        ],
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
    }


def write_new(path: Path, value: dict[str, Any]) -> None:
    require(not path.exists(), f"refusing to overwrite {path}")
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def build(output: Path) -> dict[str, Any]:
    require(not output.exists(), f"refusing to overwrite {output}")
    output.mkdir(parents=True)
    audit = compose()
    write_new(output / "audit.json", audit)
    write_new(
        output / "result-seal.json",
        {
            "schema": SEAL_SCHEMA,
            "status": "audit_frozen",
            "audit_sha256": sha256(output / "audit.json"),
        },
    )
    return audit


def verify(output: Path) -> dict[str, Any]:
    seal = load(output / "result-seal.json", "seal")
    require(seal.get("schema") == SEAL_SCHEMA, "seal schema changed")
    require(
        sha256(output / "audit.json") == seal.get("audit_sha256"),
        "audit seal changed",
    )
    audit = load(output / "audit.json", "audit")
    require(audit.get("schema") == SCHEMA, "audit schema changed")
    require(audit.get("status") == "current_seven_gate_audit_verified", "audit status changed")
    return audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    build_parser = subparsers.add_parser("build")
    build_parser.add_argument("--output", type=Path, required=True)
    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = (
            build(args.output.resolve())
            if args.command == "build"
            else verify(args.output.resolve(strict=True))
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (
        OSError,
        ValueError,
        KeyError,
        subprocess.CalledProcessError,
        Stage131Error,
    ) as error:
        raise SystemExit(f"stage131-gate-audit: {error}")


if __name__ == "__main__":
    main()

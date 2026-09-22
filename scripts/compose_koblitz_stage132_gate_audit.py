#!/usr/bin/env python3
"""Compose the current Koblitz gates with the selected n59 ell=15
cofactor-projected width optimization."""
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
STAGE131 = EVIDENCE / "stage-131-current-gate-audit-20260921"
OPT = REPO / "docs/ic/runs/koblitz-n59-cofactor-projected-l15-optimization-20260921.json"
PARAMS = REPO / "docs/ic/params/k1n59-cofactor-projected-l15-full-public-59001.json"
GATE_STATUS = EVIDENCE / "GATE_STATUS.md"
SCHEMA = "koblitz_stage132_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage132_current_gate_audit_seal.v1"
MARKER = "Current through Stage 132"


class Stage132Error(RuntimeError):
    pass


def require(value: bool, message: str) -> None:
    if not value:
        raise Stage132Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be an object")
    return value


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay_stage131() -> dict[str, Any]:
    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts/compose_koblitz_stage131_gate_audit.py"),
            "verify",
            "--output",
            str(STAGE131),
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=True,
    )
    value = json.loads(result.stdout)
    require(
        value.get("schema") == "koblitz_stage131_current_gate_audit.v1",
        "Stage-131 replay changed",
    )
    return value


def compose() -> dict[str, Any]:
    predecessor = replay_stage131()
    opt = load(OPT, "n59 ell15 optimization")
    params = load(PARAMS, "n59 ell15 params")
    require(
        opt.get("operation") == "koblitz_n59_cofactor_projected_width_optimization",
        "optimization identity changed",
    )
    require(
        params.get("name") == "k1n59-cofactor-projected-standard-l15-full-public-59001-stage132",
        "parameter identity changed",
    )
    require(params.get("targets") == [{"public_hash_seed": 59001}], "target changed")
    require(opt["params"]["sha256"] == sha256(PARAMS), "parameter pin changed")
    require(
        opt["predecessor"]["sha256"]
        == predecessor["evidence_pins"]["n59_evidence_sha256"],
        "ell14 predecessor pin changed",
    )
    factor_base = opt["factor_base"]
    require(
        factor_base["parent"] == {"kind": "standard_subspace", "ell": 15}
        and factor_base["points"] == 32934
        and factor_base["projected_columns"] == 16344,
        "ell15 factor-base shape changed",
    )
    require(
        factor_base["factor_base_logs_known_by_construction"] is False
        and factor_base["target_subgroup_enumerated"] is False
        and factor_base["scalar_preimages_retained"] is False,
        "factor-base knowledge boundary changed",
    )
    relations = opt["relation_stream"]
    require(
        relations["relations"] == 54749
        and relations["trials"] == 2600000
        and relations["canonical_relations_sha256"]
        == "e2004ac6e81979caf984e4fa745dc1a1ee99d13892a3f60615a5662179e3ad99",
        "ell15 relation stream changed",
    )
    for name in ("default_thread", "one_worker"):
        run = opt[name]
        require(run["verified"] is True, f"{name} solution is unverified")
        require(run["recovered_scalar"] == 17861472351607, f"{name} scalar changed")
        require(
            run["canonical_relations_sha256"] == relations["canonical_relations_sha256"],
            f"{name} relation hash changed",
        )
        require(run["ic_over_rho_wall_ratio"] > 1, f"{name} rho boundary changed")
    comparison = opt["comparison_to_ell14"]
    require(
        comparison["default_thread"]["ic_wall_reduction_fraction"] > 0.04
        and comparison["default_thread"]["whole_cpu_reduction_fraction"] > 0.13,
        "default-thread improvement boundary changed",
    )
    require(
        comparison["one_worker"]["ic_wall_reduction_fraction"] > 0.17
        and comparison["one_worker"]["whole_cpu_reduction_fraction"] > 0.16,
        "one-worker improvement boundary changed",
    )
    require(
        comparison["default_thread"]["rss_multiple"] > 4.8
        and comparison["one_worker"]["rss_multiple"] > 4.7,
        "memory tradeoff boundary changed",
    )
    accounting = opt["process_accounting"]
    require(
        accounting["processes"] == 4
        and accounting["total_core_seconds"] > 2681
        and accounting["peak_rss_bytes"] == 4375724032,
        "optimization process accounting changed",
    )
    require(MARKER in GATE_STATUS.read_text(), "gate-status marker changed")

    gates = dict(predecessor["gates"])
    gates["3_single_core_core_memory_conflicts_wall"] = (
        "partial_current_n59_ell15_complete_magma_missing"
    )
    gates["6_full_cost_vs_automorphism_rho"] = (
        "failed_current_n59_ell15_improved_time_at_4_8x_memory"
    )
    return {
        "schema": SCHEMA,
        "status": "current_seven_gate_audit_verified",
        "all_seven_gates_passed": False,
        "koblitz_index_calculus_sota": False,
        "claim_boundary": (
            "finite same-target n59 ell14-to-ell15 improvement with fully charged time "
            "and memory; full cost remains far behind rho, so this is not a SOTA"
        ),
        "predecessor": {
            "stage131_audit_sha256": sha256(STAGE131 / "audit.json"),
            "stage131_seal_sha256": sha256(STAGE131 / "result-seal.json"),
            "stage131_status": predecessor["status"],
        },
        "evidence_pins": {
            "optimization_sha256": sha256(OPT),
            "params_sha256": sha256(PARAMS),
            "gate_status_sha256": sha256(GATE_STATUS),
        },
        "inherited_phase_b_same_instance_matrix": predecessor[
            "inherited_phase_b_same_instance_matrix"
        ],
        "inherited_current_n53": predecessor["inherited_current_n53"],
        "inherited_current_n41": predecessor["inherited_current_n41"],
        "inherited_current_n59_standard_cap": predecessor[
            "inherited_current_n59_standard_cap"
        ],
        "inherited_current_n59_standard_frontier": predecessor[
            "inherited_current_n59_standard_frontier"
        ],
        "inherited_current_n59_cofactor_projected_ell14": predecessor[
            "current_n59_cofactor_projected"
        ],
        "current_n59_cofactor_projected_ell15": opt,
        "gates": gates,
        "next_targets": [
            "execute the frozen 160-input packet under licensed Magma F4 with complete resource and terminal receipts",
            "optimize the same n59 public target below the selected ell15 full cost without hiding its 4.8x memory tradeoff",
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
    require(sha256(output / "audit.json") == seal.get("audit_sha256"), "audit seal changed")
    audit = compose()
    require(audit == load(output / "audit.json", "audit"), "current audit changed")
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage132Error) as error:
        raise SystemExit(f"stage132-gate-audit: {error}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Compose current gates with the selected n59 witnessed compact table."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
from typing import Any

R = Path(__file__).resolve().parents[1]
E = R / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
S = E / "stage-133-current-gate-audit-20260921"
W = R / "docs/ic/runs/koblitz-n59-witnessed-compact-optimization-20260921.json"
P = R / "docs/ic/params/k1n59-cofactor-projected-l15-witnessed-full-public-59001.json"
G = E / "GATE_STATUS.md"
SC = "koblitz_stage134_current_gate_audit.v1"
SS = "koblitz_stage134_current_gate_audit_seal.v1"
MARK = "Current through Stage 134"


class Error(RuntimeError):
    pass


def req(value: bool, message: str) -> None:
    if not value:
        raise Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    req(isinstance(value, dict), f"{context} must be an object")
    return value


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay() -> dict[str, Any]:
    result = subprocess.run(
        [sys.executable, str(R / "scripts/compose_koblitz_stage133_gate_audit.py"),
         "verify", "--output", str(S)], cwd=R, text=True, capture_output=True, check=True,
    )
    value = json.loads(result.stdout)
    req(value.get("schema") == "koblitz_stage133_current_gate_audit.v1", "Stage-133 replay changed")
    return value


def compose() -> dict[str, Any]:
    predecessor = replay()
    run = load(W, "witnessed compact evidence")
    params = load(P, "witnessed compact params")
    req(run.get("operation") == "koblitz_n59_witnessed_compact_optimization", "evidence identity changed")
    req(params.get("name") == "k1n59-cofactor-projected-l15-witnessed-compact-full-public-59001-stage134", "parameter identity changed")
    req(params.get("targets") == [{"public_hash_seed": 59001}], "target changed")
    req(run["params"]["sha256"] == sha(P), "parameter pin changed")
    factor = run["factor_base"]
    req(factor["projected_columns"] == 16344 and factor["points"] == 32934, "factor base changed")
    req(factor["factor_base_logs_known_by_construction"] is False and factor["target_subgroup_enumerated"] is False and factor["scalar_preimages_retained"] is False, "knowledge boundary changed")
    table = run["pair_table"]
    req(table["tier"] == "witnessed_compact" and table["pair_entries"] == 542340645, "table shape changed")
    req(table["exact_candidate_check"].startswith("stored pair sum"), "exact check changed")
    relations = run["relation_stream"]
    req(relations["relations"] == 54749 and relations["trials"] == 2600000, "relation counts changed")
    req(relations["canonical_relations_sha256"] == "e2004ac6e81979caf984e4fa745dc1a1ee99d13892a3f60615a5662179e3ad99", "relation hash changed")
    for name in ("default_thread", "one_worker"):
        result = run[name]
        req(result["verified"] is True and result["recovered_scalar"] == 17861472351607, f"{name} solution changed")
        req(result["canonical_relations_sha256"] == relations["canonical_relations_sha256"], f"{name} relation hash changed")
        req(result["ic_over_rho_wall_ratio"] > 1, f"{name} rho boundary changed")
    comparison = run["comparison_to_ell15_compact"]
    req(comparison["default_thread"]["ic_wall_reduction_fraction"] > .32 and comparison["default_thread"]["whole_cpu_reduction_fraction"] > .28, "default improvement changed")
    req(comparison["one_worker"]["ic_wall_reduction_fraction"] > .28 and comparison["one_worker"]["whole_cpu_reduction_fraction"] > .27, "one-worker improvement changed")
    req(comparison["default_thread"]["rss_multiple"] > 1.5 and comparison["one_worker"]["rss_multiple"] > 1.5, "memory tradeoff changed")
    accounting = run["process_accounting"]
    req(accounting["processes"] == 4 and accounting["total_core_seconds"] > 2094 and accounting["peak_rss_bytes"] == 6565560320, "process accounting changed")
    req(MARK in G.read_text(), "gate marker changed")
    gates = dict(predecessor["gates"])
    gates["3_single_core_core_memory_conflicts_wall"] = "partial_current_n59_witnessed_complete_magma_missing"
    gates["6_full_cost_vs_automorphism_rho"] = "failed_current_n59_witnessed_improved_time_at_1_5x_ell15_memory"
    return {
        "schema": SC,
        "status": "current_seven_gate_audit_verified",
        "all_seven_gates_passed": False,
        "koblitz_index_calculus_sota": False,
        "claim_boundary": "finite same-target n59 witnessed-compact improvement with exact time and memory charges; full cost remains far behind rho and is not a SOTA",
        "predecessor": {"stage133_audit_sha256": sha(S / "audit.json"), "stage133_seal_sha256": sha(S / "result-seal.json"), "stage133_status": predecessor["status"]},
        "evidence_pins": {"witnessed_evidence_sha256": sha(W), "params_sha256": sha(P), "gate_status_sha256": sha(G)},
        "inherited_phase_b_same_instance_matrix": predecessor["inherited_phase_b_same_instance_matrix"],
        "inherited_current_n53": predecessor["inherited_current_n53"],
        "inherited_current_n41": predecessor["inherited_current_n41"],
        "inherited_current_n59_standard_cap": predecessor["inherited_current_n59_standard_cap"],
        "inherited_current_n59_standard_frontier": predecessor["inherited_current_n59_standard_frontier"],
        "inherited_current_n59_cofactor_projected_ell14": predecessor["inherited_current_n59_cofactor_projected_ell14"],
        "inherited_current_n59_cofactor_projected_ell15_compact": predecessor["inherited_current_n59_cofactor_projected_ell15"],
        "inherited_current_n59_collector_tuning": predecessor["current_n59_collector_tuning"],
        "current_n59_witnessed_compact": run,
        "gates": gates,
        "next_targets": [
            "execute the frozen 160-input packet under licensed Magma F4 with complete resource and terminal receipts",
            "reduce the same n59 witnessed full cost while charging its 6.6 GB peak memory",
            "obtain unaffiliated reproduction and a source-pinned novelty/correctness review",
        ],
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
    }


def write_new(path: Path, value: dict[str, Any]) -> None:
    req(not path.exists(), f"refusing overwrite {path}")
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def build(output: Path) -> dict[str, Any]:
    req(not output.exists(), f"refusing overwrite {output}")
    output.mkdir(parents=True)
    audit = compose()
    write_new(output / "audit.json", audit)
    write_new(output / "result-seal.json", {"schema": SS, "status": "audit_frozen", "audit_sha256": sha(output / "audit.json")})
    return audit


def verify(output: Path) -> dict[str, Any]:
    seal = load(output / "result-seal.json", "seal")
    req(seal.get("schema") == SS, "seal schema changed")
    req(sha(output / "audit.json") == seal.get("audit_sha256"), "audit seal changed")
    audit = compose()
    req(audit == load(output / "audit.json", "audit"), "current audit changed")
    return audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("build", "verify"):
        child = sub.add_parser(command)
        child.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = build(args.output.resolve()) if args.command == "build" else verify(args.output.resolve(strict=True))
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Error) as error:
        raise SystemExit(f"stage134-gate-audit: {error}")


if __name__ == "__main__":
    main()

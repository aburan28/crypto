#!/usr/bin/env python3
"""Compose current gates with the selected cached witnessed n59 builder."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
from typing import Any

R = Path(__file__).resolve().parents[1]
E = R / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
S = E / "stage-134-current-gate-audit-20260921"
C = R / "docs/ic/runs/koblitz-n59-witness-key-cache-optimization-20260921.json"
P = R / "docs/ic/params/k1n59-cofactor-projected-l15-witness-cached-full-public-59001.json"
G = E / "GATE_STATUS.md"
SC = "koblitz_stage135_current_gate_audit.v1"
SS = "koblitz_stage135_current_gate_audit_seal.v1"
MARK = "Current through Stage 135"


class Error(RuntimeError): pass
def req(v: bool, m: str) -> None:
    if not v: raise Error(m)
def load(p: Path, c: str) -> dict[str, Any]:
    v = json.loads(p.read_text()); req(isinstance(v, dict), f"{c} object"); return v
def sha(p: Path) -> str: return hashlib.sha256(p.read_bytes()).hexdigest()


def replay() -> dict[str, Any]:
    x = subprocess.run([sys.executable, str(R / "scripts/compose_koblitz_stage134_gate_audit.py"), "verify", "--output", str(S)], cwd=R, text=True, capture_output=True, check=True)
    v = json.loads(x.stdout); req(v.get("schema") == "koblitz_stage134_current_gate_audit.v1", "Stage-134 replay changed"); return v


def compose() -> dict[str, Any]:
    predecessor = replay(); run = load(C, "cache evidence"); params = load(P, "cache params")
    req(run.get("operation") == "koblitz_n59_witness_key_cache_optimization", "evidence identity changed")
    req(params.get("name") == "k1n59-cofactor-projected-l15-witness-cached-full-public-59001-stage135", "parameter identity changed")
    req(params.get("targets") == [{"public_hash_seed": 59001}], "target changed")
    req(run["params"]["sha256"] == sha(P), "parameter pin changed")
    req(run["factor_base"]["projected_columns"] == 16344 and run["factor_base"]["factor_base_logs_known_by_construction"] is False, "factor base changed")
    table = run["pair_table"]
    req(table["builder"] == "cached_first_pass_packed_keys" and table["cached_key_entries"] == 542340645, "cache shape changed")
    req(table["temporary_cache_bytes"] == 4338725160, "cache charge changed")
    relations = run["relation_stream"]
    req(relations["canonical_relations_sha256"] == "e2004ac6e81979caf984e4fa745dc1a1ee99d13892a3f60615a5662179e3ad99", "relation hash changed")
    for name in ("default_thread", "one_worker"):
        result = run[name]
        req(result["verified"] is True and result["recovered_scalar"] == 17861472351607, f"{name} solution changed")
        req(result["canonical_relations_sha256"] == relations["canonical_relations_sha256"], f"{name} relation hash changed")
        req(result["ic_over_rho_wall_ratio"] > 1, f"{name} rho boundary changed")
    comparison = run["comparison_to_stage134_two_pass"]
    req(comparison["default_thread"]["ic_wall_reduction_fraction"] > .10 and comparison["default_thread"]["build_wall_reduction_fraction"] > .45, "default improvement changed")
    req(comparison["one_worker"]["ic_wall_reduction_fraction"] > .10 and comparison["one_worker"]["build_wall_reduction_fraction"] > .48, "one-worker improvement changed")
    req(comparison["default_thread"]["rss_multiple"] > 1.5 and comparison["one_worker"]["rss_multiple"] > 1.6, "memory tradeoff changed")
    accounting = run["process_accounting"]
    req(accounting["processes"] == 3 and accounting["total_core_seconds"] > 1576 and accounting["peak_rss_bytes"] == 10073538560, "accounting changed")
    req(MARK in G.read_text(), "gate marker changed")
    gates = dict(predecessor["gates"])
    gates["3_single_core_core_memory_conflicts_wall"] = "partial_current_n59_cached_witness_complete_magma_missing"
    gates["6_full_cost_vs_automorphism_rho"] = "failed_current_n59_cached_witness_improved_time_at_10gb_peak"
    return {
        "schema": SC, "status": "current_seven_gate_audit_verified", "all_seven_gates_passed": False, "koblitz_index_calculus_sota": False,
        "claim_boundary": "finite same-target n59 cached witnessed-build improvement with exact time and 10 GB peak charge; full cost remains far behind rho and is not a SOTA",
        "predecessor": {"stage134_audit_sha256": sha(S / "audit.json"), "stage134_seal_sha256": sha(S / "result-seal.json"), "stage134_status": predecessor["status"]},
        "evidence_pins": {"cache_evidence_sha256": sha(C), "params_sha256": sha(P), "gate_status_sha256": sha(G)},
        "inherited_phase_b_same_instance_matrix": predecessor["inherited_phase_b_same_instance_matrix"],
        "inherited_current_n53": predecessor["inherited_current_n53"], "inherited_current_n41": predecessor["inherited_current_n41"],
        "inherited_current_n59_standard_cap": predecessor["inherited_current_n59_standard_cap"], "inherited_current_n59_standard_frontier": predecessor["inherited_current_n59_standard_frontier"],
        "inherited_current_n59_cofactor_projected_ell14": predecessor["inherited_current_n59_cofactor_projected_ell14"],
        "inherited_current_n59_cofactor_projected_ell15_compact": predecessor["inherited_current_n59_cofactor_projected_ell15_compact"],
        "inherited_current_n59_collector_tuning": predecessor["inherited_current_n59_collector_tuning"],
        "inherited_current_n59_witnessed_two_pass": predecessor["current_n59_witnessed_compact"],
        "current_n59_cached_witnessed": run, "gates": gates,
        "next_targets": ["execute the frozen 160-input packet under licensed Magma F4 with complete resource and terminal receipts", "reduce the same n59 cached witnessed full cost without hiding its 10 GB construction peak", "obtain unaffiliated reproduction and a source-pinned novelty/correctness review"],
        "licensed_magma_complete": False, "independent_external_reproduction_satisfied": False, "full_cost_gate_passed": False,
    }


def write_new(p: Path, v: dict[str, Any]) -> None: req(not p.exists(), f"overwrite {p}"); p.write_text(json.dumps(v, indent=2, sort_keys=True) + "\n")
def build(o: Path) -> dict[str, Any]:
    req(not o.exists(), f"overwrite {o}"); o.mkdir(parents=True); a = compose(); write_new(o / "audit.json", a); write_new(o / "result-seal.json", {"schema": SS, "status": "audit_frozen", "audit_sha256": sha(o / "audit.json")}); return a
def verify(o: Path) -> dict[str, Any]:
    s = load(o / "result-seal.json", "seal"); req(s.get("schema") == SS, "seal schema"); req(sha(o / "audit.json") == s.get("audit_sha256"), "seal"); a = load(o / "audit.json", "audit"); req(a.get("schema") == SC, "audit schema"); req(a.get("status") == "current_seven_gate_audit_verified", "audit status"); return a
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest="command", required=True)
    for command in ("build", "verify"): child = sub.add_parser(command); child.add_argument("--output", type=Path, required=True)
    x = p.parse_args()
    try: result = build(x.output.resolve()) if x.command == "build" else verify(x.output.resolve(strict=True)); print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Error) as error: raise SystemExit(f"stage135-gate-audit: {error}")
if __name__ == "__main__": main()

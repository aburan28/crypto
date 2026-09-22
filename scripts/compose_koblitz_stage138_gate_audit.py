#!/usr/bin/env python3
"""Compose current gates with the selected n59 sparse-rank relation tail."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
from typing import Any

R = Path(__file__).resolve().parents[1]
E = R / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
S = E / "stage-137-current-gate-audit-20260922"
T = R / "docs/ic/runs/koblitz-n59-rank-tail-20260922.json"
P = R / "docs/ic/params/k1n59-cofactor-projected-l15-rank-tail-public-59001.json"
G = E / "GATE_STATUS.md"
SC = "koblitz_stage138_current_gate_audit.v1"
SS = "koblitz_stage138_current_gate_audit_seal.v1"
MARK = "Current through Stage 138"


class Error(RuntimeError): pass
def req(v: bool, m: str) -> None:
    if not v: raise Error(m)
def load(p: Path, c: str) -> dict[str, Any]:
    v = json.loads(p.read_text()); req(isinstance(v, dict), f"{c} object"); return v
def sha(p: Path) -> str: return hashlib.sha256(p.read_bytes()).hexdigest()


def replay() -> dict[str, Any]:
    x = subprocess.run(
        [sys.executable, str(R / "scripts/compose_koblitz_stage137_gate_audit.py"), "verify", "--output", str(S)],
        cwd=R, text=True, capture_output=True, check=True,
    )
    v = json.loads(x.stdout)
    req(v.get("schema") == "koblitz_stage137_current_gate_audit.v1", "Stage-137 replay changed")
    return v


def compose() -> dict[str, Any]:
    predecessor = replay(); tail = load(T, "rank-tail evidence"); params = load(P, "rank-tail params")
    req(tail.get("operation") == "koblitz_n59_sparse_rank_targeted_tail", "evidence identity changed")
    req(params.get("name") == "k1n59-cofactor-projected-l15-rank-tail-public-59001-stage138", "parameter identity changed")
    req(params.get("targets") == [{"public_hash_seed": 59001}], "target changed")
    req(tail["params"]["sha256"] == sha(P), "parameter pin changed")
    fb = tail["factor_base"]
    req(fb["points"] == 32934 and fb["projected_columns"] == 16344, "factor-base shape changed")
    req(fb["factor_base_logs_known_by_construction"] is False and fb["target_subgroup_enumerated"] is False and fb["scalar_preimages_retained"] is False, "knowledge boundary changed")
    policy = tail["policy"]
    req(policy["name"] == "public_sparse_rank_tail_64" and policy["rank_columns_per_round"] == 64, "rank policy changed")
    req(policy["observed_attempts"] == 295 and policy["observed_pair_lookups"] == 2950000 and policy["observed_relations"] == 71, "rank-tail charge changed")
    stream = tail["relation_stream"]
    req(stream["uniform_units"] == 15 and stream["uniform_trials"] == 1500000 and stream["uniform_relations"] == 31727, "uniform stream changed")
    req(stream["combined_relations"] == 31798 and stream["combined_relations_sha256"] == "32cb1609d40513a12750254ba291fb2a8a262a6eaeb05e598b28212a53fbf3aa", "combined relation stream changed")
    for name, floor in (("selected_default_thread", 100), ("selected_one_worker", 700)):
        result = tail[name]
        req(result["verified"] is True and result["recovered_scalar"] == 17861472351607, f"{name} solution changed")
        req(result["uniform_relations_sha256"] == stream["uniform_relations_sha256"] and result["targeted_relations_sha256"] == stream["targeted_relations_sha256"] and result["combined_relations_sha256"] == stream["combined_relations_sha256"], f"{name} relation hashes changed")
        req(result["ic_over_rho_wall_ratio"] > floor, f"{name} rho boundary changed")
    reduction = tail["comparison_to_stage137"]["reduction"]
    req(reduction["default_ic_wall_fraction"] > .38 and reduction["default_whole_cpu_fraction"] > .29, "default improvement changed")
    req(reduction["one_worker_ic_wall_fraction"] > .31 and reduction["one_worker_whole_cpu_fraction"] > .31, "one-worker improvement changed")
    frontier = tail["frontier_controls"]["replays"]
    req([(x["policy"], x["uniform_units"], x["status"]) for x in frontier] == [
        ("coverage_only", 23, "solved_all_logs"), ("coverage_only", 22, "solved_all_logs"),
        ("coverage_only", 20, "solved_all_logs"), ("coverage_only", 10, "underdetermined_after_cap"),
        ("least_represented_64", 10, "underdetermined_after_cap"), ("least_represented_64", 15, "solved_all_logs"),
        ("least_represented_64", 12, "solved_all_logs"), ("least_represented_64", 14, "solved_all_logs"),
    ], "frontier changed")
    accounting = tail["process_accounting"]
    req(accounting["retained_processes"] == 19 and accounting["retained_total_core_seconds"] > 3748 and accounting["retained_peak_rss_bytes"] == 10091528192, "process accounting changed")
    req(MARK in G.read_text(), "gate marker changed")
    gates = dict(predecessor["gates"])
    gates["1_all_stage_resource_charging"] = "partial_missing_licensed_magma_and_one_preliminary_receipt"
    gates["3_single_core_core_memory_conflicts_wall"] = "partial_current_n59_rank_tail_complete_magma_missing"
    gates["6_full_cost_vs_automorphism_rho"] = "failed_current_n59_rank_tail_improved_far_behind_rho"
    return {
        "schema": SC, "status": "current_seven_gate_audit_verified",
        "all_seven_gates_passed": False, "koblitz_index_calculus_sota": False,
        "claim_boundary": "finite single-public-target sparse-rank relation-tail improvement; all policy controls are charged, but full cost remains far behind rho and this is not a SOTA",
        "predecessor": {"stage137_audit_sha256": sha(S / "audit.json"), "stage137_seal_sha256": sha(S / "result-seal.json"), "stage137_status": predecessor["status"]},
        "evidence_pins": {"rank_tail_evidence_sha256": sha(T), "params_sha256": sha(P), "gate_status_sha256": sha(G)},
        "inherited_phase_b_same_instance_matrix": predecessor["inherited_phase_b_same_instance_matrix"],
        "inherited_current_n53": predecessor["inherited_current_n53"], "inherited_current_n41": predecessor["inherited_current_n41"],
        "inherited_current_n59_standard_cap": predecessor["inherited_current_n59_standard_cap"],
        "inherited_current_n59_standard_frontier": predecessor["inherited_current_n59_standard_frontier"],
        "inherited_current_n59_cofactor_projected_ell14": predecessor["inherited_current_n59_cofactor_projected_ell14"],
        "inherited_current_n59_cofactor_projected_ell15_compact": predecessor["inherited_current_n59_cofactor_projected_ell15_compact"],
        "inherited_current_n59_collector_tuning": predecessor["inherited_current_n59_collector_tuning"],
        "inherited_current_n59_witnessed_two_pass": predecessor["inherited_current_n59_witnessed_two_pass"],
        "inherited_current_n59_cached_witnessed": predecessor["inherited_current_n59_cached_witnessed"],
        "inherited_current_n59_filter_rejection": predecessor["inherited_current_n59_filter_rejection"],
        "inherited_current_n59_targeted_tail": predecessor["current_n59_targeted_tail"],
        "current_n59_rank_tail": tail, "gates": gates,
        "next_targets": [
            "execute the frozen 160-input packet under licensed Magma F4 with complete resource and terminal receipts",
            "reduce the same n59 full cost while charging factor-base discovery, the 10 GB table and rank-tail policy search",
            "obtain unaffiliated reproduction and a source-pinned novelty/correctness review",
        ],
        "licensed_magma_complete": False, "independent_external_reproduction_satisfied": False, "full_cost_gate_passed": False,
    }


def write_new(p: Path, v: dict[str, Any]) -> None:
    req(not p.exists(), f"overwrite {p}"); p.write_text(json.dumps(v, indent=2, sort_keys=True) + "\n")
def build(o: Path) -> dict[str, Any]:
    req(not o.exists(), f"overwrite {o}"); o.mkdir(parents=True); audit = compose()
    write_new(o / "audit.json", audit)
    write_new(o / "result-seal.json", {"schema": SS, "status": "audit_frozen", "audit_sha256": sha(o / "audit.json")})
    return audit
def verify(o: Path) -> dict[str, Any]:
    seal = load(o / "result-seal.json", "seal"); req(seal.get("schema") == SS, "seal schema")
    req(sha(o / "audit.json") == seal.get("audit_sha256"), "seal changed")
    audit = compose(); req(audit == load(o / "audit.json", "audit"), "current audit changed"); return audit
def main() -> None:
    parser = argparse.ArgumentParser(); sub = parser.add_subparsers(dest="cmd", required=True)
    for command in ("build", "verify"):
        p = sub.add_parser(command); p.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = build(args.output.resolve()) if args.cmd == "build" else verify(args.output.resolve(strict=True))
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Error) as exc:
        raise SystemExit(f"stage138-gate-audit: {exc}")
if __name__ == "__main__": main()

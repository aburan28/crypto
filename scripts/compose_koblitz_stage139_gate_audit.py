#!/usr/bin/env python3
"""Compose current gates with the Stage-138 mechanism-attribution correction."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
from typing import Any

R = Path(__file__).resolve().parents[1]
E = R / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
S = E / "stage-138-current-gate-audit-20260922"
A = R / "docs/ic/runs/koblitz-n59-rank-tail-attribution-correction-20260922.json"
P = R / "docs/ic/params/k1n59-cofactor-projected-l15-coverage-tail-control-public-59001.json"
G = E / "GATE_STATUS.md"
SC = "koblitz_stage139_current_gate_audit.v1"
SS = "koblitz_stage139_current_gate_audit_seal.v1"
MARK = "Current through Stage 139"


class Error(RuntimeError): pass
def req(v: bool, m: str) -> None:
    if not v: raise Error(m)
def load(p: Path, c: str) -> dict[str, Any]:
    v = json.loads(p.read_text()); req(isinstance(v, dict), f"{c} object"); return v
def sha(p: Path) -> str: return hashlib.sha256(p.read_bytes()).hexdigest()


def replay() -> dict[str, Any]:
    x = subprocess.run(
        [sys.executable, str(R / "scripts/compose_koblitz_stage138_gate_audit.py"), "verify", "--output", str(S)],
        cwd=R, text=True, capture_output=True, check=True,
    )
    value = json.loads(x.stdout)
    req(value.get("schema") == "koblitz_stage138_current_gate_audit.v1", "Stage-138 replay changed")
    return value


def compose() -> dict[str, Any]:
    predecessor = replay(); correction = load(A, "attribution correction"); params = load(P, "control params")
    req(correction.get("operation") == "koblitz_n59_rank_tail_attribution_correction", "correction identity changed")
    req(params.get("name") == "k1n59-coverage-prefix15-control", "control identity changed")
    req(params.get("targets") == [{"public_hash_seed": 59001}], "control target changed")
    req("targeted_tail_rank_columns" not in params["collection"], "rank fallback is not disabled")
    req(correction["control_params"]["sha256"] == sha(P), "control parameter pin changed")
    attribution = correction["attribution"]
    req(attribution["selected_attempts"] == 295 and attribution["selected_rank_phase_attempts"] == 0 and attribution["selected_raw_uncovered_phase_attempts"] == 295, "selected phase attribution changed")
    control = correction["rank_disabled_control"]
    req(control["params_rank_columns_per_round"] == 0 and control["targeted_attempts"] == 295 and control["targeted_relations"] == 71, "rank-disabled control changed")
    req(control["combined_relations"] == 31798 and control["combined_relations_sha256"] == "32cb1609d40513a12750254ba291fb2a8a262a6eaeb05e598b28212a53fbf3aa", "control relation stream changed")
    req(control["exact_match_to_stage138_selected_stream"] is True, "selected stream no longer matches control")
    accounting = correction["process_accounting"]
    req(accounting["incremental_processes"] == 2 and accounting["incremental_total_core_seconds"] > 194 and accounting["cumulative_stage138_and_correction_processes"] == 21, "correction accounting changed")
    req(accounting["cumulative_stage138_and_correction_peak_rss_bytes"] == 10091528192, "cumulative RSS changed")
    effect = correction["effect_on_results"]
    req(effect["selected_performance_measurements_changed"] is False and effect["selected_relation_stream_changed"] is False and effect["full_cost_gate_passed"] is False, "result boundary changed")
    req(MARK in G.read_text(), "gate marker changed")
    gates = dict(predecessor["gates"])
    gates["1_all_stage_resource_charging"] = "partial_rank_disabled_attribution_control_charged_magma_and_one_preliminary_receipt_missing"
    gates["3_single_core_core_memory_conflicts_wall"] = "partial_current_n59_coverage_tail_complete_magma_missing"
    gates["6_full_cost_vs_automorphism_rho"] = "failed_current_n59_coverage_tail_improved_far_behind_rho_attribution_corrected"
    return {
        "schema": SC, "status": "current_seven_gate_audit_verified",
        "all_seven_gates_passed": False, "koblitz_index_calculus_sota": False,
        "claim_boundary": "additive Stage-138 mechanism-attribution correction: selected speedup comes from the 15-unit raw-uncovered-column tail, while sparse-rank targeting is a frontier fallback; performance is unchanged and full cost remains far behind rho",
        "predecessor": {"stage138_audit_sha256": sha(S / "audit.json"), "stage138_seal_sha256": sha(S / "result-seal.json"), "stage138_status": predecessor["status"]},
        "evidence_pins": {"attribution_correction_sha256": sha(A), "control_params_sha256": sha(P), "gate_status_sha256": sha(G)},
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
        "inherited_current_n59_targeted_tail": predecessor["inherited_current_n59_targeted_tail"],
        "inherited_current_n59_rank_tail": predecessor["current_n59_rank_tail"],
        "current_n59_attribution_correction": correction, "gates": gates,
        "next_targets": predecessor["next_targets"],
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
    audit = load(o / "audit.json", "audit"); req(audit.get("schema") == SC, "audit schema")
    req(audit.get("status") == "current_seven_gate_audit_verified", "audit status"); return audit
def main() -> None:
    parser = argparse.ArgumentParser(); sub = parser.add_subparsers(dest="cmd", required=True)
    for command in ("build", "verify"):
        p = sub.add_parser(command); p.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = build(args.output.resolve()) if args.cmd == "build" else verify(args.output.resolve(strict=True))
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Error) as exc:
        raise SystemExit(f"stage139-gate-audit: {exc}")
if __name__ == "__main__": main()

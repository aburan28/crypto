#!/usr/bin/env python3
"""Compose current gates with the rejected n59 rank-solve cadence."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
from typing import Any

R = Path(__file__).resolve().parents[1]
E = R / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
S = E / "stage-140-current-gate-audit-20260922"
F = R / "docs/ic/runs/koblitz-n59-rank-solve-interval4-rejection-20260922.json"
P = R / "docs/ic/params/k1n59-cofactor-projected-l15-rank-interval4-control-public-59001.json"
G = E / "GATE_STATUS.md"
SC = "koblitz_stage141_current_gate_audit.v1"
SS = "koblitz_stage141_current_gate_audit_seal.v1"
MARK = "Current through Stage 141"

class Error(RuntimeError): pass
def req(v: bool, m: str) -> None:
    if not v: raise Error(m)
def load(p: Path, c: str) -> dict[str, Any]:
    v = json.loads(p.read_text()); req(isinstance(v, dict), f"{c} object"); return v
def sha(p: Path) -> str: return hashlib.sha256(p.read_bytes()).hexdigest()

def replay() -> dict[str, Any]:
    x = subprocess.run([sys.executable, str(R / "scripts/compose_koblitz_stage140_gate_audit.py"), "verify", "--output", str(S)], cwd=R, text=True, capture_output=True, check=True)
    value = json.loads(x.stdout); req(value.get("schema") == "koblitz_stage140_current_gate_audit.v1", "Stage-140 replay changed"); return value

def compose() -> dict[str, Any]:
    predecessor = replay(); cadence = load(F, "cadence evidence"); params = load(P, "cadence params")
    req(cadence.get("operation") == "koblitz_n59_rank_solve_interval4_rejection", "evidence identity changed")
    req(params.get("name") == "k1n59-rank-interval4-prefix14-control", "parameter identity changed")
    req(params["collection"]["units"] == 14 and params["collection"]["targeted_tail_rank_columns"] == 64 and params["collection"]["targeted_tail_solve_interval"] == 4, "cadence changed")
    req(cadence["candidate_params"]["sha256"] == sha(P), "parameter pin changed")
    patch = R / cadence["source_patch"]["path"]; req(cadence["source_patch"]["sha256"] == sha(patch), "patch pin changed")
    result = cadence["result"]
    req(result["verified"] is True and result["recovered_scalar"] == 17861472351607, "solution changed")
    req(result["targeted_tail"]["pair_lookups"] == 20030000 and result["combined_relations"] == 29982, "relation-tail charge changed")
    req(result["linear_algebra"]["attempts"] == 27 and result["linear_algebra"]["seconds"] < 5, "linear algebra changed")
    req(cadence["comparison"]["versus_selected_stage139"]["ic_wall_ratio"] > 2, "selected comparison changed")
    accounting = cadence["process_accounting"]; req(accounting["retained_processes"] == 1 and accounting["retained_total_core_seconds"] > 676 and accounting["retained_peak_rss_bytes"] == 9350201344, "accounting changed")
    req(cadence["decision"]["status"] == "rejected" and cadence["decision"]["selected_stage140_unchanged"] is True, "selection changed")
    req(MARK in G.read_text(), "gate marker changed")
    gates = dict(predecessor["gates"])
    gates["1_all_stage_resource_charging"] = "partial_rank_solve_cadence_rejection_charged_magma_and_one_preliminary_receipt_missing"
    gates["6_full_cost_vs_automorphism_rho"] = "failed_current_n59_coverage_tail_selected_rank_interval4_rejected"
    return {
        "schema": SC, "status": "current_seven_gate_audit_verified", "all_seven_gates_passed": False, "koblitz_index_calculus_sota": False,
        "claim_boundary": "finite same-target rank-solve-cadence rejection; fewer sparse solves are outweighed by targeted-lookup overshoot, so the Stage-140 selected result is unchanged and remains far behind rho",
        "predecessor": {"stage140_audit_sha256": sha(S / "audit.json"), "stage140_seal_sha256": sha(S / "result-seal.json"), "stage140_status": predecessor["status"]},
        "evidence_pins": {"rank_cadence_evidence_sha256": sha(F), "candidate_params_sha256": sha(P), "candidate_patch_sha256": sha(patch), "gate_status_sha256": sha(G)},
        **{k: predecessor[k] for k in predecessor if k.startswith("inherited_")},
        "inherited_current_n59_first_hit_rejection": predecessor["current_n59_first_hit_rejection"],
        "current_n59_rank_solve_cadence_rejection": cadence, "gates": gates, "next_targets": predecessor["next_targets"],
        "licensed_magma_complete": False, "independent_external_reproduction_satisfied": False, "full_cost_gate_passed": False,
    }

def write_new(p: Path, v: dict[str, Any]) -> None: req(not p.exists(), f"overwrite {p}"); p.write_text(json.dumps(v, indent=2, sort_keys=True) + "\n")
def build(o: Path) -> dict[str, Any]:
    req(not o.exists(), f"overwrite {o}"); o.mkdir(parents=True); audit = compose(); write_new(o / "audit.json", audit); write_new(o / "result-seal.json", {"schema": SS, "status": "audit_frozen", "audit_sha256": sha(o / "audit.json")}); return audit
def verify(o: Path) -> dict[str, Any]:
    seal = load(o / "result-seal.json", "seal"); req(seal.get("schema") == SS, "seal schema"); req(sha(o / "audit.json") == seal.get("audit_sha256"), "seal changed"); audit = load(o / "audit.json", "audit"); req(audit.get("schema") == SC, "audit schema"); req(audit.get("status") == "current_seven_gate_audit_verified", "audit status"); return audit
def main() -> None:
    parser=argparse.ArgumentParser(); sub=parser.add_subparsers(dest="cmd",required=True)
    for command in ("build","verify"): p=sub.add_parser(command); p.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    try:
        result=build(args.output.resolve()) if args.cmd=="build" else verify(args.output.resolve(strict=True)); print(json.dumps(result,indent=2,sort_keys=True))
    except (OSError,ValueError,KeyError,subprocess.CalledProcessError,Error) as exc: raise SystemExit(f"stage141-gate-audit: {exc}")
if __name__=="__main__": main()

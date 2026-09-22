#!/usr/bin/env python3
"""Compose current gates with the rejected n59 first-hit targeted tail."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
from typing import Any

R = Path(__file__).resolve().parents[1]
E = R / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
S = E / "stage-139-current-gate-audit-20260922"
F = R / "docs/ic/runs/koblitz-n59-first-hit-tail-rejection-20260922.json"
P = R / "docs/ic/params/k1n59-cofactor-projected-l15-first-hit-control-public-59001.json"
G = E / "GATE_STATUS.md"
SC = "koblitz_stage140_current_gate_audit.v1"
SS = "koblitz_stage140_current_gate_audit_seal.v1"
MARK = "Current through Stage 140"


class Error(RuntimeError): pass
def req(v: bool, m: str) -> None:
    if not v: raise Error(m)
def load(p: Path, c: str) -> dict[str, Any]:
    v = json.loads(p.read_text()); req(isinstance(v, dict), f"{c} object"); return v
def sha(p: Path) -> str: return hashlib.sha256(p.read_bytes()).hexdigest()


def replay() -> dict[str, Any]:
    x = subprocess.run(
        [sys.executable, str(R / "scripts/compose_koblitz_stage139_gate_audit.py"), "verify", "--output", str(S)],
        cwd=R, text=True, capture_output=True, check=True,
    )
    value = json.loads(x.stdout)
    req(value.get("schema") == "koblitz_stage139_current_gate_audit.v1", "Stage-139 replay changed")
    return value


def compose() -> dict[str, Any]:
    predecessor = replay(); first = load(F, "first-hit evidence"); params = load(P, "first-hit params")
    req(first.get("operation") == "koblitz_n59_first_hit_tail_rejection", "evidence identity changed")
    req(params.get("name") == "k1n59-first-hit32-full-control", "parameter identity changed")
    req(params.get("targets") == [{"public_hash_seed": 59001}], "target changed")
    req(params["collection"].get("targeted_tail_first_hit") is True, "first-hit switch changed")
    req(first["candidate_params"]["sha256"] == sha(P), "parameter pin changed")
    patch = R / first["source_patch"]["path"]
    req(first["source_patch"]["sha256"] == sha(patch), "candidate patch pin changed")
    streams = first["relation_streams"]
    req(streams["uniform_relations"] == 31727 and streams["candidate"]["targeted_relations"] == 58 and streams["candidate"]["combined_relations"] == 31785, "candidate stream changed")
    req(streams["candidate"]["combined_relations_sha256"] == "dcbf0aad8eabce6f22241801f70ade6f4d948e2d1bc67d39a5dd1d2fd140e76c" and streams["candidate_repeated_hash_match"] is True, "candidate hash changed")
    abba = first["matched_default_abba"]
    req(abba["order"] == ["baseline", "candidate", "candidate", "baseline"], "ABBA order changed")
    req(abba["comparison"]["whole_wall_seconds"]["change_fraction"] > .04, "wall rejection changed")
    req(abba["comparison"]["whole_cpu_seconds"]["change_fraction"] > 0, "CPU rejection changed")
    req(abba["comparison"]["targeted_pair_lookups"]["change_fraction"] < -.08, "lookup reduction changed")
    req(all(x["verified"] is True and x["recovered_scalar"] == 17861472351607 for x in abba["runs"]), "solution changed")
    accounting = first["process_accounting"]
    req(accounting["retained_processes"] == 9 and accounting["retained_total_core_seconds"] > 3806 and accounting["retained_peak_rss_bytes"] == 10080059392, "process accounting changed")
    req(first["decision"]["status"] == "rejected" and first["decision"]["selected_stage139_unchanged"] is True, "selection changed")
    req(MARK in G.read_text(), "gate marker changed")
    gates = dict(predecessor["gates"])
    gates["1_all_stage_resource_charging"] = "partial_first_hit_rejection_charged_magma_and_one_preliminary_receipt_missing"
    gates["6_full_cost_vs_automorphism_rho"] = "failed_current_n59_coverage_tail_selected_first_hit_rejected"
    return {
        "schema": SC, "status": "current_seven_gate_audit_verified",
        "all_seven_gates_passed": False, "koblitz_index_calculus_sota": False,
        "claim_boundary": "finite same-target first-hit optimization rejection; pair lookups fall but matched default full wall and CPU regress, so the Stage-139 selected result is unchanged and remains far behind rho",
        "predecessor": {"stage139_audit_sha256": sha(S / "audit.json"), "stage139_seal_sha256": sha(S / "result-seal.json"), "stage139_status": predecessor["status"]},
        "evidence_pins": {"first_hit_evidence_sha256": sha(F), "candidate_params_sha256": sha(P), "candidate_patch_sha256": sha(patch), "gate_status_sha256": sha(G)},
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
        "inherited_current_n59_rank_tail": predecessor["inherited_current_n59_rank_tail"],
        "inherited_current_n59_attribution_correction": predecessor["current_n59_attribution_correction"],
        "current_n59_first_hit_rejection": first, "gates": gates,
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
        raise SystemExit(f"stage140-gate-audit: {exc}")
if __name__ == "__main__": main()

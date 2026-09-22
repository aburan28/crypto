#!/usr/bin/env python3
"""Compose current gates with the charged n59 ell15 collector tuning."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
from typing import Any

R = Path(__file__).resolve().parents[1]
E = R / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
S = E / "stage-132-current-gate-audit-20260921"
T = R / "docs/ic/runs/koblitz-n59-l15-collector-tuning-20260921.json"
G = E / "GATE_STATUS.md"
SC = "koblitz_stage133_current_gate_audit.v1"
SS = "koblitz_stage133_current_gate_audit_seal.v1"
MARK = "Current through Stage 133"


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
        [sys.executable, str(R / "scripts/compose_koblitz_stage132_gate_audit.py"),
         "verify", "--output", str(S)],
        cwd=R, text=True, capture_output=True, check=True,
    )
    value = json.loads(result.stdout)
    req(value.get("schema") == "koblitz_stage132_current_gate_audit.v1", "Stage-132 replay changed")
    return value


def compose() -> dict[str, Any]:
    predecessor = replay()
    tuning = load(T, "collector tuning")
    req(tuning.get("operation") == "koblitz_n59_l15_collector_tuning", "tuning identity changed")
    req(tuning["instance"]["public_target_seed"] == 59001, "target changed")
    sweep = tuning["equal_scan_window_sweep"]
    req(sweep["summands_scanned_per_arm"] == 102400000, "sweep work changed")
    req([x["collection_window"] for x in sweep["arms"]] == [512, 1024, 2048], "window panel changed")
    req(sweep["selected_window"] == 1024, "window selection changed")
    req(sweep["throughput_ratios_vs_selected"]["window512"] < 1, "window512 boundary changed")
    req(sweep["throughput_ratios_vs_selected"]["window2048"] < 1, "window2048 boundary changed")
    rejected = tuning["rejected_scratch_reuse"]
    patch = R / rejected["source_patch"]
    req(sha(patch) == rejected["source_patch_sha256"], "rejected patch pin changed")
    req(rejected["relation_hash_preserved"] is True, "scratch correctness changed")
    req(all(x > 1 for x in rejected["relation_unit_wall_ratios_vs_baseline"]), "scratch timing boundary changed")
    req(rejected["decision"] == "rejected_wall_regression", "scratch decision changed")
    accounting = tuning["new_process_accounting"]
    req(accounting["processes"] == 4 and accounting["total_core_seconds"] > 1078, "process accounting changed")
    req(MARK in G.read_text(), "gate marker changed")
    return {
        "schema": SC,
        "status": "current_seven_gate_audit_verified",
        "all_seven_gates_passed": False,
        "koblitz_index_calculus_sota": False,
        "claim_boundary": "finite n59 ell15 collector tuning with two rejected windows and rejected scratch reuse; selected full-cost result unchanged and not a SOTA",
        "predecessor": {"stage132_audit_sha256": sha(S / "audit.json"), "stage132_seal_sha256": sha(S / "result-seal.json"), "stage132_status": predecessor["status"]},
        "evidence_pins": {"tuning_sha256": sha(T), "rejected_patch_sha256": sha(patch), "gate_status_sha256": sha(G)},
        "inherited_phase_b_same_instance_matrix": predecessor["inherited_phase_b_same_instance_matrix"],
        "inherited_current_n53": predecessor["inherited_current_n53"],
        "inherited_current_n41": predecessor["inherited_current_n41"],
        "inherited_current_n59_standard_cap": predecessor["inherited_current_n59_standard_cap"],
        "inherited_current_n59_standard_frontier": predecessor["inherited_current_n59_standard_frontier"],
        "inherited_current_n59_cofactor_projected_ell14": predecessor["inherited_current_n59_cofactor_projected_ell14"],
        "inherited_current_n59_cofactor_projected_ell15": predecessor["current_n59_cofactor_projected_ell15"],
        "current_n59_collector_tuning": tuning,
        "gates": predecessor["gates"],
        "next_targets": [
            "execute the frozen 160-input packet under licensed Magma F4 with complete resource and terminal receipts",
            "test a directly witnessed compact table on the same n59 target while charging its extra memory",
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
        raise SystemExit(f"stage133-gate-audit: {error}")


if __name__ == "__main__":
    main()

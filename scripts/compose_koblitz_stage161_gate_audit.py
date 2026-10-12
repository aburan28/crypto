#!/usr/bin/env python3
"""Compose or historically verify the Stage-161 seven-gate audit."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE160_AUDIT = GATES / "stage-160-current-gate-audit-20260923"
STAGE161 = GATES / "stage-161-fast-mask-hash-20260923"
STAGE161_VERIFICATION = GATES / "stage-161-fast-hash-verification.json"
GATE_STATUS = GATES / "GATE_STATUS.md"
SCOREBOARD = REPO / "docs/index-calculus-scoreboard.html"

SCHEMA = "koblitz_stage161_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage161_current_gate_audit_seal.v1"


class Stage161AuditError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage161AuditError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be an object")
    return value


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay(script: str, *args: str) -> dict[str, Any]:
    completed = subprocess.run(
        [sys.executable, str(REPO / "scripts" / script), *args],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )
    require(
        completed.returncode == 0,
        f"{script} replay failed ({completed.returncode}): {completed.stderr.strip()}",
    )
    value = json.loads(completed.stdout)
    require(isinstance(value, dict), f"{script} replay must be an object")
    return value


def compose() -> dict[str, Any]:
    stage160 = replay(
        "compose_koblitz_stage160_gate_audit.py",
        "verify",
        "--output",
        str(STAGE160_AUDIT),
    )
    stage161 = replay("verify_koblitz_stage161_fast_hash.py")
    committed_verification = load(STAGE161_VERIFICATION, "Stage-161 verification")
    require(stage161 == committed_verification, "Stage-161 verification output changed")
    result = load(STAGE161 / "result.json", "Stage-161 result")
    seal = load(STAGE161 / "result-seal.json", "Stage-161 result seal")
    require(
        stage160.get("schema") == "koblitz_stage160_current_gate_audit.v1"
        and stage160.get("status") == "current_seven_gate_audit_verified",
        "Stage-160 historical audit changed",
    )
    require(
        stage161.get("status") == "verified"
        and stage161.get("result_sha256") == seal.get("result_sha256")
        and stage161.get("inventory_sha256") == seal.get("inventory_sha256"),
        "Stage-161 sealed result changed",
    )
    require(
        result.get("status") == "complete_fast_mask_hash_true_positive"
        and result["selected_native_f4"]["source_witness_valid"] is True
        and result["comparisons"]["selected_over_stage160"]["wall_ratio"] < 0.823
        and result["comparisons"]["selected_over_same_host_direct_mitm"]["wall_ratio"] > 20,
        "Stage-161 scientific boundary changed",
    )
    require(
        result["campaign_accounting"]["complete_campaign_cost"] is None
        and result["licensed_magma_f4_complete"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["full_cost_gate_passed"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "Stage-161 claim boundary widened",
    )
    gate_text = GATE_STATUS.read_text()
    scoreboard_text = SCOREBOARD.read_text()
    require(
        "Stage 161 trusted-mask hashing" in gate_text
        and "Stage 161 trusted-mask hashing optimization" in gate_text,
        "current gate status lacks Stage 161",
    )
    require(
        'id="phase-b-native-f4-single-target"' in scoreboard_text
        and "still 20.22&times; direct MITM" in scoreboard_text
        and "60.645599" in scoreboard_text,
        "canonical scoreboard lacks Stage 161",
    )
    return {
        "schema": SCHEMA,
        "status": "current_seven_gate_audit_verified",
        "predecessor": {
            "stage160_audit_sha256": sha256(STAGE160_AUDIT / "audit.json"),
            "stage160_seal_sha256": sha256(STAGE160_AUDIT / "result-seal.json"),
            "stage160_status": stage160["status"],
        },
        "stage161": {
            "result_sha256": sha256(STAGE161 / "result.json"),
            "result_seal_sha256": sha256(STAGE161 / "result-seal.json"),
            "verification_sha256": sha256(STAGE161_VERIFICATION),
            "blind_instance_id": stage161["blind_instance_id"],
            "source_commit": stage161["source_commit"],
            "selected_wall_seconds": stage161["selected_wall_seconds"],
            "selected_core_seconds": stage161["selected_core_seconds"],
            "selected_peak_rss_bytes": stage161["selected_peak_rss_bytes"],
            "build_phase_speedup": stage161["build_phase_speedup"],
            "whole_target_speedup": stage161["whole_target_speedup"],
            "same_host_f4_over_mitm_wall": stage161["same_host_f4_over_mitm_wall"],
            "complete_campaign_cost": None,
        },
        "gates": result["gates"],
        "next_targets": [
            "execute the full frozen 160-input packet under licensed Magma F4",
            "run a complete native-F4 panel rather than one selected n59 target",
            "reduce selected native-F4 cost below same-target direct MITM or retain it only as a solver control",
            "preserve complete outer metering for every future development command",
            "obtain unaffiliated reproduction and a source-pinned novelty/correctness review",
        ],
        "gate_status_sha256": sha256(GATE_STATUS),
        "scoreboard_sha256": sha256(SCOREBOARD),
        "all_seven_gates_passed": False,
        "licensed_magma_f4_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
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
    seal = load(output / "result-seal.json", "Stage-161 audit seal")
    require(seal.get("schema") == SEAL_SCHEMA, "Stage-161 audit seal schema changed")
    require(sha256(output / "audit.json") == seal.get("audit_sha256"), "Stage-161 audit seal changed")
    stored = load(output / "audit.json", "Stage-161 audit")
    stage160 = replay(
        "compose_koblitz_stage160_gate_audit.py",
        "verify",
        "--output",
        str(STAGE160_AUDIT),
    )
    stage161 = replay("verify_koblitz_stage161_fast_hash.py")
    committed_verification = load(STAGE161_VERIFICATION, "Stage-161 verification")
    require(stage161 == committed_verification, "Stage-161 verification output changed")
    require(
        stored.get("schema") == SCHEMA
        and stored.get("status") == "current_seven_gate_audit_verified"
        and stored["predecessor"]["stage160_audit_sha256"]
        == sha256(STAGE160_AUDIT / "audit.json")
        and stored["predecessor"]["stage160_seal_sha256"]
        == sha256(STAGE160_AUDIT / "result-seal.json")
        and stored["stage161"]["result_sha256"] == sha256(STAGE161 / "result.json")
        and stored["stage161"]["result_seal_sha256"]
        == sha256(STAGE161 / "result-seal.json")
        and stored["stage161"]["verification_sha256"] == sha256(STAGE161_VERIFICATION)
        and stored["stage161"]["selected_wall_seconds"] == stage161["selected_wall_seconds"]
        and stored["all_seven_gates_passed"] is False
        and stored["koblitz_index_calculus_sota"] is False,
        "historical Stage-161 audit dependencies changed",
    )
    require(
        stage160.get("status") == "current_seven_gate_audit_verified",
        "Stage-160 historical audit changed",
    )
    return stored


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    build_parser = sub.add_parser("build")
    build_parser.add_argument("--output", type=Path, required=True)
    verify_parser = sub.add_parser("verify")
    verify_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        value = (
            build(args.output.resolve())
            if args.command == "build"
            else verify(args.output.resolve(strict=True))
        )
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, Stage161AuditError) as error:
        raise SystemExit(f"stage161-gate-audit: {error}")


if __name__ == "__main__":
    main()

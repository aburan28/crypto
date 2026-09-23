#!/usr/bin/env python3
"""Compose the current seven-gate audit with the Stage-159 F4 supplement."""

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
STAGE128 = GATES / "stage-128-current-gate-audit-20260922"
STAGE159 = GATES / "stage-159-native-f4-single-target-20260922"
STAGE159_VERIFICATION = GATES / "stage-159-native-f4-verification.json"
GATE_STATUS = GATES / "GATE_STATUS.md"
SCOREBOARD = REPO / "docs/index-calculus-scoreboard.html"

SCHEMA = "koblitz_stage159_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage159_current_gate_audit_seal.v1"


class Stage159AuditError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage159AuditError(message)


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
    stage128 = replay(
        "compose_koblitz_stage128_gate_audit.py",
        "verify",
        "--output",
        str(STAGE128),
    )
    stage159 = replay("verify_koblitz_stage159_native_f4.py")
    committed_verification = load(STAGE159_VERIFICATION, "Stage-159 verification")
    require(stage159 == committed_verification, "Stage-159 verification output changed")
    result = load(STAGE159 / "result.json", "Stage-159 result")
    seal = load(STAGE159 / "result-seal.json", "Stage-159 result seal")
    require(
        stage128.get("schema") == "koblitz_stage128_current_gate_audit.v1"
        and stage128.get("status") == "current_seven_gate_audit_verified",
        "Stage-128 historical audit changed",
    )
    require(
        stage159.get("status") == "verified"
        and stage159.get("result_sha256") == seal.get("result_sha256")
        and stage159.get("inventory_sha256") == seal.get("inventory_sha256"),
        "Stage-159 sealed result changed",
    )
    require(
        result.get("status") == "complete_single_target_native_f4_true_positive"
        and result["selected_native_f4"]["source_witness_valid"] is True
        and result["comparisons"]["selected_over_full_fixed_x1"]["wall_ratio"] < 0.488
        and result["comparisons"]["selected_over_same_host_direct_mitm"]["wall_ratio"] > 33,
        "Stage-159 scientific boundary changed",
    )
    require(
        result["campaign_accounting"]["complete_campaign_cost"] is None
        and result["licensed_magma_f4_complete"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["full_cost_gate_passed"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "Stage-159 claim boundary widened",
    )
    gate_text = GATE_STATUS.read_text()
    scoreboard_text = SCOREBOARD.read_text()
    require(
        "Current through Stage 128, plus the additive Stage 159 native-F4" in gate_text
        and "Stage 159 native-F4 single-target supplement" in gate_text,
        "current gate status lacks Stage 159",
    )
    require(
        'id="phase-b-native-f4-single-target"' in scoreboard_text
        and "33.40&times; direct MITM" in scoreboard_text,
        "canonical scoreboard lacks Stage 159",
    )
    return {
        "schema": SCHEMA,
        "status": "current_seven_gate_audit_verified",
        "predecessor": {
            "stage128_audit_sha256": sha256(STAGE128 / "audit.json"),
            "stage128_seal_sha256": sha256(STAGE128 / "result-seal.json"),
            "stage128_status": stage128["status"],
        },
        "stage159": {
            "result_sha256": sha256(STAGE159 / "result.json"),
            "result_seal_sha256": sha256(STAGE159 / "result-seal.json"),
            "verification_sha256": sha256(STAGE159_VERIFICATION),
            "blind_instance_id": stage159["blind_instance_id"],
            "selected_wall_seconds": stage159["selected_wall_seconds"],
            "selected_core_seconds": stage159["selected_core_seconds"],
            "selected_peak_rss_bytes": stage159["selected_peak_rss_bytes"],
            "selected_word_xors": stage159["selected_word_xors"],
            "same_host_f4_over_mitm_wall": stage159["same_host_f4_over_mitm_wall"],
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
    seal = load(output / "result-seal.json", "Stage-159 audit seal")
    require(seal.get("schema") == SEAL_SCHEMA, "Stage-159 audit seal schema changed")
    require(sha256(output / "audit.json") == seal.get("audit_sha256"), "Stage-159 audit seal changed")
    current = compose()
    require(current == load(output / "audit.json", "Stage-159 audit"), "current Stage-159 audit changed")
    return current


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    build_parser = sub.add_parser("build")
    build_parser.add_argument("--output", type=Path, required=True)
    verify_parser = sub.add_parser("verify")
    verify_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        value = build(args.output.resolve()) if args.command == "build" else verify(args.output.resolve(strict=True))
        print(json.dumps(value, indent=2, sort_keys=True))
    except (
        OSError,
        ValueError,
        KeyError,
        subprocess.CalledProcessError,
        Stage159AuditError,
    ) as error:
        raise SystemExit(f"stage159-gate-audit: {error}")


if __name__ == "__main__":
    main()

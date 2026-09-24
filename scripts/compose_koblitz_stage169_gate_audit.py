#!/usr/bin/env python3
"""Compose the current seven-gate audit with the Stage-169 fresh holdout."""

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
PREDECESSOR = GATES / "stage-168-current-gate-audit-20260923"
STAGE = GATES / "stage-169-fresh-target-holdout-20260923"
VERIFICATION = GATES / "stage-169-fresh-target-holdout-verification.json"
GATE_STATUS = GATES / "GATE_STATUS.md"
SCOREBOARD = REPO / "docs/index-calculus-scoreboard.html"
SCHEMA = "koblitz_stage169_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage169_current_gate_audit_seal.v1"


class AuditError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


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
    require(completed.returncode == 0, f"{script} failed: {completed.stderr.strip()}")
    return json.loads(completed.stdout)


def compose() -> dict[str, Any]:
    predecessor = replay(
        "compose_koblitz_stage168_gate_audit.py",
        "verify",
        "--output",
        str(PREDECESSOR),
    )
    verification = replay("verify_koblitz_stage169_fresh_holdout.py")
    require(
        verification == load(VERIFICATION, "Stage-169 verification"),
        "verification changed",
    )
    result = load(STAGE / "result.json", "Stage-169 result")
    seal = load(STAGE / "result-seal.json", "Stage-169 seal")
    require(
        predecessor.get("status") == "current_seven_gate_audit_verified"
        and verification.get("status") == "verified"
        and verification["result_sha256"] == seal["result_sha256"]
        and verification["fresh_target_holdout_complete"] is True
        and result["comparisons"]["f4_over_direct_mitm"]["wall_ratio"] > 6.9,
        "scientific boundary changed",
    )
    require(
        result["preregistration"]["truth_consulted_for_selection"] is False
        and result["campaign_accounting"]["complete_campaign_cost"] is None
        and result["licensed_magma_f4_complete"] is False
        and result["independent_external_reproduction_satisfied"] is False
        and result["full_cost_gate_passed"] is False
        and result["koblitz_index_calculus_sota"] is False,
        "claim boundary widened",
    )
    gate_text = GATE_STATUS.read_text()
    scoreboard = SCOREBOARD.read_text()
    require(
        "Stage 169 preregistered fresh-target holdout" in gate_text
        and "8,261.489665 sequential wall-seconds" in gate_text,
        "gate status lacks Stage 169",
    )
    require(
        "preregistered fresh n = 59 holdout" in scoreboard
        and "20.539819" in scoreboard,
        "scoreboard lacks Stage 169",
    )
    return {
        "schema": SCHEMA,
        "status": "current_seven_gate_audit_verified",
        "predecessor": {
            "stage168_audit_sha256": sha256(PREDECESSOR / "audit.json"),
            "stage168_seal_sha256": sha256(PREDECESSOR / "result-seal.json"),
            "stage168_status": predecessor["status"],
        },
        "stage169": {
            "result_sha256": sha256(STAGE / "result.json"),
            "result_seal_sha256": sha256(STAGE / "result-seal.json"),
            "verification_sha256": sha256(VERIFICATION),
            "blind_instance_id": verification["blind_instance_id"],
            "source_commit": verification["source_commit"],
            "classification": verification["classification"],
            "f4_wall_seconds": verification["f4_wall_seconds"],
            "f4_total_core_seconds": verification["f4_total_core_seconds"],
            "f4_single_core_seconds": None,
            "f4_peak_rss_bytes": verification["f4_peak_rss_bytes"],
            "f4_systems_completed": verification["f4_systems_completed"],
            "f4_word_xors": verification["f4_word_xors"],
            "direct_mitm_wall_seconds": verification["direct_mitm_wall_seconds"],
            "f4_over_direct_mitm_wall": verification[
                "f4_over_direct_mitm_wall"
            ],
            "fresh_target_holdout_complete": True,
            "complete_campaign_cost": None,
        },
        "gates": result["gates"],
        "next_targets": [
            "execute the full frozen 160-input packet under licensed Magma F4",
            "expand preregistered native-F4 holdouts beyond one positive and one negative target",
            "reduce fresh-target F4 full cost relative to direct MITM and automorphism rho",
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
    seal = load(output / "result-seal.json", "Stage-169 audit seal")
    require(seal.get("schema") == SEAL_SCHEMA, "audit seal schema changed")
    require(
        sha256(output / "audit.json") == seal.get("audit_sha256"),
        "audit seal changed",
    )
    current = compose()
    stored = load(output / "audit.json", "Stage-169 audit")
    current_scientific = dict(current)
    stored_scientific = dict(stored)
    for value in (current_scientific, stored_scientific):
        value.pop("gate_status_sha256", None)
        value.pop("scoreboard_sha256", None)
    require(current_scientific == stored_scientific, "current scientific audit changed")
    return current


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("build").add_argument("--output", type=Path, required=True)
    sub.add_parser("verify").add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        value = (
            build(args.output.resolve())
            if args.command == "build"
            else verify(args.output.resolve(strict=True))
        )
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, AuditError) as error:
        raise SystemExit(f"stage169-gate-audit: {error}")


if __name__ == "__main__":
    main()

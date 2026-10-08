#!/usr/bin/env python3
"""Compose the current seven-gate audit with the Stage-166 dense-product result."""

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
PREDECESSOR = GATES / "stage-165-current-gate-audit-20260923"
STAGE = GATES / "stage-166-dense-product-cancellation-20260923"
VERIFICATION = GATES / "stage-166-dense-product-cancellation-verification.json"
GATE_STATUS = GATES / "GATE_STATUS.md"
SCOREBOARD = REPO / "docs/index-calculus-scoreboard.html"
SCHEMA = "koblitz_stage166_current_gate_audit.v1"
SEAL_SCHEMA = "koblitz_stage166_current_gate_audit_seal.v1"


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
        "compose_koblitz_stage165_gate_audit.py",
        "verify",
        "--output",
        str(PREDECESSOR),
    )
    verification = replay("verify_koblitz_stage166_dense_mul.py")
    require(
        verification == load(VERIFICATION, "Stage-166 verification"),
        "verification changed",
    )
    result = load(STAGE / "result.json", "Stage-166 result")
    seal = load(STAGE / "result-seal.json", "Stage-166 seal")
    require(
        predecessor.get("status") == "current_seven_gate_audit_verified"
        and verification.get("status") == "verified"
        and verification["result_sha256"] == seal["result_sha256"]
        and result["comparisons"]["same_binary_sorted_product_control"][
            "selected_wall_speedup"
        ]
        > 1.08
        and result["comparisons"]["selected_over_same_host_direct_mitm"][
            "wall_ratio"
        ]
        > 6.6,
        "scientific boundary changed",
    )
    require(
        result["mechanism"]["uses_target_labels"] is False
        and result["mechanism"]["uses_discrete_log_labels"] is False
        and result["mechanism"]["uses_known_witness"] is False
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
        "Stage 166 dense F4 product cancellation" in gate_text
        and "7,400.453979 sequential wall-seconds" in gate_text,
        "gate status lacks Stage 166",
    )
    require(
        "still 6.63&times; direct MITM" in scoreboard and "19.887753" in scoreboard,
        "scoreboard lacks Stage 166",
    )
    return {
        "schema": SCHEMA,
        "status": "current_seven_gate_audit_verified",
        "predecessor": {
            "stage165_audit_sha256": sha256(PREDECESSOR / "audit.json"),
            "stage165_seal_sha256": sha256(PREDECESSOR / "result-seal.json"),
            "stage165_status": predecessor["status"],
        },
        "stage166": {
            "result_sha256": sha256(STAGE / "result.json"),
            "result_seal_sha256": sha256(STAGE / "result-seal.json"),
            "verification_sha256": sha256(VERIFICATION),
            "blind_instance_id": verification["blind_instance_id"],
            "source_commit": verification["source_commit"],
            "selected_wall_seconds": verification["selected_wall_seconds"],
            "selected_core_seconds": verification["selected_core_seconds"],
            "selected_peak_rss_bytes": verification["selected_peak_rss_bytes"],
            "selected_word_xors": verification["selected_word_xors"],
            "same_binary_control_speedup": verification[
                "same_binary_control_speedup"
            ],
            "same_host_f4_over_mitm_wall": verification[
                "same_host_f4_over_mitm_wall"
            ],
            "complete_campaign_cost": None,
        },
        "gates": result["gates"],
        "next_targets": [
            "execute the full frozen 160-input packet under licensed Magma F4",
            "run a fresh-target native-F4 schedule holdout rather than infer generality from the selected target",
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
    seal = load(output / "result-seal.json", "Stage-166 audit seal")
    require(seal.get("schema") == SEAL_SCHEMA, "audit seal schema changed")
    require(
        sha256(output / "audit.json") == seal.get("audit_sha256"),
        "audit seal changed",
    )
    current = compose()
    stored = load(output / "audit.json", "Stage-166 audit")
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
        raise SystemExit(f"stage166-gate-audit: {error}")


if __name__ == "__main__":
    main()

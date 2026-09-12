#!/usr/bin/env python3
"""Compose and verify the hosted Stage 58 n=53 unknown-scalar result."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage57_n53_unknown as stage57
import verify_koblitz_stage49_width_results as common


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
DEFAULT_OUTPUT = STAGE / "stage-58-n53-unknown-result-20260912"
SCHEMA = "koblitz_stage58_unknown_result.v1"
SEAL_SCHEMA = "koblitz_stage58_unknown_result_seal.v1"
CONFIG = {
    "run_id": 34714165792,
    "commit": "40557b15a174d2d3ac8ecba4b4e108856d9eb2ea",
    "artifact_id": 10304845829,
    "artifact_digest": "sha256:78b2ef8c4147ad09638655c6b22c59aa26f2eda7159a41dd4348150c641803ed",
    "archive": "koblitz-stage58-n53-unknown-evidence-34714165792.tar.gz",
    "archive_root": "koblitz-stage58-n53-unknown-34714165792",
    "prefix": "stage57",
    "job": "n53-unknown-scalar",
    "runner": stage57,
}


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def result() -> dict[str, Any]:
    with tempfile.TemporaryDirectory() as directory:
        root = common.extract(58, CONFIG, Path(directory))
        verification = stage57.verify(root / "stage57-build", root / "stage57-run")
        embedded = load(root / "stage57-verification.json", "embedded Stage-58 verification")
        require(verification == embedded, "Stage-58 embedded verification changed")
        run = load(root / "stage57-run/result.json", "Stage-58 run")
        build_seal = custody.sha256_file(root / "stage57-build/result-seal.json", "Stage-58 build seal")
        run_seal = custody.sha256_file(root / "stage57-run/result-seal.json", "Stage-58 run seal")
    require(verification.get("status") == "n53_public_unknown_scalar_verified", "Stage-58 verification status changed")
    require(verification.get("target_scalar_constructed_or_supplied") is False, "Stage-58 target scalar was constructed")
    require(verification.get("factor_base_logs_known_by_construction") is False, "Stage-58 factor-base logs were constructed")
    require(run.get("unknown_scalar_recovery_verified") is True, "Stage-58 recovery is unverified")
    require(run.get("direct_summary", {}).get("published_fixture_scalar") is None and run.get("rho", {}).get("published_fixture_scalar") is None, "Stage-58 retained an expected scalar")
    return {
        "schema": SCHEMA,
        "status": "hosted_n53_public_unknown_scalar_verified",
        "source_commit": CONFIG["commit"],
        "archive": common.archive_identity(58, CONFIG),
        "workflow": common.workflow_accounting(58, CONFIG),
        "artifact": common.artifact_accounting(58, CONFIG),
        "build_seal_sha256": build_seal,
        "run_seal_sha256": run_seal,
        "verification": verification,
        "comparison": run["comparison"],
        "direct": run["direct_summary"],
        "rho": run["rho"],
        "processes": {"direct": run["direct_process"]["metrics"], "rho": run["rho_process"]["metrics"]},
        "resources": {
            "build_outer": run["build_binding"]["outer_resources"],
            "science_outer": run["outer_resources"],
            "charged_total_core_seconds_available": run["charged_total_core_seconds_available"],
            "charged_sequential_wall_seconds_available": run["charged_sequential_wall_seconds_available"],
            "maximum_sampled_process_tree_rss_bytes": run["maximum_sampled_process_tree_rss_bytes"],
        },
        "factor_base_algebraically_point_defined": True,
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "target_scalar_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
        "unknown_scalar_recovery_verified": True,
        "sat_conflicts": None,
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "remaining_blockers": [
            "n=53 unknown-scalar direct remains 4.428199 times same-target rho wall",
            "fresh build plus direct remains 25.722944 times rho wall",
            "licensed same-instance Magma and unaffiliated review remain absent",
            "finite n=53 recovery does not establish an asymptotic improvement",
        ],
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-58 output must be new")
    output.mkdir()
    custody.write_json_new(output / "verification.json", result())
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "hosted_result_frozen", "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "result-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "verification.json", "committed Stage-58 result")
    require(committed == result(), "committed Stage-58 result differs from source evidence")
    seal = load(output / "result-seal.json", "Stage-58 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-58 result seal changed")
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-58 result inventory changed")
    return committed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    create = sub.add_parser("compose")
    create.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    check = sub.add_parser("verify")
    check.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    try:
        value = freeze(args.output.resolve()) if args.command == "compose" else verify(args.output.resolve())
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, VerificationError, common.VerificationError, stage57.Stage57Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage58-result: {error}")


if __name__ == "__main__":
    main()

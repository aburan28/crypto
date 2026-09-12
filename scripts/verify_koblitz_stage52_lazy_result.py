#!/usr/bin/env python3
"""Compose and verify the hosted Stage 52 lazy-label n=53 result."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage44_n53_parallel as stage44
import run_koblitz_stage48_n53_width as stage48
import verify_koblitz_stage49_width_results as common


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
PREDECESSOR = STAGE / "stage-49-n53-width-results-20260912/verification.json"
DEFAULT_OUTPUT = STAGE / "stage-52-n53-lazy-result-20260912"
SCHEMA = "koblitz_stage52_lazy_result.v1"
SEAL_SCHEMA = "koblitz_stage52_lazy_result_seal.v1"
CONFIG = {
    "run_id": 34711534312,
    "commit": "5e3997e9751cd0b370d482875cced796e682887f",
    "artifact_id": 10303597397,
    "artifact_digest": "sha256:a948751620c24650231d6b9813a00136e530a461ef9b4aef6cf340a52e10a620",
    "archive": "koblitz-stage52-n53-lazy-evidence-34711534312.tar.gz",
    "archive_root": "koblitz-stage52-n53-lazy-34711534312",
    "prefix": "stage48",
    "job": "n53-parallel-width",
    "runner": stage48,
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
    predecessor = load(PREDECESSOR, "Stage-49 predecessor")
    require(predecessor.get("status") == "hosted_n53_scratch_and_width_results_verified", "Stage-49 predecessor changed")
    hosted = common.one_result(52, CONFIG)
    verification = hosted["verification"]
    require(verification.get("source_commit") == CONFIG["commit"], "Stage-52 source changed")
    require(verification.get("relations") == 189, "Stage-52 relation count changed")
    require(verification.get("candidate_over_rho_wall_ratio") == 5.16632272222771, "Stage-52 rho ratio changed")
    require(verification.get("candidate_whole_process_crossover") is False, "Stage-52 boundary changed")
    with tempfile.TemporaryDirectory() as directory:
        root = common.extract(52, CONFIG, Path(directory))
        rows = stage44.parse_lines(root / "stage48-run/evidence/parallel-4096.stdout")
        base = next(row for row in rows if row.get("kind") == "point_defined_factor_base")
    support_bytes = base["support_table_allocated_bytes"]
    require(support_bytes == 1_275_068_416, "Stage-52 support allocation changed")
    return {
        "schema": SCHEMA,
        "status": "hosted_n53_lazy_label_result_verified",
        "source_commit": CONFIG["commit"],
        "predecessor_stage49": {"path": str(PREDECESSOR.relative_to(REPO)), "sha256": custody.sha256_file(PREDECESSOR, "Stage-49 result")},
        "hosted": hosted,
        "support_table_allocated_bytes": support_bytes,
        "prior_support_table_allocated_bytes": 1_409_286_144,
        "support_table_bytes_saved": 134_217_728,
        "preferred_query_mode": "pair_pair_parallel_4096",
        "same_public_target": True,
        "exact_relation_hash_equality": True,
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
        "fixture_scalar_used_by_collector": False,
        "sat_conflicts": None,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "remaining_blockers": [
            "preferred n=53 direct remains 5.166323 times same-target rho wall",
            "fresh build plus preferred direct remains 22.668683 times rho wall",
            "unknown-scalar n=53, licensed same-instance Magma, and unaffiliated review remain absent",
            "lazy label loading changes finite constants, not the asymptotic exponent",
        ],
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-52 output must be new")
    output.mkdir()
    custody.write_json_new(output / "verification.json", result())
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "hosted_result_frozen", "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "result-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "verification.json", "committed Stage-52 result")
    require(committed == result(), "committed Stage-52 result differs from source evidence")
    seal = load(output / "result-seal.json", "Stage-52 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-52 result seal changed")
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-52 result inventory changed")
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
    except (OSError, ValueError, KeyError, VerificationError, common.VerificationError, stage48.Stage48Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage52-result: {error}")


if __name__ == "__main__":
    main()

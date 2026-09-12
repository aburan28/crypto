#!/usr/bin/env python3
"""Compose and verify the hosted Stage 56 PCLMUL n=53 result."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage55_n53_pclmul as stage55
import verify_koblitz_stage49_width_results as common


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
PREDECESSOR = STAGE / "stage-52-n53-lazy-result-20260912/verification.json"
DEFAULT_OUTPUT = STAGE / "stage-56-n53-pclmul-result-20260912"
SCHEMA = "koblitz_stage56_pclmul_result.v1"
SEAL_SCHEMA = "koblitz_stage56_pclmul_result_seal.v1"
CONFIG = {
    "run_id": 34712955835,
    "commit": "e0daf815f4c82440c53ac4a2b53c79566f1084aa",
    "artifact_id": 10304460074,
    "artifact_digest": "sha256:372494e8db1ca6aec025e9edb52afd5e0847fed94c414b488b58d54d538583f4",
    "archive": "koblitz-stage56-n53-pclmul-evidence-34712955835.tar.gz",
    "archive_root": "koblitz-stage56-n53-pclmul-34712955835",
    "prefix": "stage55",
    "job": "n53-pclmul-ab",
    "runner": stage55,
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
    predecessor = load(PREDECESSOR, "Stage-52 predecessor")
    require(predecessor.get("status") == "hosted_n53_lazy_label_result_verified", "Stage-52 predecessor changed")
    hosted = common.one_result(56, CONFIG)
    verification = hosted["verification"]
    require(verification.get("source_commit") == CONFIG["commit"], "Stage-56 source changed")
    require(verification.get("relations") == 189, "Stage-56 relation count changed")
    require(verification.get("candidate_direct_wall_speedup") == 1.4186949503361426, "Stage-56 wall speedup changed")
    require(verification.get("candidate_direct_core_seconds_ratio") == 0.6143520974888967, "Stage-56 CPU ratio changed")
    require(verification.get("candidate_over_rho_wall_ratio") == 3.688473082311116, "Stage-56 rho ratio changed")
    require(verification.get("candidate_whole_process_crossover") is False, "Stage-56 boundary changed")
    with tempfile.TemporaryDirectory() as directory:
        root = common.extract(56, CONFIG, Path(directory))
        run = load(root / "stage55-run/result.json", "Stage-56 run")
    require(run["baseline_summary"].get("field_mul_backend") == "portable_sparse_carryless", "Stage-56 baseline backend changed")
    require(run["candidate_summary"].get("field_mul_backend") == "x86_64_pclmulqdq", "Stage-56 candidate backend changed")
    return {
        "schema": SCHEMA,
        "status": "hosted_n53_pclmul_result_verified",
        "source_commit": CONFIG["commit"],
        "predecessor_stage52": {"path": str(PREDECESSOR.relative_to(REPO)), "sha256": custody.sha256_file(PREDECESSOR, "Stage-52 result")},
        "hosted": hosted,
        "preferred_query_mode": "pair_pair_parallel_4096",
        "preferred_field_mul_backend": "x86_64_pclmulqdq",
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
            "preferred n=53 direct remains 3.688473 times same-target rho wall",
            "fresh build plus preferred direct remains 21.182086 times rho wall",
            "unknown-scalar n=53, licensed same-instance Magma, and unaffiliated review remain absent",
            "PCLMUL changes finite constants, not the asymptotic exponent",
        ],
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-56 output must be new")
    output.mkdir()
    custody.write_json_new(output / "verification.json", result())
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "hosted_result_frozen", "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "result-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "verification.json", "committed Stage-56 result")
    require(committed == result(), "committed Stage-56 result differs from source evidence")
    seal = load(output / "result-seal.json", "Stage-56 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-56 result seal changed")
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-56 result inventory changed")
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
    except (OSError, ValueError, KeyError, VerificationError, common.VerificationError, stage55.Stage55Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage56-result: {error}")


if __name__ == "__main__":
    main()

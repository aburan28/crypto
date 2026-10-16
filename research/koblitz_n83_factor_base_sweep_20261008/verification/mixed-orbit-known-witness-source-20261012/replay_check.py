#!/usr/bin/env python3
"""Exercise control-receipt replay with constructed data; never invoke SAT."""

import argparse
import copy
import json
from pathlib import Path

import control_supervisor as gate


HERE = Path(__file__).resolve().parent


def check(small: Path, envelope: Path) -> dict:
    if (gate.sha256_file(small) != gate.SMALL_SHA
            or gate.sha256_file(envelope) != gate.ENVELOPE_SHA):
        raise ValueError("retained archive SHA-256 mismatch")
    fixture = json.loads((HERE / "fixture.json").read_text())
    config = {
        "schema": "n83.mixed-orbit-known-witness-config/v1",
        "scope": "known_answer_control_only", "study": gate.STUDY,
        "fixture_blake3": gate.FIXTURE_BLAKE3,
        "small_object_blake3": gate.SMALL_BLAKE3,
        "envelope_object_blake3": gate.ENVELOPE_BLAKE3,
        "small_summands": 3, "residual_summands": 2,
        "conflict_budget": gate.CONFLICTS, "max_models": gate.MODELS,
        "budget_seconds": gate.WORKER_WALL,
        "memory_cgroup_limit_bytes": 4 * 1024**3,
        "memory_cgroup_swap_limit_bytes": 0,
    }
    candidate = {
        "a": "0", "b": "1", "small_point_indices": [0, 166, 332],
        "residual_points": [fixture["terms"][3]["point"], fixture["terms"][4]["point"]],
    }
    result = {
        "schema": "n83.mixed-orbit-known-witness-result/v1",
        "scope": "known_answer_control_only", "status": "PASS_known_witness_solver_access",
        "target": fixture["equation"]["target"],
        "natural_relations": None, "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None, "candidate": candidate,
        "search": {"solver_calls": 1, "conflicts": 1, "models": 1,
                   "invalid_models": 0, "patterns_attempted": 1,
                   "exhausted": False, "refuted": False},
        "base_import_elapsed_ms_L0": 8.0, "worker_wall_ms_L0": 434.0,
        "memory_cgroup_peak_bytes": 468664320,
    }
    positive = gate.control_receipt_error(config, result, fixture, small, envelope)
    wrong_small = copy.deepcopy(result)
    wrong_small["candidate"]["small_point_indices"][0] = 1
    small_error = gate.control_receipt_error(config, wrong_small, fixture, small, envelope)
    wrong_residual = copy.deepcopy(result)
    wrong_residual["candidate"]["residual_points"][0] = fixture["terms"][0]["point"]
    residual_error = gate.control_receipt_error(config, wrong_residual,
                                                 fixture, small, envelope)
    wrong_fixture = {**config, "fixture_blake3": "0" * 64}
    hash_error = gate.control_receipt_error(wrong_fixture, result, fixture, small, envelope)
    if (positive is not None
            or small_error != "control candidate failed independent archive replay"
            or residual_error != "control candidate failed independent archive replay"
            or hash_error != "control config fixture_blake3 mismatch"):
        raise ValueError("constructed control replay did not separate good and bad receipts")
    return {
        "schema": "n83.known-witness-constructed-replay/v1",
        "scope": "source_only_audit_test",
        "constructed_receipt_only": True,
        "solver_search_executed": False,
        "small_sha256": gate.SMALL_SHA,
        "envelope_sha256": gate.ENVELOPE_SHA,
        "fixture_sha256": gate.sha256_file(HERE / "fixture.json"),
        "valid_constructed_witness_replay": "PASS",
        "wrong_small_index": small_error,
        "wrong_residual_domain": residual_error,
        "wrong_fixture_hash": hash_error,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small", type=Path, required=True)
    parser.add_argument("--envelope", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "replay-check.json")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(check(args.small, args.envelope), indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text() != rendered:
            raise ValueError("stored constructed-replay record differs")
        print("PASS: constructed-replay record matches exact archives")
    else:
        args.output.write_text(rendered)
        print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()

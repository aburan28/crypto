#!/usr/bin/env python3
"""Exercise the K600 partial-relation auditor on exact retained archives.

This constructs an auditor input from a public known-answer group equation.
It does not invoke the relation solver or measure natural-target yield.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile

STUDY_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STUDY_DIR))
import k600_mixed_orbit_supervisor as gate  # noqa: E402

FIXTURE = (STUDY_DIR / "verification/mixed-orbit-known-witness-source-20261012/fixture.json")
FIXTURE_SHA256 = "fc203db28eed0164e47b0c755c3d9a09f80df1acd5e4c7e09dfe172fa5f0a538"


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def replay(small: Path, envelope: Path, output: Path) -> dict:
    small = small.resolve(strict=True)
    envelope = envelope.resolve(strict=True)
    if digest(small) != gate.SMALL_SHA256 or digest(envelope) != gate.ENVELOPE_SHA256:
        raise ValueError("retained archive SHA-256 mismatch")
    if digest(FIXTURE) != FIXTURE_SHA256:
        raise ValueError("known-answer fixture SHA-256 mismatch")
    fixture = json.loads(FIXTURE.read_text())
    if (fixture.get("scope") != "known_answer_control_only"
            or fixture.get("schema") != "n83.mixed-orbit-known-witness/v1"):
        raise ValueError("fixture provenance mismatch")
    terms = fixture["terms"]
    small_terms = [row for row in terms if row["domain"] == "small"]
    residual_terms = [row for row in terms if row["domain"] == "residual"]
    indices = [row["global_point_index"] for row in small_terms]
    if (indices != [0, 166, 332] or len(residual_terms) != 2
            or any(row["column"] < gate.SMALL_ORBITS for row in residual_terms)):
        raise ValueError("fixture does not define the K600 split")
    expected_sum = gate.parse_point(fixture["equation"]["target"])
    points = [gate.parse_point(row["point"]) for row in terms]
    if expected_sum is None or any(point is None for point in points):
        raise ValueError("invalid fixture point")
    if gate._group_sum(points) != expected_sum:
        raise ValueError("fixture group equation mismatch")

    order = gate.PRIMARY_ORDER
    a, b = int(gate.FIRST_TRIAL_A), int(gate.FIRST_TRIAL_B)
    neg_a_generator = gate.point_mul(gate.PRIMARY_GENERATOR, (order - a) % order)
    difference = gate.point_add(expected_sum, neg_a_generator)
    if difference is None:
        raise ValueError("constructed target is identity")
    target = gate.point_mul(difference, pow(b, -1, order))
    if (target is None or gate.point_mul(target, order) is not None
            or not gate.replay_equation(gate.FIRST_TRIAL_A, gate.FIRST_TRIAL_B,
                                        points[:3], points[3:], target)):
        raise ValueError("constructed target failed group replay")

    event = {
        "trial": 0, "a": gate.FIRST_TRIAL_A, "b": gate.FIRST_TRIAL_B,
        "elapsed_process_wall_ms": 0.0,
        "search": {"solver_calls": 1, "conflicts": 20, "models": 1,
                   "domain_clauses": gate.DOMAIN_CLAUSES_3_PLUS_2,
                   "invalid_models": 0, "group_rejected_models": 0,
                   "exhausted": False, "refuted": False},
        "candidate": {"source_id": 1, "a": gate.FIRST_TRIAL_A,
                      "b": gate.FIRST_TRIAL_B, "small_point_indices": indices,
                      "residual_points": [row["point"] for row in residual_terms]},
        "graph": {"kind": "pivot", "merge_depth": 0},
    }
    result = {
        "trials_attempted": 1, "solver_calls": 1, "solver_conflicts": 20,
        "solver_models": 1, "invalid_models": 0, "group_rejected_models": 0,
        "solver_exhausted_trials": 0, "solver_refuted_trials": 0,
        "graph_accepted_candidates": 1, "graph_pivot_rows": 1,
        "graph_capped_rows": 0, "completed_rows": 0,
        "graph_full_from_partials": 0, "graph_pivot_reductions": 0,
    }
    with tempfile.TemporaryDirectory(prefix="n83-k600-archive-replay-") as folder:
        event_file = Path(folder) / "events.jsonl"

        def check(candidate: dict) -> str | None:
            event_file.write_text(json.dumps(candidate, sort_keys=True) + "\n")
            return gate.replay_events(event_file, small, envelope, target, result)

        if (error := check(event)) is not None:
            raise ValueError(f"correct constructed relation rejected: {error}")
        wrong_small = copy.deepcopy(event)
        wrong_small["candidate"]["small_point_indices"][1] += 1
        wrong_small_error = check(wrong_small)
        wrong_residual = copy.deepcopy(event)
        wrong_residual["candidate"]["residual_points"][0] = small_terms[0]["point"]
        wrong_residual_error = check(wrong_residual)
        wrong_domain = copy.deepcopy(event)
        wrong_domain["search"]["domain_clauses"] -= 1
        wrong_domain_error = check(wrong_domain)
    if (wrong_small_error != "partial point equation failed independent replay"
            or wrong_residual_error != "residual point outside pinned envelope suffix"
            or wrong_domain_error != "event K600 domain-clause count mismatch"):
        raise ValueError("negative replay controls did not fail at their expected gates")
    receipt = {
        "schema": "n83.k600-real-archive-replay/v1",
        "scope": "constructed_known_answer_auditor_control_only",
        "status": "PASS_constructed_replay",
        "replay_source_sha256": digest(Path(__file__)),
        "supervisor_source_sha256": digest(Path(gate.__file__)),
        "fixture_sha256": FIXTURE_SHA256,
        "small_object_sha256": gate.SMALL_SHA256,
        "envelope_object_sha256": gate.ENVELOPE_SHA256,
        "small_indices": indices,
        "residual_columns": [row["column"] for row in residual_terms],
        "constructed_target": [str(target[0]), str(target[1])],
        "trial_coefficients": {"a": gate.FIRST_TRIAL_A, "b": gate.FIRST_TRIAL_B},
        "domain_clauses": gate.DOMAIN_CLAUSES_3_PLUS_2,
        "controls": {
            "correct_relation": "PASS",
            "wrong_small_index": wrong_small_error,
            "residual_outside_suffix": wrong_residual_error,
            "wrong_domain_clause_count": wrong_domain_error,
        },
        "solver_search_executed": False,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    with output.open("x", encoding="utf-8") as stream:
        json.dump(receipt, stream, sort_keys=True, indent=2)
        stream.write("\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small", type=Path, required=True)
    parser.add_argument("--envelope", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = replay(args.small, args.envelope, args.output)
    print(json.dumps({"status": result["status"], "receipt": str(args.output)}))


if __name__ == "__main__":
    main()

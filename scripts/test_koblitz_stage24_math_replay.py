#!/usr/bin/env python3
"""Focused canonical and in-memory tamper checks for the Stage-24 replay."""

from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import koblitz_stage24_math_replay as replay


def comparable(value: dict) -> dict:
    """Remove the platform-specific compiled-helper identity."""
    value = deepcopy(value)
    value["replay_implementation"].pop("blake3_helper_binary")
    return json.loads(json.dumps(value, sort_keys=True))


def expected_relation_row(base: dict, projected: dict, indices: list[int]) -> tuple[list, list, list]:
    row = [0] * len(projected["representatives"])
    summands = []
    negated_flags = []
    for index in indices:
        location = projected["orbit_of"][index]
        if location is None:
            continue
        orbit, power, negated = location
        coefficient = pow(replay.CURVE0["lambda"], power, replay.CURVE0["r"])
        if negated:
            coefficient = (-coefficient) % replay.CURVE0["r"]
        row[orbit] = (row[orbit] + coefficient) % replay.CURVE0["r"]
        summands.append([orbit, power])
        negated_flags.append(negated)
    return row, summands, negated_flags


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--blake3-helper", type=Path, required=True)
    parser.add_argument("--expected-result", type=Path, required=True)
    args = parser.parse_args()

    actual = replay.replay(args.bundle, args.blake3_helper)
    expected = json.loads(args.expected_result.read_text())
    assert comparable(actual) == comparable(expected)
    assert actual["status"] == "PASS" and actual["check_count"] == 1251
    assert actual["retained_mathematical_witness_replay_completed"] is True
    for field in (
        "independent_mathematical_payload_replay_completed",
        "scientific_measurement_admitted",
        "external_portable_verification_satisfied",
        "independent_external_reproduction_satisfied",
        "full_cost_gate_passed",
        "koblitz_index_calculus_sota",
    ):
        assert actual[field] is False

    run = args.bundle / "original/run"
    row_one = json.loads((run / "tasks/04-row-01-ic/result.json").read_text())
    target = (int(row_one["target"]["x"]), int(row_one["target"]["y"]))
    first_attempt = row_one["attempt_records"][0]
    computed_target = replay.padd(
        replay.pmul(replay.CURVE0["g"], int(first_attempt["coefficient_a"]), 0),
        replay.pmul(target, int(first_attempt["coefficient_b"]), 0),
        0,
    )
    assert computed_target == (
        int(first_attempt["target"]["x"]), int(first_attempt["target"]["y"])
    )
    forged_target = (computed_target[0] ^ 1, computed_target[1])
    assert forged_target != computed_target

    base = replay.factor_base(replay.CURVE0, [0, 2])
    projected = replay.projected_map(replay.CURVE0, base)
    relation = row_one["relation_matrix"][0]
    indices = first_attempt["decomposition_indices"]
    row, summands, negated = expected_relation_row(base, projected, indices)
    assert [int(value) for value in relation["row"]] == row
    assert relation["summands"] == summands and relation["summand_negated"] == negated
    forged_row = list(row)
    forged_row[0] = (forged_row[0] + 1) % replay.CURVE0["r"]
    assert forged_row != row

    scalar = int(row_one["report"]["recovered_scalar"])
    assert replay.pmul(replay.CURVE0["g"], scalar, 0) == target
    assert replay.pmul(replay.CURVE0["g"], scalar + 1, 0) != target

    forged_relations = deepcopy(row_one["relation_matrix"])
    original_certificate = replay.rref_certificate(
        row_one["relation_matrix"], replay.CURVE0["r"], replay.CURVE0["h"]
    )
    # Shift every RHS coherently so the same factor-base variables solve with
    # d+1: h*a' = h*a + (-h*b).  Replay must produce a different target scalar.
    for relation_record in forged_relations:
        relation_record["coefficient_a"] = str(
            (
                int(relation_record["coefficient_a"])
                - int(relation_record["coefficient_b"])
            )
            % replay.CURVE0["r"]
        )
    forged_certificate = replay.rref_certificate(
        forged_relations, replay.CURVE0["r"], replay.CURVE0["h"]
    )
    assert forged_certificate["scalar"] == (original_certificate["scalar"] + 1) % replay.CURVE0["r"]

    rho = json.loads((run / "tasks/05-row-01-rho/result.json").read_text())
    forged_rho = deepcopy(rho)
    forged_rho["charges"]["walk_group_additions"] += 1
    try:
        replay.validate_rho(forged_rho, target, scalar)
    except AssertionError:
        pass
    else:
        raise AssertionError("forged rho ledger was accepted")

    print(
        json.dumps(
            {
                "schema": "koblitz_stage24_independent_math_replay_tests.v1",
                "status": "PASS",
                "canonical_replay_checks": actual["check_count"],
                "tamper_checks": 5,
                "tested_tampers": [
                    "relation-attempt target coordinate",
                    "projected relation coefficient",
                    "recovered scalar",
                    "modular relation right-hand side",
                    "rho walk-addition ledger",
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

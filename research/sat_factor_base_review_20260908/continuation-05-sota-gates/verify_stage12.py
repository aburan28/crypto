#!/usr/bin/env python3
"""Verify the indexed projected-orbit map against the frozen reference."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


SEMANTIC_FIELDS = [
    "divisor_indices",
    "divisor_polynomial",
    "linearised_exponents",
    "dimension",
    "abscissae",
    "rational_points",
    "signed_frobenius_orbits_before_projection",
    "projected_points",
    "projected_signed_frobenius_orbits",
    "m_cofactor_admissible",
]


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def semantics(candidate: dict) -> dict:
    return {field: candidate[field] for field in SEMANTIC_FIELDS}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("new_directory", type=Path)
    parser.add_argument("reference_directory", type=Path)
    parser.add_argument("stage8_directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    rows = []
    for curve_a in [0, 1]:
        new = load(args.new_directory / f"a{curve_a}" / "discovery.json")
        old = load(args.reference_directory / f"a{curve_a}" / "discovery.json")
        new_meter = load(args.new_directory / f"a{curve_a}" / "discovery.metrics.json")
        old_meter = load(args.reference_directory / f"a{curve_a}" / "discovery.metrics.json")
        assert new_meter["returncode"] == old_meter["returncode"] == 0
        assert not new_meter["timed_out"] and not old_meter["timed_out"]
        assert [semantics(candidate) for candidate in new["candidates"]] == [
            semantics(candidate) for candidate in old["candidates"]
        ]
        assert semantics(new["selected"]) == semantics(old["selected"])
        candidate_projection_rows = []
        for new_candidate, old_candidate in zip(new["candidates"], old["candidates"]):
            new_ns = new_candidate["timing_ns"]["projection_census"]
            old_ns = old_candidate["timing_ns"]["projection_census"]
            candidate_projection_rows.append(
                {
                    "divisor_indices": new_candidate["divisor_indices"],
                    "reference_projection_ns": old_ns,
                    "indexed_projection_ns": new_ns,
                    "speedup": old_ns / new_ns,
                }
            )
        rows.append(
            {
                "curve_a": curve_a,
                "candidate_projection_rows": candidate_projection_rows,
                "reference_wall_seconds": old_meter["metrics"]["wall_seconds"],
                "indexed_wall_seconds": new_meter["metrics"]["wall_seconds"],
                "reference_core_seconds": old_meter["metrics"]["total_core_seconds"],
                "indexed_core_seconds": new_meter["metrics"]["total_core_seconds"],
                "process_core_speedup": old_meter["metrics"]["total_core_seconds"]
                / new_meter["metrics"]["total_core_seconds"],
            }
        )

    new_ic = load(args.new_directory / "n15-ic-53.json")
    new_ic_meter = load(args.new_directory / "n15-ic-53.metrics.json")
    old_ic = load(args.stage8_directory / "ic-53.json")
    old_ic_meter = load(args.stage8_directory / "ic-53.metrics.json")
    assert new_ic_meter["returncode"] == old_ic_meter["returncode"] == 0
    assert new_ic["target"] == old_ic["target"]
    assert new_ic["factor_base_predicate"] == old_ic["factor_base_predicate"]
    for field in [
        "relations",
        "trials",
        "sat_calls",
        "sat_models",
        "sat_refutations",
        "sat_unknowns",
        "sat_invalid_models",
        "sat_conflicts",
        "direct_relation",
        "recovered_scalar",
        "verified_unknown_scalar_recovery",
    ]:
        assert new_ic["report"][field] == old_ic["report"][field]

    old_total_core = sum(row["reference_core_seconds"] for row in rows)
    new_total_core = sum(row["indexed_core_seconds"] for row in rows)
    result = {
        "schema": "koblitz_projected_orbit_index_result.v1",
        "semantic_equivalence": True,
        "degree23_discovery_rows": rows,
        "degree23_combined": {
            "reference_core_seconds": old_total_core,
            "indexed_core_seconds": new_total_core,
            "process_core_speedup": old_total_core / new_total_core,
        },
        "degree15_scalar_blind_smoke": {
            "secret": 53,
            "relations": new_ic["report"]["relations"],
            "trials": new_ic["report"]["trials"],
            "conflicts": new_ic["report"]["sat_conflicts"],
            "recovered_and_verified": new_ic["report"]["verified_unknown_scalar_recovery"],
            "reference_core_seconds": old_ic_meter["metrics"]["total_core_seconds"],
            "indexed_core_seconds": new_ic_meter["metrics"]["total_core_seconds"],
        },
        "claim_boundary": "Setup-cost engineering with identical algebraic outputs and recovery counters; not a relation-solving or SOTA result.",
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verify orbit-derived cofactor admission against the pointwise reference."""

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
    parser.add_argument("original_directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    rows = []
    for curve_a in [0, 1]:
        new = load(args.new_directory / f"a{curve_a}" / "discovery.json")
        reference = load(args.reference_directory / f"a{curve_a}" / "discovery.json")
        new_meter = load(args.new_directory / f"a{curve_a}" / "discovery.metrics.json")
        reference_meter = load(args.reference_directory / f"a{curve_a}" / "discovery.metrics.json")
        assert new_meter["returncode"] == reference_meter["returncode"] == 0
        assert not new_meter["timed_out"] and not reference_meter["timed_out"]
        assert [semantics(candidate) for candidate in new["candidates"]] == [
            semantics(candidate) for candidate in reference["candidates"]
        ]
        assert semantics(new["selected"]) == semantics(reference["selected"])

        candidates = []
        for new_candidate, reference_candidate in zip(new["candidates"], reference["candidates"]):
            new_ns = new_candidate["timing_ns"]["cofactor_admission"]
            reference_ns = reference_candidate["timing_ns"]["cofactor_admission"]
            candidates.append(
                {
                    "divisor_indices": new_candidate["divisor_indices"],
                    "reference_admission_ns": reference_ns,
                    "orbit_admission_ns": new_ns,
                    "speedup": reference_ns / new_ns,
                }
            )
        rows.append(
            {
                "curve_a": curve_a,
                "candidates": candidates,
                "reference_wall_seconds": reference_meter["metrics"]["wall_seconds"],
                "orbit_wall_seconds": new_meter["metrics"]["wall_seconds"],
                "reference_core_seconds": reference_meter["metrics"]["total_core_seconds"],
                "orbit_core_seconds": new_meter["metrics"]["total_core_seconds"],
                "process_core_speedup": reference_meter["metrics"]["total_core_seconds"]
                / new_meter["metrics"]["total_core_seconds"],
            }
        )

    original_core = 0.0
    for curve_a in [0, 1]:
        original_meter = load(args.original_directory / f"a{curve_a}" / "discovery.metrics.json")
        assert original_meter["returncode"] == 0 and not original_meter["timed_out"]
        original_core += original_meter["metrics"]["total_core_seconds"]
    reference_core = sum(row["reference_core_seconds"] for row in rows)
    orbit_core = sum(row["orbit_core_seconds"] for row in rows)
    result = {
        "schema": "koblitz_orbit_cofactor_admission_result.v1",
        "semantic_equivalence": True,
        "degree23_rows": rows,
        "combined": {
            "original_stage10_core_seconds": original_core,
            "indexed_projection_stage12_core_seconds": reference_core,
            "orbit_admission_stage14_core_seconds": orbit_core,
            "incremental_process_core_speedup": reference_core / orbit_core,
            "cumulative_process_core_speedup": original_core / orbit_core,
        },
        "claim_boundary": "Public setup-cost engineering with unchanged algebraic outputs; not a relation-solving or SOTA result.",
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

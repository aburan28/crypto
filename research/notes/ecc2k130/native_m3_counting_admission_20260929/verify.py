#!/usr/bin/env python3
"""Independent exact-integer and parent-admission replay of the m3 envelope."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import prod
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[4]
ADMISSION = ROOT / "research/ecc2k130_leaf_native_m3_gate_20260926/admission.py"
INPUT = ROOT / "research/notes/ecc2k130/m10_export_capacity_20260925/inputs/balanced_m10_result.json"
FROZEN = Path(__file__).with_name("FROZEN.json")
EXPECTED_INPUT_HASH = "1bf1dc09dec347dcd421266bc4ece4b529931f21874fc3e6a26b5620f7a2cc25"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_balanced_product(total: int) -> int:
    quotient, remainder = divmod(total, 3)
    return quotient ** (3 - remainder) * (quotient + 1) ** remainder


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    row = json.loads(args.result.read_text())
    lock = json.loads(FROZEN.read_text())
    assert row["frozen_manifest_sha256"] == sha(FROZEN)
    for relative, digest in lock["sha256"].items():
        assert sha(ROOT / relative) == digest, relative
    assert row["protocol_commit"] == lock["protocol_commit"]
    assert row["domain"] == "ECC2K130-NATIVE-M3-COUNTING-ADMISSION-20260929-v1"
    assert row["classification"] == "PASS_BOUND"
    assert sha(INPUT) == EXPECTED_INPUT_HASH
    assert row["input_sha256"][str(INPUT.relative_to(ROOT))] == EXPECTED_INPUT_HASH
    parent = json.loads(INPUT.read_text())
    q = parent["q"]
    assert q == int(row["subgroup_order"])
    assert parent["physical_f0_points"] == row["m10_balanced_physical_points_per_slot"]
    total = 10 * parent["physical_f0_points"]
    assert total == row["fixed_total_physical_slot_choice_budget"]

    # An exhaustive three-way partition checks the balancing identity on a
    # small independent domain before using it for the large integer brackets.
    for budget in range(0, 61):
        exhaustive = max(first * second * (budget - first - second)
                         for first in range(budget + 1)
                         for second in range(budget - first + 1))
        assert exhaustive == independent_balanced_product(budget)

    slots = row["fixed_budget_balanced_m3_slots"]
    assert len(slots) == 3 and sum(slots) == total
    product = prod(slots)
    assert product == independent_balanced_product(total)
    assert product == int(row["fixed_budget_ordered_tuple_ceiling"])
    hit = min(Fraction(1), Fraction(product, q))
    probes = row["illustrative_target_probes"]
    assert probes == 1_000_000 and row["rows_per_probe_max"] == 1
    expected = probes * hit
    assert row["fixed_budget_hit_probability_upper"] == f"{hit.numerator}/{hit.denominator}"
    assert row["fixed_budget_expected_supported_probes_upper"] == f"{expected.numerator}/{expected.denominator}"
    any_hit = min(Fraction(1), expected)
    assert row["fixed_budget_probability_any_hit_upper"] == f"{any_hit.numerator}/{any_hit.denominator}"

    assert row["native_m3_implicit_slot_dimensions"] == [44, 44, 43]
    maxima = [(1 << (dimension + 1)) - 1 for dimension in (44, 44, 43)]
    assert row["native_m3_slot_physical_maxima"] == maxima
    assert row["native_m3_combined_physical_maximum"] == sum(maxima)
    assert [item["support_ceiling_threshold"] for item in row["thresholds"]] == ["1/100", "1/2"]
    for item in row["thresholds"]:
        numerator, denominator = map(int, item["support_ceiling_threshold"].split("/"))
        threshold = -(-(numerator * q) // denominator)
        budget = item["minimum_total_physical_slot_choices"]
        before = independent_balanced_product(budget - 1)
        at = independent_balanced_product(budget)
        assert before < threshold <= at
        assert item["tuple_threshold"] == str(threshold)
        assert item["maximum_tuples_at_previous_total"] == str(before)
        assert item["maximum_tuples_at_minimum_total"] == str(at)
        assert sum(item["balanced_slot_sizes"]) == budget
        assert prod(item["balanced_slot_sizes"]) == at
        floor = (max(0, budget - 6) + 4 * 262 - 1) // (4 * 262)
        assert item["conditional_minimum_projected_signed_orbit_log_classes"] == floor

    command = [sys.executable, str(ADMISSION), "--ordered-slot-sizes",
               *map(str, slots), "--target-space-size", str(q),
               "--targets", str(probes), "--required-rank", "1",
               "--rows-per-target-max", "1"]
    child = subprocess.run(command, check=True, text=True, capture_output=True)
    admitted = json.loads(child.stdout)
    assert admitted["classification"] == "NECESSARY_BOUND_ONLY"
    assert admitted["support_cardinality_upper"] == str(product)
    assert admitted["hit_probability_upper_fraction"] == row["fixed_budget_hit_probability_upper"]
    assert admitted["expected_supported_probes_upper_fraction"] == row["fixed_budget_expected_supported_probes_upper"]
    assert admitted["full_rank_probability_upper_fraction"] == row["fixed_budget_probability_any_hit_upper"]
    receipt = {
        "classification": "PASS_INDEPENDENT_REPLAY",
        "result_sha256": sha(args.result),
        "verifier_sha256": sha(Path(__file__)),
        "parent_admission_sha256": sha(ADMISSION),
        "parent_admission_command": command,
        "parent_admission_stdout_sha256": hashlib.sha256(child.stdout.encode()).hexdigest(),
        "checks": ["frozen_input", "small_budget_exhaustive_partition",
                   "fixed_budget_fractions", "threshold_minimality",
                   "physical_slot_maxima", "conditional_orbit_floor",
                   "parent_admission_rank_one_parity"],
    }
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"classification": receipt["classification"],
                      "checks": len(receipt["checks"])}, sort_keys=True))


if __name__ == "__main__":
    main()

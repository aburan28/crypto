#!/usr/bin/env python3
"""Exact necessary support envelope for ordered native m3 factor slots."""
from __future__ import annotations

import argparse
import ast
from fractions import Fraction
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
INPUT = ROOT / "research/notes/ecc2k130/m10_export_capacity_20260925/inputs/balanced_m10_result.json"
CHAIN = ROOT / "research/notes/ecc2k130/native_m3_chain_20260929/chain.py"
ADMISSION = ROOT / "research/ecc2k130_leaf_native_m3_gate_20260926/admission.py"
FROZEN = Path(__file__).with_name("FROZEN.json")
HASHES = {
    INPUT: "1bf1dc09dec347dcd421266bc4ece4b529931f21874fc3e6a26b5620f7a2cc25",
    CHAIN: "03b3e97b74e35a100ab60a7c5be0cea09108405e185cccda74c6942a94b21944",
    ADMISSION: "2a03a30c337d5dfb4d0d30ddff102293db6c8461ea1520138ef268ea5c8f0ece",
}
Q = 680564733841876926932320129493409985129
TARGET_PROBES = 1_000_000
THRESHOLDS = ((1, 100), (1, 2))
PARENT_COMMIT = "ce79b74645c8e564ec361918dfb72bd990cc56d6"
PROTOCOL_COMMIT = "bfa61af45cab9c05e10c85117369e9fb3005fb22"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def balanced_slots(total: int) -> tuple[int, int, int]:
    if total < 0:
        raise ValueError("negative total")
    each, remainder = divmod(total, 3)
    return tuple(each + int(i < remainder) for i in range(3))  # type: ignore[return-value]


def maximum_product(total: int) -> int:
    a, b, c = balanced_slots(total)
    return a * b * c


def minimum_budget(threshold: int) -> int:
    if threshold <= 0:
        raise ValueError("nonpositive tuple threshold")
    low, high = 0, 1
    while maximum_product(high) < threshold:
        high *= 2
    while low + 1 < high:
        middle = (low + high) // 2
        if maximum_product(middle) < threshold:
            low = middle
        else:
            high = middle
    return high


def fraction_string(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def frozen_dims() -> tuple[int, int, int]:
    tree = ast.parse(CHAIN.read_text())
    matches = [node for node in tree.body if isinstance(node, ast.Assign)
               and any(isinstance(target, ast.Name) and target.id == "LEAF_DIMS"
                       for target in node.targets)]
    if len(matches) != 1:
        raise AssertionError("missing or repeated LEAF_DIMS")
    value = ast.literal_eval(matches[0].value)
    if value != (44, 44, 43):
        raise AssertionError("frozen native m3 dimensions changed")
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    lock = json.loads(FROZEN.read_text())
    if lock["protocol_commit"] != PROTOCOL_COMMIT:
        raise AssertionError("protocol commit lock drift")
    for relative, expected in lock["sha256"].items():
        if sha256(ROOT / relative) != expected:
            raise AssertionError(f"source/input lock drift: {relative}")

    actual_hashes = {}
    for path, expected in HASHES.items():
        actual = sha256(path)
        if actual != expected:
            raise AssertionError(f"frozen input changed: {path.relative_to(ROOT)}")
        actual_hashes[str(path.relative_to(ROOT))] = actual
    source = json.loads(INPUT.read_text())
    if (source["q"] != Q or source["physical_f0_points"] != 7977
            or source["arm"] != {"m": 10, "d": 13}
            or source["ordered_physical_tuples"] != 7977 ** 10):
        raise AssertionError("m10 input fields changed")
    dims = frozen_dims()
    physical_budget = 10 * source["physical_f0_points"]
    slots = balanced_slots(physical_budget)
    tuple_ceiling = maximum_product(physical_budget)
    hit = min(Fraction(1), Fraction(tuple_ceiling, Q))
    expected = TARGET_PROBES * hit
    thresholds = []
    for numerator, denominator in THRESHOLDS:
        required_tuples = (numerator * Q + denominator - 1) // denominator
        total = minimum_budget(required_tuples)
        log_classes = (max(0, total - 6) + 1047) // 1048
        thresholds.append({
            "support_ceiling_threshold": f"{numerator}/{denominator}",
            "tuple_threshold": str(required_tuples),
            "minimum_total_physical_slot_choices": total,
            "balanced_slot_sizes": list(balanced_slots(total)),
            "maximum_tuples_at_previous_total": str(maximum_product(total - 1)),
            "maximum_tuples_at_minimum_total": str(maximum_product(total)),
            "conditional_minimum_projected_signed_orbit_log_classes": log_classes,
        })
    maximum_slot_physical = [(1 << (dimension + 1)) - 1 for dimension in dims]
    result = {
        "domain": "ECC2K130-NATIVE-M3-COUNTING-ADMISSION-20260929-v1",
        "classification": "PASS_BOUND",
        "parent_main_commit": PARENT_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
        "frozen_manifest_sha256": sha256(FROZEN),
        "input_sha256": actual_hashes,
        "subgroup_order": str(Q),
        "target_distribution": "uniform_marginal_one_q_coset",
        "illustrative_target_probes": TARGET_PROBES,
        "rows_per_probe_max": 1,
        "m10_balanced_physical_points_per_slot": source["physical_f0_points"],
        "fixed_total_physical_slot_choice_budget": physical_budget,
        "fixed_budget_balanced_m3_slots": list(slots),
        "fixed_budget_ordered_tuple_ceiling": str(tuple_ceiling),
        "fixed_budget_hit_probability_upper": fraction_string(hit),
        "fixed_budget_expected_supported_probes_upper": fraction_string(expected),
        "fixed_budget_probability_any_hit_upper": fraction_string(min(Fraction(1), expected)),
        "native_m3_implicit_slot_dimensions": list(dims),
        "native_m3_slot_physical_maxima": maximum_slot_physical,
        "native_m3_combined_physical_maximum": sum(maximum_slot_physical),
        "thresholds": thresholds,
        "scope": "exact necessary counts; no observed factor counts, PDP yield, rank or cost",
    }
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"classification": result["classification"],
                      "fixed_budget_hit_probability_upper": result["fixed_budget_hit_probability_upper"],
                      "threshold_totals": [row["minimum_total_physical_slot_choices"]
                                           for row in thresholds]}, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Independently replay the n37 native-base counting proof by dynamic programming."""

import argparse
from decimal import Decimal, localcontext
import gzip
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
INPUTS = {
    "frozen_sha256": ROOT / "research/notes/ecc2k130/disjoint_cold_v2_20261001/FROZEN.json",
    "raw_sha256": ROOT / "research/koblitz_isogeny_descent_37_results_20260925/raw.json.gz",
    "contract_sha256": ROOT / "research/koblitz_isogeny_descent_37_20260925/contract.json",
}
EXPECTED_HASHES = {
    "frozen_sha256": "da958a3f1117dd1b88703fed2055c5c9f64255516e1e920f35d193bf53c6a88d",
    "raw_sha256": "eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90",
    "contract_sha256": "ba5e9844517ca00d1a798b2caf4552ac37ffab73e875dfcba1e04d5fa89ed8c0",
}


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def ball_by_dp(dimensions, max_norm):
    """Count integer vectors by exact L1 norm, adding one dimension at a time."""
    exact = [1] + [0] * max_norm
    for _ in range(dimensions):
        updated = exact.copy()
        for norm in range(1, max_norm + 1):
            updated[norm] += 2 * sum(exact[norm - magnitude] for magnitude in range(1, norm + 1))
        exact = updated
    return sum(exact)


def replay():
    bound_bytes = (HERE / "BOUND.json").read_bytes()
    bound = json.loads(bound_bytes)
    if bound["schema"] != "ecc2k130-n37-native-column-capacity-v1" or bound["status"] != "PASS":
        raise ValueError("wrong bound schema or status")
    if bound["producer_sha256"] != sha256((HERE / "bound.py").read_bytes()):
        raise ValueError("producer hash changed")
    actual_hashes = {label: sha256(path.read_bytes()) for label, path in INPUTS.items()}
    if bound["input_sha256"] != actual_hashes or actual_hashes != EXPECTED_HASHES:
        raise ValueError("source input hashes changed")
    frozen = json.loads(INPUTS["frozen_sha256"].read_bytes())
    raw = json.loads(gzip.decompress(INPUTS["raw_sha256"].read_bytes()))
    contract = json.loads(INPUTS["contract_sha256"].read_bytes())
    spec = frozen["specs"]["n37_L1024"]
    k, order = spec["K"], spec["subgroup_order"]
    if (k, order, spec["n"], spec["L"], spec["a"], spec["automorphism_size"]) != (
            42, 230603167, 37, 1024, 0, 74):
        raise ValueError("source control changed")
    if (raw["contract"] != contract or contract["factor_base_subgroup_order"] != order or
            raw["isogeny_certificate"]["source_group_order"] !=
            raw["isogeny_certificate"]["target_group_order"] or
            raw["isogeny_certificate"]["degree"] != 73):
        raise ValueError("descent certificate changed")
    group_order = raw["isogeny_certificate"]["source_group_order"]
    if group_order != 137439487532 or group_order % 73 == 0:
        raise ValueError("rational isogeny is not certified bijective")
    if (bound["source_and_leaf_group_order"], bound["rational_isogeny_bijective"],
            bound["isogeny_degree"], bound["field_degree"], bound["source_curve_a"]) != (
            group_order, True, 73, 37, 0):
        raise ValueError("source/leaf group interpretation changed")
    if (bound["equal_log_columns"], bound["subgroup_order"],
            bound["source_signed_frobenius_orbit_size"],
            bound["source_physical_points_if_all_orbits_full"],
            bound["native_signed_physical_points_without_extra_action"]) != (
            k, order, 74, 74 * k, 2 * k):
        raise ValueError("column or physical-point interpretation changed")
    rows = bound["rows"]
    if [row["summands_at_most"] for row in rows] != list(range(2, 8)):
        raise ValueError("incomplete arity grid")
    for row in rows:
        m = row["summands_at_most"]
        count = ball_by_dp(k, m)
        if count != row["formal_coefficient_vectors"]:
            raise ValueError(f"coefficient count mismatch at m={m}")
        with localcontext() as ctx:
            ctx.prec = 40
            ratio = format(Decimal(count) / Decimal(order), ".12f")
            ceiling = format(Decimal(min(count, order)) / Decimal(order), ".12f")
        if (row["formal_count_over_subgroup_decimal"],
                row["uniform_target_support_ceiling_decimal"]) != (ratio, ceiling):
            raise ValueError(f"decimal ratio mismatch at m={m}")
        if (row["uniform_target_support_ceiling_numerator"],
                row["uniform_target_support_ceiling_denominator"]) != (min(count, order), order):
            raise ValueError(f"uniform support ceiling mismatch at m={m}")
        threshold = row["first_columns_with_counting_ceiling_one"]
        if not (ball_by_dp(threshold - 1, m) < order <= ball_by_dp(threshold, m)):
            raise ValueError(f"first-column threshold mismatch at m={m}")
    if (bound["first_arity_with_counting_ceiling_one_at_k42"], bound["decision"],
            bound["full_dlp_speedup"]) != (
            6, "NATIVE_K42_UP_TO_FIVE_SUMMANDS_CAPACITY_NO_GO", None):
        raise ValueError("decision or scope changed")
    return {
        "schema": "ecc2k130-n37-native-column-capacity-replay-v1",
        "status": "PASS",
        "verifier_sha256": sha256(Path(__file__).read_bytes()),
        "bound_sha256": sha256(bound_bytes),
        "independent_method": "dynamic_programming_over_exact_l1_norm",
        "checked_arities": [row["summands_at_most"] for row in rows],
        "first_arity_with_counting_ceiling_one_at_k42": 6,
        "decision": bound["decision"],
        "full_dlp_speedup": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    with args.out.open("x") as output:
        json.dump(replay(), output, indent=2)
        output.write("\n")


if __name__ == "__main__":
    main()

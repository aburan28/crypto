#!/usr/bin/env python3
"""Exact necessary support ceilings for a fixed n37 native 42-column base."""

import argparse
from decimal import Decimal, localcontext
import gzip
import hashlib
import json
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
FROZEN = ROOT / "research/notes/ecc2k130/disjoint_cold_v2_20261001/FROZEN.json"
RAW = ROOT / "research/koblitz_isogeny_descent_37_results_20260925/raw.json.gz"
CONTRACT = ROOT / "research/koblitz_isogeny_descent_37_20260925/contract.json"
EXPECTED = {
    "frozen_sha256": "da958a3f1117dd1b88703fed2055c5c9f64255516e1e920f35d193bf53c6a88d",
    "raw_sha256": "eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90",
    "contract_sha256": "ba5e9844517ca00d1a798b2caf4552ac37ffab73e875dfcba1e04d5fa89ed8c0",
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def formal_vectors(k, m):
    """Number of integer coefficient vectors with L1 norm at most m."""
    return sum(2**s * comb(k, s) * comb(m, s) for s in range(min(k, m) + 1))


def first_capacity_k(group_order, m):
    upper = 1
    while formal_vectors(upper, m) < group_order:
        upper *= 2
    lower = 0
    while lower + 1 < upper:
        middle = (lower + upper) // 2
        if formal_vectors(middle, m) >= group_order:
            upper = middle
        else:
            lower = middle
    return upper


def produce():
    frozen_bytes, raw_bytes, contract_bytes = FROZEN.read_bytes(), RAW.read_bytes(), CONTRACT.read_bytes()
    actual = {
        "frozen_sha256": digest(frozen_bytes),
        "raw_sha256": digest(raw_bytes),
        "contract_sha256": digest(contract_bytes),
    }
    if actual != EXPECTED:
        raise ValueError(f"frozen input hash mismatch: {actual}")
    frozen = json.loads(frozen_bytes)
    raw = json.loads(gzip.decompress(raw_bytes))
    contract = json.loads(contract_bytes)
    spec = frozen["specs"]["n37_L1024"]
    certificate = raw["isogeny_certificate"]
    k, r = spec["K"], spec["subgroup_order"]
    if (spec["n"], spec["L"], spec["a"], spec["automorphism_size"], k) != (37, 1024, 0, 74, 42):
        raise ValueError("the selected source control changed")
    if (contract["field_exponent"], contract["factor_base_subgroup_order"],
            contract["source_group_order"], contract["prime_degree"]) != (37, r, 137439487532, 73):
        raise ValueError("the archived subgroup or descent changed")
    if raw["contract"] != contract or raw["contract_sha256"] != actual["contract_sha256"]:
        raise ValueError("the archived run does not match its frozen contract")
    if (certificate["degree"], certificate["source_group_order"],
            certificate["target_group_order"]) != (73, contract["source_group_order"],
                                                   contract["source_group_order"]):
        raise ValueError("the archived isogeny certificate changed")
    if certificate["source_group_order"] % 73 == 0 or r % 73 == 0:
        raise ValueError("the isogeny need not be injective on rational points")

    rows = []
    with localcontext() as ctx:
        ctx.prec = 40
        for m in range(2, 8):
            count = formal_vectors(k, m)
            rows.append({
                "summands_at_most": m,
                "formal_coefficient_vectors": count,
                "formal_count_over_subgroup_decimal": format(Decimal(count) / Decimal(r), ".12f"),
                "uniform_target_support_ceiling_numerator": min(count, r),
                "uniform_target_support_ceiling_denominator": r,
                "uniform_target_support_ceiling_decimal": format(Decimal(min(count, r)) / Decimal(r), ".12f"),
                "first_columns_with_counting_ceiling_one": first_capacity_k(r, m),
            })
    return {
        "schema": "ecc2k130-n37-native-column-capacity-v1",
        "status": "PASS",
        "producer_sha256": digest(Path(__file__).read_bytes()),
        "input_sha256": actual,
        "source_control": "n37_L1024",
        "field_degree": 37,
        "source_curve_a": 0,
        "isogeny_degree": 73,
        "source_and_leaf_group_order": contract["source_group_order"],
        "rational_isogeny_bijective": True,
        "subgroup_order": r,
        "equal_log_columns": k,
        "source_signed_frobenius_orbit_size": spec["automorphism_size"],
        "source_physical_points_if_all_orbits_full": k * spec["automorphism_size"],
        "native_signed_physical_points_without_extra_action": 2 * k,
        "rows": rows,
        "first_arity_with_counting_ceiling_one_at_k42": next(
            row["summands_at_most"] for row in rows
            if row["formal_coefficient_vectors"] >= r),
        "decision": "NATIVE_K42_UP_TO_FIVE_SUMMANDS_CAPACITY_NO_GO",
        "full_dlp_speedup": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    with args.out.open("x") as output:
        json.dump(produce(), output, indent=2)
        output.write("\n")


if __name__ == "__main__":
    main()

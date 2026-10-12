#!/usr/bin/env python3
"""Independent exact-integer replay of the m83 relation-capacity screen."""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import hashlib
from itertools import combinations_with_replacement
import json
from math import comb
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
FILES = {
    "m83_summary": ROOT / "research/notes/ecc2k130/m83_solver_matrix_20260927/results/validated_summary.json",
    "m83_results": ROOT / "research/notes/ecc2k130/m83_solver_matrix_20260927/RESULTS.md",
    "challenge_source": ROOT / "research/ecc2k130_relations/relations.py",
}
FROZEN = {
    "m83_summary": "ff5b525fd188fb30f13b1959eec823ae3a239424c8c493503382ccc751c00261",
    "m83_results": "00a1e65642713dc04a71e158386e6d0048806fe87d73233766ea041afcc0e687",
    "challenge_source": "0180501ef8fe00f5c54f910822a28cd86dc3c2c1af164348976c1c2e25cc763f",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_scientific(numerator: int, denominator: int) -> str:
    with localcontext() as context:
        context.prec = 40
        return format(Decimal(numerator) / Decimal(denominator), ".8E")


def read_order(path: Path, pattern: str) -> int:
    values = re.findall(pattern, path.read_text(), flags=re.MULTILINE)
    assert len(values) == 1, (path, values)
    return int(values[0])


def replay(output: dict) -> dict:
    actual_hashes = {name: digest(path) for name, path in FILES.items()}
    assert actual_hashes == output["input_sha256"] == FROZEN
    assert digest(HERE / "derive.py") == output["derive_source_sha256"]
    order83 = read_order(FILES["m83_results"], r"r=(\d+)")
    order131 = read_order(FILES["challenge_source"], r"^CHALLENGE_ELL = (\d+)$")
    assert output["subgroup_orders"] == {"m83": order83, "m131": order131}
    source = json.loads(FILES["m83_summary"].read_text())["full_orbit_relation_control"]
    assert len(source) == 2 and sorted((x["seed"], x["signed_orbit_size"]) for x in source) == [
        (260938, 332), (260939, 498)
    ]
    assert all(x["pair_attempts"] == comb(x["signed_orbit_size"] + 1, 2)
               and x["natural_rank"] == x["natural_verified_queries"] == 0 for x in source)

    observed = output["archived_m83_bases"]
    assert len(observed) == 6
    expected_keys = {(x["seed"], x["signed_orbit_size"], m)
                     for x in source for m in (3, 4, 10)}
    assert {(x["seed"], x["base_points"], x["summands"]) for x in observed} == expected_keys
    for row in observed:
        b, m = row["base_points"], row["summands"]
        tuples = comb(b + m - 1, m)
        one = min(order83, tuples)
        four = min(order83, 4 * one)
        assert row["unordered_tuples"] == tuples
        assert row["one_uniform_target_ceiling_numerator"] == one
        assert row["four_uniform_target_union_ceiling_numerator"] == four
        assert row["denominator"] == order83
        assert row["one_uniform_target_ceiling"] == exact_scientific(one, order83)
        assert row["four_uniform_target_union_ceiling"] == exact_scientific(four, order83)
        assert row["observed_natural_hits"] == (0 if m == 3 else None)
        assert row["observed_natural_rank"] == (0 if m == 3 else None)

    thresholds = output["planning_thresholds"]
    assert len(thresholds) == 8
    assert {(x["field_degree"], x["summands"]) for x in thresholds} == {
        (degree, m) for degree in (83, 131) for m in (3, 4, 7, 10)
    }
    for row in thresholds:
        order = {83: order83, 131: order131}[row["field_degree"]]
        b, m = row["smallest_base_for_1pct_counting_ceiling"], row["summands"]
        assert b >= 1 and 100 * comb(b + m - 1, m) >= order
        assert b == 1 or 100 * comb(b + m - 2, m) < order
        tuples, pairs = comb(b + m - 1, m), comb(b + 1, 2)
        assert row["unordered_tuples_at_threshold"] == tuples
        assert row["one_target_counting_ceiling"] == exact_scientific(min(order, tuples), order)
        assert row["unordered_pair_entries"] == pairs
        assert row["hypothetical_pair_bytes_at_16_each"] == 16 * pairs
        assert row["hypothetical_pair_gib_at_16_each"] == exact_scientific(16 * pairs, 1024**3)

    m3 = [row for row in observed if row["summands"] == 3]
    joint = min(order83, sum(row["four_uniform_target_union_ceiling_numerator"] for row in m3))
    assert output["joint_eight_target_m3_union_ceiling_numerator"] == joint
    assert output["joint_eight_target_m3_union_ceiling_denominator"] == order83
    assert output["joint_eight_target_m3_union_ceiling"] == exact_scientific(joint, order83)
    assert all(100 * row["four_uniform_target_union_ceiling_numerator"] < order83 for row in m3)
    assert output["decision"] == "M83_M3_FIXED_BASES_UNDERPOWERED"
    assert all(output[key] is None for key in (
        "full_dlp_cost", "matched_rho_cost", "common_operation_S", "ecc2k130_speedup"))

    # Exhaustive small-group controls test the combinatorial inequality itself,
    # including collisions and repeated points, without importing derive.py.
    small_controls = []
    for modulus, base, m in ((37, (1, 3, 7, 9), 3),
                             (19, (0, 1, 4), 4),
                             (11, (0, 1), 2)):
        support = {sum(points) % modulus
                   for points in combinations_with_replacement(base, m)}
        ceiling = min(modulus, comb(len(base) + m - 1, m))
        assert len(support) <= ceiling
        small_controls.append({"modulus": modulus, "base_points": list(base),
                               "summands": m, "actual_support": len(support),
                               "counting_ceiling": ceiling})
    return {
        "verified": True, "decision": output["decision"],
        "input_sha256": actual_hashes,
        "derive_source_sha256": digest(HERE / "derive.py"),
        "verifier_source_sha256": digest(HERE / "verify.py"),
        "archived_rows_checked": len(observed),
        "thresholds_checked": len(thresholds),
        "small_group_controls": small_controls,
        "m83_m3_joint_eight_target_upper": output["joint_eight_target_m3_union_ceiling"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = replay(json.loads(args.input.read_text()))
    result["derived_json_sha256"] = digest(args.input)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

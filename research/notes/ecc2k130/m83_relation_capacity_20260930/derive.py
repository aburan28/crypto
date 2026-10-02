#!/usr/bin/env python3
"""Exact support-count screen for the frozen m83 full-orbit controls."""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
SUMMARY = ROOT / "research/notes/ecc2k130/m83_solver_matrix_20260927/results/validated_summary.json"
RESULTS = ROOT / "research/notes/ecc2k130/m83_solver_matrix_20260927/RESULTS.md"
CHALLENGE = ROOT / "research/ecc2k130_relations/relations.py"
EXPECTED_SHA256 = {
    "m83_summary": "ff5b525fd188fb30f13b1959eec823ae3a239424c8c493503382ccc751c00261",
    "m83_results": "00a1e65642713dc04a71e158386e6d0048806fe87d73233766ea041afcc0e687",
    "challenge_source": "0180501ef8fe00f5c54f910822a28cd86dc3c2c1af164348976c1c2e25cc763f",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def count_multisets(base_points: int, summands: int) -> int:
    """C(B+m-1,m), evaluated by an exact multiplicative recurrence."""
    assert base_points >= 1 and summands >= 1
    answer = 1
    for step in range(1, summands + 1):
        product = answer * (base_points + step - 1)
        assert product % step == 0
        answer = product // step
    return answer


def sci(numerator: int, denominator: int) -> str:
    with localcontext() as context:
        context.prec = 40
        return format(Decimal(numerator) / Decimal(denominator), ".8E")


def minimum_base_for_fraction(order: int, summands: int, numerator: int, denominator: int) -> int:
    low, high = 0, 1
    while denominator * count_multisets(high, summands) < numerator * order:
        high *= 2
    while high - low > 1:
        middle = (low + high) // 2
        if denominator * count_multisets(middle, summands) >= numerator * order:
            high = middle
        else:
            low = middle
    return high


def recorded_order(path: Path, expression: str) -> int:
    matches = re.findall(expression, path.read_text(), flags=re.MULTILINE)
    if len(matches) != 1:
        raise ValueError(f"expected one subgroup order in {path}, found {matches}")
    return int(matches[0])


def derive() -> dict:
    inputs = {name: sha256(path) for name, path in (
        ("m83_summary", SUMMARY), ("m83_results", RESULTS),
        ("challenge_source", CHALLENGE))}
    if inputs != EXPECTED_SHA256:
        raise ValueError("frozen capacity input SHA-256 changed")
    order83 = recorded_order(RESULTS, r"r=(\d+)")
    order131 = recorded_order(CHALLENGE, r"^CHALLENGE_ELL = (\d+)$")
    if (order83, order131) != (
        2417851639230796216685689,
        680564733841876926932320129493409985129,
    ):
        raise ValueError("frozen subgroup order changed")
    summary = json.loads(SUMMARY.read_text())
    controls = summary["full_orbit_relation_control"]
    if sorted((row["seed"], row["signed_orbit_size"]) for row in controls) != [
        (260938, 332), (260939, 498)
    ]:
        raise ValueError("frozen full-orbit control bases changed")

    archived = []
    for row in controls:
        base = row["signed_orbit_size"]
        if (row["pair_attempts"] != base * (base + 1) // 2
                or row["natural_rank"] != 0 or row["natural_verified_queries"] != 0):
            raise ValueError("frozen m83 control pair count or outcome changed")
        for summands in (3, 4, 10):
            tuples = count_multisets(base, summands)
            supported = min(order83, tuples)
            four = min(order83, 4 * supported)
            archived.append({
                "seed": row["seed"], "base_points": base, "summands": summands,
                "unordered_tuples": tuples,
                "one_uniform_target_ceiling_numerator": supported,
                "four_uniform_target_union_ceiling_numerator": four,
                "denominator": order83,
                "one_uniform_target_ceiling": sci(supported, order83),
                "four_uniform_target_union_ceiling": sci(four, order83),
                "observed_natural_hits": 0 if summands == 3 else None,
                "observed_natural_rank": row["natural_rank"] if summands == 3 else None,
            })

    thresholds = []
    for degree, order in ((83, order83), (131, order131)):
        for summands in (3, 4, 7, 10):
            base = minimum_base_for_fraction(order, summands, 1, 100)
            tuples = count_multisets(base, summands)
            pairs = base * (base + 1) // 2
            thresholds.append({
                "field_degree": degree, "summands": summands,
                "smallest_base_for_1pct_counting_ceiling": base,
                "unordered_tuples_at_threshold": tuples,
                "one_target_counting_ceiling": sci(min(order, tuples), order),
                "unordered_pair_entries": pairs,
                "hypothetical_pair_bytes_at_16_each": 16 * pairs,
                "hypothetical_pair_gib_at_16_each": sci(16 * pairs, 1024**3),
            })

    m3 = [row for row in archived if row["summands"] == 3]
    joint = min(order83, sum(row["four_uniform_target_union_ceiling_numerator"] for row in m3))
    underpowered = all(
        100 * row["four_uniform_target_union_ceiling_numerator"] < order83
        for row in m3
    )
    return {
        "schema_version": 1,
        "protocol": "m83_relation_capacity_20260930/PROTOCOL.md",
        "input_sha256": inputs,
        "derive_source_sha256": sha256(HERE / "derive.py"),
        "subgroup_orders": {"m83": order83, "m131": order131},
        "target_law": "uniform subgroup target; four natural targets per archived seed",
        "bound": "supported m-sums <= min(r, C(B+m-1,m)); four-target probability <= 4*support/r",
        "archived_m83_bases": archived,
        "joint_eight_target_m3_union_ceiling_numerator": joint,
        "joint_eight_target_m3_union_ceiling_denominator": order83,
        "joint_eight_target_m3_union_ceiling": sci(joint, order83),
        "planning_thresholds": thresholds,
        "decision": "M83_M3_FIXED_BASES_UNDERPOWERED" if underpowered else "HYPOTHESIS_FALSIFIED",
        "full_dlp_cost": None,
        "matched_rho_cost": None,
        "common_operation_S": None,
        "ecc2k130_speedup": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.write_text(json.dumps(derive(), indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

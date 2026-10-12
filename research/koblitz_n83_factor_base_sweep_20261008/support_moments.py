#!/usr/bin/env python3
"""Exact support-count screen for the retained degree-83 factor bases.

For a uniform nonidentity subgroup target Q, each unordered m-multiset of
factor-base points has one fixed sum. If M is the number of such multisets and
Z of them sum to identity, the expected number of decompositions is
(M-Z)/(r-1). Therefore P(any decomposition) <= min(1, M/(r-1)). This is a
first-moment ceiling, not a measured yield or a solver-runtime prediction.
"""

import gzip
import hashlib
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PANEL = ROOT / "pilot-01"
ORDERS = {
    0: 2417851639230796216685689,
    1: 8569786107849059,
}
SUMMANDS = range(2, 7)
ORBIT_POINTS_PER_COLUMN = 166


def support_count(points: int, summands: int) -> int:
    """Unordered sums with repeated factor-base points allowed."""
    assert points > 0 and summands > 0
    return math.comb(points + summands - 1, summands)


def vacuous_ceiling_column_threshold(order: int, summands: int) -> int:
    """First K for which M/(r-1) reaches one; this is not sufficiency."""
    low, high = 0, 1
    while support_count(ORBIT_POINTS_PER_COLUMN * high, summands) < order - 1:
        high *= 2
    while high - low > 1:
        middle = (low + high) // 2
        if support_count(ORBIT_POINTS_PER_COLUMN * middle, summands) >= order - 1:
            high = middle
        else:
            low = middle
    return high


def compact_ratio(numerator: int, denominator: int) -> str:
    with localcontext() as context:
        context.prec = 24
        return f"{Decimal(numerator) / Decimal(denominator):.6E}"


def assert_header_binding(row: dict, order: int) -> str:
    object_path = PANEL / row["object"]
    with gzip.open(object_path, "rt") as stream:
        header = json.loads(stream.readline())
    assert int(header["subgroup_order"]) == order
    assert header["curve_a"] == row["a"]
    assert header["point_count"] == row["points"]
    assert header["orbit_columns"] == row["columns"]
    return header["curve_slug"]


def summarize() -> dict:
    manifest_bytes = (PANEL / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    replay = json.loads((PANEL / "replay.json").read_bytes())
    assert manifest["status"] == "completed_factor_base_panel"
    assert replay["status"] == "PASS"
    assert len(manifest["bases"]) == manifest["completed_base_count"] == 54

    by_arm_and_size: dict[tuple[int, int], list[dict]] = {}
    for row in manifest["bases"]:
        arm, columns, points = row["a"], row["columns"], row["points"]
        assert arm in ORDERS and points == ORBIT_POINTS_PER_COLUMN * columns
        by_arm_and_size.setdefault((arm, columns), []).append(row)

    cases = []
    thresholds = []
    for arm, columns in sorted(by_arm_and_size):
        rows = by_arm_and_size[(arm, columns)]
        order = ORDERS[arm]
        points = rows[0]["points"]
        assert all(row["points"] == points for row in rows)
        slug = assert_header_binding(rows[0], order)
        for summands in SUMMANDS:
            multisets = support_count(points, summands)
            ceiling_numerator = min(multisets, order - 1)
            cases.append({
                "curve_a": arm,
                "curve_slug": slug,
                "subgroup_order": str(order),
                "orbit_columns": columns,
                "point_count": points,
                "base_variants_at_size": len(rows),
                "distinct_point_sets_at_size": len({row["point_set_blake3"] for row in rows}),
                "summands": summands,
                "unordered_multisets_with_repetition": str(multisets),
                "uniform_all_target_expected_count": {
                    "numerator": str(multisets), "denominator": str(order)
                },
                "uniform_nonidentity_target_hit_ceiling": {
                    "numerator": str(ceiling_numerator),
                    "denominator": str(order - 1),
                    "decimal": compact_ratio(ceiling_numerator, order - 1),
                },
            })
    for arm, order in ORDERS.items():
        for summands in SUMMANDS:
            threshold = vacuous_ceiling_column_threshold(order, summands)
            if threshold > 1:
                assert support_count(ORBIT_POINTS_PER_COLUMN * (threshold - 1), summands) < order - 1
            assert support_count(ORBIT_POINTS_PER_COLUMN * threshold, summands) >= order - 1
            thresholds.append({
                "curve_a": arm,
                "summands": summands,
                "first_columns_with_vacuous_markov_ceiling": threshold,
                "interpretation": "necessary threshold for this upper bound to reach one; no sufficiency or runtime claim",
            })

    assert len(cases) == 2 * 3 * len(SUMMANDS)
    return {
        "schema": "n83.support-moment-screen/v1",
        "study": manifest["study"],
        "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "replay_panel_manifest_blake3": replay["panel_manifest_blake3"],
        "target_model": "uniform nonidentity point of the stated prime-order subgroup",
        "counting_model": "unordered m-multisets of distinct factor-base points; repetitions allowed",
        "identity_sum_multisets": "unknown; using zero gives a conservative upper bound",
        "scope_limit": "fixed public fixtures, solver discovery probability, runtime, rank, and large-prime partials are not inferred",
        "cases": cases,
        "thresholds": thresholds,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    output = PANEL / "support-moments.json"
    encoded = (json.dumps(summarize(), indent=2, sort_keys=True) + "\n").encode()
    if output.exists():
        assert output.read_bytes() == encoded, "existing immutable screen disagrees"
    else:
        with output.open("xb") as stream:
            stream.write(encoded)
    print(output)


if __name__ == "__main__":
    main()

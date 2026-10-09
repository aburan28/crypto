#!/usr/bin/env python3
"""Exact rank-query necessary conditions for one-row full-smooth oracles.

The input is the immutable support-moments receipt.  A query target must be
uniform over the *whole* subgroup and the oracle may return at most one row
per query.  Query outcomes need not be independent.  This is a probability
ceiling, not a runtime estimate or a model of the existing capped driver.
"""

import hashlib
import json
import math
from pathlib import Path

from support_moments import ORDERS, PANEL, support_count


def ceil_div(numerator: int, denominator: int) -> int:
    assert numerator >= 0 and denominator > 0
    return (numerator + denominator - 1) // denominator


def success_numerator(multisets: int, order: int) -> int:
    """Numerator of the capped one-query success ceiling over denominator r."""
    assert multisets > 0 and order > 1
    return min(multisets, order)


def rank_probability_ceiling(
    queries: int, required_rank: int, multisets: int, order: int
) -> tuple[int, int]:
    """Exact Markov ceiling for rank K from at most one row per query."""
    assert queries >= 0 and required_rank > 0
    if queries < required_rank:
        return (0, 1)
    numerator = min(queries * success_numerator(multisets, order), order * required_rank)
    denominator = order * required_rank
    divisor = math.gcd(numerator, denominator)
    return (numerator // divisor, denominator // divisor)


def half_probability_query_floor(required_rank: int, multisets: int, order: int) -> int:
    """Necessary query count before the ceiling can reach one half."""
    assert required_rank > 0
    return max(
        required_rank,
        ceil_div(required_rank * order, 2 * success_numerator(multisets, order)),
    )


def summarize() -> dict:
    support_bytes = (PANEL / "support-moments.json").read_bytes()
    support = json.loads(support_bytes)
    manifest_bytes = (PANEL / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    replay = json.loads((PANEL / "replay.json").read_bytes())
    if support["schema"] != "n83.support-moment-screen/v1":
        raise ValueError("unexpected support-moments schema")
    if support["manifest_sha256"] != hashlib.sha256(manifest_bytes).hexdigest():
        raise ValueError("support screen does not bind this panel manifest")
    if support["replay_panel_manifest_blake3"] != replay["panel_manifest_blake3"]:
        raise ValueError("support screen does not bind this replay")
    if replay["status"] != "PASS" or manifest["status"] != "completed_factor_base_panel":
        raise ValueError("panel or replay is incomplete")
    if support["study"] != manifest["study"]:
        raise ValueError("study mismatch")
    expected = {
        (row["a"], row["columns"], summands)
        for row in manifest["bases"]
        for summands in range(2, 7)
    }
    if len(expected) != 30 or len(support["cases"]) != 30:
        raise ValueError("unexpected panel or support-screen case count")

    cases = []
    seen = set()
    for case in support["cases"]:
        arm = case["curve_a"]
        order = int(case["subgroup_order"])
        columns = case["orbit_columns"]
        summands = case["summands"]
        points = case["point_count"]
        multisets = int(case["unordered_multisets_with_repetition"])
        key = (arm, columns, summands)
        if key in seen or key not in expected:
            raise ValueError("duplicate or unexpected support-screen case")
        seen.add(key)
        if order != ORDERS[arm] or points != 166 * columns:
            raise ValueError("wrong subgroup order or base size")
        if multisets != support_count(points, summands):
            raise ValueError("wrong multiset count")
        one_query = success_numerator(multisets, order)
        floor = half_probability_query_floor(columns, multisets, order)
        below = rank_probability_ceiling(floor - 1, columns, multisets, order)
        at = rank_probability_ceiling(floor, columns, multisets, order)
        if not (2 * below[0] < below[1] and 2 * at[0] >= at[1]):
            raise ValueError("half-probability threshold is not minimal")
        cases.append({
            "curve_a": arm,
            "curve_slug": case["curve_slug"],
            "orbit_columns": columns,
            "point_count": points,
            "summands": summands,
            "required_rank_rows": columns,
            "uniform_all_target_one_query_hit_ceiling": {
                "numerator": str(one_query), "denominator": str(order)
            },
            "minimum_queries_before_half_rank_probability_is_possible": str(floor),
            "probability_ceiling_at_preceding_query_count": {
                "numerator": str(below[0]), "denominator": str(below[1])
            },
        })

    if seen != expected:
        raise ValueError("support-screen case is missing")

    return {
        "schema": "n83.rank-query-screen/v1",
        "study": support["study"],
        "support_moments_sha256": hashlib.sha256(support_bytes).hexdigest(),
        "target_model": "each query target is uniform over the whole prime-order subgroup; independence is not required",
        "oracle_model": "at most one full-smooth relation row per query; K independent rows required for rank K",
        "bound": "P(rank >= K after q queries) <= 0 if q < K; otherwise min(1, q*min(M,r)/(r*K)), where M=binomial(B+m-1,m)",
        "scope_limit": "does not apply to biased query targets, multiple rows returned per query, large-prime partials, or solvers whose accepted event differs; no runtime inference",
        "cases": cases,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    output = PANEL / "rank-query-screen.json"
    encoded = (json.dumps(summarize(), indent=2, sort_keys=True) + "\n").encode()
    if output.exists():
        assert output.read_bytes() == encoded, "existing immutable screen disagrees"
    else:
        with output.open("xb") as stream:
            stream.write(encoded)
    print(output)


if __name__ == "__main__":
    main()

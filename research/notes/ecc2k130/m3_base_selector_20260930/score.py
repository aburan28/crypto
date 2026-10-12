#!/usr/bin/env python3
"""Charged, target-independent support and first-witness row score for m3."""
from __future__ import annotations

import hashlib
import json


def score_candidate(curve, base: list[tuple[int, int]], universe: set,
                    meter, rank_tracker_type, modulus: int = 421) -> dict:
    """Score one ordered eight-point base before a challenge Q is available.

    The traversal uses the actual solver's first-third/first-pair order.
    It charges all 36 pair additions and 288 pair-plus-third additions.
    `universe` is the complete set of nonzero subgroup points; it is not a
    sample of the later relation-target stream.
    """
    assert len(base) == len(set(base)) == 8
    assert len(universe) == modulus - 1 == 420
    before = meter.snapshot()
    pairs = [(i, j, curve.add(base[i], base[j]))
             for i in range(8) for j in range(i, 8)]
    assert len(pairs) == 36
    first = {}
    for k in range(8):
        for i, j, pair in pairs:
            target = curve.add(pair, base[k])
            if target is None:
                continue
            assert target in universe, "base summand left the prime subgroup"
            first.setdefault(target, (k, i, j))
    after = meter.snapshot()
    # The enclosing phase owns elapsed CPU. Two nested snapshots can never
    # have identical cpu_ns even when their arithmetic counters agree.
    cost = {key: value for key, value in meter.delta(before, after).items()
            if key != "cpu_ns"}
    assert cost["group_add"] == 324, cost
    assert len(first) <= 120

    rank = rank_tracker_type(width=8, modulus=modulus)
    rows = set()
    digest_rows = []
    for target in sorted(first):
        k, i, j = first[target]
        coefficients = [0] * 8
        for index in (i, j, k):
            coefficients[index] += 1
        rows.add(tuple(coefficients))
        rank.add(coefficients, 0)
        digest_rows.append([target[0], target[1], k, i, j])
    digest = hashlib.sha256(json.dumps(
        digest_rows, separators=(",", ":")).encode()).hexdigest()
    return {"distinct_support": len(first),
            "first_witness_base_row_rank": len(rank.pivots),
            "distinct_first_witness_rows": len(rows),
            "first_witness_sha256": digest,
            "score_cost": cost,
            "score_mod_r_ops": dict(sorted(rank.ops.items()))}


def choose(scores: list[dict]) -> int | None:
    """Return the protocol's candidate index, or None if rank is deficient."""
    assert len(scores) == 2
    eligible = [index for index, score in enumerate(scores)
                if score["first_witness_base_row_rank"] == 8]
    if not eligible:
        return None
    return max(eligible, key=lambda index: (
        scores[index]["distinct_support"],
        scores[index]["distinct_first_witness_rows"],
        -index))

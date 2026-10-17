#!/usr/bin/env python3
"""Independent Python replay of the frozen P-192 adaptive-base audit."""

import hashlib
import json
import sys
from pathlib import Path

GX = int("188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012", 16)
GY = int("07192B95FFC8DA78631011ED6B24CDD573F977A11E794811", 16)
BB = int("64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1", 16)


def add(p1, p2, modulus):
    if p1 is None:
        return p2
    if p2 is None:
        return p1
    x1, y1 = p1
    x2, y2 = p2
    if x1 == x2 and (y1 + y2) % modulus == 0:
        return None
    if p1 == p2:
        slope = ((3 * x1 * x1 - 3) * pow(2 * y1, -1, modulus)) % modulus
    else:
        slope = ((y2 - y1) * pow(x2 - x1, -1, modulus)) % modulus
    x3 = (slope * slope - x1 - x2) % modulus
    return x3, (slope * (x1 - x3) - y1) % modulus


def neg(point, modulus):
    return None if point is None else (point[0], -point[1] % modulus)


def mul(scalar, point, modulus):
    acc = None
    while scalar:
        if scalar & 1:
            acc = add(acc, point, modulus)
        point = add(point, point, modulus)
        scalar >>= 1
    return acc


def key(point):
    return "O" if point is None else f"{point[0]}:{point[1]}"


def on_curve(point, modulus):
    if point is None:
        return True
    x, y = point
    return (y * y - x * x * x + 3 * x - BB) % modulus == 0


def cell(pre, coefficients, modulus):
    def chi(point):
        if point is None:
            return "I"
        x = point[0]
        value = 0
        for c in reversed(coefficients):
            value = (value * x + c) % modulus
        if value == 0:
            return "0"
        return "+" if pow(value, (modulus - 1) // 2, modulus) == 1 else "-"

    return chi(pre) + chi(add(pre, pre, modulus))


def main():
    if len(sys.argv) != 6:
        raise SystemExit("usage: verify_adaptive_private_columns.py BASE DENSITY GUIDE MIXED OUTPUT")
    raw = [Path(p).read_bytes() for p in sys.argv[1:]]
    base, density, guide, mixed, result = [json.loads(b) for b in raw]
    assert result["input_sha256"] == [hashlib.sha256(b).hexdigest() for b in raw[:4]]
    modulus = int(base["field_prime"])
    order = int(result["subgroup_order"])
    assert modulus == int(density["p"]) == int(result["field_prime"])
    assert order == int(density["r"])
    generator = GX, GY
    assert on_curve(generator, modulus) and mul(order, generator, modulus) is None
    original = [(int(v["x"]), int(v["y"])) for v in base["points"]]
    assert len(original) == 16
    for point in original:
        assert on_curve(point, modulus) and mul(order, point, modulus) is None
    guided = [neg(mul(2 * int(s), generator, modulus), modulus)
              for s in guide["accepted_scalars_decimal"]]
    assert len(guided) == 64
    signed = original + [neg(p, modulus) for p in original] + guided + [neg(p, modulus) for p in guided]
    assert len(signed) == len({key(p) for p in signed}) == 160
    coefficients = [int(v) for v in density["sextic_coefficients_decimal"]]
    candidates = result["candidates"]
    assert len(candidates) == 64
    keys = set()
    cells = {}
    query_points = []
    old_keys = {key(point) for point in signed}
    for j, record in enumerate(candidates):
        assert record["index"] == j
        scalar = int(record["query_scalar"])
        assert scalar == int(mixed["independent_accepted_queries"][j]["scalar"])
        pre = mul(scalar, generator, modulus)
        assert cell(pre, coefficients, modulus) == "++"
        query = neg(add(pre, pre, modulus), modulus)
        query_points.append(query)
        assert key(query) == record["query_point_key"]
        i, k = record["base_indices"]
        assert i == j % 16 and k == (7 * j + 3) % 16 and i != k
        new = add(add(query, neg(original[i], modulus), modulus), neg(original[k], modulus), modulus)
        assert on_curve(new, modulus) and mul(order, new, modulus) is None
        assert add(add(original[i], original[k], modulus), new, modulus) == query
        assert key(new) == record["new_point_key"]
        assert new is not None
        assert key(new) not in old_keys
        assert key(new) not in keys
        keys.add(key(new))
        pre_new = mul((order - 1) // 2, new, modulus)
        label = cell(pre_new, coefficients, modulus)
        assert label == record["new_point_character_cell"]
        cells[label] = cells.get(label, 0) + 1
    assert cells == result["new_point_character_cells"]
    full_keys = {}
    full_duplicates = full_old = full_identity = 0
    base_pairs = [(i, k, add(original[i], original[k], modulus))
                  for i in range(16) for k in range(i, 16)]
    assert len(base_pairs) == 136
    for j, query in enumerate(query_points):
        for i, k, base_pair in base_pairs:
            new = add(query, neg(base_pair, modulus), modulus)
            new_key = key(new)
            full_identity += new is None
            full_old += new_key in old_keys
            full_duplicates += new_key in full_keys
            full_keys[new_key] = (j, i, k)
    assert keys <= set(full_keys)
    assert result["full_pair_candidate_rows"] == 64 * 136 == 8704
    assert len(full_keys) == result["full_pair_distinct_points"]
    assert (full_duplicates, full_old, full_identity) == (
        result["full_pair_duplicate_points"],
        result["full_pair_points_already_in_old_base"],
        result["full_pair_identity_points"]) == (0, 0, 0)
    assert hashlib.sha256("\n".join(sorted(full_keys)).encode()).hexdigest() == result["full_pair_keys_sha256"]
    pair = [[add(signed[i], signed[j], modulus) for j in range(i, 160)] for i in range(160)]
    support_hits = 0
    full_support_hits = 0
    for i in range(160):
        for j in range(i, 160):
            partial = pair[i][j - i]
            for k in range(j, 160):
                target_key = key(add(partial, signed[k], modulus))
                support_hits += target_key in keys
                full_support_hits += target_key in full_keys
    assert support_hits == result["old_three_sum_witnesses"] == 0
    assert full_support_hits == result["full_pair_old_three_sum_witnesses"] == 0
    assert result["full_pair_private_identity_block"] is True
    assert (result["full_pair_rank_if_private"],
            result["full_pair_unknown_columns_if_private"],
            result["full_pair_old_unknown_nullity_if_private"]) == (8704, 8720, 16)
    assert (result["new_point_duplicates"], result["new_points_already_in_old_base"],
            result["identity_new_points"]) == (0, 0, 0)
    assert result["private_identity_block"] is True
    assert (result["augmented_relation_rank_if_private"],
            result["augmented_unknown_columns_if_private"],
            result["old_unknown_nullity_if_private"]) == (64, 80, 16)
    print(json.dumps({"status": "ADAPTIVE_PRIVATE_COLUMNS_INDEPENDENT_PASS",
                      "checked_queries": 64, "distinct_new_points": 64,
                      "old_support_hits": support_hits, "character_cells": cells,
                      "full_pair_candidates": len(full_keys),
                      "full_pair_old_support_hits": full_support_hits}, sort_keys=True))


if __name__ == "__main__":
    main()

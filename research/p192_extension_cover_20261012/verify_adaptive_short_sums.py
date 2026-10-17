#!/usr/bin/env python3
"""Independent affine-curve replay of adaptive two- and four-sum membership."""

import hashlib
import json
import sys
import time
from pathlib import Path

from verify_adaptive_private_columns import GX, GY, add, key, mul, neg, on_curve


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def parse_key(value):
    return None if value == "O" else tuple(map(int, value.split(":")))


def main():
    if len(sys.argv) != 5:
        raise SystemExit("usage: verify_adaptive_short_sums.py BASE GUIDE ADAPTIVE RESULT")
    start = time.monotonic()
    raw = [Path(p).read_bytes() for p in sys.argv[1:]]
    base, guide, adaptive, result = map(json.loads, raw)
    assert result["input_sha256"] == list(map(digest, raw[:3]))
    p, r = int(result["field_prime"]), int(result["subgroup_order"])
    assert p == int(base["field_prime"]) == int(adaptive["field_prime"])
    assert r == int(adaptive["subgroup_order"])
    assert adaptive["input_sha256"][0] == digest(raw[0])
    assert adaptive["input_sha256"][2] == digest(raw[1])
    unknown = [(int(row["x"]), int(row["y"])) for row in base["points"]]
    assert len(unknown) == 16 and all(on_curve(v, p) for v in unknown)
    guided = [neg(mul(2 * int(s), (GX, GY), p), p)
              for s in guide["accepted_scalars_decimal"]]
    assert len(guided) == 64
    old = unknown + [neg(v, p) for v in unknown] + guided + [neg(v, p) for v in guided]
    assert len(old) == len({key(v) for v in old}) == 160
    pairs = [(i, j, add(old[i], old[j], p)) for i in range(160) for j in range(i, 160)]
    index = {}
    for i, j, value in pairs:
        index.setdefault(key(value), (i, j))
    assert len(pairs) == result["old_unordered_pairs_with_repetition"] == 12880
    assert len(index) == result["distinct_pair_sum_points"]
    two_hits, four_hits = [], []
    for k, row in enumerate(adaptive["candidates"]):
        assert row["index"] == k
        c = parse_key(row["new_point_key"])
        assert on_curve(c, p) and mul(r, c, p) is None
        assert add(add(unknown[row["base_indices"][0]],
                       unknown[row["base_indices"][1]], p), c, p) == parse_key(row["query_point_key"])
        if key(c) in index:
            two_hits.append({"candidate_index": k, "old_indices": list(index[key(c)])})
        for i, j, value in pairs:
            complement = add(c, neg(value, p), p)
            if key(complement) in index:
                u, v = index[key(complement)]
                assert add(value, add(old[u], old[v], p), p) == c
                four_hits.append({"candidate_index": k, "old_indices": [i, j, u, v]})
                break
    assert len(adaptive["candidates"]) == result["selected_new_points"] == 64
    assert result["two_sum_membership_tests"] == 64
    assert result["four_sum_complement_tests_upper_bound"] == 64 * 12880
    assert two_hits == result["two_sum_hits"]
    assert four_hits == result["four_sum_hits"]
    print(json.dumps({"status": "ADAPTIVE_SHORT_SUMS_INDEPENDENT_PASS",
                      "selected": 64, "pair_sums": len(pairs), "distinct_pair_sums": len(index),
                      "two_hits": len(two_hits), "four_hits": len(four_hits),
                      "elapsed_ms": round(1000 * (time.monotonic() - start))}, sort_keys=True))


if __name__ == "__main__":
    main()

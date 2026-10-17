#!/usr/bin/env python3
"""Independent signed-multiset count for the frozen 16+64 Type I base.

This uses no elliptic-curve arithmetic. Equality with the native target count
certifies that its repeated targets are exhausted by formal cancellation.
"""

import json
import sys
from collections import Counter
from itertools import combinations_with_replacement
from pathlib import Path


def term(index: int) -> tuple[int, int]:
    if index < 16:
        return index, 1
    if index < 32:
        return index - 16, -1
    if index < 96:
        return index - 16, 1
    return index - 80, -1


def formal(indices: tuple[int, int, int]) -> tuple[tuple[int, int], ...]:
    coefficients: Counter[int] = Counter()
    for index in indices:
        column, sign = term(index)
        coefficients[column] += sign
    return tuple(sorted((column, value) for column, value in coefficients.items() if value))


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: verify_formal_collision_counts.py NATIVE.json")
    native_bytes = Path(sys.argv[1]).read_bytes()
    native = json.loads(native_bytes)
    vectors = {"bearing": Counter(), "anchored": Counter()}
    known_only = set()
    zero_unknown = set()
    counts = Counter()
    for indices in combinations_with_replacement(range(160), 3):
        vector = formal(indices)
        bearing = any(column < 16 for column, _ in vector)
        kind = "bearing" if bearing else "anchored"
        vectors[kind][vector] += 1
        if bearing:
            counts["bearing_triples"] += 1
        elif indices[0] < 32:
            counts["zero_unknown_triples"] += 1
            zero_unknown.add(vector)
        else:
            counts["known_only_triples"] += 1
            known_only.add(vector)
    assert sum(counts.values()) == native["all_unordered_triples"] == 695520
    for name in ("bearing_triples", "known_only_triples", "zero_unknown_triples"):
        assert counts[name] == native[name], name
    assert len(vectors["bearing"]) == native["bearing_targets"] == 333184
    assert len(vectors["anchored"]) == native["anchored_targets"] == 349696
    assert len(known_only) == native["known_only_targets"] == 349696
    assert len(zero_unknown) == native["zero_unknown_distinct_targets"] == 128
    assert zero_unknown <= known_only
    assert sum(n - 1 for n in vectors["bearing"].values()) == native["repeated_bearing_triples"] == 2528
    assert sum(n - 1 for n in vectors["anchored"].values()) == native["repeated_anchored_triples"] == 10112
    print(json.dumps({
        "status": "FORMAL_COLLISION_COMBINATORIAL_PASS",
        "bearing_formal_vectors": len(vectors["bearing"]),
        "anchored_formal_vectors": len(vectors["anchored"]),
        "zero_unknown_vectors": len(zero_unknown),
        "zero_unknown_subset_known_only": True,
        "native_target_counts_equal_formal_vector_counts": True,
    }, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Independent bit-serial/Fermat factor-only verification; no six-sum oracle."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "rotated_pdp_corpus_20260925" / "verify.py"
spec = importlib.util.spec_from_file_location("pdp_preflight_verify", SOURCE)
assert spec and spec.loader
pdp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pdp)


def targets() -> list[int]:
    seen, values, counter = set(), [], 0
    while len(values) < 64:
        payload = f"ECC2K130-ROTATED-SLOT-PLACEMENT-20260925-v1/target/{counter}".encode("ascii")
        k = 1 + int.from_bytes(hashlib.sha256(payload).digest(), "big") % 130872
        if k not in seen:
            seen.add(k)
            values.append(k)
        counter += 1
    return values


def verify(expected: dict) -> None:
    assert expected["no_tuple_sums_computed"] is True
    assert expected["positions"] == list(range(6)) and expected["normal_coordinate_index"] == 18
    assert expected["shared_targets"] == targets()
    field = pdp.parent_verify.GF(19, 0x80027)
    pdp.rabin_prime_degree(field)
    curve = pdp.parent_verify.E(field)
    assert pdp.source_order(19) == 4 * 130873
    h = pdp.generator(curve, 130873)
    assert h == (385982, 301867)
    normal = []
    element = 3
    for _ in range(19):
        normal.append(element)
        element = field.square(element)
    assert element == 3 and pdp.rank(normal) == 19 and field.trace(3) == 1
    assert len(expected["arms"]) == 6
    for position, arm in enumerate(expected["arms"]):
        assert arm["position"] == position
        bases = []
        for i in range(6):
            b = [normal[i], normal[i + 6]]
            if i == position:
                b.append(normal[18])
            bases.append(b)
        assert [pdp.rank(b) for b in bases] == arm["slot_ranks"]
        assert pdp.rank([v for b in bases for v in b]) == arm["joint_rank"] == 13
        factors = []
        for b in bases:
            xs = {0}
            for v in b:
                xs.update({x ^ v for x in tuple(xs)})
            factors.append(sorted(p for x in xs for p in pdp.lifts(curve, x)))
        assert [[list(point) for point in factor] for factor in factors] == arm["factor_points"]
        assert [len(factor) for factor in factors] == arm["factor_sizes"]
        projected = [[curve.scalar(point, 4) for point in factor] for factor in factors]
        signed_counts = [len(set(values)) for values in projected]
        assert signed_counts == arm["slot_signed_projected_counts"]
        pair_sets = [{min(point, pdp.neg(point)) for point in values if point is not None}
                     for values in projected]
        pair_counts = [len(values) for values in pair_sets]
        assert pair_counts == arm["slot_pair_column_counts"]
        assert sum(pair_counts) == arm["column_budget"]
        assert len(set().union(*pair_sets)) == arm["global_column_count"]
        assert math.prod(map(len, factors)) == arm["tuple_count"]
    print("All six factor-only placements, ranks, columns, point lists and targets independently verified; no tuple sums read.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected", type=Path, required=True)
    args = parser.parse_args()
    verify(json.loads(args.expected.read_text()))

#!/usr/bin/env python3
"""Factor-only freeze for six n19 rank-two-plus-b18 slots; no tuple sums."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import produce


def build() -> dict:
    f = produce.pdp.parent.Field(19, [0, 1, 2, 5])
    f.rabin_prime_degree()
    curve = produce.pdp.parent.Curve(f)
    assert produce.pdp.parent.source_group_order(19) == 4 * 130873
    h = produce.pdp.generator(curve, 130873)
    assert h == (385982, 301867)
    normal = produce.pdp.parent.normal_conjugates(f, 3)
    assert len(normal) == 19 and produce.pdp.parent.rank(normal) == 19
    arms = []
    for position in range(6):
        bases = [[normal[i], normal[i + 6]] + ([normal[18]] if i == position else [])
                 for i in range(6)]
        ranks = [produce.pdp.parent.rank(b) for b in bases]
        assert ranks == [2 + (i == position) for i in range(6)]
        assert produce.pdp.parent.rank([v for b in bases for v in b]) == 13
        factors = [sorted(point for x in produce.pdp.parent.all_x(b)
                          for point in produce.pdp.lifts(curve, x)) for b in bases]
        projected = [[curve.scalar(point, 4) for point in factor] for factor in factors]
        pairs = [len({min(point, curve.neg(point)) for point in values if point is not None})
                 for values in projected]
        union = {min(point, curve.neg(point)) for values in projected
                 for point in values if point is not None}
        arms.append({"position": position, "slot_ranks": ranks, "joint_rank": 13,
                     "factor_sizes": list(map(len, factors)),
                     "factor_points": [[[x, y] for x, y in factor] for factor in factors],
                     "slot_signed_projected_counts": [len(set(values)) for values in projected],
                     "slot_pair_column_counts": pairs, "column_budget": sum(pairs),
                     "global_column_count": len(union),
                     "tuple_count": math.prod(map(len, factors))})
    return {"model": "n19-m6-beta3-rank2-plus-b18-factor-only-v1",
            "n": 19, "q": 130873, "beta": 3,
            "normal_coordinate_index": 18,
            "positions": list(range(6)),
            "shared_targets": produce.selected_targets(),
            "arms": arms,
            "no_tuple_sums_computed": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.write_text(json.dumps(build(), sort_keys=True, separators=(",", ":")) + "\n")

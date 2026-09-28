#!/usr/bin/env python3
"""Post-outcome exact count of producer integer convolution updates.

This does not alter the frozen six-arm census or its decision. It reconstructs
the support-only prefix sets and counts each labelled-factor update executed
by the frozen producer's sparse recurrence.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import produce

Q = 130873


def reconstruct(raw: Path) -> dict:
    f = produce.pdp.parent.Field(19, [0, 1, 2, 5])
    f.rabin_prime_degree()
    curve = produce.pdp.parent.Curve(f)
    h = produce.pdp.generator(curve, Q)
    assert h == (385982, 301867)
    conjugates = produce.pdp.parent.normal_conjugates(f, 3)
    lookup = {}
    current = None
    for k in range(Q):
        lookup[current] = k
        current = curve.add(current, h)
    assert current is None and len(lookup) == Q
    arms = []
    for position in range(6):
        result = json.loads((raw / f"p{position}" / "result.json").read_text())
        factors = []
        for i in range(6):
            basis = [conjugates[i + 6 * j] for j in range(3)]
            if i == position:
                basis.append(conjugates[18])
            factors.append(sorted(point for x in produce.pdp.parent.all_x(basis)
                                  for point in produce.pdp.lifts(curve, x)))
        assert [len(group) for group in factors] == result["factor_sizes"]
        prefix = {0}
        supports = [1]
        updates = []
        for factor in factors:
            shifts = [lookup[curve.scalar(point, 4)] for point in factor]
            updates.append(len(prefix) * len(shifts))
            prefix = {(k + shift) % Q for k in prefix for shift in shifts}
            supports.append(len(prefix))
        assert len(prefix) == result["support"]
        arms.append({"position": position, "prefix_supports": supports,
                     "integer_convolution_updates_by_stage": updates,
                     "integer_convolution_updates_total": sum(updates)})
    return {"status": "post-outcome cost reconstruction", "arms": arms}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.write_text(json.dumps(reconstruct(args.raw), sort_keys=True, separators=(",", ":")) + "\n")

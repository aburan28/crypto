#!/usr/bin/env python3
"""Workflow parameter files for ledger §20, by PROTOCOL.md's rules.

One file per (curve, column count, descent summands, seed set):

- base: `subgroup_orbits`, `2n · columns` points, seeded by the set;
- window `max(1, round(|F|/32))`;
- unit `clamp(round(2r / (n |F|^2)), 16, 65536)` probes, one planned unit,
  at most 100,000;
- descent cap `max(10^5, 64 · ceil(2.02 r / |F|^2))` trials;
- `m = 3` pair-table collection aimed at the least-mentioned columns;
- 32 targets, `random_seed` = 100 · seed + i.

    python3 make_params.py <a> <n> <columns> <descent m> <seed> <out.json>
"""
from __future__ import annotations

import json
import math
import sys

# r of every curve the protocol names, and of the smoke-test curve.
R = {
    (1, 17): 65587,
    (1, 19): 262543,
    (1, 23): 4196903,
    (1, 45): 29264761,
    (0, 37): 230603167,
    (1, 43): 4644189029,
    (1, 47): 106781081677,
    (0, 41): 549756390943,
    (0, 53): 21044858204113,
    (0, 61): 162888033982417,
}
TARGETS = 32


def params(a: int, n: int, columns: int, m: int, seed: int) -> dict:
    r = R[(a, n)]
    points = 2 * n * columns
    window = max(1, round(points / 32))
    unit = min(65536, max(16, round(2 * r / (n * points * points))))
    cap = max(100_000, 64 * math.ceil(2.02 * r / (points * points)))
    return {
        "schema_version": 1,
        "name": f"k{a}n{n}-c{columns}-m{m}-s{seed}",
        "curve": {"degree": n, "curve_a": a},
        "summands": 3,
        "descent_summands": m,
        "collection_window": window,
        "collection_aim": True,
        "solver": "pair_table",
        "seed": seed,
        "max_trials": min(cap, 2**32 - 1),
        "collection": {"unit_trials": unit, "units": 1, "max_units": 100_000},
        "factor_base": {"mode": "spec", "spec": {"kind": "subgroup_orbits", "points": points, "seed": seed}},
        "targets": [{"random_seed": 100 * seed + i} for i in range(TARGETS)],
    }


def main() -> None:
    a, n, columns, m, seed = (int(x) for x in sys.argv[1:6])
    with open(sys.argv[6], "x") as f:
        json.dump(params(a, n, columns, m, seed), f, indent=1)
        f.write("\n")


if __name__ == "__main__":
    main()

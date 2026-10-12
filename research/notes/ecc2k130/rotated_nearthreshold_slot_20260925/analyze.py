#!/usr/bin/env python3
"""Exact equal-budget same-q comparison of six frozen position arrays."""
from __future__ import annotations

import argparse
import itertools
import json
import struct
from pathlib import Path

Q = 130873


def compare(root: Path) -> dict:
    rows = [json.loads((root / f"p{position}" / "result.json").read_text())
            for position in range(6)]
    assert [row["position"] for row in rows] == list(range(6))
    arrays = []
    for position, row in enumerate(rows):
        raw = (root / f"p{position}" / "counts.u32le").read_bytes()
        assert len(raw) == 4 * Q
        values = struct.unpack(f"<{Q}I", raw)
        assert sum(values) == row["tuple_count"]
        assert sum(bool(v) for v in values) == row["support"]
        arrays.append(values)
    pairs = []
    for left, right in itertools.combinations(range(6), 2):
        a, b = rows[left], rows[right]
        equal = (a["tuple_count"], a["column_budget"]) == (b["tuple_count"], b["column_budget"])
        intersection = sum(bool(x) and bool(y) for x, y in zip(arrays[left], arrays[right]))
        symmetric = a["support"] + b["support"] - 2 * intersection
        membership_flips = sum(bool(arrays[left][k]) != bool(arrays[right][k])
                               for k in [r["k"] for r in a["sample_targets"]])
        assert [r["k"] for r in a["sample_targets"]] == [r["k"] for r in b["sample_targets"]]
        pairs.append({"positions": [left, right], "equal_tuple_and_column_budget": equal,
                      "tuple_counts": [a["tuple_count"], b["tuple_count"]],
                      "column_budgets": [a["column_budget"], b["column_budget"]],
                      "support_counts": [a["support"], b["support"]],
                      "intersection": intersection,
                      "support_symmetric_difference": symmetric,
                      "absolute_support_difference": abs(a["support"] - b["support"]),
                      "sample_membership_flips": membership_flips,
                      "sample_ranks": [a["sample_rank"], b["sample_rank"]],
                      "second_moments": [a["second_moment"], b["second_moment"]],
                      "passes_position_gate": equal and abs(a["support"] - b["support"]) >= 1309
                      and symmetric >= 2618})
    equal_pairs = [p for p in pairs if p["equal_tuple_and_column_budget"]]
    return {"q": Q, "positions": list(range(6)), "arms": [
        {k: row[k] for k in ("position", "slot_ranks", "factor_sizes", "slot_signed_projected_counts",
                             "slot_pair_column_counts", "column_budget", "global_column_count",
                             "tuple_count", "support", "misses", "second_moment", "collision_pairs",
                             "sample_support", "sample_rank", "counts_sha256")}
        for row in rows], "pairs": pairs,
        "equal_budget_pair_count": len(equal_pairs),
        "position_gate": "PASS" if any(p["passes_position_gate"] for p in equal_pairs)
        else "NO_EQUAL_BUDGET_PAIR" if not equal_pairs else "FAIL"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.write_text(json.dumps(compare(args.raw), sort_keys=True, separators=(",", ":")) + "\n")

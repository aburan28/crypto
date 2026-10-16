#!/usr/bin/env python3
"""Exact residual-column closure allowance for archived N83 pilot pairs."""

from __future__ import annotations

import argparse
import hashlib
import json
from math import comb
from pathlib import Path


HERE = Path(__file__).resolve().parent
SUPPORT = HERE.parent / "pilot-pair-support-screen-20261014/screen.json"
OUTPUT = HERE / "closure.json"
SUPPORT_SHA256 = "2c2229b9bfc78935cbc6c3ca7c6612f8fade2a40a0e0592f93a99b855d41dd45"
ORDER = 2417851639230796216685689
SHAPES = ((5, 1), (5, 2), (6, 1), (6, 2))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def ceil_div(numerator: int, denominator: int) -> int:
    require(numerator > 0 and denominator > 0, "positive fraction required")
    return (numerator + denominator - 1) // denominator


def half_count_floor(multisets: int, threshold: int) -> int:
    """Queries before Markov's bound permits half-chance of threshold hits."""
    return ceil_div(threshold * (ORDER - 1), 2 * multisets)


def canonical(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def calculate(support_path: Path = SUPPORT) -> dict:
    require(hashlib.sha256(support_path.read_bytes()).hexdigest() == SUPPORT_SHA256,
            "frozen archive-backed support receipt changed")
    source = json.loads(support_path.read_text())
    require(source.get("schema") == "n83.pilot-pair-uniform-query-support/v1"
            and source.get("status") == "EXACT_CONDITIONAL_COUNTING_BOUND"
            and source.get("subgroup_order") == str(ORDER)
            and source.get("points_per_signed_frobenius_orbit") == 166
            and source.get("ordered_pair_count") == 648
            and source.get("worker_admissible_pair_count") == 558
            and source.get("new_n83_solver_search_executed") is False
            and source.get("selected_best_total_runtime") is None,
            "parent support scope changed")
    entries = []
    for pair in source["ordered_pairs"]:
        base = pair["base_orbits"]
        pool = pair["pool_orbits"]
        residual = pair["residual_signed_orbits"]
        require(base in (64, 256, 600) and pool in (64, 256, 600)
                and 0 <= residual <= pool, "orbit counts changed")
        if residual < 2:
            require(pair["worker_admission"] == "UNAVAILABLE_RESIDUAL_LT_TWO"
                    and pair["shape_bounds"] == [], "unavailable pair changed")
            continue
        require(pair["worker_admission"] == "ADMISSIBLE_MIN_TWO_RESIDUAL_ORBITS"
                and len(pair["shape_bounds"]) == len(SHAPES),
                "admissible shapes changed")
        # For any t valid partial rows, the residual coefficient matrix has
        # at most L= residual orbit columns. Its left nullity is >= t-L.
        # These combinations may be dependent after projection to base columns.
        threshold = base + residual
        for bound, (summands, residual_summands) in zip(pair["shape_bounds"], SHAPES):
            small_summands = summands - residual_summands
            multiset_count = (comb(166 * base + small_summands - 1, small_summands)
                              * comb(166 * residual + residual_summands - 1,
                                     residual_summands))
            require((bound["summands"], bound["residual_summands"])
                    == (summands, residual_summands)
                    and bound["small_points"] == 166 * base
                    and bound["residual_points"] == 166 * residual
                    and int(bound["point_multisets_including_repetition"])
                    == multiset_count, "support count changed")
            floor = half_count_floor(multiset_count, threshold)
            require(2 * (floor - 1) * multiset_count < threshold * (ORDER - 1)
                    <= 2 * floor * multiset_count,
                    "Markov threshold arithmetic changed")
            entries.append({
                "base_object_blake3": pair["base_object_blake3"],
                "pool_object_blake3": pair["pool_object_blake3"],
                "base_policy": pair["base_policy"], "base_seed": pair["base_seed"],
                "base_orbits": base,
                "pool_policy": pair["pool_policy"], "pool_seed": pair["pool_seed"],
                "pool_orbits": pool, "residual_orbits": residual,
                "summands": summands, "residual_summands": residual_summands,
                "residual_nullity_threshold_partials": threshold,
                "half_threshold_decomposition_query_floor": floor,
            })
    require(len(entries) == 2232, "admissible pair/shape population changed")
    entries.sort(key=lambda row: (row["base_object_blake3"],
                                  row["pool_object_blake3"], row["summands"],
                                  row["residual_summands"]))
    selected = [row for row in entries
                if row["base_policy"] == "public_x_sequential"
                and row["base_seed"] == 2026100802
                and row["base_orbits"] == 600
                and row["pool_policy"] in ("public_x_gray_prefix", "public_x_hash")
                and row["pool_seed"] == 2026100801
                and row["pool_orbits"] == 600
                and (row["summands"], row["residual_summands"]) == (5, 2)]
    selected.sort(key=lambda row: row["pool_policy"])
    require(len(selected) == 2 and [row["pool_policy"] for row in selected]
            == ["public_x_gray_prefix", "public_x_hash"],
            "matched Gray/hash pair changed")
    summaries = []
    for summands, residual_summands in SHAPES:
        for base in (64, 256, 600):
            subset = [row for row in entries if row["base_orbits"] == base
                      and (row["summands"], row["residual_summands"])
                      == (summands, residual_summands)]
            floors = sorted(row["half_threshold_decomposition_query_floor"]
                            for row in subset)
            summaries.append({"base_orbits": base, "summands": summands,
                              "residual_summands": residual_summands,
                              "pairs": len(subset),
                              "lowest_floor": floors[0],
                              "lower_median_floor": floors[(len(floors) - 1) // 2],
                              "highest_floor": floors[-1]})
    return {
        "schema": "n83.pilot-pair-residual-nullity-allowance/v1",
        "status": "EXACT_LINEAR_ALGEBRA_AND_CONDITIONAL_COUNTING",
        "parent_support_sha256": SUPPORT_SHA256,
        "subgroup_order": str(ORDER),
        "admissible_pair_shape_count": len(entries),
        "full_table_sha256": hashlib.sha256(canonical(entries)).hexdigest(),
        "selected_gray_hash_3_plus_2": selected,
        "summaries": summaries,
        "new_n83_solver_search_executed": False,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    expected = canonical(calculate())
    if args.write:
        OUTPUT.write_bytes(expected)
    else:
        require(OUTPUT.read_bytes() == expected, "saved closure receipt changed")
    print(json.dumps({"status": "PASS", "receipt_sha256":
                      hashlib.sha256(expected).hexdigest()}))


if __name__ == "__main__":
    main()

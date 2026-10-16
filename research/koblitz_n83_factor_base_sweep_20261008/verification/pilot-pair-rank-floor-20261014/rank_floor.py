#!/usr/bin/env python3
"""Conditional uniform-query floors for K relation decompositions per pilot pair."""

import argparse
import hashlib
import json
from math import comb
from pathlib import Path


HERE = Path(__file__).resolve().parent
SUPPORT = HERE.parent / "pilot-pair-support-screen-20261014/screen.json"
OUTPUT = HERE / "rank-floor.json"
SUPPORT_SHA256 = "2c2229b9bfc78935cbc6c3ca7c6612f8fade2a40a0e0592f93a99b855d41dd45"
ORDER = 2417851639230796216685689
SHAPES = ((5, 1), (5, 2), (6, 1), (6, 2))
SIZES = (64, 256, 600)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ceil_div(a, b):
    return (a + b - 1) // b


def query_floor(order, multiset_count, desired_rows):
    """First A for which the Markov ceiling on >= K decompositions reaches 1/2."""
    require(order > 1 and multiset_count > 0 and desired_rows > 0,
            "positive order, multiset count and row goal required")
    return ceil_div(desired_rows * (order - 1), 2 * multiset_count)


def calculate():
    require(sha256(SUPPORT) == SUPPORT_SHA256,
            "frozen pilot-pair support receipt changed")
    source = json.loads(SUPPORT.read_text())
    require(source.get("schema") == "n83.pilot-pair-uniform-query-support/v1"
            and source.get("status") == "EXACT_CONDITIONAL_COUNTING_BOUND"
            and source.get("subgroup_order") == str(ORDER)
            and source.get("points_per_signed_frobenius_orbit") == 166
            and source.get("ordered_pair_count") == 648
            and source.get("worker_admissible_pair_count") == 558
            and source.get("new_n83_solver_search_executed") is False
            and source.get("selected_best_total_runtime") is None,
            "source support scope or population changed")
    entries = []
    for pair in source["ordered_pairs"]:
        base_orbits = pair["base_orbits"]
        residual_orbits = pair["residual_signed_orbits"]
        require(base_orbits in SIZES and 0 <= residual_orbits <=
                pair["pool_orbits"], "pair size changed")
        if residual_orbits < 2:
            require(pair["worker_admission"] ==
                    "UNAVAILABLE_RESIDUAL_LT_TWO" and
                    pair["shape_bounds"] == [],
                    "unavailable pair has support bounds")
            continue
        require(pair["worker_admission"] ==
                "ADMISSIBLE_MIN_TWO_RESIDUAL_ORBITS" and
                len(pair["shape_bounds"]) == 4,
                "admissible pair shape count changed")
        for bound, (summands, residual_summands) in zip(
                pair["shape_bounds"], SHAPES):
            small_summands = summands - residual_summands
            B = 166 * base_orbits
            L = 166 * residual_orbits
            M = (comb(B + small_summands - 1, small_summands)
                 * comb(L + residual_summands - 1, residual_summands))
            require((bound["summands"], bound["residual_summands"]) ==
                    (summands, residual_summands)
                    and bound["small_points"] == B
                    and bound["residual_points"] == L
                    and int(bound["point_multisets_including_repetition"]) == M
                    and bound["uniform_nonidentity_hit_ceiling"] == {
                        "numerator": str(min(M, ORDER - 1)),
                        "denominator": str(ORDER - 1)},
                    "support multiset count or query denominator differs")
            floor = query_floor(ORDER, M, base_orbits)
            require((floor - 1) * 2 * M < base_orbits * (ORDER - 1)
                    and floor * 2 * M >= base_orbits * (ORDER - 1),
                    "query floor arithmetic differs")
            entries.append({
                "base_object_blake3": pair["base_object_blake3"],
                "pool_object_blake3": pair["pool_object_blake3"],
                "base_policy": pair["base_policy"],
                "base_seed": pair["base_seed"],
                "base_orbits": base_orbits,
                "pool_policy": pair["pool_policy"],
                "pool_seed": pair["pool_seed"],
                "pool_orbits": pair["pool_orbits"],
                "residual_orbits": residual_orbits,
                "summands": summands,
                "residual_summands": residual_summands,
                "half_K_decomposition_query_floor": floor,
            })
    require(len(entries) == 558 * 4, "admissible pair/shape count differs")
    summaries = []
    for summands, residual_summands in SHAPES:
        for K in SIZES:
            values = sorted(entry["half_K_decomposition_query_floor"]
                            for entry in entries
                            if entry["base_orbits"] == K and
                            (entry["summands"], entry["residual_summands"])
                            == (summands, residual_summands))
            require(values, "empty rank-goal panel")
            middle = len(values) // 2
            summaries.append({
                "summands": summands,
                "residual_summands": residual_summands,
                "desired_decompositions_K": K,
                "pair_count": len(values),
                "minimum_query_floor": values[0],
                "lower_median_query_floor": values[middle - 1],
                "upper_median_query_floor": values[middle],
                "maximum_query_floor": values[-1],
                "pairs_with_floor_at_least_100": sum(v >= 100 for v in values),
                "pairs_with_floor_at_least_1000": sum(v >= 1000 for v in values),
            })
    require([s["pair_count"] for s in summaries[:3]] == [206, 186, 166],
            "base-size pair population differs")
    selected = [entry for entry in entries
                if entry["base_policy"] == "public_x_sequential"
                and entry["base_seed"] == 2026100802
                and entry["base_orbits"] == entry["pool_orbits"] == 600
                and entry["pool_seed"] == 2026100801
                and (entry["summands"], entry["residual_summands"])
                == (5, 2)
                and entry["pool_policy"] in
                    ("public_x_gray_prefix", "public_x_hash")]
    require(len(selected) == 2 and
            {(row["pool_policy"], row["residual_orbits"],
              row["half_K_decomposition_query_floor"]) for row in selected}
            == {("public_x_gray_prefix", 97, 33975),
                ("public_x_hash", 600, 889)},
            "selected same-base support contrast differs")
    canonical = sorted((row["base_object_blake3"],
                        row["pool_object_blake3"], row["summands"],
                        row["residual_summands"],
                        row["half_K_decomposition_query_floor"])
                       for row in entries)
    digest = hashlib.sha256(json.dumps(canonical, separators=(",", ":"))
                            .encode()).hexdigest()
    return {
        "schema": "n83.pilot-pair-K-decomposition-query-floor/v1",
        "status": "EXACT_CONDITIONAL_MARKOV_FLOOR",
        "parent_support_sha256": SUPPORT_SHA256,
        "source_sha256": sha256(Path(__file__)),
        "subgroup_order": str(ORDER),
        "admissible_pair_count": 558,
        "pair_shape_count": len(entries),
        "all_pair_shape_floors_sha256": digest,
        "shape_and_K_summaries": summaries,
        "selected_same_base_contrast": sorted(selected,
                                              key=lambda row: row["pool_policy"]),
        "proof_scope": {
            "query_law": "each query has a uniform nonidentity marginal conditional on fixed archived domains; independence between queries is unnecessary",
            "decomposition_count": "each point multiset contributes at most one folded relation row; rank cannot exceed decomposition count",
            "markov": "Pr(at least K decompositions after A queries) <= min(1,A*M/((r-1)*K))",
            "floor": "ceil(K*(r-1)/(2*M)) queries are necessary before this upper bound permits probability one half",
            "K": "hypothetical independent-row goal equal to small-base orbit-column count, not a proof that every one-target algorithm needs K rows",
            "excluded": ["nonuniform or target-adaptive query schedules",
                         "additional relation sources or combined search shapes",
                         "observed yield, useful rank or solver throughput",
                         "a cold-runtime or factor-base winner"],
        },
        "new_n83_solver_search_executed": False,
        "selected_best_total_runtime": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true",
                        help="create the immutable aggregate receipt")
    args = parser.parse_args()
    row = calculate()
    rendered = json.dumps(row, sort_keys=True, indent=2) + "\n"
    if args.write:
        with OUTPUT.open("x") as target:
            target.write(rendered)
    else:
        require(OUTPUT.read_text() == rendered,
                "frozen K-decomposition floor receipt differs")
    print(json.dumps({"status": row["status"],
                      "pair_shapes": row["pair_shape_count"],
                      "K600_five_two_best_floor":
                          next(item["minimum_query_floor"] for item in
                               row["shape_and_K_summaries"]
                               if (item["summands"],
                                   item["residual_summands"],
                                   item["desired_decompositions_K"]) ==
                               (5, 2, 600))}, sort_keys=True))


if __name__ == "__main__":
    main()

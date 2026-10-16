#!/usr/bin/env python3
"""Exact uniform-query support ceilings for archived N83 pilot pairs."""

import argparse
import hashlib
import json
from math import comb
from pathlib import Path


HERE = Path(__file__).resolve().parent
OVERLAP = HERE.parent / "pilot-overlap-census-20261013/overlap.json"
OUTPUT = HERE / "screen.json"
ORDER = 2417851639230796216685689
POINTS_PER_SIGNED_ORBIT = 166
SHAPES = ((5, 1), (5, 2), (6, 1), (6, 2))
CURVE = "icv1-f2m83-tm6151469093347-debefd74"
S3_PREFIX = (
    "s3://crypto-autoresearcher/factor-bases/icv1/etc/"
    "koblitz_n83_factor_base_sweep_20261008/a0/objects/"
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ceil_div(a, b):
    return (a + b - 1) // b


def multiset_count(base_points, residual_points, summands, residual_summands):
    """Count small/residual point multisets, including repeated points."""
    small_summands = summands - residual_summands
    require(base_points > 0 and residual_points > 0 and
            small_summands > 0 and residual_summands > 0,
            "positive disjoint-domain shape required")
    return (comb(base_points + small_summands - 1, small_summands)
            * comb(residual_points + residual_summands - 1,
                   residual_summands))


def support_bound(base_orbits, residual_orbits, summands,
                  residual_summands, order=ORDER):
    require(order > 1 and base_orbits > 0 and residual_orbits > 0,
            "invalid support input")
    base_points = POINTS_PER_SIGNED_ORBIT * base_orbits
    residual_points = POINTS_PER_SIGNED_ORBIT * residual_orbits
    count = multiset_count(base_points, residual_points, summands,
                           residual_summands)
    denominator = order - 1
    return {
        "summands": summands,
        "residual_summands": residual_summands,
        "small_points": base_points,
        "residual_points": residual_points,
        "point_multisets_including_repetition": str(count),
        "uniform_nonidentity_hit_ceiling": {
            "numerator": str(min(count, denominator)),
            "denominator": str(denominator),
        },
        "ceiling_nonvacuous": count < denominator,
        "minimum_queries_before_half_hit_ceiling_is_possible":
            ceil_div(denominator, 2 * count),
    }


def calculate():
    source = json.loads(OVERLAP.read_text())
    require(source.get("schema") == "n83.primary-pilot-exact-overlap-census/v1"
            and source.get("status") == "PASS_archive_backed_exact_orbit_overlap"
            and source.get("curve_slug") == CURVE
            and source.get("primary_archive_count") == 27
            and source.get("ordered_cross_stream_pair_count") == 648
            and source.get("passes_current_min_two_residual_screen_count") == 558
            and source.get("selected_best_total_runtime") is None,
            "archived overlap receipt does not match this screen")
    objects = source["archive_objects"]
    require(len(objects) == 27, "archived object count differs")
    by_digest = {}
    for obj in objects:
        digest = obj["compressed_blake3"]
        require(digest not in by_digest and obj["independent_full_point_row_replay"]
                and obj["s3_uri"] == S3_PREFIX + digest + ".jsonl.gz"
                and obj["orbit_columns"] in (64, 256, 600),
                "archived object identity differs")
        by_digest[digest] = obj
    pairs = source["ordered_pairs"]
    require(len(pairs) == 648, "ordered pair count differs")
    by_direction = {}
    for pair in pairs:
        base_hash = pair["base_object_blake3"]
        pool_hash = pair["pool_object_blake3"]
        key = (base_hash, pool_hash)
        require(base_hash in by_digest and pool_hash in by_digest
                and key not in by_direction and base_hash != pool_hash,
                "ordered pair identity differs")
        base, pool = by_digest[base_hash], by_digest[pool_hash]
        require((pair["base_policy"], pair["base_seed"],
                 pair["base_orbits"]) ==
                (base["policy"], base["seed"], base["orbit_columns"])
                and (pair["pool_policy"], pair["pool_seed"],
                     pair["pool_orbits"]) ==
                (pool["policy"], pool["seed"], pool["orbit_columns"])
                and (base["policy"], base["seed"]) !=
                (pool["policy"], pool["seed"]),
                "pair archive metadata differs")
        shared = pair["shared_signed_orbits"]
        residual = pair["residual_signed_orbits"]
        require(0 <= shared <= min(pair["base_orbits"], pair["pool_orbits"])
                and residual == pair["pool_orbits"] - shared
                and pair["shared_x_values"] == 83 * shared
                and pair["passes_current_min_two_residual_screen"] ==
                    (residual >= 2)
                and len(pair["source_derived_exact_model_shapes"]) ==
                    (4 if residual >= 2 else 0),
                "exact overlap or worker-admission invariant differs")
        by_direction[key] = pair
    require(len(by_direction) == 648, "ordered pair population differs")
    for (base_hash, pool_hash), pair in by_direction.items():
        reverse = by_direction.get((pool_hash, base_hash))
        require(reverse is not None and
                reverse["shared_signed_orbits"] == pair["shared_signed_orbits"],
                "reciprocal overlap differs")

    rows = []
    for pair in pairs:
        residual = pair["residual_signed_orbits"]
        row = {key: pair[key] for key in (
            "base_object_blake3", "pool_object_blake3", "base_policy",
            "base_seed", "base_orbits", "pool_policy", "pool_seed",
            "pool_orbits", "shared_signed_orbits", "residual_signed_orbits")}
        row["worker_admission"] = ("ADMISSIBLE_MIN_TWO_RESIDUAL_ORBITS"
                                   if residual >= 2 else
                                   "UNAVAILABLE_RESIDUAL_LT_TWO")
        row["shape_bounds"] = ([support_bound(pair["base_orbits"], residual,
                                               summands, large)
                                for summands, large in SHAPES]
                               if residual >= 2 else [])
        rows.append(row)
    summaries = []
    for summands, large in SHAPES:
        bounds = [bound for row in rows for bound in row["shape_bounds"]
                  if (bound["summands"], bound["residual_summands"]) ==
                  (summands, large)]
        counts = [int(bound["point_multisets_including_repetition"])
                  for bound in bounds]
        require(len(bounds) == 558, "shape population differs")
        summaries.append({
            "summands": summands,
            "residual_summands": large,
            "admissible_pair_count": len(bounds),
            "nonvacuous_hit_ceiling_count": sum(bound["ceiling_nonvacuous"]
                                                for bound in bounds),
            "minimum_point_multisets": str(min(counts)),
            "maximum_point_multisets": str(max(counts)),
            "maximum_hit_ceiling": {
                "numerator": str(min(max(counts), ORDER - 1)),
                "denominator": str(ORDER - 1),
            },
        })
    require(sum(row["worker_admission"] ==
                "UNAVAILABLE_RESIDUAL_LT_TWO" for row in rows) == 90,
            "nonadmitted pair count differs")
    return {
        "schema": "n83.pilot-pair-uniform-query-support/v1",
        "status": "EXACT_CONDITIONAL_COUNTING_BOUND",
        "curve_slug": CURVE,
        "subgroup_order": str(ORDER),
        "points_per_signed_frobenius_orbit": POINTS_PER_SIGNED_ORBIT,
        "overlap_receipt_sha256": sha256(OVERLAP),
        "source_sha256": sha256(Path(__file__)),
        "ordered_pair_count": len(rows),
        "worker_admissible_pair_count": 558,
        "shape_summaries": summaries,
        "ordered_pairs": rows,
        "proof_scope": {
            "query_distribution": "uniform nonidentity subgroup point, fixed archived domains",
            "point_multisets": "C(B+u-1,u)*C(L+v-1,v), repetitions allowed",
            "single_query_hit": "Pr(any decomposition) <= min(1,M/(r-1))",
            "multiple_query_floor": "ceil((r-1)/(2*M)) is necessary before the union-bound hit ceiling can reach one half; each query must have a uniform marginal",
            "excluded": ["frozen target or deterministic nonuniform schedule",
                         "observed solver success or relation yield",
                         "linear algebra, target logarithm or total runtime ranking"],
        },
        "new_n83_solver_search_executed": False,
        "selected_best_total_runtime": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true",
                        help="create the immutable support receipt")
    args = parser.parse_args()
    result = calculate()
    rendered = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if args.write:
        with OUTPUT.open("x") as target:
            target.write(rendered)
    else:
        require(OUTPUT.read_text() == rendered,
                "frozen support receipt differs from source or input")
    print(json.dumps({"status": result["status"],
                      "admissible_pairs": result["worker_admissible_pair_count"],
                      "shape_nonvacuous_counts":
                          [row["nonvacuous_hit_ceiling_count"]
                           for row in result["shape_summaries"]]},
                     sort_keys=True))


if __name__ == "__main__":
    main()

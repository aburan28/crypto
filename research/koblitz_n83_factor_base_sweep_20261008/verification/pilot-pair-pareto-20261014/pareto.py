#!/usr/bin/env python3
"""Recount the source-only 3+2 frontier of every retained ordered pilot pair."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
SCREEN = HERE.parent / "pilot-pair-support-screen-20261014" / "screen.json"
CAPACITY = HERE.parent / "paired-capacity-screen-20261013" / "screen.py"
MANIFEST = STUDY / "pilot-01" / "manifest.json"
RECEIPT = HERE / "receipt.json"
PINS = {
    SCREEN: "2c2229b9bfc78935cbc6c3ca7c6612f8fade2a40a0e0592f93a99b855d41dd45",
    CAPACITY: "116617e9de7552e48b548e632e80121508515e7079878305de0b975c8a924c5d",
    MANIFEST: "5bd862987844f709e03beeab9fd265e3d2024011bcadcb82f5af9c72018df6a6",
}
BASE = "3888e2597fd3519807c9d722f49f83329c93ce5cb0929501cf7cb7cbc3cab1a0"
POOLS = {
    "gray_pool": "966a0493c96674f0af6d184615ed05309d08fba5c7342957b4186cd50366dc47",
    "hash_pool": "b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e",
}


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_capacity():
    spec = importlib.util.spec_from_file_location("n83_pareto_capacity", CAPACITY)
    require(spec is not None and spec.loader is not None, "capacity source missing")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def dominates(left: tuple[int, int, int], right: tuple[int, int, int]) -> bool:
    """Strict product-order dominance; equal metric triples remain ties."""
    return all(a <= b for a, b in zip(left, right)) and any(
        a < b for a, b in zip(left, right))


def pair_metrics(pair: dict, order: int, capacity) -> dict:
    k, residual = pair["base_orbits"], pair["residual_signed_orbits"]
    require(type(k) is int and type(residual) is int and k > 0 and residual >= 2,
            "admitted pair shape malformed")
    shapes = [shape for shape in pair["shape_bounds"]
              if shape["summands"] == 5 and shape["residual_summands"] == 2]
    require(len(shapes) == 1, "five-summand 3+2 shape absent or duplicated")
    shape = shapes[0]
    base_points, residual_points = 166 * k, 166 * residual
    require(shape["small_points"] == base_points
            and shape["residual_points"] == residual_points,
            "support domain point counts changed")
    multisets = comb(base_points + 2, 3) * comb(residual_points + 1, 2)
    require(shape["point_multisets_including_repetition"] == str(multisets),
            "support multiset formula changed")
    query_floor = (k * (order - 1) + 2 * multisets - 1) // (2 * multisets)
    source = capacity.requirements(k, k + residual, 5, 2)
    require(source["current_worker_cap_pass"] is True,
            "admitted pair exceeds current worker cap")
    return {
        "base_object_blake3": pair["base_object_blake3"],
        "pool_object_blake3": pair["pool_object_blake3"],
        "base_orbits": k,
        "residual_signed_orbits": residual,
        "conditional_half_K_row_query_floor": query_floor,
        "source_domain_clauses": source["domain_clauses"],
        "source_total_variables": source["total_variables"],
    }


def metrics(row: dict) -> tuple[int, int, int]:
    return (row["conditional_half_K_row_query_floor"],
            row["source_domain_clauses"], row["source_total_variables"])


def make_receipt() -> dict:
    for path, expected in PINS.items():
        require(sha256(path) == expected, f"pinned source bytes changed: {path.name}")
    screen = json.loads(SCREEN.read_text())
    manifest = json.loads(MANIFEST.read_text())
    capacity = load_capacity()
    require(screen["schema"] == "n83.pilot-pair-uniform-query-support/v1"
            and screen["status"] == "EXACT_CONDITIONAL_COUNTING_BOUND"
            and screen["ordered_pair_count"] == 648
            and screen["worker_admissible_pair_count"] == 558,
            "parent pair census changed")
    order = int(screen["subgroup_order"])
    require(order == 2417851639230796216685689,
            "primary subgroup order changed")
    objects = {row["compressed_blake3"]: row for row in manifest["bases"]}
    require(len(objects) == len(manifest["bases"]), "manifest duplicate object identity")
    seen = set()
    admitted = []
    unavailable = 0
    for pair in screen["ordered_pairs"]:
        key = pair["base_object_blake3"], pair["pool_object_blake3"]
        require(key not in seen and key[0] != key[1], "duplicate or self pair")
        seen.add(key)
        require(all(digest in objects for digest in key), "pair object missing from manifest")
        base, pool = (objects[digest] for digest in key)
        for prefix, row in (("base", base), ("pool", pool)):
            require(pair[f"{prefix}_orbits"] == row["columns"]
                    and pair[f"{prefix}_policy"] == row["policy"]
                    and pair[f"{prefix}_seed"] == row["seed"],
                    "pair identity differs from stored object")
        if pair["worker_admission"] == "UNAVAILABLE_RESIDUAL_LT_TWO":
            require(pair["residual_signed_orbits"] < 2 and not pair["shape_bounds"],
                    "unavailable pair carries a shape")
            unavailable += 1
            continue
        require(pair["worker_admission"] == "ADMISSIBLE_MIN_TWO_RESIDUAL_ORBITS",
                "unknown worker admission")
        admitted.append(pair_metrics(pair, order, capacity))
    require(len(seen) == 648 and len(admitted) == 558 and unavailable == 90,
            "pair partition changed")
    admitted.sort(key=lambda row: (row["base_object_blake3"], row["pool_object_blake3"]))
    shapes = [metrics(row) for row in admitted]
    frontier_flags = [not any(dominates(other, shape) for other in shapes)
                      for shape in shapes]
    frontier = defaultdict(list)
    for row, keep in zip(admitted, frontier_flags):
        if keep:
            frontier[metrics(row)].append(row)
    frontier_shapes = []
    for shape, group in sorted(frontier.items()):
        example = group[0]
        frontier_shapes.append({
            "conditional_half_K_row_query_floor": shape[0],
            "source_domain_clauses": shape[1],
            "source_total_variables": shape[2],
            "tied_pair_count": len(group),
            "example_base_object_blake3": example["base_object_blake3"],
            "example_pool_object_blake3": example["pool_object_blake3"],
            "example_base_orbits": example["base_orbits"],
            "example_residual_signed_orbits": example["residual_signed_orbits"],
        })
    selected = {}
    for role, pool_hash in POOLS.items():
        row = next((row for row in admitted if row["base_object_blake3"] == BASE
                    and row["pool_object_blake3"] == pool_hash), None)
        require(row is not None, "selected paired arm absent")
        dominators = [other for other in admitted if dominates(metrics(other), metrics(row))]
        dominators.sort(key=lambda other: (metrics(other), other["base_object_blake3"],
                                           other["pool_object_blake3"]))
        selected[role] = {
            **row, "source_pareto": not dominators,
            "strict_dominator_pair_count": len(dominators),
            "best_query_floor_dominator": dominators[0] if dominators else None,
        }
    require(len(frontier_shapes) == 24
            and sum(row["tied_pair_count"] for row in frontier_shapes) == 311,
            "source frontier census changed")
    return {
        "schema": "n83.pilot-pair-source-pareto/v1",
        "status": "PASS_SOURCE_ONLY_NO_RUNTIME_WINNER",
        "curve_slug": screen["curve_slug"],
        "subgroup_order": screen["subgroup_order"],
        "screen_sha256": PINS[SCREEN],
        "capacity_source_sha256": PINS[CAPACITY],
        "manifest_sha256": PINS[MANIFEST],
        "source_sha256": sha256(Path(__file__)),
        "ordered_pair_count": len(seen),
        "worker_admissible_pair_count": len(admitted),
        "unavailable_pair_count": unavailable,
        "source_pareto_pair_count": sum(frontier_flags),
        "source_pareto_unique_shape_count": len(frontier_shapes),
        "frontier_shapes": frontier_shapes,
        "selected_paired_arms": selected,
        "new_n83_solver_search_executed": False,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = make_receipt()
    encoded = (json.dumps(receipt, sort_keys=True, indent=2) + "\n").encode()
    if args.write:
        RECEIPT.write_bytes(encoded)
    else:
        require(RECEIPT.read_bytes() == encoded, "saved Pareto receipt differs from replay")
    print(json.dumps({key: receipt[key] for key in (
        "status", "worker_admissible_pair_count", "source_pareto_pair_count",
        "source_pareto_unique_shape_count")}, sort_keys=True))


if __name__ == "__main__":
    main()

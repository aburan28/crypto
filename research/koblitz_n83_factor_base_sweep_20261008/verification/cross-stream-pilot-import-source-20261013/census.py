#!/usr/bin/env python3
"""Reproduce the pinned primary-arm pilot cross-stream addressability census."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
PANEL = HERE.parent.parent / "pilot-01" / "manifest.json"
RECEIPT = HERE / "census.json"
PANEL_SHA256 = "5bd862987844f709e03beeab9fd265e3d2024011bcadcb82f5af9c72018df6a6"
POLICIES = {"public_x_sequential", "public_x_hash", "public_x_gray_prefix"}
SIZES = {64, 256, 600}
SEEDS = {2026100801, 2026100802, 2026100803}


def generate(panel_bytes: bytes) -> dict:
    if hashlib.sha256(panel_bytes).hexdigest() != PANEL_SHA256:
        raise ValueError("pilot manifest SHA-256 changed")
    panel = json.loads(panel_bytes)
    rows = panel["bases"]
    if (
        panel["schema"] != "n83.factor-base-panel/v1"
        or panel["status"] != "completed_factor_base_panel"
        or panel["completed_base_count"] != 54
        or panel["expected_base_count"] != 54
        or len(rows) != 54
    ):
        raise ValueError("pilot panel metadata changed")
    primary = [row for row in rows if row["a"] == 0]
    if len(primary) != 27:
        raise ValueError("primary-arm row count changed")
    hashes = set()
    for row in primary:
        digest = row["compressed_blake3"]
        if (
            row["policy"] not in POLICIES
            or row["seed"] not in SEEDS
            or row["columns"] not in SIZES
            or row["points"] != 166 * row["columns"]
            or row["object"] != f"objects/{digest}.jsonl.gz"
            or digest in hashes
        ):
            raise ValueError("primary-arm row metadata changed")
        hashes.add(digest)
    primary.sort(key=lambda row: (row["policy"], row["seed"], row["columns"]))
    relation = {"small_lt_pool": 0, "small_eq_pool": 0, "small_gt_pool": 0}
    stream = {"policy_only": 0, "seed_only": 0, "policy_and_seed": 0}
    for small, pool in itertools.permutations(primary, 2):
        policy_differs = small["policy"] != pool["policy"]
        seed_differs = small["seed"] != pool["seed"]
        if not (policy_differs or seed_differs):
            continue
        relation[
            "small_lt_pool"
            if small["columns"] < pool["columns"]
            else "small_gt_pool"
            if small["columns"] > pool["columns"]
            else "small_eq_pool"
        ] += 1
        stream[
            "policy_and_seed"
            if policy_differs and seed_differs
            else "policy_only"
            if policy_differs
            else "seed_only"
        ] += 1
    if relation != {key: 216 for key in relation} or stream != {
        "policy_only": 162,
        "seed_only": 162,
        "policy_and_seed": 324,
    }:
        raise ValueError("cross-stream Cartesian count changed")
    return {
        "schema": "n83.cross-stream-pilot-source-census/v1",
        "status": "source_addressable_unmeasured",
        "pilot_manifest_sha256": PANEL_SHA256,
        "primary_rows": len(primary),
        "other_arm_rows_excluded": len(rows) - len(primary),
        "ordered_cross_stream_pairs": sum(relation.values()),
        "size_relation_counts": relation,
        "stream_difference_counts": stream,
        "primary_objects": [
            {
                key: row[key]
                for key in ("policy", "seed", "columns", "compressed_blake3", "s3_uri")
            }
            for row in primary
        ],
        "archive_import_executed": False,
        "solver_search_executed": False,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = (json.dumps(generate(PANEL.read_bytes()), indent=2, sort_keys=True) + "\n").encode()
    if args.check:
        if RECEIPT.read_bytes() != output:
            raise SystemExit("FAIL: cross-stream census receipt differs")
        print("PASS: pinned pilot cross-stream source census")
    else:
        RECEIPT.write_bytes(output)


if __name__ == "__main__":
    main()

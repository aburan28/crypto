#!/usr/bin/env python3
"""Admit and materialize a freshly constructed native n53 base in exact order."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

BASE_ARRAYS_SHA256 = "2de8ec46916999dfe5e1ac95e68ea15a3f001eec27e013fa96d9525719b4f684"
POINT_KEYS = ("factor_base_point_coordinates", "factor_base_point_labels", "factor_base_representatives")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode() + b"\n"


def materialize(training_raw: Path, header_path: Path, receipt_path: Path) -> dict:
    raw_bytes = training_raw.read_bytes()
    raw = json.loads(raw_bytes)
    assert all(raw.get(key) is None for key in (
        "factor_base_input_path", "factor_base_input_hash", "factor_base_input_blake3"
    )), "training used a pre-existing base"
    generated = raw["compact_orbit_base_header"]
    assert generated["kind"] == "point_defined_factor_base"
    assert generated["selection_mode"] == "ascending_x_v1"
    assert (generated["n"], generated["a"], generated["cofactor"],
            generated["orbit_columns"], generated["factor_base_points"],
            generated["signed_automorphism_size"], generated["scanned_x"]) == (
                53, 0, 428, 220, 23320, 106, 465
            )
    assert generated["subgroup_order"] == 21044858204113
    assert generated["generator"] == [198217578752339, 7929897206038174]
    arrays = {key: generated[key] for key in POINT_KEYS}
    array_hash = sha(canonical(arrays))
    assert array_hash == BASE_ARRAYS_SHA256, "fresh native base differs from #823 exact ordered arrays"
    points = [tuple(point) for point in arrays["factor_base_point_coordinates"]]
    labels = [tuple(label) for label in arrays["factor_base_point_labels"]]
    reps = [tuple(point) for point in arrays["factor_base_representatives"]]
    assert len(points) == len(labels) == len(set(points)) == 23320
    assert len(reps) == len(set(reps)) == 220
    assert all(len(point) == 2 for point in points + reps)
    assert all(0 <= col < 220 and 0 < coeff < generated["subgroup_order"] for col, coeff in labels)
    fresh = dict(generated)
    # The loader reads this field but does not calculate it. Give it a pinned
    # SHA-256 of the exact ordered native arrays, with the algorithm explicit.
    fresh["base_hash"] = array_hash
    fresh["base_hash_algorithm"] = "sha256_canonical_ordered_point_arrays_v1"
    fresh["field_x_values_scanned"] = generated["scanned_x"]
    header_bytes = canonical(fresh)
    assert not header_path.exists() and not receipt_path.exists()
    header_path.write_bytes(header_bytes)
    receipt = {
        "classification": "FRESH_NATIVE_BASE_MATCHES_PINNED_ORDER",
        "training_raw_sha256": sha(raw_bytes),
        "fresh_header_sha256": sha(header_bytes),
        "ordered_point_arrays_sha256": array_hash,
        "base_hash_algorithm": fresh["base_hash_algorithm"],
        "points": 23320, "orbit_columns": 220, "scanned_x": 465,
        "selection_mode": generated["selection_mode"],
    }
    receipt_path.write_bytes(json.dumps(receipt, indent=2, sort_keys=True).encode() + b"\n")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-raw", type=Path, required=True)
    parser.add_argument("--header", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(materialize(args.training_raw, args.header, args.receipt), sort_keys=True))


if __name__ == "__main__":
    main()

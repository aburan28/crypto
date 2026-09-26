#!/usr/bin/env python3
"""Materialize a newly constructed n53 base and compare it to the certified arm.

The certified header is a read-only validator.  All points and labels written
for the point-query child come from the fresh training producer's output.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

BASE_HASH = "d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71"
FIXTURE_GZIP_HASH = "23397af2ef668aed0775bcb409e1ae19555357ded635452818c9a3812f679d08"
MATCH_KEYS = (
    "n", "a", "subgroup_order", "field_modulus_low_terms",
    "generator", "orbit_columns", "signed_automorphism_size",
    "factor_base_points", "factor_base_point_coordinates",
    "factor_base_point_labels", "factor_base_representatives",
)
POINT_KEYS = (
    "factor_base_point_coordinates", "factor_base_point_labels",
    "factor_base_representatives",
)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True).encode("ascii") + b"\n"


def materialize(training_raw: Path, fixture_gzip: Path, header_path: Path,
                receipt_path: Path) -> dict:
    raw = json.loads(training_raw.read_bytes())
    if any(raw.get(key) is not None for key in (
        "factor_base_input_path", "factor_base_input_hash", "factor_base_input_blake3"
    )):
        raise AssertionError("training loaded a pre-existing factor base")
    generated = raw.get("compact_orbit_base_header")
    if not isinstance(generated, dict) or generated.get("kind") != "point_defined_factor_base":
        raise AssertionError("fresh producer did not emit a complete base header")
    fixture_gz = fixture_gzip.read_bytes()
    if sha(fixture_gz) != FIXTURE_GZIP_HASH:
        raise AssertionError("certified fixture gzip SHA drift")
    if generated.get("selection_mode") != "legacy_rank_fixture_lcg_v1":
        raise AssertionError("training did not use the frozen rank-fixture LCG rule")
    fixture_bytes = gzip.decompress(fixture_gz)
    fixture = json.loads(fixture_bytes)
    if fixture.get("base_hash") != BASE_HASH:
        raise AssertionError("certified point-set hash drift")
    if generated.get("cofactor") != int(fixture["cofactor"]):
        raise AssertionError("constructed cofactor differs from certified base")
    for key in MATCH_KEYS:
        if generated.get(key) != fixture.get(key):
            raise AssertionError(f"constructed base differs from certified base: {key}")
    if generated.get("scanned_x") != fixture.get("field_x_values_scanned"):
        raise AssertionError("fresh base scan bound differs from certified arm")
    if (generated["n"], generated["a"], generated["orbit_columns"],
            generated["factor_base_points"], generated["signed_automorphism_size"]) != (53, 0, 220, 23320, 106):
        raise AssertionError("wrong n53 base geometry")
    point_serialization = canonical({key: generated[key] for key in POINT_KEYS})
    if point_serialization != canonical({key: fixture[key] for key in POINT_KEYS}):
        raise AssertionError("canonical points/labels/representatives bytes differ")
    pairs = list(zip(generated["factor_base_point_coordinates"],
                     generated["factor_base_point_labels"]))
    assert len(pairs) == len({tuple(point) for point, _ in pairs}) == 23320
    semantic_rows = sorted([*point, *label] for point, label in pairs)
    semantic_sha = sha(canonical(semantic_rows))
    assert semantic_sha == sha(canonical(sorted(
        [*point, *label] for point, label in zip(
            fixture["factor_base_point_coordinates"],
            fixture["factor_base_point_labels"]))))
    fresh = dict(generated)
    # The producer owns the value; the certified header serialized it as text.
    fresh["cofactor"] = str(generated["cofactor"])
    fresh["base_hash"] = BASE_HASH
    fresh["field_x_values_scanned"] = generated["scanned_x"]
    header_bytes = canonical(fresh)
    header_path.write_bytes(header_bytes)
    result = {
        "classification": "FRESH_BASE_MATCHES_CERTIFIED_ARM",
        "training_raw_sha256": sha(training_raw.read_bytes()),
        "certified_fixture_gzip_sha256": sha(fixture_gz),
        "certified_fixture_header_sha256": sha(fixture_bytes),
        "canonical_point_data_sha256": sha(point_serialization),
        "semantic_point_label_map_sha256": semantic_sha,
        "point_order_identical": True,
        "selection_mode": generated["selection_mode"],
        "fresh_header_sha256": sha(header_bytes),
        "point_set_sha256": BASE_HASH,
        "points": 23320, "orbit_columns": 220,
        "scanned_x": generated["scanned_x"],
        "match_keys": list(MATCH_KEYS),
    }
    receipt_path.write_bytes(json.dumps(result, indent=2, sort_keys=True).encode() + b"\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-raw", type=Path, required=True)
    parser.add_argument("--fixture-gzip", type=Path, required=True)
    parser.add_argument("--header", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    result = materialize(args.training_raw, args.fixture_gzip, args.header, args.receipt)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact uniform-query capacity and rank floors for audited N83 standard bases."""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
ORDER = 2417851639230796216685689
CURVE = "icv1-f2m83-tm6151469093347-debefd74"
INPUTS = {
    "construction": ("standard-d1-d17-s3-20261015/construction.json",
                     "f601a75a66d0cf30bd33ee574e47880b842f4f82174c5e1c974147f3d2621784"),
    "frontier_s3": ("standard-d1-d17-s3-20261015/s3-replay.json",
                    "fba81340cde66a0df6622fd6ce773f89a4c13524caa5350300602ffe13f24e32"),
    "frontier_audit": ("standard-d1-d17-s3-20261015/audit-s3.json",
                       "c1e9933eb900f7587ecf7d480c8645c94c5d2ee2f24fa566eaff932cd14e0eef"),
    "d18_manifest": ("standard-d18-s3-20261015/manifest.json",
                     "ce635f6f7e2d2cb8ac2be7272fb91276784ce81eaad05927fc209b8035bbc419"),
    "d18_s3": ("standard-d18-s3-20261015/s3-replay.json",
               "bebaccfa61510703a7a95b35b938e8275997a36b3ef7e6a876dfaf7c5ce5e695"),
    "d18_audit": ("standard-d18-s3-20261015/audit.json",
                  "e3ddfa8c3f0e75426ac4ea4d18be5e6a0099fc1cbb946b9d8e25feeec3a4fef3"),
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_pinned(path: Path, digest: str) -> dict:
    data = path.read_bytes()
    if sha256(data) != digest:
        raise ValueError(f"pinned input changed: {path}")
    return json.loads(data)


def admitted_manifests() -> list[dict]:
    records = {name: load_pinned(PARENT / path, digest)
               for name, (path, digest) in INPUTS.items()}
    constructed = records["construction"]
    s3 = records["frontier_s3"]
    audited = records["frontier_audit"]
    if (constructed.get("schema") != "n83.standard-subspace-frontier-construction/v1"
            or s3.get("status") != "PASS_ALL_16_WRITE_ONCE_ARCHIVES"
            or audited.get("status") != "PASS_ALL_VALID_DIMENSIONS"):
        raise ValueError("dimension-frontier evidence is not admitted")
    c_rows = {row["dimension"]: row for row in constructed["rows"]}
    s_rows = {row["dimension"]: row for row in s3["rows"]}
    a_rows = {row["dimension"]: row for row in audited["rows"]}
    if (set(c_rows) != set(range(1, 18))
            or set(s_rows) != set(range(2, 18))
            or set(a_rows) != set(range(2, 18))
            or c_rows[1]["status"] != "CONSTRUCT_FAILURE"):
        raise ValueError("standard frontier dimensions changed")
    manifests = []
    for dimension in range(2, 18):
        build, remote, check = c_rows[dimension], s_rows[dimension], a_rows[dimension]
        manifest = build["manifest"]
        if (build["status"] != "PASS_NATIVE_CONSTRUCT_AND_EXACT_REPLAY"
                or remote["status"] != "PASS_S3_BYTES_NATIVE_REPLAY"
                or check["status"] != "PASS_BYTES_LABELS_SAMPLED_GROUP"
                or remote["object_bytes"] != manifest["compressed_bytes"]
                or not remote["object_key"].endswith(manifest["object"])
                or remote["object_sha256"] != check["compressed_sha256"]
                or remote["manifest_sha256"] != check["manifest_sha256"]
                or check["compressed_blake3"] != manifest["compressed_blake3"]
                or check["projected_point_count"] != manifest["projected_point_count"]
                or check["signed_columns"] != manifest["signed_columns"]):
            raise ValueError(f"dimension {dimension} archive evidence disagrees")
        manifests.append(manifest)
    m18 = records["d18_manifest"]
    s18 = records["d18_s3"]
    a18 = records["d18_audit"]
    if (s18.get("status") != "PASS_FULL_OBJECT_AND_MANIFEST_BYTES"
            or a18.get("status") != "PASS_BYTES_ALL_SIGN_LABELS_SAMPLED_GROUP"
            or s18["items"]["manifest"]["sha256"] != INPUTS["d18_manifest"][1]
            or s18["items"]["object"]["bytes"] != m18["compressed_bytes"]
            or s18["items"]["object"]["sha256"] != a18["compressed_sha256"]
            or a18["compressed_blake3"] != m18["compressed_blake3"]
            or a18["point_rows"] != m18["projected_point_count"]
            or a18["sign_columns_checked"] != m18["signed_columns"]):
        raise ValueError("dimension 18 archive evidence disagrees")
    manifests.append(m18)
    for dimension, manifest in enumerate(manifests, 2):
        points = manifest["projected_point_count"]
        columns = manifest["signed_columns"]
        if (manifest["schema"] != "n83.standard-subspace-factor-base-manifest/v1"
                or manifest["curve_slug"] != CURVE
                or manifest["source_dimension"] != dimension
                or manifest["selected_best_total_runtime"] is not None
                or type(points) is not int or points <= 0
                or type(columns) is not int or 2 * columns != points):
            raise ValueError(f"dimension {dimension} manifest identity changed")
    return manifests


def rank_floor(columns: int, multisets: int, order: int = ORDER) -> int:
    """First query count where Markov permits probability 1/2 of `columns` rows."""
    if columns <= 0 or multisets <= 0 or order <= 1:
        raise ValueError("rank-floor arguments must be positive")
    numerator = columns * (order - 1)
    denominator = 2 * multisets
    return (numerator + denominator - 1) // denominator


def decimal_ratio(numerator: int, denominator: int) -> str:
    with localcontext() as context:
        context.prec = 40
        return format(Decimal(numerator) / Decimal(denominator), ".8g")


def generate() -> dict:
    manifests = admitted_manifests()
    rows = []
    for manifest in manifests:
        dimension = manifest["source_dimension"]
        points = manifest["projected_point_count"]
        columns = manifest["signed_columns"]
        for arity in range(2, 7):
            multisets = math.comb(points + arity - 1, arity)
            mean = decimal_ratio(multisets, ORDER - 1)
            rows.append({
                "dimension": dimension,
                "arity": arity,
                "projected_points": points,
                "sign_only_columns": columns,
                "base_compressed_blake3": manifest["compressed_blake3"],
                "unordered_multisets": str(multisets),
                "uniform_nonidentity_mean_display": mean,
                "single_query_hit_ceiling_numerator": str(min(multisets, ORDER - 1)),
                "single_query_hit_ceiling_denominator": str(ORDER - 1),
                "queries_before_half_hit_ceiling": rank_floor(1, multisets),
                "queries_before_half_column_row_ceiling": rank_floor(columns, multisets),
                "nominal_source_coordinate_bits": arity * dimension if arity in (5, 6) else None,
                "nominal_all_affine_chain_variables":
                    arity * dimension + 83 * (arity - 2) if arity in (5, 6) else None,
                "nominal_all_affine_chain_equations":
                    83 * (arity - 1) if arity in (5, 6) else None,
            })
    canonical = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode()
    return {
        "schema": "n83.standard-subspace-capacity-rank/v1",
        "status": "PASS_EXACT_CONDITIONAL_CAPACITY_SCREEN",
        "curve_slug": CURVE,
        "subgroup_order": str(ORDER),
        "uniform_query_domain": "nonidentity_primary_subgroup",
        "factor_base_domain": "cofactor_projected_standard_polynomial_subspace",
        "row_goal": "one_independent_row_per_sign_only_column",
        "input_sha256": {name: digest for name, (_, digest) in INPUTS.items()},
        "generator_sha256": sha256(Path(__file__).read_bytes()),
        "full_table_sha256": sha256(canonical),
        "rows": rows,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = generate()
    content = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode()
    receipt = HERE / "screen.json"
    if args.write:
        with receipt.open("xb") as stream:
            stream.write(content)
    elif receipt.read_bytes() != content:
        raise ValueError("screen receipt differs from exact replay")
    print(json.dumps({"status": "PASS_EXACT_CONDITIONAL_CAPACITY_SCREEN",
                      "cells": len(result["rows"]), "receipt_sha256": sha256(content)}))


if __name__ == "__main__":
    main()

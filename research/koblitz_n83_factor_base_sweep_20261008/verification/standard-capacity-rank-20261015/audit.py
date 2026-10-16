#!/usr/bin/env python3
"""Independently replay all 85 archived-base capacity cells with integer inequalities."""

from __future__ import annotations

from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
ORDER = 2417851639230796216685689
FILES = {
    "construction": "standard-d1-d17-s3-20261015/construction.json",
    "frontier_s3": "standard-d1-d17-s3-20261015/s3-replay.json",
    "frontier_audit": "standard-d1-d17-s3-20261015/audit-s3.json",
    "d18_manifest": "standard-d18-s3-20261015/manifest.json",
    "d18_s3": "standard-d18-s3-20261015/s3-replay.json",
    "d18_audit": "standard-d18-s3-20261015/audit.json",
}


def insist(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def independent_multisets(points: int, arity: int) -> int:
    numerator = 1
    for term in range(points, points + arity):
        numerator *= term
    quotient, remainder = divmod(numerator, math.factorial(arity))
    insist(remainder == 0, "multiset count is not integral")
    return quotient


def audit(receipt: dict) -> dict:
    insist(receipt.get("schema") == "n83.standard-subspace-capacity-rank/v1"
           and receipt.get("status") == "PASS_EXACT_CONDITIONAL_CAPACITY_SCREEN"
           and receipt.get("subgroup_order") == str(ORDER)
           and receipt.get("row_goal") == "one_independent_row_per_sign_only_column"
           and receipt.get("uniform_query_domain") == "nonidentity_primary_subgroup"
           and receipt.get("selected_best_total_runtime") is None,
           "screen scope or interpretation changed")
    sources = {}
    for name, relative in FILES.items():
        raw = (PARENT / relative).read_bytes()
        insist(sha256(raw) == receipt["input_sha256"][name], f"{name} changed")
        sources[name] = json.loads(raw)
    insist(sha256((HERE / "screen.py").read_bytes()) == receipt["generator_sha256"],
           "generator source changed")
    frontier = sources["construction"]
    built = {row["dimension"]: row["manifest"] for row in frontier["rows"]
             if row["status"] == "PASS_NATIVE_CONSTRUCT_AND_EXACT_REPLAY"}
    built[18] = sources["d18_manifest"]
    remote = {row["dimension"]: row for row in sources["frontier_s3"]["rows"]}
    checked = {row["dimension"]: row for row in sources["frontier_audit"]["rows"]}
    insist(set(built) == set(range(2, 19))
           and set(remote) == set(range(2, 18))
           and set(checked) == set(range(2, 18))
           and sources["frontier_s3"]["status"] == "PASS_ALL_16_WRITE_ONCE_ARCHIVES"
           and sources["frontier_audit"]["status"] == "PASS_ALL_VALID_DIMENSIONS"
           and sources["d18_s3"]["status"] == "PASS_FULL_OBJECT_AND_MANIFEST_BYTES"
           and sources["d18_audit"]["status"] == "PASS_BYTES_ALL_SIGN_LABELS_SAMPLED_GROUP",
           "audited archive population changed")
    rows = receipt.get("rows")
    insist(isinstance(rows, list) and len(rows) == 85, "screen cell count differs")
    seen = set()
    for row in rows:
        dimension, arity = row["dimension"], row["arity"]
        insist(type(dimension) is int and 2 <= dimension <= 18
               and type(arity) is int and 2 <= arity <= 6
               and (dimension, arity) not in seen,
               "duplicate or invalid dimension/arity cell")
        seen.add((dimension, arity))
        manifest = built[dimension]
        points, columns = manifest["projected_point_count"], manifest["signed_columns"]
        if dimension < 18:
            insist(remote[dimension]["object_sha256"] == checked[dimension]["compressed_sha256"]
                   and checked[dimension]["compressed_blake3"] == manifest["compressed_blake3"],
                   "S3/base audit binding differs")
        else:
            insist(sources["d18_s3"]["items"]["object"]["sha256"]
                   == sources["d18_audit"]["compressed_sha256"]
                   and sources["d18_audit"]["compressed_blake3"] == manifest["compressed_blake3"],
                   "dimension 18 S3/base binding differs")
        multisets = independent_multisets(points, arity)
        insist(row["projected_points"] == points
               and row["sign_only_columns"] == columns
               and row["base_compressed_blake3"] == manifest["compressed_blake3"]
               and row["unordered_multisets"] == str(multisets)
               and row["single_query_hit_ceiling_numerator"] == str(min(multisets, ORDER - 1))
               and row["single_query_hit_ceiling_denominator"] == str(ORDER - 1),
               "point count, multiset count, or exact ceiling differs")
        for goal, field in ((1, "queries_before_half_hit_ceiling"),
                            (columns, "queries_before_half_column_row_ceiling")):
            count = row[field]
            insist(type(count) is int and count >= 1
                   and (count - 1) * 2 * multisets < goal * (ORDER - 1)
                   and count * 2 * multisets >= goal * (ORDER - 1),
                   f"dimension {dimension} arity {arity} has a non-tight floor")
        with localcontext() as context:
            context.prec = 50
            exact = Decimal(multisets) / Decimal(ORDER - 1)
            display = Decimal(row["uniform_nonidentity_mean_display"])
            insist(abs(display - exact) <= exact * Decimal("0.0000001"),
                   "mean display disagrees with exact ratio")
        if arity in (5, 6):
            insist(row["nominal_source_coordinate_bits"] == arity * dimension
                   and row["nominal_all_affine_chain_variables"]
                   == arity * dimension + 83 * (arity - 2)
                   and row["nominal_all_affine_chain_equations"] == 83 * (arity - 1),
                   "nominal source chain size differs")
        else:
            insist(all(row[key] is None for key in (
                "nominal_source_coordinate_bits", "nominal_all_affine_chain_variables",
                "nominal_all_affine_chain_equations")),
                "unadmitted small-arity model size was inferred")
    insist(seen == {(dimension, arity) for dimension in range(2, 19)
                    for arity in range(2, 7)}, "cell grid is incomplete")
    canonical = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode()
    insist(sha256(canonical) == receipt["full_table_sha256"], "full table hash differs")
    return {"schema": "n83.standard-subspace-capacity-rank-independent-audit/v1",
            "status": "PASS_ALL_85_CELLS_AND_TIGHT_INTEGER_FLOORS",
            "cells_checked": len(rows),
            "screen_sha256": sha256((HERE / "screen.json").read_bytes()),
            "full_table_sha256": sha256(canonical),
            "selected_best_total_runtime": None}


def main() -> None:
    if len(sys.argv) != 2 or sys.argv[1] not in ("--write", "--check"):
        raise SystemExit("usage: audit.py (--write|--check)")
    result = audit(json.loads((HERE / "screen.json").read_text()))
    content = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode()
    path = HERE / "audit.json"
    if sys.argv[1] == "--write":
        with path.open("xb") as stream:
            stream.write(content)
    elif path.read_bytes() != content:
        raise ValueError("independent audit receipt differs on replay")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

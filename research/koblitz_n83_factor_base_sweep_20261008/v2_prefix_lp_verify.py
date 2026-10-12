#!/usr/bin/env python3
"""Check the exact retained K=1,182 prefix of the K=2,048 N83 object.

This is an object-identity and residual-pool check, not a solver or runtime
measurement. Both full gzip streams are read through EOF, including their CRCs.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from math import comb
from pathlib import Path


STUDY = "koblitz_n83_factor_base_sweep_20261008"
PREFIX = f"s3://crypto-autoresearcher/factor-bases/icv1/etc/{STUDY}/v2-size-frontier/a0/"
ORDER = 2417851639230796216685689
POINTS_PER_ORBIT = 166
FROZEN = {
    1182: {
        "compressed_blake3": "76e1cbe4ea5a9d26e5b24a2d27f9401271dc471846529190a0416d8e09837a92",
        "compressed_sha256": "eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380",
        "bytes": 5956114,
        "point_set_blake3": "429cde4bc514bf1026e4993e269192f858cb7410ade777178a1c0db7dbf0001a",
    },
    2048: {
        "compressed_blake3": "82990cda5d57b1c4c4b635673a881ad233c24f707e7ac644a11478aeddddf851",
        "compressed_sha256": "aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a",
        "bytes": 10313647,
        "point_set_blake3": "d528f5828a05a00add4252212bc909ece481f33093e4074ff7dd8c6a51d3c97e",
    },
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def manifest_row(path: Path, columns: int, object_path: Path) -> dict:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    rows = manifest.get("bases")
    if (manifest.get("study") != STUDY
            or manifest.get("schema") != "n83.factor-base-panel/v2-size-frontier"
            or not isinstance(rows, list) or len(rows) != 1):
        raise ValueError("unexpected retained panel manifest")
    row = rows[0]
    frozen = FROZEN[columns]
    expected_object = f"objects/{frozen['compressed_blake3']}.jsonl.gz"
    expected = {
        "a": 0,
        "columns": columns,
        "points": columns * POINTS_PER_ORBIT,
        "policy": "public_x_hash",
        "seed": 2026100801,
        "bytes": frozen["bytes"],
        "compressed_blake3": frozen["compressed_blake3"],
        "point_set_blake3": frozen["point_set_blake3"],
        "object": expected_object,
        "s3_uri": PREFIX + expected_object,
    }
    if any(type(row.get(key)) is not type(value) or row[key] != value
           for key, value in expected.items()):
        raise ValueError("retained manifest row disagrees with the frozen object")
    if object_path.stat().st_size != frozen["bytes"]:
        raise ValueError("compressed object length differs from manifest")
    if sha256_file(object_path) != frozen["compressed_sha256"]:
        raise ValueError("compressed object SHA-256 differs from frozen download")
    return row


def check_prefix(
    small_object: Path,
    large_object: Path,
    small_columns: int,
    large_columns: int,
    points_per_orbit: int,
) -> dict:
    if small_columns <= 0 or large_columns <= small_columns or points_per_orbit <= 0:
        raise ValueError("invalid orbit geometry")
    with gzip.open(small_object, "rb") as small, gzip.open(large_object, "rb") as large:
        sh = json.loads(small.readline())
        lh = json.loads(large.readline())
        common = (
            "study", "schema", "curve_a", "curve_slug", "n", "policy", "seed",
            "closure", "modulus_low_terms", "subgroup_order", "cofactor",
            "frobenius_lambda", "generator",
        )
        if any(sh.get(key) != lh.get(key) for key in common):
            raise ValueError("object headers disagree on curve or construction")
        if (sh.get("orbit_columns") != small_columns
                or lh.get("orbit_columns") != large_columns
                or sh.get("point_count") != small_columns * points_per_orbit
                or lh.get("point_count") != large_columns * points_per_orbit):
            raise ValueError("object header orbit geometry is inconsistent")
        sr, lr = sh.get("representatives"), lh.get("representatives")
        if (not isinstance(sr, list) or not isinstance(lr, list)
                or len(sr) != small_columns or len(lr) != large_columns
                or sr != lr[:small_columns]):
            raise ValueError("smaller representative list is not an exact prefix")
        shared_digest = hashlib.sha256()
        residual_digest = hashlib.sha256()
        for index in range(small_columns * points_per_orbit):
            row = small.readline()
            if not row or row != large.readline():
                raise ValueError(f"point-row prefix differs at index {index}")
            if json.loads(row).get("column") != index // points_per_orbit:
                raise ValueError(f"small point-row label differs at index {index}")
            shared_digest.update(row)
        if small.readline():
            raise ValueError("smaller object has trailing point rows")
        residual_points = (large_columns - small_columns) * points_per_orbit
        for index in range(residual_points):
            row = large.readline()
            if not row:
                raise ValueError("larger object ends before residual pool")
            if json.loads(row).get("column") != small_columns + index // points_per_orbit:
                raise ValueError(f"residual point-row label differs at index {index}")
            residual_digest.update(row)
        if large.readline():
            raise ValueError("larger object has trailing point rows")
    return {
        "small_columns": small_columns,
        "large_columns": large_columns,
        "small_point_rows": small_columns * points_per_orbit,
        "residual_columns": large_columns - small_columns,
        "residual_point_rows": residual_points,
        "representative_prefix": True,
        "point_row_prefix": True,
        "shared_rows_sha256": shared_digest.hexdigest(),
        "residual_rows_sha256": residual_digest.hexdigest(),
    }


def verify(
    small_manifest: Path,
    small_object: Path,
    large_manifest: Path,
    large_object: Path,
) -> dict:
    small_row = manifest_row(small_manifest, 1182, small_object)
    large_row = manifest_row(large_manifest, 2048, large_object)
    geometry = check_prefix(small_object, large_object, 1182, 2048, POINTS_PER_ORBIT)
    small_points = geometry["small_point_rows"]
    residual_points = geometry["residual_point_rows"]
    counts = {}
    for m in (5, 6):
        counts[str(m)] = {
            str(j): str(
                comb(small_points + m - j - 1, m - j)
                * comb(residual_points + j - 1, j)
            )
            for j in (0, 1, 2)
        }
    return {
        "schema": "n83.v2-prefix-large-prime-pool/v1",
        "study": STUDY,
        "status": "PASS_exact_prefix_only",
        "small_s3_uri": small_row["s3_uri"],
        "large_s3_uri": large_row["s3_uri"],
        "small_compressed_blake3": small_row["compressed_blake3"],
        "large_compressed_blake3": large_row["compressed_blake3"],
        "small_compressed_sha256": FROZEN[1182]["compressed_sha256"],
        "large_compressed_sha256": FROZEN[2048]["compressed_sha256"],
        **geometry,
        "nominal_unordered_multiset_counts_by_m_and_residual_count": counts,
        "uniform_nonidentity_group_denominator": str(ORDER - 1),
        "relation_yield": None,
        "rank_gain": None,
        "total_index_calculus_runtime_ms": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("small_manifest", type=Path)
    parser.add_argument("small_object", type=Path)
    parser.add_argument("large_manifest", type=Path)
    parser.add_argument("large_object", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = verify(
        args.small_manifest, args.small_object, args.large_manifest, args.large_object
    )
    with args.output.open("x", encoding="utf-8") as target:
        json.dump(result, target, sort_keys=True, indent=2)
        target.write("\n")
    print(json.dumps({"status": result["status"], "residual_columns": result["residual_columns"]}))


if __name__ == "__main__":
    main()

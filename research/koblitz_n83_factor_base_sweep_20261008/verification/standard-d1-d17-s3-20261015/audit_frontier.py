#!/usr/bin/env python3
"""Independently check N83 standard-subspace archive bytes, labels, and nesting."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import blake3


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
STUDY = HERE.parents[1]
sys.path.insert(0, str(STUDY))
SPEC = importlib.util.spec_from_file_location("n83_frontier_group", STUDY / "chain_search_supervisor.py")
assert SPEC is not None and SPEC.loader is not None
GROUP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GROUP)
ORDER = 2417851639230796216685689
LIMIT = 1 << 83
CURVE = "icv1-f2m83-tm6151469093347-debefd74"
SOURCE = "examples/koblitz_n83_standard_subspace_export.rs"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def committed_bytes(commit: str, path: str) -> bytes:
    return subprocess.run(
        ["git", "show", f"{commit}:{path}"], cwd=ROOT, check=True, capture_output=True
    ).stdout


def canonical(value: object) -> bool:
    return (isinstance(value, str) and value.isascii() and value.isdecimal()
            and str(int(value)) == value and int(value) < LIMIT)


def audit_one(directory: Path, dimension: int, source_commit: str) -> tuple[dict, set[tuple[int, int]]]:
    manifest_bytes = (directory / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    require(manifest.get("schema") == "n83.standard-subspace-factor-base-manifest/v1"
            and manifest.get("curve_slug") == CURVE
            and manifest.get("source_dimension") == dimension
            and manifest.get("source_commit") == source_commit
            and manifest.get("selected_best_total_runtime") is None,
            f"dimension {dimension}: wrong manifest identity")
    require(blake3.blake3(committed_bytes(source_commit, SOURCE)).hexdigest()
            == manifest["source_blake3"], f"dimension {dimension}: committed source mismatch")
    tree = subprocess.run(
        ["git", "rev-parse", f"{source_commit}:src"], cwd=ROOT,
        check=True, capture_output=True, text=True
    ).stdout.strip()
    require(tree == manifest["src_git_tree"], f"dimension {dimension}: src tree mismatch")
    for name, field in (("Cargo.toml", "cargo_manifest_blake3"),
                        ("Cargo.lock", "cargo_lock_blake3")):
        require(blake3.blake3(committed_bytes(source_commit, name)).hexdigest() == manifest[field],
                f"dimension {dimension}: {name} mismatch")
    object_hash = manifest["compressed_blake3"]
    key = f"objects/{object_hash}.jsonl.gz"
    prefix = "s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/standard-subspace/a0/"
    require(manifest["object"] == key and manifest["s3_uri"] == prefix + key,
            f"dimension {dimension}: content-addressed key mismatch")
    compressed = (directory / key).read_bytes()
    require(len(compressed) == manifest["compressed_bytes"]
            and blake3.blake3(compressed).hexdigest() == object_hash,
            f"dimension {dimension}: compressed bytes mismatch")
    plain = gzip.decompress(compressed)
    require(len(plain) == manifest["plain_bytes"] and plain.endswith(b"\n")
            and blake3.blake3(plain).hexdigest() == manifest["plain_blake3"],
            f"dimension {dimension}: plain bytes mismatch")
    lines = plain.splitlines()
    count = manifest["projected_point_count"]
    columns = manifest["signed_columns"]
    require(type(count) is int and type(columns) is int and count == 2 * columns
            and len(lines) == count + 1, f"dimension {dimension}: row count mismatch")
    header = json.loads(lines[0])
    require(header.get("schema") == "n83.standard-subspace-factor-base/v1"
            and header.get("curve_slug") == CURVE
            and header.get("source_dimension") == dimension
            and header.get("source_point_count") == manifest["source_point_count"]
            and header.get("projected_point_count") == count
            and header.get("signed_columns") == columns
            and header.get("subgroup_order") == str(ORDER)
            and header.get("cofactor") == "4"
            and header.get("modulus_low_terms") == [0, 1, 2, 45]
            and header.get("closure") == "sign_only_after_public_cofactor_projection"
            and header.get("subspace_membership_after_projection") is False
            and header.get("source_blake3") == manifest["source_blake3"]
            and header.get("selected_best_total_runtime") is None,
            f"dimension {dimension}: header mismatch")
    require(header.get("source_x_rule") == "polynomial_basis_low_bits_0_through_dimension_minus_1",
            f"dimension {dimension}: x domain mismatch")
    pairs: list[list[tuple[int, int, bool]]] = [[] for _ in range(columns)]
    points: set[tuple[int, int]] = set()
    digest = blake3.blake3()
    previous = None
    on_curve = subgroup = 0
    for index, line in enumerate(lines[1:]):
        row = json.loads(line)
        require(set(row) == {"index", "point", "signed_column", "negated"}
                and type(row["index"]) is int and row["index"] == index
                and set(row["point"]) == {"x", "y"}
                and canonical(row["point"]["x"]) and canonical(row["point"]["y"]),
                f"dimension {dimension}: malformed row {index}")
        x, y = int(row["point"]["x"]), int(row["point"]["y"])
        require(previous is None or previous < (x, y), f"dimension {dimension}: unsorted row {index}")
        previous = (x, y)
        points.add((x, y))
        column, negated = row["signed_column"], row["negated"]
        require(type(column) is int and 0 <= column < columns and type(negated) is bool,
                f"dimension {dimension}: bad label {index}")
        pairs[column].append((x, y, negated))
        digest.update(f"{x}:{y}\n".encode())
        if index % 256 == 0:
            require(GROUP.field_mul(y, y) ^ GROUP.field_mul(x, y)
                    == GROUP.field_mul(GROUP.field_mul(x, x), x) ^ 1,
                    f"dimension {dimension}: off-curve sample {index}")
            on_curve += 1
        if index % 4096 == 0:
            require(GROUP.point_mul((x, y), ORDER) is None,
                    f"dimension {dimension}: outside subgroup sample {index}")
            subgroup += 1
    require(len(points) == count and digest.hexdigest() == manifest["projected_point_set_blake3"]
            and header["projected_point_set_blake3"] == manifest["projected_point_set_blake3"],
            f"dimension {dimension}: projected point set mismatch")
    for column, pair in enumerate(pairs):
        require(len(pair) == 2 and {pair[0][2], pair[1][2]} == {False, True},
                f"dimension {dimension}: incomplete sign pair {column}")
        positive = next((x, y) for x, y, negated in pair if not negated)
        negative = next((x, y) for x, y, negated in pair if negated)
        require(positive[0] == negative[0] and positive[1] ^ positive[0] == negative[1],
                f"dimension {dimension}: damaged sign pair {column}")
    return ({"dimension": dimension, "status": "PASS_BYTES_LABELS_SAMPLED_GROUP",
             "source_point_count": manifest["source_point_count"], "projected_point_count": count,
             "signed_columns": columns, "compressed_bytes": len(compressed),
             "compressed_blake3": object_hash,
             "compressed_sha256": hashlib.sha256(compressed).hexdigest(),
             "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
             "on_curve_samples": on_curve, "subgroup_samples": subgroup}, points)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive_root", type=Path)
    parser.add_argument("--d18-dir", type=Path)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    construction = json.loads((HERE / "construction.json").read_text())
    source_commit = construction["source_commit"]
    results = []
    previous: set[tuple[int, int]] = set()
    for row in construction["rows"]:
        dimension = row["dimension"]
        if dimension == 1:
            require(row["status"] == "CONSTRUCT_FAILURE"
                    and "cofactor projection contains no nonidentity point" in row["stderr_text"],
                    "dimension 1 disposition changed")
            continue
        require(row["status"] == "PASS_NATIVE_CONSTRUCT_AND_EXACT_REPLAY",
                f"dimension {dimension} native construction did not pass")
        result, points = audit_one(args.archive_root / f"d{dimension}", dimension, source_commit)
        require(previous <= points, f"dimension {dimension}: previous projected set is not nested")
        result["new_projected_points_vs_previous"] = len(points - previous)
        previous = points
        results.append(result)
    d18_nested = None
    if args.d18_dir is not None:
        manifest = json.loads((args.d18_dir / "manifest.json").read_text())
        compressed = (args.d18_dir / manifest["object"]).read_bytes()
        require(blake3.blake3(compressed).hexdigest() == manifest["compressed_blake3"],
                "dimension 18 comparator bytes changed")
        lines = gzip.decompress(compressed).splitlines()[1:]
        points = {(int(row["point"]["x"]), int(row["point"]["y"]))
                  for line in lines for row in [json.loads(line)]}
        require(len(points) == manifest["projected_point_count"] and previous <= points,
                "dimension 17 projected points are not nested in dimension 18")
        d18_nested = {"status": "PASS_D17_SUBSET_D18", "d18_compressed_blake3": manifest["compressed_blake3"],
                      "d18_projected_point_count": len(points), "new_projected_points": len(points - previous)}
    receipt = {"schema": "n83.standard-subspace-frontier-independent-audit/v1",
               "status": "PASS_ALL_VALID_DIMENSIONS", "source_commit": source_commit,
               "dimension_1_status": "NO_NONIDENTITY_COFACTOR_PROJECTION",
               "rows": results, "d18_nesting": d18_nested,
               "selected_best_total_runtime": None}
    if args.write or args.receipt is not None:
        path = args.receipt if args.receipt is not None else HERE / "audit.json"
        with path.open("x") as stream:
            json.dump(receipt, stream, indent=2, sort_keys=True)
            stream.write("\n")
    print(json.dumps({"status": receipt["status"], "valid_dimensions": len(results),
                      "d18_nesting": d18_nested}, sort_keys=True))


if __name__ == "__main__":
    main()

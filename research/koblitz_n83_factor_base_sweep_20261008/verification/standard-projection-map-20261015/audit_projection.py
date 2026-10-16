#!/usr/bin/env python3
"""Independent N83 source/projected mapping audit against an archived base."""

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
STUDY = HERE.parents[1]
ROOT = HERE.parents[3]
sys.path.insert(0, str(STUDY))
SPEC = importlib.util.spec_from_file_location("n83_map_group", STUDY / "chain_search_supervisor.py")
assert SPEC is not None and SPEC.loader is not None
GROUP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GROUP)
ORDER = 2417851639230796216685689
SOURCE = "examples/koblitz_n83_standard_projection_map.rs"


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def point(value: object) -> tuple[int, int]:
    require(isinstance(value, dict) and set(value) == {"x", "y"}, "bad point object")
    out = []
    for name in ("x", "y"):
        text = value[name]
        require(isinstance(text, str) and text.isascii() and text.isdecimal()
                and str(int(text)) == text and int(text) < (1 << 83), "bad point coordinate")
        out.append(int(text))
    return (out[0], out[1])


def on_curve(p: tuple[int, int]) -> bool:
    x, y = p
    return (GROUP.field_mul(y, y) ^ GROUP.field_mul(x, y)
            == GROUP.field_mul(GROUP.field_mul(x, x), x) ^ 1)


def frobenius(point_: tuple[int, int], phase: int, negated: bool) -> tuple[int, int]:
    x, y = point_
    for _ in range(phase):
        x, y = GROUP.field_mul(x, x), GROUP.field_mul(y, y)
    return (x, y ^ x if negated else y)


def point_digest(points: list[tuple[int, int]]) -> str:
    digest = blake3.blake3()
    for x, y in sorted(points):
        digest.update(f"{x}:{y}\n".encode())
    return digest.hexdigest()


def audit(map_file: Path, base_directory: Path, commit: str) -> dict:
    data = map_file.read_bytes()
    lines = data.splitlines()
    header = json.loads(lines[0])
    dimension = header.get("source_dimension")
    require(header.get("schema") == "n83.standard-subspace-projection-map/v1"
            and header.get("curve_slug") == "icv1-f2m83-tm6151469093347-debefd74"
            and type(dimension) is int and 2 <= dimension <= 18
            and header.get("source_x_rule") == "polynomial_basis_low_bits_0_through_dimension_minus_1"
            and header.get("subgroup_order") == str(ORDER)
            and header.get("cofactor") == "4"
            and header.get("frobenius_lambda") == str(GROUP.PRIMARY_FROBENIUS_LAMBDA)
            and header.get("selected_best_total_runtime") is None,
            "projection-map header identity changed")
    committed = subprocess.run(["git", "show", f"{commit}:{SOURCE}"], cwd=ROOT,
                               check=True, capture_output=True).stdout
    require(blake3.blake3(committed).hexdigest() == header["source_blake3"],
            "projection-map committed source differs")
    source_count = header["source_point_count"]
    projected_count = header["projected_point_count"]
    columns_count = header["projected_signed_frobenius_columns"]
    require(all(type(v) is int and v > 0 for v in (source_count, projected_count, columns_count))
            and len(lines) == 1 + columns_count + projected_count + source_count,
            "projection-map row count mismatch")
    rows = iter(json.loads(line) for line in lines[1:])
    columns = []
    for index in range(columns_count):
        row = next(rows)
        require(row.get("type") == "column" and row.get("index") == index,
                f"bad column row {index}")
        columns.append(point(row["representative"]))
    require(columns == sorted(set(columns)), "canonical columns unsorted or duplicated")
    base_manifest = json.loads((base_directory / "manifest.json").read_text())
    require(base_manifest["source_dimension"] == dimension
            and base_manifest["projected_point_count"] == projected_count
            and base_manifest["projected_point_set_blake3"] == header["projected_point_set_blake3"]
            and base_manifest["signed_columns"] == header["sign_only_columns"],
            "archived base does not match map header")
    base_compressed = (base_directory / base_manifest["object"]).read_bytes()
    require(blake3.blake3(base_compressed).hexdigest() == base_manifest["compressed_blake3"],
            "base archive compressed bytes changed")
    base_rows = [json.loads(line) for line in gzip.decompress(base_compressed).splitlines()[1:]]
    require(len(base_rows) == projected_count, "base archive row count differs")
    projected = []
    sampled_frobenius = sampled_subgroup = 0
    for index in range(projected_count):
        row = next(rows)
        base = base_rows[index]
        require(row.get("type") == "projected" and row.get("index") == index
                and row.get("point") == base["point"]
                and row.get("sign_only_column") == base["signed_column"]
                and row.get("sign_only_negated") == base["negated"],
                f"projected row {index} differs from archived base")
        col, phase, neg = row["column"], row["phase"], row["negated"]
        require(type(col) is int and 0 <= col < columns_count
                and type(phase) is int and 0 <= phase < 83 and type(neg) is bool,
                f"bad projected label {index}")
        p = point(row["point"])
        projected.append((p, col, phase, neg))
        if index % 64 == 0:
            require(frobenius(columns[col], phase, neg) == p,
                    f"projected Frobenius equation fails at {index}")
            sampled_frobenius += 1
        if index % 1024 == 0:
            require(on_curve(p) and GROUP.point_mul(p, ORDER) is None,
                    f"projected subgroup sample fails at {index}")
            sampled_subgroup += 1
    require(point_digest([entry[0] for entry in projected]) == header["projected_point_set_blake3"],
            "projected point-set digest differs")
    observed_columns = {col for _, col, _, _ in projected}
    require(len(observed_columns) == columns_count, "some projected columns have no point")
    base_pairs: dict[int, list[tuple[tuple[int, int], int, int, bool]]] = {}
    for index, base in enumerate(base_rows):
        base_pairs.setdefault(base["signed_column"], []).append(projected[index])
    require(len(base_pairs) == header["sign_only_columns"], "sign-only column coverage differs")
    for base_col, members in base_pairs.items():
        require(len(members) == 2 and members[0][1] == members[1][1]
                and members[0][2] == members[1][2]
                and members[0][3] != members[1][3],
                f"sign-only pair {base_col} has inconsistent Frobenius labels")
    source_points = []
    seen_projected = set()
    excluded = sampled_cofactor = 0
    for index in range(source_count):
        row = next(rows)
        require(row.get("type") == "source" and row.get("index") == index,
                f"bad source row {index}")
        p = point(row["point"])
        source_points.append(p)
        require(row.get("source_x_code") == str(p[0]) and p[0] < (1 << dimension),
                f"source x code outside subspace at {index}")
        projected_index = row.get("projected_index")
        if projected_index is None:
            require(row.get("column") is None and row.get("phase") is None
                    and row.get("negated") is None, f"excluded source {index} has label")
            excluded += 1
        else:
            require(type(projected_index) is int and 0 <= projected_index < projected_count
                    and projected_index not in seen_projected,
                    f"source projection {index} missing or duplicated")
            seen_projected.add(projected_index)
            _, col, phase, neg = projected[projected_index]
            require((row.get("column"), row.get("phase"), row.get("negated"))
                    == (col, phase, neg), f"source row {index} label differs")
        if index % 512 == 0:
            projected_point = GROUP.point_mul(p, 4)
            require(on_curve(p) and (projected_point is None if projected_index is None
                                     else projected_point == projected[projected_index][0]),
                    f"source cofactor projection sample fails at {index}")
            sampled_cofactor += 1
    require(len(seen_projected) == projected_count
            and excluded == header["excluded_source_points"]
            and point_digest(source_points) == header["source_point_set_blake3"],
            "source projection coverage or digest differs")
    return {"schema": "n83.standard-projection-map-independent-audit/v1",
            "status": "PASS_ROWS_LABELS_SAMPLED_GROUP",
            "source_commit": commit, "source_dimension": dimension,
            "map_plain_bytes": len(data), "map_plain_sha256": hashlib.sha256(data).hexdigest(),
            "map_plain_blake3": blake3.blake3(data).hexdigest(),
            "base_compressed_blake3": base_manifest["compressed_blake3"],
            "source_point_count": source_count, "projected_point_count": projected_count,
            "sign_only_columns": header["sign_only_columns"],
            "projected_signed_frobenius_columns": columns_count,
            "excluded_source_points": excluded,
            "sampled_frobenius_rows": sampled_frobenius,
            "sampled_subgroup_rows": sampled_subgroup,
            "sampled_source_cofactor_rows": sampled_cofactor,
            "selected_best_total_runtime": None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("map_file", type=Path)
    parser.add_argument("base_directory", type=Path)
    parser.add_argument("source_commit")
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    result = audit(args.map_file, args.base_directory, args.source_commit)
    if args.receipt:
        with args.receipt.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

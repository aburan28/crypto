#!/usr/bin/env python3
"""Independent byte, sign-label, and sampled group replay for the N83 d18 base."""

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
SPEC = importlib.util.spec_from_file_location("n83_standard_group", STUDY / "chain_search_supervisor.py")
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


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_decimal(value: object) -> bool:
    return isinstance(value, str) and value.isdecimal() and str(int(value)) == value and int(value) < LIMIT


def check_sign_pair(rows: list[tuple[int, int, bool]], column: int) -> None:
    require(len(rows) == 2, f"column {column} is not an exact sign pair")
    positive = next((row for row in rows if row[2] is False), None)
    negative = next((row for row in rows if row[2] is True), None)
    require(positive is not None and negative is not None
            and positive[0] == negative[0]
            and positive[1] ^ positive[0] == negative[1],
            f"column {column} has a damaged negation label")


def audit(directory: Path) -> dict:
    manifest_bytes = (directory / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    require(manifest.get("schema") == "n83.standard-subspace-factor-base-manifest/v1"
            and manifest.get("curve_slug") == CURVE
            and manifest.get("source_dimension") == 18
            and manifest.get("source_point_count") == 261447
            and manifest.get("projected_point_count") == 261444
            and manifest.get("signed_columns") == 130722
            and manifest.get("selected_best_total_runtime") is None,
            "manifest scope or exact prior-art counts changed")
    commit = manifest["source_commit"]
    require(isinstance(commit, str) and len(commit) == 40, "source commit malformed")
    source = subprocess.run(["git", "show", f"{commit}:{SOURCE}"], cwd=ROOT,
                            check=True, capture_output=True).stdout
    require(blake3.blake3(source).hexdigest() == manifest["source_blake3"],
            "committed exporter source does not match manifest")
    tree = subprocess.run(["git", "rev-parse", f"{commit}:src"], cwd=ROOT,
                          check=True, capture_output=True, text=True).stdout.strip()
    require(tree == manifest["src_git_tree"], "compiled library tree changed")
    for name, field in (("Cargo.toml", "cargo_manifest_blake3"),
                        ("Cargo.lock", "cargo_lock_blake3")):
        data = subprocess.run(["git", "show", f"{commit}:{name}"], cwd=ROOT,
                              check=True, capture_output=True).stdout
        require(blake3.blake3(data).hexdigest() == manifest[field],
                f"{name} differs from frozen construction")
    object_hash = manifest["compressed_blake3"]
    key = f"objects/{object_hash}.jsonl.gz"
    require(manifest["object"] == key and
            manifest["s3_uri"] ==
            "s3://crypto-autoresearcher/factor-bases/icv1/etc/"
            "koblitz_n83_factor_base_sweep_20261008/standard-subspace/a0/" + key,
            "content-addressed object identity changed")
    compressed = (directory / key).read_bytes()
    require(len(compressed) == manifest["compressed_bytes"]
            and blake3.blake3(compressed).hexdigest() == object_hash,
            "compressed bytes do not match manifest")
    plain = gzip.decompress(compressed)
    require(len(plain) == manifest["plain_bytes"]
            and blake3.blake3(plain).hexdigest() == manifest["plain_blake3"],
            "gzip or plain bytes do not match manifest")
    local_plain = directory / "factor-base.jsonl"
    if local_plain.exists():
        require(local_plain.read_bytes() == plain, "optional local plain copy differs")
    lines = plain.splitlines()
    require(len(lines) == 261445, "archive row count changed")
    header = json.loads(lines[0])
    require(header.get("schema") == "n83.standard-subspace-factor-base/v1"
            and header.get("curve_slug") == CURVE
            and header.get("source_dimension") == 18
            and header.get("source_x_rule") ==
            "polynomial_basis_low_bits_0_through_dimension_minus_1"
            and header.get("closure") == "sign_only_after_public_cofactor_projection"
            and header.get("subspace_membership_after_projection") is False
            and header.get("source_point_count") == 261447
            and header.get("projected_point_count") == 261444
            and header.get("signed_columns") == 130722
            and header.get("source_blake3") == manifest["source_blake3"]
            and header.get("projected_point_set_blake3") ==
            manifest["projected_point_set_blake3"]
            and header.get("subgroup_order") == str(ORDER)
            and header.get("cofactor") == "4"
            and header.get("modulus_low_terms") == [0, 1, 2, 45]
            and header.get("selected_best_total_runtime") is None,
            "archive header changed")
    labels: list[list[tuple[int, int, bool]]] = [[] for _ in range(130722)]
    digest = blake3.blake3()
    previous = None
    on_curve_samples = 0
    subgroup_samples = 0
    for index, line in enumerate(lines[1:]):
        row = json.loads(line)
        point = row["point"]
        require(row["index"] == index
                and canonical_decimal(point.get("x"))
                and canonical_decimal(point.get("y")),
                f"malformed point row {index}")
        x, y = int(point["x"]), int(point["y"])
        require(previous is None or previous < (x, y),
                f"point order or uniqueness failed at {index}")
        previous = (x, y)
        column = row["signed_column"]
        negated = row["negated"]
        require(type(column) is int and 0 <= column < len(labels)
                and type(negated) is bool, f"bad sign label at {index}")
        labels[column].append((x, y, negated))
        digest.update(f"{x}:{y}\n".encode())
        if index % 256 == 0:
            require(GROUP.field_mul(y, y) ^ GROUP.field_mul(x, y)
                    == GROUP.field_mul(GROUP.field_mul(x, x), x) ^ 1,
                    f"sampled point {index} is off curve")
            on_curve_samples += 1
        if index % 4096 == 0:
            require(GROUP.point_mul((x, y), ORDER) is None,
                    f"sampled point {index} is outside prime subgroup")
            subgroup_samples += 1
    require(digest.hexdigest() == manifest["projected_point_set_blake3"],
            "projected point-set digest changed")
    for column, pair in enumerate(labels):
        check_sign_pair(pair, column)
    return {
        "schema": "n83.standard-subspace-d18-independent-audit/v1",
        "status": "PASS_BYTES_ALL_SIGN_LABELS_SAMPLED_GROUP",
        "source_commit": commit,
        "source_blake3": manifest["source_blake3"],
        "object": key,
        "compressed_blake3": object_hash,
        "compressed_sha256": sha256(compressed),
        "plain_sha256": sha256(plain),
        "manifest_sha256": sha256(manifest_bytes),
        "point_rows": 261444,
        "sign_columns_checked": len(labels),
        "on_curve_samples": on_curve_samples,
        "subgroup_samples": subgroup_samples,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive_dir", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = audit(args.archive_dir)
    if args.write:
        with (HERE / "audit.json").open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

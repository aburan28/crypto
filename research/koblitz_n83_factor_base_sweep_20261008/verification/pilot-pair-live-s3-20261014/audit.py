#!/usr/bin/env python3
"""Replay three retained N83 S3 archives and their selected-pair overlap."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import blake3


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
MANIFEST = STUDY / "pilot-01" / "manifest.json"
SCREEN = HERE.parent / "pilot-pair-support-screen-20261014" / "screen.json"
HEADS = HERE / "heads.json"
RECEIPT = HERE / "audit.json"
ROLES = {
    "base": "3888e2597fd3519807c9d722f49f83329c93ce5cb0929501cf7cb7cbc3cab1a0",
    "gray_pool": "966a0493c96674f0af6d184615ed05309d08fba5c7342957b4186cd50366dc47",
    "hash_pool": "b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e",
}
FIELD_LIMIT = 1 << 83


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_archive(path: Path, manifest_row: dict) -> tuple[dict, dict]:
    compressed = path.read_bytes()
    digest = blake3.blake3(compressed).hexdigest()
    require(digest == manifest_row["compressed_blake3"], "compressed BLAKE3 mismatch")
    require(len(compressed) == manifest_row["bytes"], "compressed length mismatch")

    columns: dict[int, set[tuple[int, int]]] = {}
    labels: dict[int, set[tuple[int, bool]]] = {}
    all_points: set[tuple[int, int]] = set()
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        header = json.loads(next(stream))
        require(header["schema"] == "n83.factor-base/v1", "archive schema mismatch")
        require(header["closure"] == "signed_frobenius", "archive closure mismatch")
        require(header["curve_slug"] == "icv1-f2m83-tm6151469093347-debefd74", "curve mismatch")
        require(header["subgroup_order"] == "2417851639230796216685689", "order mismatch")
        require(header["n"] == 83, "field degree mismatch")
        require(header["modulus_low_terms"] == [0, 1, 2, 45], "field modulus mismatch")
        require(header["orbit_columns"] == manifest_row["columns"] == 600, "column count mismatch")
        require(header["point_count"] == manifest_row["points"] == 99_600, "point count mismatch")
        require(header["policy"] == manifest_row["policy"], "policy mismatch")
        require(header["seed"] == manifest_row["seed"], "seed mismatch")
        require(len(header["representatives"]) == 600, "representative count mismatch")
        for line in stream:
            row = json.loads(line)
            column = row["column"]
            phase = row["phase"]
            negative = row["negative"]
            require(type(column) is int and 0 <= column < 600, "invalid column")
            require(type(phase) is int and 0 <= phase < 83, "invalid phase")
            require(type(negative) is bool, "invalid sign")
            require(0 < int(row["coefficient"]) < int(header["subgroup_order"]), "invalid coefficient")
            point = tuple(map(int, row["point"]))
            require(len(point) == 2 and all(0 <= value < FIELD_LIMIT for value in point), "invalid point")
            require(point not in all_points, "duplicate point across archived columns")
            all_points.add(point)
            columns.setdefault(column, set()).add(point)
            label = phase, negative
            bucket = labels.setdefault(column, set())
            require(label not in bucket, "duplicate signed Frobenius label")
            bucket.add(label)

    require(len(all_points) == 99_600 and len(columns) == 600, "archive size mismatch")
    require(set(columns) == set(range(600)), "missing archive column")
    require(all(len(points) == 166 for points in columns.values()), "orbit size mismatch")
    require(all(len(bucket) == 166 for bucket in labels.values()), "orbit label mismatch")
    for column, representative in enumerate(header["representatives"]):
        require(tuple(map(int, representative)) in columns[column], "representative outside orbit")

    canonical = {min(points): frozenset(points) for points in columns.values()}
    require(len(canonical) == 600, "duplicate signed Frobenius orbit")
    canonical_bytes = "".join(
        f"{x},{y}\n" for x, y in sorted(canonical)
    ).encode("ascii")
    details = {
        "compressed_blake3": digest,
        "compressed_bytes": len(compressed),
        "compressed_sha256": sha256(compressed),
        "canonical_orbit_key_sha256": sha256(canonical_bytes),
        "columns": len(columns),
        "points": len(all_points),
        "policy": manifest_row["policy"],
        "seed": manifest_row["seed"],
        "s3_uri": manifest_row["s3_uri"],
    }
    return details, {"orbits": canonical, "representatives": {
        tuple(map(int, point)) for point in header["representatives"]
    }}


def audit(archive_dir: Path) -> dict:
    manifest_bytes = MANIFEST.read_bytes()
    screen_bytes = SCREEN.read_bytes()
    head_bytes = HEADS.read_bytes()
    manifest = json.loads(manifest_bytes)
    screen = json.loads(screen_bytes)
    heads = json.loads(head_bytes)
    rows = {row["compressed_blake3"]: row for row in manifest["bases"]}
    objects = {}
    domains = {}
    for role, digest in ROLES.items():
        row = rows[digest]
        path = archive_dir / f"{digest}.jsonl.gz"
        objects[role], domains[role] = load_archive(path, row)
        require(heads["objects"][role]["compressed_blake3"] == digest, "S3 HEAD identity mismatch")
        require(heads["objects"][role]["content_length"] == objects[role]["compressed_bytes"],
                "S3 HEAD length mismatch")

    pairs = {}
    base_orbits = domains["base"]["orbits"]
    for role in ("gray_pool", "hash_pool"):
        pool_orbits = domains[role]["orbits"]
        common = base_orbits.keys() & pool_orbits.keys()
        for key in common:
            require(base_orbits[key] == pool_orbits[key], "canonical-key collision")
        point_overlap = 166 * len(common)
        require(
            point_overlap == len(set().union(*base_orbits.values()) & set().union(*pool_orbits.values())),
            "partial orbit overlap",
        )
        screen_row = next(
            row for row in screen["ordered_pairs"]
            if row["base_object_blake3"] == ROLES["base"]
            and row["pool_object_blake3"] == ROLES[role]
        )
        require(screen_row["worker_admission"] == "ADMISSIBLE_MIN_TWO_RESIDUAL_ORBITS", "worker admission mismatch")
        require(screen_row["shared_signed_orbits"] == len(common), "archived shared count mismatch")
        require(screen_row["residual_signed_orbits"] == 600 - len(common), "archived residual count mismatch")
        pairs[role] = {
            "shared_signed_orbits": len(common),
            "residual_signed_orbits": 600 - len(common),
            "residual_points": 166 * (600 - len(common)),
            "shared_points": point_overlap,
            "literal_header_representatives_shared": len(
                domains["base"]["representatives"] & domains[role]["representatives"]
            ),
            "parent_support_row_matched": True,
        }

    return {
        "schema": "n83.live-s3-selected-pair-audit/v1",
        "status": "PASS_ARCHIVE_BYTES_AND_ORBIT_OVERLAP",
        "manifest_sha256": sha256(manifest_bytes),
        "support_screen_sha256": sha256(screen_bytes),
        "s3_head_sha256": sha256(head_bytes),
        "source_sha256": sha256(Path(__file__).read_bytes()),
        "objects": objects,
        "pairs": pairs,
        "new_n83_solver_search_executed": False,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-dir", type=Path, required=True)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    value = audit(args.archive_dir)
    encoded = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()
    if args.write:
        RECEIPT.write_bytes(encoded)
    else:
        require(RECEIPT.read_bytes() == encoded, "saved receipt differs from replay")
    print(json.dumps({"status": value["status"], "pairs": value["pairs"]}, sort_keys=True))


if __name__ == "__main__":
    main()

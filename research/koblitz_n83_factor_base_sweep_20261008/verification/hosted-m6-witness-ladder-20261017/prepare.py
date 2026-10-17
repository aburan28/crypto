#!/usr/bin/env python3
"""Rebuild two planted 4+2 N83 witnesses from the exact S3 factor bases."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import blake3

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
PAIRED_AUDIT = STUDY / "verification/k600-gray-hash-paired-protocol-20261014"
sys.path.insert(0, str(PAIRED_AUDIT))
sys.path.insert(0, str(STUDY))
import receipt_audit as prior  # noqa: E402
from chain_search_supervisor import PRIMARY_ORDER, parse_point, point_add, point_mul  # noqa: E402

SOURCE_RECEIPT = STUDY / "verification/hosted-m6-paired-run-39054543374/factor-bases-s3-check.json"
BASE = "3888e2597fd3519807c9d722f49f83329c93ce5cb0929501cf7cb7cbc3cab1a0"
POOLS = (
    ("gray_pool", "966a0493c96674f0af6d184615ed05309d08fba5c7342957b4186cd50366dc47", 503),
    ("hash_pool", "b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e", 0),
)
ORBIT = 166


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decimal_point(point: tuple[int, int]) -> list[str]:
    return [str(point[0]), str(point[1])]


def checked_partial(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
    for point in points:
        require(parse_point(decimal_point(point)) == point, "selected point is off curve")
        require(point_mul(point, PRIMARY_ORDER) is None, "selected point outside subgroup")
    partial = None
    intermediates = []
    for point in points:
        partial = point_add(partial, point)
        require(partial is not None, "witness has identity partial")
        intermediates.append(partial)
    require(point_mul(intermediates[-1], PRIMARY_ORDER) is None,
            "target outside subgroup")
    return intermediates


def build(directory: Path) -> dict:
    receipt = json.loads(SOURCE_RECEIPT.read_text())
    require(receipt["status"] == "PASS_3_FROZEN_S3_OBJECTS", "S3 receipt status differs")
    rows = {row["compressed_blake3"]: row for row in receipt["objects"]}
    require(set(rows) == {BASE, *(pool for _, pool, _ in POOLS)}, "S3 identities differ")
    archives = {}
    for digest, row in rows.items():
        path = directory / f"{digest}.jsonl.gz"
        require(path.is_file() and path.stat().st_size == row["bytes"]
                and sha256(path) == row["compressed_sha256"]
                and blake3.blake3(path.read_bytes()).hexdigest() == digest,
                f"S3 archive bytes differ: {digest}")
        archives[digest] = prior.archive_points(path)
    base = archives[BASE]
    base_set = set(base)
    base_x = {point[0] for point in base}
    base_reps = [base[i] for i in range(0, len(base), ORBIT)]
    require(len(base_reps) == 600 and len({point[0] for point in base_reps}) == 600,
            "base representative count differs")
    arms = []
    for role, pool_digest, shared_expected in POOLS:
        pool = archives[pool_digest]
        residual_reps = []
        shared = 0
        for offset in range(0, len(pool), ORBIT):
            orbit = pool[offset:offset + ORBIT]
            overlap = sum(point in base_set for point in orbit)
            require(overlap in (0, ORBIT), "partial signed-orbit overlap")
            if overlap:
                shared += 1
            else:
                require(all(point[0] not in base_x for point in orbit),
                        "residual x-domain overlap")
                residual_reps.append(orbit[0])
        require(shared == shared_expected
                and len(residual_reps) == 600 - shared_expected
                and len({point[0] for point in residual_reps}) == len(residual_reps),
                "residual representative count differs")
        selected = sorted(base_reps[:2], key=lambda point: point[0])
        selected += base_reps[2:4] + residual_reps[:2]
        partials = checked_partial(selected)
        require(selected[0][0] < selected[1][0], "first-pair order differs")
        arms.append({
            "role": role, "base_object_blake3": BASE, "pool_object_blake3": pool_digest,
            "shared_signed_orbits": shared, "residual_signed_orbits": len(residual_reps),
            "base_reps": [str(point[0]) for point in base_reps],
            "residual_reps": [str(point[0]) for point in residual_reps],
            "selected_points": [decimal_point(point) for point in selected],
            "ordered_xs": [str(point[0]) for point in selected],
            "partial_xs": [str(point[0]) for point in partials[1:-1]],
            "target": decimal_point(partials[-1]), "target_x": str(partials[-1][0]),
            "identity_mask": 0,
        })
    return {
        "schema": "n83.m6-planted-witness-ladder/v1",
        "status": "PASS_S3_BOUND_GROUP_WITNESSES",
        "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
        "s3_receipt_sha256": sha256(SOURCE_RECEIPT),
        "archive_sha256": {digest: row["compressed_sha256"] for digest, row in rows.items()},
        "summands": 6, "base_summands": 4, "residual_summands": 2,
        "arms": arms,
        "natural_relations": None, "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bases", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    saved = HERE / "fixture.json"
    result = build(args.bases)
    encoded = (json.dumps(result, sort_keys=True, indent=2) + "\n").encode()
    if args.check:
        require(saved.read_bytes() == encoded, "saved six-summand fixture differs")
    else:
        with saved.open("xb") as stream:
            stream.write(encoded)
    print(json.dumps({"status": result["status"], "fixture_sha256": hashlib.sha256(encoded).hexdigest(),
                      "arms": [row["role"] for row in result["arms"]]}, sort_keys=True))


if __name__ == "__main__":
    main()

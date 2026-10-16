#!/usr/bin/env python3
"""Independent minimum-point replay of every N83 cross-stream orbit intersection."""

from __future__ import annotations

import argparse
import blake3
import gzip
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
PARENTS = {
    "census": (STUDY / "verification/cross-stream-pilot-import-source-20261013/census.json",
               "d1a490e2e73552d96a5ee327e34367c5bb32fbc6f6d591eb9c0e5f029118a285"),
    "manifest": (STUDY / "pilot-01/manifest.json",
                 "5bd862987844f709e03beeab9fd265e3d2024011bcadcb82f5af9c72018df6a6"),
    "s3_audit": (STUDY / "verification/s3-stream-integrity-20261014/audit.json",
                 "589d1c2d8577d596624e75df1ea15dfa09b52d61d4b1081eb2e01093331da843"),
}


def insist(ok: bool, why: str) -> None:
    if not ok:
        raise ValueError(why)


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def minimum_point_orbits(path: Path, expected: dict) -> tuple[set[tuple[int, int]], str]:
    raw = path.read_bytes()
    insist(blake3.blake3(raw).hexdigest() == expected["compressed_blake3"],
           "local archive does not match S3 manifest")
    plain = gzip.decompress(raw)
    insist(blake3.blake3(plain).hexdigest() == expected["plain_blake3"],
           "plain archive does not match S3 manifest")
    iterator = iter(plain.splitlines())
    header = json.loads(next(iterator))
    k = expected["columns"]
    insist(header["orbit_columns"] == k and header["point_count"] == 166 * k
           and header["closure"] == "signed_frobenius" and header["curve_a"] == 0
           and header["policy"] == expected["policy"] and header["seed"] == expected["seed"]
           and len(header["representatives"]) == k, "archive identity or width changed")
    minima = set()
    all_points = set()
    for column in range(k):
        labels = set()
        points = set()
        positive_zero = None
        for _ in range(166):
            item = json.loads(next(iterator))
            label = item["phase"], item["negative"]
            point = tuple(int(v) for v in item["point"])
            insist(item["column"] == column and 0 <= label[0] < 83
                   and type(label[1]) is bool and len(point) == 2
                   and label not in labels and point not in points,
                   "incomplete or repeated signed orbit row")
            labels.add(label)
            points.add(point)
            if label == (0, False):
                positive_zero = point
        insist(len(labels) == len(points) == 166
               and positive_zero == tuple(map(int, header["representatives"][column]))
               and all_points.isdisjoint(points), "orbit phase, representative, or disjointness differs")
        all_points.update(points)
        minimum = min(points)
        insist(minimum not in minima, "duplicate orbit minimum")
        minima.add(minimum)
    insist(next(iterator, None) is None and len(all_points) == 166 * k,
           "archive has extra or missing rows")
    return minima, sha(raw)


def check_pairs(pairs: list[dict], populations: dict[tuple, set]) -> dict:
    expected_keys = {(base, pool) for base in populations for pool in populations
                     if base[:2] != pool[:2]}
    insist(len(expected_keys) == 648 and len(pairs) == 648, "ordered pair grid differs")
    seen = set()
    overlap_histogram: dict[str, int] = {}
    for row in pairs:
        base, pool = tuple(row["base"]), tuple(row["pool"])
        insist((base, pool) in expected_keys and (base, pool) not in seen,
               "unexpected or duplicate pair")
        seen.add((base, pool))
        common = len(populations[base].intersection(populations[pool]))
        insist(row["shared_orbits"] == common
               and row["base_only_orbits"] == base[2] - common
               and row["residual_orbits"] == pool[2] - common
               and row["union_orbits"] == base[2] + pool[2] - common
               and row["residual_points"] == 166 * (pool[2] - common),
               "pair overlap or residual count changed")
        label = str(common)
        overlap_histogram[label] = overlap_histogram.get(label, 0) + 1
    insist(seen == expected_keys, "ordered pair grid is incomplete")
    return overlap_histogram


def audit(objects_dir: Path, receipt: dict, receipt_sha256: str) -> dict:
    insist(receipt.get("schema") == "n83.primary-cross-stream-orbit-overlap/v1"
           and receipt.get("status") == "PASS_27_ARCHIVES_648_EXACT_PAIRS"
           and receipt.get("selected_best_total_runtime") is None
           and receipt.get("archive_count") == 27 and receipt.get("pair_count") == 648,
           "receipt claim differs")
    sources = {}
    for name, (path, digest) in PARENTS.items():
        data = path.read_bytes()
        insist(sha(data) == digest == receipt["input_sha256"][name],
               f"{name} provenance differs")
        sources[name] = json.loads(data)
    insist(sha((HERE / "overlap.py").read_bytes()) == receipt["generator_sha256"],
           "generator source changed")
    census = sources["census"]
    manifest = sources["manifest"]
    s3_audit = sources["s3_audit"]
    insist(census["pilot_manifest_sha256"] == PARENTS["manifest"][1]
           and census["ordered_cross_stream_pairs"] == 648
           and s3_audit["status"] == "PASS" and s3_audit["object_count"] == 57,
           "parent admission differs")
    manifest_by_uri = {r["s3_uri"]: r for r in manifest["bases"] if r["a"] == 0}
    s3_by_uri = {r["uri"]: r for r in s3_audit["objects"]}
    populations = {}
    rows = receipt["archives"]
    insist(len(rows) == 27 and len(census["primary_objects"]) == 27,
           "archive census count differs")
    for census_row, saved in zip(sorted(census["primary_objects"],
                                        key=lambda r: (r["policy"], r["seed"], r["columns"])), rows):
        key = (census_row["policy"], census_row["seed"], census_row["columns"])
        uri = census_row["s3_uri"]
        parent = manifest_by_uri[uri]
        observed = s3_by_uri[uri]
        insist((saved["policy"], saved["seed"], saved["columns"]) == key
               and saved["s3_uri"] == uri
               and saved["compressed_blake3"] == parent["compressed_blake3"]
               == census_row["compressed_blake3"] == observed["observed_blake3"]
               and observed["status"] == "PASS" and observed["observed_bytes"] == parent["bytes"],
               "archive S3 binding changed")
        path = objects_dir / (parent["compressed_blake3"] + ".jsonl.gz")
        minima, compressed_sha = minimum_point_orbits(path, parent)
        insist(saved["compressed_sha256"] == compressed_sha and len(minima) == key[2],
               "archive compressed bytes or orbit count changed")
        populations[key] = minima
    histogram = check_pairs(receipt["pairs"], populations)
    return {"schema": "n83.primary-cross-stream-orbit-overlap-independent-audit/v1",
            "status": "PASS_27_ARCHIVES_648_PAIRS_MINIMUM_POINT_REPLAY",
            "receipt_sha256": receipt_sha256,
            "archives_checked": 27, "pairs_checked": 648,
            "overlap_histogram": histogram,
            "selected_best_total_runtime": None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--objects-dir", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    original = (HERE / "overlap.json").read_bytes()
    result = audit(args.objects_dir, json.loads(original), sha(original))
    data = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode()
    target = HERE / "audit.json"
    if args.write:
        with target.open("xb") as stream:
            stream.write(data)
    else:
        insist(target.read_bytes() == data, "audit receipt differs on independent replay")
    print(json.dumps({"status": result["status"], "archives": 27, "pairs": 648,
                      "receipt_sha256": sha(data)}))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact signed-orbit intersections of the 27 archived primary pilot bases."""

from __future__ import annotations

import argparse
import blake3
import gzip
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
ORDER = 2417851639230796216685689
SOURCES = {
    "census": (STUDY / "verification/cross-stream-pilot-import-source-20261013/census.json",
               "d1a490e2e73552d96a5ee327e34367c5bb32fbc6f6d591eb9c0e5f029118a285"),
    "manifest": (STUDY / "pilot-01/manifest.json",
                 "5bd862987844f709e03beeab9fd265e3d2024011bcadcb82f5af9c72018df6a6"),
    "s3_audit": (STUDY / "verification/s3-stream-integrity-20261014/audit.json",
                 "589d1c2d8577d596624e75df1ea15dfa09b52d61d4b1081eb2e01093331da843"),
}


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def pinned_inputs() -> tuple[dict, dict, dict]:
    loaded = {}
    for name, (path, digest) in SOURCES.items():
        raw = path.read_bytes()
        require(sha256(raw) == digest, f"pinned {name} changed")
        loaded[name] = json.loads(raw)
    census, manifest, audit = (loaded[name] for name in SOURCES)
    require(census["status"] == "source_addressable_unmeasured"
            and census["ordered_cross_stream_pairs"] == 648
            and census["pilot_manifest_sha256"] == SOURCES["manifest"][1]
            and manifest["status"] == "completed_factor_base_panel"
            and audit["status"] == "PASS"
            and audit["object_count"] == 57,
            "parent census, manifest, or S3 audit differs")
    return census, manifest, audit


def identity(row: dict) -> tuple[str, int, int]:
    return row["policy"], row["seed"], row["columns"]


def archive_orbits(path: Path, row: dict) -> tuple[list[str], str]:
    compressed = path.read_bytes()
    require(blake3.blake3(compressed).hexdigest() == row["compressed_blake3"],
            f"compressed archive changed: {identity(row)}")
    plain = gzip.decompress(compressed)
    require(blake3.blake3(plain).hexdigest() == row["plain_blake3"],
            f"plain archive changed: {identity(row)}")
    lines = plain.splitlines()
    header = json.loads(lines[0])
    k = row["columns"]
    require(header["curve_a"] == 0 and header["closure"] == "signed_frobenius"
            and header["orbit_columns"] == k and header["point_count"] == 166 * k
            and header["policy"] == row["policy"] and header["seed"] == row["seed"]
            and len(lines) == 1 + 166 * k,
            f"archive header/count differs: {identity(row)}")
    powers = [pow(int(header["frobenius_lambda"]), phase, ORDER) for phase in range(83)]
    orbit_hashes = []
    all_points = set()
    for col in range(k):
        points = set()
        labels = set()
        for line in lines[1 + col * 166:1 + (col + 1) * 166]:
            item = json.loads(line)
            phase, negative = item["phase"], item["negative"]
            require(item["column"] == col and type(phase) is int and 0 <= phase < 83
                    and type(negative) is bool,
                    f"orbit row label differs: {identity(row)}, column {col}")
            coefficient = -powers[phase] if negative else powers[phase]
            require(int(item["coefficient"]) == coefficient % ORDER,
                    f"Frobenius coefficient differs: {identity(row)}, column {col}")
            point = tuple(map(int, item["point"]))
            require(len(point) == 2 and point not in points,
                    f"duplicate/malformed orbit point: {identity(row)}, column {col}")
            points.add(point)
            labels.add((phase, negative))
        require(len(points) == 166 and len(labels) == 166 and all_points.isdisjoint(points),
                f"incomplete or overlapping orbit: {identity(row)}, column {col}")
        all_points.update(points)
        canonical = json.dumps(sorted(points), separators=(",", ":")).encode()
        orbit_hashes.append(sha256(canonical))
    require(len(set(orbit_hashes)) == k and len(all_points) == 166 * k,
            f"archive has duplicate orbit: {identity(row)}")
    return orbit_hashes, sha256(compressed)


def make(objects_dir: Path) -> dict:
    census, manifest, audit = pinned_inputs()
    primary = sorted(census["primary_objects"], key=identity)
    require(len(primary) == 27 and len({identity(row) for row in primary}) == 27,
            "primary pilot archive census differs")
    manifest_rows = {(r["policy"], r["seed"], r["columns"]): r for r in manifest["bases"] if r["a"] == 0}
    s3_rows = {r["uri"]: r for r in audit["objects"]}
    entries = []
    sets = {}
    for row in primary:
        key = identity(row)
        source = manifest_rows[key]
        remote = s3_rows[row["s3_uri"]]
        require(source["compressed_blake3"] == row["compressed_blake3"]
                and source["s3_uri"] == row["s3_uri"]
                and remote["status"] == "PASS"
                and remote["observed_blake3"] == row["compressed_blake3"]
                and remote["observed_bytes"] == source["bytes"],
                f"S3 custody differs: {key}")
        path = objects_dir / (row["compressed_blake3"] + ".jsonl.gz")
        hashes, compressed_sha = archive_orbits(path, source)
        sets[key] = set(hashes)
        entries.append({"policy": key[0], "seed": key[1], "columns": key[2],
                        "s3_uri": row["s3_uri"],
                        "compressed_blake3": row["compressed_blake3"],
                        "compressed_sha256": compressed_sha,
                        "orbit_set_sha256": sha256("\n".join(sorted(hashes)).encode())})
    pairs = []
    for base in primary:
        bk = identity(base)
        for pool in primary:
            pk = identity(pool)
            if bk[:2] == pk[:2]:
                continue
            shared = len(sets[bk] & sets[pk])
            pairs.append({"base": list(bk), "pool": list(pk), "shared_orbits": shared,
                          "base_only_orbits": bk[2] - shared,
                          "residual_orbits": pk[2] - shared,
                          "union_orbits": bk[2] + pk[2] - shared,
                          "residual_points": 166 * (pk[2] - shared)})
    require(len(pairs) == 648, "ordered cross-stream pair grid differs")
    return {"schema": "n83.primary-cross-stream-orbit-overlap/v1",
            "status": "PASS_27_ARCHIVES_648_EXACT_PAIRS",
            "input_sha256": {name: digest for name, (_, digest) in SOURCES.items()},
            "generator_sha256": sha256(Path(__file__).read_bytes()),
            "archive_count": 27, "pair_count": 648,
            "archives": entries, "pairs": pairs,
            "selected_best_total_runtime": None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--objects-dir", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = make(args.objects_dir)
    content = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode()
    output = HERE / "overlap.json"
    if args.write:
        with output.open("xb") as stream:
            stream.write(content)
    else:
        require(output.read_bytes() == content, "overlap receipt differs on replay")
    print(json.dumps({"status": result["status"], "archives": 27, "pairs": 648,
                      "receipt_sha256": sha256(content)}))


if __name__ == "__main__":
    main()

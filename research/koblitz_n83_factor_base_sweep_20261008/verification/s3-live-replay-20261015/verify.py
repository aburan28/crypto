#!/usr/bin/env python3
"""Replay the live N83 S3 archive inventory against downloaded bytes."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path, PurePosixPath

from blake3 import blake3


HERE = Path(__file__).resolve().parent
PILOT = HERE.parents[1] / "pilot-01" / "manifest.json"
BUCKET = "crypto-autoresearcher"
PREFIX = "factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/"
SUFFIX = ".jsonl.gz"
CHUNK = 1 << 20


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def hashes(path: Path) -> tuple[str, str, str, int]:
    compressed_blake3 = blake3()
    compressed_sha256 = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(CHUNK), b""):
            compressed_blake3.update(chunk)
            compressed_sha256.update(chunk)
    plain_blake3 = blake3()
    plain_bytes = 0
    # Reading through EOF also verifies the gzip stream and CRC.
    with gzip.open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(CHUNK), b""):
            plain_blake3.update(chunk)
            plain_bytes += len(chunk)
    return (compressed_blake3.hexdigest(), compressed_sha256.hexdigest(),
            plain_blake3.hexdigest(), plain_bytes)


def replay(inventory_path: Path, download_root: Path) -> dict:
    inventory_bytes = inventory_path.read_bytes()
    inventory = json.loads(inventory_bytes)
    listed = inventory["Objects"]
    require(isinstance(listed, list), "missing S3 listing")
    require(len({item["Key"] for item in listed}) == len(listed),
            "duplicate S3 listing key")
    objects = {item["Key"]: item for item in listed
               if item["Key"].endswith(SUFFIX)}
    require(len(objects) == 57, "expected 57 stored N83 archives")
    pilot_bytes = PILOT.read_bytes()
    pilot = json.loads(pilot_bytes)
    bases = pilot["bases"]
    require(len(bases) == 54, "expected 54 pilot bases")
    expected = {}
    for base in bases:
        uri = base["s3_uri"]
        prefix = f"s3://{BUCKET}/"
        require(uri.startswith(prefix), "pilot archive bucket differs")
        key = uri[len(prefix):]
        require(key not in expected, "duplicate pilot archive")
        expected[key] = base
    require(not (expected.keys() - objects.keys()), "pilot archive absent from S3 listing")

    rows = []
    local_names = set()
    for key, item in sorted(objects.items()):
        require(key.startswith(PREFIX), f"unexpected archive prefix: {key}")
        relative = PurePosixPath(key[len(PREFIX):])
        require(not relative.is_absolute() and ".." not in relative.parts,
                f"unsafe archive path: {key}")
        local_names.add(relative.as_posix())
        path = download_root.joinpath(*relative.parts)
        require(path.is_file(), f"download absent: {key}")
        require(path.stat().st_size == item["Size"], f"S3 size mismatch: {key}")
        b3, sha, plain_b3, plain_size = hashes(path)
        require(path.name == b3 + SUFFIX, f"content address mismatch: {key}")
        base = expected.get(key)
        if base is not None:
            require(base["bytes"] == item["Size"], f"pilot size mismatch: {key}")
            require(base["compressed_blake3"] == b3,
                    f"pilot compressed digest mismatch: {key}")
            require(base["plain_blake3"] == plain_b3,
                    f"pilot plain digest mismatch: {key}")
        rows.append({"s3_uri": f"s3://{BUCKET}/{key}",
                     "compressed_bytes": item["Size"],
                     "compressed_blake3": b3,
                     "compressed_sha256": sha,
                     "plain_blake3": plain_b3,
                     "plain_bytes": plain_size,
                     "pilot_manifest_member": base is not None})
    actual_names = {p.relative_to(download_root).as_posix()
                    for p in download_root.rglob(f"*{SUFFIX}")}
    require(local_names == actual_names, "downloaded archive set differs from S3 listing")
    require(sum(row["pilot_manifest_member"] for row in rows) == 54,
            "pilot manifest count mismatch")
    return {"schema": "n83.live-s3-factor-base-replay/v1",
            "status": "PASS_57_LIVE_ARCHIVES",
            "bucket": BUCKET,
            "prefix": PREFIX,
            "listed_object_count": len(listed),
            "archive_count": len(rows),
            "pilot_archive_count": 54,
            "additional_archive_count": len(rows) - 54,
            "compressed_bytes": sum(row["compressed_bytes"] for row in rows),
            "plain_bytes": sum(row["plain_bytes"] for row in rows),
            "inventory_sha256": hashlib.sha256(inventory_bytes).hexdigest(),
            "pilot_manifest_sha256": hashlib.sha256(pilot_bytes).hexdigest(),
            "archives": rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=HERE / "list.json")
    parser.add_argument("--download-root", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, default=HERE / "audit.json")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    receipt = replay(args.inventory, args.download_root)
    if args.check:
        require(json.loads(args.receipt.read_text()) == receipt,
                "saved receipt differs from fresh replay")
    else:
        with args.receipt.open("x", encoding="utf-8") as stream:
            json.dump(receipt, stream, indent=2, sort_keys=True)
            stream.write("\n")
    print(json.dumps({key: value for key, value in receipt.items()
                      if key != "archives"}, sort_keys=True))


if __name__ == "__main__":
    main()

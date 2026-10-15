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
SOURCE_MANIFESTS = {
    "pilot": PILOT,
    "v2_k1182": HERE.parent / "matched-k1182-k2048-900s-20261011"
    / "evidence-01" / "panel-k1182" / "manifest.json",
    "v2_k2048": HERE.parent / "matched-k1182-k2048-900s-20261011"
    / "evidence-01" / "panel-k2048" / "manifest.json",
    "v3_structured_d12": HERE.parent / "structured-d12-orbit-20261010"
    / "manifest.json",
}
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
    expected = {}
    source_hashes = {}
    for source, manifest_path in SOURCE_MANIFESTS.items():
        raw = manifest_path.read_bytes()
        source_hashes[source] = hashlib.sha256(raw).hexdigest()
        parsed = json.loads(raw)
        bases = parsed["bases"] if "bases" in parsed else [parsed]
        require(len(bases) == (54 if source == "pilot" else 1),
                f"unexpected source-manifest base count: {source}")
        for base in bases:
            uri = base["s3_uri"]
            prefix = f"s3://{BUCKET}/"
            require(uri.startswith(prefix), f"archive bucket differs: {source}")
            key = uri[len(prefix):]
            require(key not in expected, f"duplicate source archive: {key}")
            expected[key] = (source, base)
    require(len(expected) == 57 and expected.keys() == objects.keys(),
            "live archive keys differ from all source manifests")

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
        source, base = expected[key]
        recorded_size = base["bytes"] if "bytes" in base else base["compressed_bytes"]
        require(recorded_size == item["Size"], f"source size mismatch: {key}")
        require(base["compressed_blake3"] == b3,
                f"source compressed digest mismatch: {key}")
        require(base["plain_blake3"] == plain_b3,
                f"source plain digest mismatch: {key}")
        rows.append({"s3_uri": f"s3://{BUCKET}/{key}",
                     "compressed_bytes": item["Size"],
                     "compressed_blake3": b3,
                     "compressed_sha256": sha,
                     "plain_blake3": plain_b3,
                     "plain_bytes": plain_size,
                     "source_manifest": source,
                     "pilot_manifest_member": source == "pilot"})
    actual_names = {p.relative_to(download_root).as_posix()
                    for p in download_root.rglob(f"*{SUFFIX}")}
    require(local_names == actual_names, "downloaded archive set differs from S3 listing")
    require(sum(row["pilot_manifest_member"] for row in rows) == 54,
            "pilot manifest count mismatch")
    return {"schema": "n83.live-s3-factor-base-replay/v2",
            "status": "PASS_57_LIVE_ARCHIVES_ALL_SOURCE_MANIFESTS",
            "bucket": BUCKET,
            "prefix": PREFIX,
            "listed_object_count": len(listed),
            "archive_count": len(rows),
            "pilot_archive_count": 54,
            "additional_archive_count": len(rows) - 54,
            "compressed_bytes": sum(row["compressed_bytes"] for row in rows),
            "plain_bytes": sum(row["plain_bytes"] for row in rows),
            "inventory_sha256": hashlib.sha256(inventory_bytes).hexdigest(),
            "pilot_manifest_sha256": source_hashes["pilot"],
            "source_manifest_sha256": source_hashes,
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

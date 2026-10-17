#!/usr/bin/env python3
"""Reconstruct and hash-check every file in the N83 finalist raw archive."""

import argparse
import hashlib
import json
import os
import shutil
import tarfile
from pathlib import Path, PurePosixPath

HERE = Path(__file__).resolve().parent


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe(name):
    parts = PurePosixPath(name).parts
    if (not parts or any(part in ("", ".", "..") for part in parts)
            or PurePosixPath(name).is_absolute()):
        raise ValueError(f"unsafe archive path: {name}")
    return Path(*parts)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError(f"output already exists: {output}")
    manifest_path = HERE / "raw_dedup_manifest.json"
    receipt = json.loads((HERE / "raw_dedup_receipt.json").read_text())
    manifest = json.loads(manifest_path.read_text())
    archive = HERE / "raw_dedup.tar.gz"
    if (manifest["schema"] != "n83.finalist-dedup-archive/v1"
            or manifest["run_id"] != 39005111289
            or sha256(archive) != receipt["archive_sha256"]
            or archive.stat().st_size != receipt["archive_bytes"]
            or sha256(manifest_path) != receipt["manifest_sha256"]):
        raise ValueError("archive identity differs")
    entries = manifest["entries"]
    blobs = {item["blob"] for item in entries.values()}
    for name, item in entries.items():
        safe(name)
        safe(item["blob"])
        if item["blob"] not in entries:
            raise ValueError(f"missing canonical blob for {name}")
        canonical = entries[item["blob"]]
        if (canonical["blob"] != item["blob"]
                or canonical["sha256"] != item["sha256"]
                or canonical["bytes"] != item["bytes"]):
            raise ValueError(f"duplicate binding differs: {name}")
    output.mkdir(parents=True)
    with tarfile.open(archive, "r:gz") as tar:
        all_members = tar.getmembers()
        members = {member.name: member for member in all_members}
        if (len(members) != len(all_members)
                or set(members) != blobs | {"MANIFEST.json"}):
            raise ValueError("stored blobs or manifest differ")
        manifest_member = tar.extractfile(members["MANIFEST.json"])
        if manifest_member is None or manifest_member.read() != manifest_path.read_bytes():
            raise ValueError("embedded manifest differs")
        for name in sorted(blobs):
            member = members[name]
            if not member.isfile() or member.size != entries[name]["bytes"]:
                raise ValueError(f"stored member differs: {name}")
            stream = tar.extractfile(member)
            if stream is None:
                raise ValueError(f"stored member unavailable: {name}")
            dest = output / safe(name)
            dest.parent.mkdir(parents=True, exist_ok=True)
            with stream, dest.open("wb") as sink:
                shutil.copyfileobj(stream, sink, 1024 * 1024)
            if sha256(dest) != entries[name]["sha256"]:
                raise ValueError(f"stored member hash differs: {name}")
    for name, item in sorted(entries.items()):
        if name != item["blob"]:
            dest = output / safe(name)
            dest.parent.mkdir(parents=True, exist_ok=True)
            os.link(output / safe(item["blob"]), dest)
    if (len(entries) != receipt["file_count"]
            or len(blobs) != receipt["stored_blob_count"]
            or len({item["sha256"] for name, item in entries.items()
                    if name.endswith("/base_dump.jsonl")}) !=
            receipt["unique_base_dump_count"]):
        raise ValueError("archive counts differ")
    print(json.dumps({"status": "PASS", "file_count": len(entries),
                      "stored_blob_count": len(blobs),
                      "unique_base_dump_count": receipt["unique_base_dump_count"]},
                     sort_keys=True))


if __name__ == "__main__":
    main()

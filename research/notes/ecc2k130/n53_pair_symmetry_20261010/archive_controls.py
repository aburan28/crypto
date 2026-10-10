#!/usr/bin/env python3
"""Losslessly archive and independently hash the development controls."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
import tarfile


HERE = Path(__file__).resolve().parent
CONTROLS = HERE / "controls"
ARCHIVE = HERE / "controls.tar.gz"
MANIFEST = HERE / "controls_manifest.json"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def verify() -> dict:
    manifest = json.loads(MANIFEST.read_text())
    assert digest(ARCHIVE.read_bytes()) == manifest["archive_sha256"]
    with tarfile.open(ARCHIVE, "r:gz") as archive:
        members = archive.getmembers()
        assert [member.name for member in members] == sorted(manifest["files"])
        for member in members:
            assert member.isfile()
            data = archive.extractfile(member).read()
            assert digest(data) == manifest["files"][member.name]["sha256"]
            assert len(data) == manifest["files"][member.name]["bytes"]
    return manifest


def create() -> dict:
    assert CONTROLS.is_dir() and not ARCHIVE.exists() and not MANIFEST.exists()
    paths = sorted(path for path in CONTROLS.rglob("*") if path.is_file())
    assert paths and all(not path.is_symlink() for path in paths)
    files = {}
    with ARCHIVE.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w|") as archive:
                for path in paths:
                    name = str(path.relative_to(HERE))
                    data = path.read_bytes()
                    files[name] = {"sha256": digest(data), "bytes": len(data)}
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    info.mtime = info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    info.mode = 0o644
                    archive.addfile(info, io.BytesIO(data))
    manifest = {"schema": "n53-pair-symmetry-controls-v1", "files": files,
                "archive_sha256": digest(ARCHIVE.read_bytes())}
    MANIFEST.write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
    verify()
    shutil.rmtree(CONTROLS)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--create", action="store_true")
    args = parser.parse_args()
    manifest = create() if args.create else verify()
    print(json.dumps({"status": "PASS", "files": len(manifest["files"]),
                      "archive_sha256": manifest["archive_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()

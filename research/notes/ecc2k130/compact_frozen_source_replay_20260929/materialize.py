#!/usr/bin/env python3
"""Materialize the frozen source tree of a compact-orbit evaluation.

A compact-orbit FROZEN.json pins the SHA-256 of every source file its
evaluation builds and runs, including shared files (the index-calculus crate,
Cargo.toml, the workflow itself) that later PRs legitimately edit.  Rather
than re-pinning the freeze, the evaluation is replayed from a copy of the
checkout in which each drifted pinned file is replaced by a committed gzip
snapshot of its frozen bytes.  Every pinned file in the copy must hash to the
frozen digest; a drifted file without a matching snapshot fails closed.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
SNAPSHOT_DIR = HERE / "historical_snapshots"
# koblitz_index_calculus.rs, the largest snapshot, is 669,217 bytes; the cap
# guards against a corrupted gzip stream, not a size any source approaches.
MAX_SNAPSHOT_BYTES = 8 * 1024 * 1024


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_manifest(snapshot_dir: Path = SNAPSHOT_DIR) -> dict[str, dict]:
    manifest = json.loads((snapshot_dir.parent / "SNAPSHOTS.json").read_text())
    assert manifest["schema"] == "compact-frozen-source-snapshots-v1"
    return {entry["sha256"]: entry for entry in manifest["snapshots"]}


def snapshot_bytes(path: str, digest: str, manifest: dict[str, dict],
                   snapshot_dir: Path = SNAPSHOT_DIR) -> bytes:
    entry = manifest.get(digest)
    assert entry is not None, f"{path}: drifted and no snapshot of {digest}"
    assert entry["path"] == path, f"{path}: snapshot {digest} is {entry['path']}"
    packed = (snapshot_dir / f"{digest}.gz").read_bytes()
    assert sha(packed) == entry["gzip_sha256"], f"{path}: snapshot gzip hash"
    with gzip.open(snapshot_dir / f"{digest}.gz", "rb") as stream:
        raw = stream.read(MAX_SNAPSHOT_BYTES + 1)
    assert len(raw) <= MAX_SNAPSHOT_BYTES, f"{path}: snapshot too large"
    assert len(raw) == entry["bytes"], f"{path}: snapshot length"
    assert sha(raw) == digest, f"{path}: snapshot content hash"
    return raw


def materialize(root: Path, freezes: list[Path], dest: Path,
                snapshot_dir: Path = SNAPSHOT_DIR) -> dict:
    assert not dest.exists(), "never overwrite a materialized tree"
    pins: dict[str, str] = {}
    for freeze in freezes:
        for path, digest in json.loads(freeze.read_text())["source_sha256"].items():
            assert pins.setdefault(path, digest) == digest, f"{path}: conflicting pins"
    manifest = load_manifest(snapshot_dir)
    shutil.copytree(root, dest, symlinks=True,
                    ignore=shutil.ignore_patterns(".git", "target"))
    replayed, live = {}, []
    for path, digest in sorted(pins.items()):
        target = dest / path
        if target.is_file() and sha(target.read_bytes()) == digest:
            live.append(path)
            continue
        current = sha(target.read_bytes()) if target.is_file() else None
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(snapshot_bytes(path, digest, manifest, snapshot_dir))
        replayed[path] = {"frozen_sha256": digest, "live_sha256": current}
    for path, digest in pins.items():
        assert sha((dest / path).read_bytes()) == digest, path
    return {"schema": "compact-frozen-source-materialization-v1",
            "freezes": [str(freeze.relative_to(root)) for freeze in freezes],
            "pinned_files": len(pins), "live": live, "replayed_from_snapshot": replayed}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", type=Path, action="append", required=True)
    parser.add_argument("--dest", type=Path, required=True)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    root = HERE.parents[3]
    receipt = materialize(root, [freeze.resolve() for freeze in args.freeze],
                          args.dest.resolve())
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.receipt:
        args.receipt.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()

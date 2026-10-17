#!/usr/bin/env python3
"""Select exact frozen N83 pilot objects from a fresh hosted regeneration."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import time

import blake3


EXPECTED = {
    "3888e2597fd3519807c9d722f49f83329c93ce5cb0929501cf7cb7cbc3cab1a0":
        ("public_x_sequential", 2026100802,
         "8463aaf30365d6f6f5d6b189aff4e8afcb6c8305ec27c216700776cbdccba29d"),
    "966a0493c96674f0af6d184615ed05309d08fba5c7342957b4186cd50366dc47":
        ("public_x_gray_prefix", 2026100801,
         "df02963d1706fa7a64fd7558d88db79b23ea34801ac58062c4c500113e8c57fe"),
    "b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e":
        ("public_x_hash", 2026100801,
         "3e85c9e419fcddd3c9faed3c84470bf96fcad9178e48e1d73629369335841ecf"),
}


def digest(path: Path, constructor) -> str:
    state = constructor()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            state.update(block)
    return state.hexdigest()


def plain_blake3(path: Path) -> str:
    state = blake3.blake3()
    with gzip.open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            state.update(block)
    return state.hexdigest()


def prepare(panel: Path, archives: Path, source_commit: str,
            workflow_commit: str, started_unix_seconds: int) -> dict:
    manifest_path = panel / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if (manifest.get("schema") != "n83.factor-base-panel/v1"
            or manifest.get("status") != "completed_factor_base_panel"
            or manifest.get("completed_base_count") != 54
            or manifest.get("expected_base_count") != 54
            or manifest.get("source_commit") != source_commit
            or not isinstance(manifest.get("bases"), list)
            or len(manifest["bases"]) != 54):
        raise ValueError("fresh pilot manifest is incomplete or from another source")
    rows = {row.get("compressed_blake3"): row for row in manifest["bases"]}
    if len(rows) != 54:
        raise ValueError("duplicate pilot object hashes")
    if archives.exists():
        raise ValueError("archive output must be fresh")
    checked = []
    for object_hash, (policy, seed, expected_sha256) in EXPECTED.items():
        row = rows.get(object_hash)
        object_name = f"{object_hash}.jsonl.gz"
        if (row is None or row.get("a") != 0 or row.get("columns") != 600
                or row.get("policy") != policy or row.get("seed") != seed
                or row.get("object") != f"objects/{object_name}"):
            raise ValueError(f"frozen construction identity changed: {object_hash}")
        path = panel / "objects" / object_name
        if (path.stat().st_size != row.get("bytes")
                or digest(path, blake3.blake3) != object_hash
                or digest(path, hashlib.sha256) != expected_sha256
                or plain_blake3(path) != row.get("plain_blake3")):
            raise ValueError(f"frozen archive bytes changed: {object_hash}")
        checked.append((object_name, path, expected_sha256))
    archives.mkdir(parents=True)
    for object_name, path, expected_sha256 in checked:
        with path.open("rb") as source, (archives / object_name).open("xb") as target:
            shutil.copyfileobj(source, target, 1 << 20)
        if digest(archives / object_name, hashlib.sha256) != expected_sha256:
            raise ValueError(f"copied archive differs: {object_name}")
    return {
        "schema": "n83.hosted-regenerated-archives/v1",
        "status": "PASS_3_FROZEN_OBJECTS",
        "transport": "fresh_producer_regeneration_from_pinned_checkout",
        "source_commit": source_commit,
        "workflow_dispatch_head_sha": workflow_commit,
        "pilot_manifest_sha256": digest(manifest_path, hashlib.sha256),
        "pilot_completed_base_count": 54,
        "pilot_construction_seconds_L0": manifest.get("elapsed_seconds_informational"),
        "hosted_preparation_wall_seconds_L0": time.time() - started_unix_seconds,
        "archive_sha256": {name: sha for name, _, sha in checked},
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--archives", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--workflow-commit", required=True)
    parser.add_argument("--started-unix-seconds", type=int, required=True)
    args = parser.parse_args()
    result = prepare(args.panel, args.archives, args.source_commit,
                     args.workflow_commit, args.started_unix_seconds)
    with args.receipt.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write("\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

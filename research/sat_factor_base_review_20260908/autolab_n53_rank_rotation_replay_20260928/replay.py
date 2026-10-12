#!/usr/bin/env python3
"""Replay the released n53 rotation freeze with pinned historical source bytes.

The original `check_protocol.py` and `verify_archive.py` run unmodified. Only
the frozen sources listed in SNAPSHOTS are read from committed gzip snapshots
instead of the live tree; each snapshot must decompress to the exact SHA-256
recorded in the original FROZEN.json. It never builds or runs an n53 producer.
"""
from __future__ import annotations

import argparse
import contextlib
import gzip
import hashlib
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ORIGINAL = HERE.parent / "autolab_n53_rank_rotation_20260925"
ORIGINAL_FROZEN_SHA256 = "e787ea56ed6948dcd9dcbe52266aff76d94c4ee9412ebd7a5bd6dd9d8561fa92"
RELEASE_SOURCE_HEAD = "9759251e79797478386e0dad503da4286574a0b8"

# Original `source_sha256` keys whose bytes are pinned by a snapshot taken at
# RELEASE_SOURCE_HEAD. `rust_producer` and `cargo_toml` are shared top-level
# files that unrelated work on main keeps editing; `workflow` is the PR
# workflow, whose hash-only job now routes through this replay. Every other
# frozen source and input lives in this research thread's own directories and
# stays verified against the live tree by the original checker.
SNAPSHOTS = {
    "rust_producer": "rust_producer.rs.gz",
    "cargo_toml": "cargo_toml.toml.gz",
    "workflow": "workflow.yml.gz",
}
SNAPSHOT_DIR = HERE / "historical_snapshots"
# Largest snapshot (rust_producer) decompresses to 292,406 bytes; the cap only
# bounds a truncated or corrupted stream.
SNAPSHOT_MAX_BYTES = 1024 * 1024
WORKFLOW = ".github/workflows/n53-target-cyclic-rank.yml"

NEW_FILES = (
    ("PROTOCOL.md", "replay.py", "test_replay.py")
    + tuple(f"historical_snapshots/{name}" for name in SNAPSHOTS.values())
    + (WORKFLOW,)
)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def original_modules():
    if str(ORIGINAL) not in sys.path:
        sys.path.insert(0, str(ORIGINAL))
    import check_protocol
    import verify_archive
    assert Path(check_protocol.__file__).resolve().parent == ORIGINAL
    assert Path(verify_archive.__file__).resolve().parent == ORIGINAL
    return check_protocol, verify_archive


def verify_freeze() -> dict:
    own = json.loads((HERE / "FROZEN.json").read_text())
    assert set(own["files"]) == set(NEW_FILES)
    for name in NEW_FILES:
        path = ROOT / name if name.startswith(".github/") else HERE / name
        assert sha(path.read_bytes()) == own["files"][name], name
    raw = (ORIGINAL / "FROZEN.json").read_bytes()
    assert sha(raw) == ORIGINAL_FROZEN_SHA256, "original FROZEN.json changed"
    frozen = json.loads(raw)
    assert frozen["schema"] == "n53_target_cyclic_rank_factorial_freeze_v1"
    assert frozen["status"] == "released_for_one_outcome"
    assert set(SNAPSHOTS) <= set(frozen["source_sha256"])
    return frozen


def read_snapshot(key: str, frozen: dict, snapshot_dir: Path = SNAPSHOT_DIR) -> bytes:
    with gzip.open(snapshot_dir / SNAPSHOTS[key], "rb") as stream:
        raw = stream.read(SNAPSHOT_MAX_BYTES + 1)
    assert len(raw) <= SNAPSHOT_MAX_BYTES, key
    assert sha(raw) == frozen["source_sha256"][key], key
    return raw


@contextlib.contextmanager
def pinned_sources(frozen: dict, snapshot_dir: Path = SNAPSHOT_DIR):
    check_protocol, _ = original_modules()
    live = check_protocol.sources
    with tempfile.TemporaryDirectory(prefix="n53-rotation-snapshots-") as temporary:
        pinned = {}
        for key, name in SNAPSHOTS.items():
            path = Path(temporary) / name[: -len(".gz")]
            path.write_bytes(read_snapshot(key, frozen, snapshot_dir))
            pinned[key] = path

        def sources() -> dict[str, Path]:
            paths = live()
            assert set(pinned) <= set(paths)
            paths.update(pinned)
            return paths

        check_protocol.sources = sources
        try:
            yield pinned
        finally:
            check_protocol.sources = live


def live_digests() -> dict[str, str]:
    check_protocol, _ = original_modules()
    paths = check_protocol.sources()
    return {key: sha(paths[key].read_bytes()) for key in SNAPSHOTS}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("freeze", "preflight", "archive"), required=True)
    parser.add_argument("--bundle", type=Path)
    args = parser.parse_args()
    frozen = verify_freeze()
    if args.mode == "freeze":
        print("REPLAY_AND_ORIGINAL_FREEZE_PASS")
        return
    check_protocol, verify_archive = original_modules()
    live = live_digests()
    with pinned_sources(frozen):
        if args.mode == "preflight":
            value = check_protocol.preflight()
            result = {
                "status": "PINNED_SNAPSHOT_PREFLIGHT_PASS",
                "release_main_head": value["release_main_head"],
                "targets": value["schedule_proof"]["target_count"],
            }
        else:
            assert args.bundle is not None, "--bundle is required for archive mode"
            result = verify_archive.verify(args.bundle)
    result["release_source_head"] = RELEASE_SOURCE_HEAD
    result["pinned_sha256"] = {key: frozen["source_sha256"][key] for key in SNAPSHOTS}
    result["live_sha256"] = live
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

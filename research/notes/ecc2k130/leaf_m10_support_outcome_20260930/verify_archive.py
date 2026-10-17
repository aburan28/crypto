#!/usr/bin/env python3
"""Verify the complete immutable hosted artifact and separate-host receipt."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parent / "leaf_m10_support_20260930/evidence"
RAW_RUN = "run_36707051160"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    manifest = json.loads((HERE / "MANIFEST.json").read_text())
    run = json.loads((HERE / "HOSTED_RUN.json").read_text())
    result = json.loads((EVIDENCE / "result.json").read_text())
    hosted = json.loads((EVIDENCE / "replay.json").read_text())
    mac = json.loads((EVIDENCE / "mac_replay.json").read_text())
    assert manifest["schema"] == "ecc2k130-leaf-m10-support-archive-v1"
    assert run["status"] == "completed" and run["conclusion"] == "success"
    assert run["headSha"] == result["source_head"] == manifest["source_head"]
    assert run["url"] == manifest["hosted_run_url"]
    assert result["status"] == "PASS_CENSUS"
    assert hosted == mac and hosted["status"] == "PASS"
    assert hosted["result_sha256"] == sha(EVIDENCE / "result.json")
    assert len(manifest["files"]) == 102
    paths = set()
    hosted_bytes = 0
    hosted_files = 0
    for item in manifest["files"]:
        relative = Path(item["path"])
        assert not relative.is_absolute() and ".." not in relative.parts
        assert item["path"] not in paths
        paths.add(item["path"])
        path = EVIDENCE / relative
        assert path.is_file() and path.stat().st_size == item["bytes"]
        assert sha(path) == item["sha256"], path
        if item["provenance"] == "hosted_artifact":
            hosted_files += 1
            hosted_bytes += item["bytes"]
        else:
            assert item["provenance"] == "mac_replay"
            assert item["path"] == "mac_replay.json"
    # The sibling raw-run directory is a separately checksummed receipt copy
    # (its own SHA256SUMS); the manifest covers everything outside it.
    present = {path.relative_to(EVIDENCE).as_posix()
               for path in EVIDENCE.rglob("*")
               if path.is_file() and RAW_RUN not in path.relative_to(EVIDENCE).parts}
    assert paths == present
    raw = EVIDENCE / RAW_RUN
    if raw.is_dir():
        listed = set()
        for line in (raw / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split("  ", 1)
            assert sha(raw / name) == digest, name
            listed.add(name)
        assert listed == {path.relative_to(raw).as_posix()
                          for path in raw.rglob("*")
                          if path.is_file() and path.name != "SHA256SUMS"}
    assert hosted_files == manifest["hosted_artifact_files"] == 101
    assert hosted_bytes == manifest["hosted_artifact_bytes"]
    print(f"PASS {len(paths)} files, {hosted_bytes} hosted bytes")


if __name__ == "__main__":
    main()

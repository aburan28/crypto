#!/usr/bin/env python3
"""Freeze check and archive-only independent complete-tuple replay."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

from run import HERE, check_freeze, sha, tree_hashes


def substantive(value):
    if isinstance(value, dict):
        return {k: substantive(v) for k, v in value.items()
                if not k.endswith("wall_seconds") and not k.endswith("cpu_seconds")
                and not k.endswith("peak_rss_bytes")}
    if isinstance(value, list):
        return [substantive(v) for v in value]
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    freeze = check_freeze()
    if args.evidence is None:
        print("Frozen corpus, parent source and merged #762 archive hashes verified; no outcome read.")
        return
    receipt = json.loads((args.evidence / "receipt.json").read_text())
    assert receipt["status"] == "success" and receipt["freeze"] == freeze
    assert [item["name"] for item in receipt["commands"]] == [
        "n13-m5-producer", "n13-m5-verify", "n19-m6-producer", "n19-m6-verify"]
    assert all(item["exit_code"] == 0 for item in receipt["commands"])
    tar_path = args.evidence / "raw.tar.gz"
    assert sha(tar_path) == (args.evidence / "raw.tar.gz.sha256").read_text().strip()
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        with tarfile.open(tar_path, "r:gz") as stream:
            members = stream.getmembers()
            assert members and all(member.isfile() or member.isdir() for member in members)
            assert all(Path(member.name).parts[0] == "raw" and ".." not in Path(member.name).parts
                       for member in members)
            stream.extractall(root, filter="data")
        raw = root / "raw"
        assert tree_hashes(raw) == receipt["raw_sha256"]
        for arm in ("n13-m5", "n19-m6"):
            replay = root / f"{arm}-fresh-verify.json"
            subprocess.run([sys.executable, str(HERE / "verify.py"), "--arm", arm,
                            "--archive", str(raw / arm), "--out", str(replay)], check=True)
            archived = json.loads((raw / f"{arm}-verify.json").read_text())
            fresh = json.loads(replay.read_text())
            assert substantive(archived) == substantive(fresh), arm
    print("Both full point histograms, projected misses and SHA-selected target streams replayed.")


if __name__ == "__main__":
    main()

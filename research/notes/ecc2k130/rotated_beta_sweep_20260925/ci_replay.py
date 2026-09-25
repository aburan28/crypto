#!/usr/bin/env python3
"""Hash-only preflight or archive-only independent complete-tuple replay."""
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
    subprocess.run([sys.executable, str(HERE / "verify.py"), "--selection-only"],
                   check=True, stdout=subprocess.PIPE, timeout=90)
    if args.evidence is None:
        print("Frozen hash stream, selection, source/input and #767 archive verified; no support outcome read.")
        return
    selection = json.loads((HERE / "selection.json").read_text())
    betas = [row["beta"] for row in selection["selected"]]
    receipt = json.loads((args.evidence / "receipt.json").read_text())
    assert receipt["status"] == "success" and receipt["freeze"] == freeze
    expected_commands = ["selection-verify"]
    for beta in betas:
        expected_commands.extend([f"beta-{beta}-producer", f"beta-{beta}-verify"])
    assert [item["name"] for item in receipt["commands"]] == expected_commands
    assert all(item["exit_code"] == 0 for item in receipt["commands"])
    expected_files = {f"beta-{beta}/{name}" for beta in betas
                      for name in ("factors.json", "fixed_targets.json", "target_counts.u32le",
                                   "full_histogram.jsonl", "projected_histogram.jsonl", "summary.json")}
    expected_files.update(f"beta-{beta}-verify.json" for beta in betas)
    tar_path = args.evidence / "raw.tar.gz"
    assert sha(tar_path) == (args.evidence / "raw.tar.gz.sha256").read_text().strip()
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        with tarfile.open(tar_path, "r:gz") as stream:
            members = stream.getmembers()
            assert members and all(member.isfile() or member.isdir() for member in members)
            assert all(not Path(member.name).is_absolute() and
                       Path(member.name).parts[0] == "raw" and
                       ".." not in Path(member.name).parts for member in members)
            names = {str(Path(member.name).relative_to("raw")) for member in members if member.isfile()}
            assert names == expected_files
            stream.extractall(root, filter="data")
        raw = root / "raw"
        assert tree_hashes(raw) == receipt["raw_sha256"]
        for beta in betas:
            replay = root / f"beta-{beta}-fresh-verify.json"
            subprocess.run([sys.executable, str(HERE / "verify.py"),
                            "--beta", str(beta), "--archive", str(raw / f"beta-{beta}"),
                            "--out", str(replay)], check=True, timeout=630)
            archived = json.loads((raw / f"beta-{beta}-verify.json").read_text())
            fresh = json.loads(replay.read_text())
            assert substantive(archived) == substantive(fresh), beta
    print("All four full histograms, all q subgroup targets, and all #767 fixed point labels replayed.")


if __name__ == "__main__":
    main()

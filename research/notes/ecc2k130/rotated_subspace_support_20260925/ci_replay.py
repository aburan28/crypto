#!/usr/bin/env python3
"""Verify frozen sources/inputs and independently replay an archived outcome."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCES = ["gate.py", "verify.py", "run.py", "ci_replay.py"]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tree_hashes(root: Path) -> dict[str, str]:
    return {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob("*")) if p.is_file()}


def substantive(obj):
    if isinstance(obj, dict):
        return {key: substantive(value) for key, value in obj.items()
                if not key.endswith("wall_seconds") and not key.endswith("cpu_seconds")
                and not key.endswith("peak_rss_bytes")}
    if isinstance(obj, list):
        return [substantive(value) for value in obj]
    return obj


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    sources = {name: sha(HERE / name) for name in SOURCES}
    inputs = tree_hashes(HERE / "inputs")
    assert sources == frozen["source_sha256"], "source hash drift"
    assert inputs == frozen["input_sha256"], "input hash drift"
    manifest = json.loads((HERE / "inputs" / "input_manifest.json").read_text())
    assert manifest["domain"] == "ECC2K130-ROTATED-SUBSPACE-20260925-v1"
    assert len(manifest["cells"]) == 5
    for cell in manifest["cells"]:
        for kind in ("cov", "density"):
            entry = cell[kind]
            assert sha(HERE / "inputs" / entry["path"]) == entry["sha256"]
    if args.evidence is None:
        print("Frozen source and input hashes verified; no outcome archive yet.")
        return
    receipt = json.loads((args.evidence / "receipt.json").read_text())
    assert receipt["status"] == "success"
    assert receipt["source_sha256"] == sources
    assert receipt["input_sha256"] == inputs
    assert [item["name"] for item in receipt["commands"]] == ["toy", "density", "verify"]
    assert all(item["exit_code"] == 0 for item in receipt["commands"])
    tar_path = args.evidence / "raw.tar.gz"
    assert sha(tar_path) == (args.evidence / "raw.tar.gz.sha256").read_text().strip()
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        with tarfile.open(tar_path, "r:gz") as stream:
            members = stream.getmembers()
            assert all(Path(member.name).parts[0] == "raw" for member in members)
            assert all(member.isfile() or member.isdir() for member in members)
            stream.extractall(root, filter="data")
        raw = root / "raw"
        assert tree_hashes(raw) == receipt["raw_sha256"]
        replay = root / "independent_replay.json"
        subprocess.run([sys.executable, str(HERE / "verify.py"),
                        "--archive", str(raw), "--inputs", str(HERE / "inputs"),
                        "--output", str(replay)], check=True)
        archived = json.loads((raw / "verify_report.json").read_text())
        fresh = json.loads(replay.read_text())
        assert substantive(archived) == substantive(fresh), "independent replay changed"
    print("Frozen archive hashes, all exact small-curve decisions, witnesses, and density rows replayed.")


if __name__ == "__main__":
    main()

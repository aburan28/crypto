#!/usr/bin/env python3
"""Fail-closed source freeze; archive-only exact independent replay."""
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
REPO = HERE.parents[3]
PDP = HERE.parent / "rotated_pdp_corpus_20260925"
SUBSPACE = HERE.parent / "rotated_subspace_support_20260925"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_freeze() -> dict:
    frozen_path = HERE / "FROZEN.json"
    assert sha(frozen_path) == (HERE / "FREEZE.sha256").read_text().strip()
    frozen = json.loads(frozen_path.read_text())
    for filename, expected in frozen["local_sha256"].items():
        assert sha(HERE / filename) == expected, filename
    for name, expected in frozen["parent_sha256"].items():
        assert sha((PDP if name.startswith("pdp/") else SUBSPACE) / name.split("/", 1)[1]) == expected, name
    assert sha(PDP / "evidence" / "raw.tar.gz") == frozen["parent_archive_sha256"]
    assert sha(HERE.parent / "rotated_slot_placement_20260925" / "FROZEN.json") == frozen["fullrank_parent_freeze_sha256"]
    assert sha(HERE.parents[3] / ".github/workflows/ecc2k130-rotated-nearthreshold-slot.yml") == frozen["workflow_sha256"]
    subprocess.run([sys.executable, str(HERE.parent / "rotated_slot_placement_20260925" / "ci_replay.py")],
                   check=True, timeout=60)
    subprocess.run([sys.executable, str(HERE / "verify_preflight.py"),
                    "--expected", str(HERE / "PREFLIGHT.json")], check=True, timeout=60)
    return frozen


def hashes(root: Path) -> dict[str, str]:
    return {str(path.relative_to(root)): sha(path) for path in sorted(root.rglob("*")) if path.is_file()}


def stable(value):
    if isinstance(value, dict):
        return {key: stable(item) for key, item in value.items()
                if key not in ("wall_seconds", "cpu_seconds", "peak_rss_bytes")}
    if isinstance(value, list):
        return [stable(item) for item in value]
    return value


def replay(evidence: Path) -> None:
    freeze = check_freeze()
    receipt = json.loads((evidence / "receipt.json").read_text())
    assert receipt["status"] == "success"
    assert receipt["freeze_sha256"] == sha(HERE / "FROZEN.json")
    assert all(item["exit_code"] == 0 for item in receipt["children"])
    archive = evidence / "raw.tar.gz"
    assert sha(archive) == (evidence / "raw.tar.gz.sha256").read_text().strip()
    with tempfile.TemporaryDirectory() as name:
        root = Path(name)
        with tarfile.open(archive, "r:gz") as tar:
            members = tar.getmembers()
            assert members and all((m.isfile() or m.isdir()) and Path(m.name).parts[0] == "raw"
                                   and ".." not in Path(m.name).parts for m in members)
            tar.extractall(root, filter="data")
        raw = root / "raw"
        assert hashes(raw) == receipt["raw_sha256"]
        for position in freeze["positions"]:
            folder = raw / f"p{position}"
            verified = root / f"p{position}-fresh-verify.json"
            subprocess.run([sys.executable, str(HERE / "verify.py"), "--position", str(position),
                            "--data", str(folder), "--out", str(verified)], check=True, timeout=630)
            assert stable(json.loads(verified.read_text())) == stable(json.loads((folder / "verify.json").read_text()))
        fresh = root / "analysis-fresh.json"
        subprocess.run([sys.executable, str(HERE / "analyze.py"), "--raw", str(raw),
                        "--out", str(fresh)], check=True, timeout=30)
        assert json.loads(fresh.read_text()) == json.loads((raw / "analysis.json").read_text())
    print("All six complete subgroup arrays, 64 shared targets per arm, ranks and pair decisions replayed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    check_freeze()
    if args.evidence:
        replay(args.evidence)
    else:
        print("Frozen source and independently replayed factor-only preflight verified; no six-sum outcome read.")

#!/usr/bin/env python3
"""Recheck a committed clean panel, including fresh independent replay."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


HERE = Path(__file__).resolve().parent
# The frozen bytes of pinned sources that later commits changed, kept in the
# repository by the same snapshot store the compact-orbit replays use.
SNAPSHOT_DIR = HERE / "source_snapshots" / "historical_snapshots"
_MATERIALIZE = HERE.parent / "compact_frozen_source_replay_20260929" / "materialize.py"


def _snapshot_loader():
    spec = importlib.util.spec_from_file_location("frozen_source_materialize", _MATERIALIZE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_sources(root: Path, pins: dict[str, str], snapshot_dir: Path = SNAPSHOT_DIR) -> dict:
    """Every pinned source must still exist as the exact bytes the panel ran.

    A file still matching its pin in the tree passes as it is.  A file a
    later commit changed passes only if a committed snapshot holds exactly
    the pinned bytes; a drifted file without one fails closed.  This checks
    that the panel's source is durably available, not that today's tree is
    the tree it ran on, which every later edit to a shared file would break.
    """
    loader = _snapshot_loader()
    manifest = loader.load_manifest(snapshot_dir)
    live, snapshot = [], []
    for path, checksum in sorted(pins.items()):
        target = root / path
        if target.is_file() and digest(target.read_bytes()) == checksum:
            live.append(path)
            continue
        loader.snapshot_bytes(path, checksum, manifest, snapshot_dir)
        snapshot.append(path)
    return {"live": live, "snapshot": snapshot}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--bundle", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    bundle = args.bundle.resolve()
    archive = bundle / "clean_evidence.tar.gz"
    external = json.loads((bundle / "clean_archive_manifest.json").read_text())
    summary_bytes = (bundle / "clean_summary.json").read_bytes()
    assert digest(archive.read_bytes()) == external["archive_sha256"]
    assert archive.stat().st_size == external["archive_bytes"]
    assert digest(summary_bytes) == external["summary_sha256"]
    summary = json.loads(summary_bytes)
    expected = {f"n{n}_seed{i}" for n in (37, 41) for i in range(3)}
    assert len(summary["runs"]) == len(expected)
    assert {f"n{r['n']}_seed{r['seed_index']}" for r in summary["runs"]} == expected

    with tarfile.open(archive, "r:gz") as tar:
        members = {m.name: m for m in tar if m.isfile()}
        payload = {name: tar.extractfile(member).read() for name, member in members.items()}
    assert payload["clean_summary.json"] == summary_bytes
    lines = payload["SHA256SUMS"].decode().splitlines()
    hashed = {}
    for line in lines:
        checksum, name = line.split("  ", 1)
        assert name not in hashed
        hashed[name] = checksum
        assert digest(payload[name]) == checksum
    assert set(hashed) == set(payload) - {"SHA256SUMS"}
    assert {name.split("/", 1)[0] for name in payload if "/" in name} == expected

    verifier = Path(__file__).with_name("verify.py")
    with tempfile.TemporaryDirectory() as temporary:
        tmp = Path(temporary)
        for row in summary["runs"]:
            name = f"n{row['n']}_seed{row['seed_index']}"
            manifest = json.loads(payload[f"{name}/manifest.json"])
            receipt = json.loads(payload[f"{name}/replay.json"])
            assert manifest["clean_checkout"] and manifest["includes_pinned_main_ref"]
            assert manifest["replay_verdict"] == receipt["verdict"] == "PASS"
            sources = check_sources(root, manifest["source_sha256"])
            assert receipt["q"] == row["q"]
            assert receipt["recovered_scalar"] == row["recovered_scalar"]
            for arm in ("ic", "rho"):
                data = payload[f"{name}/{arm}.jsonl"]
                assert digest(data) == receipt[f"{arm}_stdout_sha256"]
                (tmp / f"{arm}.jsonl").write_bytes(data)
            result = subprocess.run(
                [sys.executable, str(verifier), str(tmp / "ic.jsonl"),
                 str(tmp / "rho.jsonl"), "--out", str(tmp / "replay.json")],
                cwd=root, capture_output=True, text=True,
            )
            assert result.returncode == 0, (name, result.stderr)
            assert json.loads((tmp / "replay.json").read_text()) == receipt, name
            print(f"{name}: fresh independent replay PASS ({len(sources['live'])} pinned sources "
                  f"live, {len(sources['snapshot'])} from committed snapshots)", flush=True)
    print("Committed clean-source archive PASS", flush=True)


if __name__ == "__main__":
    main()

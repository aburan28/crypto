#!/usr/bin/env python3
"""Rehash the v2 archive and independently replay every successful cell."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from archive import equivalent  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/disjoint_cold_v2_20261001"))
import verify_cold  # noqa: E402
from analyze_cold import analyze  # noqa: E402


def extract(raw: Path, destination: Path) -> None:
    with tarfile.open(raw, "r:gz") as bundle:
        for member in bundle:
            assert member.isfile() and len(Path(member.name).parts) >= 2
            assert not Path(member.name).is_absolute()
            assert ".." not in Path(member.name).parts
            stream = bundle.extractfile(member)
            assert stream is not None
            path = destination / member.name
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("wb") as output:
                while chunk := stream.read(1 << 20):
                    output.write(chunk)


def replay(archive: Path, analysis_path: Path) -> dict:
    expected = json.loads(analysis_path.read_text())
    actual = analyze(archive)
    assert equivalent(actual, expected)
    manifest = json.loads((archive / "MANIFEST.json").read_text())
    passed = []
    with tempfile.TemporaryDirectory(prefix="disjoint-v2-archive-") as temporary:
        destination = Path(temporary)
        for cell, entry in manifest["cases"].items():
            if entry["status"] != "PASS":
                continue
            extract(archive / entry["raw_path"], destination)
            run_dir = destination / f"disjoint-cold-v2-{cell}" / cell
            replayed = verify_cold.verify(cell, run_dir, relocated=True)
            hosted = json.loads((archive / entry["receipt_path"]).read_text())
            assert equivalent(replayed, hosted), cell
            if entry["second_host_matches_hosted"] is True:
                second = json.loads((archive / entry["second_host_replay_path"]).read_text())
                assert equivalent(replayed, second), cell
            passed.append(cell)
    return {"status": "PASS", "replayed_successful_cells": passed,
            "censored_cells": [cell for cell, entry in manifest["cases"].items()
                               if entry["status"] != "PASS"],
            "manifest_sha256": actual["manifest_sha256"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(replay(args.archive.resolve(), args.analysis.resolve()),
                     sort_keys=True))

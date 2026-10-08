#!/usr/bin/env python3
"""Rehash the sealed run and independently replay every archived cold child."""
from __future__ import annotations

import json
from pathlib import Path
import sys
import tarfile
import tempfile
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from archive import equivalent  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/disjoint_cold_q_20260930"))
import verify_cold  # noqa: E402
from analyze_cold import analyze  # noqa: E402
from pairing_addendum import replay as pairing_replay  # noqa: E402

ARCHIVE = HERE / "evidence_run_36794339148"


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


def replay() -> dict:
    expected = json.loads((HERE / "ANALYSIS.json").read_text())
    actual = analyze(ARCHIVE)
    assert equivalent(actual, expected)
    manifest = json.loads((ARCHIVE / "MANIFEST.json").read_text())
    passed = []
    with tempfile.TemporaryDirectory(prefix="disjoint-cold-archive-") as temporary:
        dest = Path(temporary)
        # Replay successful hosted receipts with the unchanged frozen verifier.
        # The pairing addendum changes the verifier module's target hook, so
        # deliberately run the failed-cell diagnostic only after these checks.
        for cell, entry in manifest["cases"].items():
            if entry["status"] != "PASS":
                continue
            extract(ARCHIVE / entry["raw_path"], dest)
            run_dir = dest / f"disjoint-cold-{cell}" / cell
            replayed = verify_cold.verify(cell, run_dir, relocated=True)
            hosted = json.loads((ARCHIVE / entry["receipt_path"]).read_text())
            second = json.loads((ARCHIVE / entry["second_host_replay_path"]).read_text())
            assert equivalent(replayed, hosted) and equivalent(replayed, second), cell
            passed.append(cell)
        assert len(passed) == 5

        cell = "n37_L1024"
        entry = manifest["cases"][cell]
        assert entry["status"] == "FAIL"
        extract(ARCHIVE / entry["raw_path"], dest)
        run_dir = dest / f"disjoint-cold-{cell}" / cell
        try:
            verify_cold.verify(cell, run_dir, relocated=True)
        except AssertionError:
            failure = traceback.format_exc()
            assert 'assert set(record["pinned_intermediates"])' in failure
        else:
            raise AssertionError("frozen verifier unexpectedly accepted failed cell")
        addendum = pairing_replay(run_dir)
        archived = json.loads((ARCHIVE / manifest["diagnostics"][
            "n37_L1024_pairing_addendum.json"]["path"]).read_text())
        assert equivalent(addendum, archived)
    return {"status": "PASS", "frozen_pass_cells": passed,
            "frozen_failed_cell": "n37_L1024",
            "pairing_addendum_verified_target_logs": 15360,
            "manifest_sha256": actual["manifest_sha256"]}


if __name__ == "__main__":
    print(json.dumps(replay(), sort_keys=True))

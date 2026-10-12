#!/usr/bin/env python3
"""Fail-closed checks for the compact-orbit frozen-source materializer."""
from __future__ import annotations

import gzip
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import materialize as m


class MaterializeTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = Path(tempfile.mkdtemp())
        self.root = self.tmp / "root"
        (self.root / "src").mkdir(parents=True)
        self.frozen = b"frozen bytes\n"
        (self.root / "src/a.rs").write_bytes(b"drifted bytes\n")
        (self.root / "src/b.rs").write_bytes(b"unchanged\n")
        self.freeze = self.root / "FROZEN.json"
        self.freeze.write_text(json.dumps({"source_sha256": {
            "src/a.rs": m.sha(self.frozen), "src/b.rs": m.sha(b"unchanged\n")}}))
        self.snapshots = self.tmp / "replay/historical_snapshots"
        self.snapshots.mkdir(parents=True)
        self.write_snapshot(self.frozen)

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp)

    def write_snapshot(self, raw: bytes, path: str = "src/a.rs",
                       digest: str | None = None) -> None:
        digest = digest or m.sha(raw)
        packed = gzip.compress(raw, mtime=0)
        (self.snapshots / f"{digest}.gz").write_bytes(packed)
        (self.snapshots.parent / "SNAPSHOTS.json").write_text(json.dumps({
            "schema": "compact-frozen-source-snapshots-v1",
            "snapshots": [{"path": path, "sha256": digest, "bytes": len(raw),
                           "gzip_sha256": m.sha(packed)}]}))

    def run_materialize(self) -> dict:
        return m.materialize(self.root, [self.freeze], self.tmp / "out", self.snapshots)

    def test_drifted_file_replayed_from_snapshot(self) -> None:
        receipt = self.run_materialize()
        self.assertEqual((self.tmp / "out/src/a.rs").read_bytes(), self.frozen)
        self.assertEqual(list(receipt["replayed_from_snapshot"]), ["src/a.rs"])
        self.assertEqual(receipt["live"], ["src/b.rs"])
        self.assertEqual((self.root / "src/a.rs").read_bytes(), b"drifted bytes\n")

    def test_missing_snapshot_fails_closed(self) -> None:
        (self.root / "src/b.rs").write_bytes(b"also drifted\n")
        with self.assertRaises(AssertionError):
            self.run_materialize()

    def test_tampered_snapshot_fails_closed(self) -> None:
        self.write_snapshot(b"not the frozen bytes\n", digest=m.sha(self.frozen))
        with self.assertRaises(AssertionError):
            self.run_materialize()

    def test_snapshot_for_other_path_rejected(self) -> None:
        self.write_snapshot(self.frozen, path="src/other.rs")
        with self.assertRaises(AssertionError):
            self.run_materialize()

    def test_evaluation_commit_restores_complete_build_tree(self) -> None:
        for path in ("examples", "experiments", "benches", "docs"):
            (self.root / path).mkdir()
            (self.root / path / "context.txt").write_text("frozen context\n")
        (self.root / "Cargo.toml").write_text("[package]\n")
        (self.root / "src/a.rs").write_bytes(b"base caller\n")
        (self.root / "src/b.rs").write_bytes(b"historical dependency\n")
        subprocess.run(["git", "init", "-q"], cwd=self.root, check=True)
        subprocess.run(["git", "config", "user.email", "ci@example.invalid"],
                       cwd=self.root, check=True)
        subprocess.run(["git", "config", "user.name", "CI test"],
                       cwd=self.root, check=True)
        subprocess.run(["git", "add", "Cargo.toml", "src", "examples",
                        "experiments", "benches", "docs"], cwd=self.root, check=True)
        subprocess.run(["git", "commit", "-qm", "evaluation source"],
                       cwd=self.root, check=True)
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=self.root, text=True).strip()
        (self.root / "src/a.rs").write_bytes(self.frozen)
        (self.root / "src/b.rs").write_bytes(b"incompatible current dependency\n")
        (self.root / "src/current_only.rs").write_bytes(b"untracked module\n")
        self.freeze.write_text(json.dumps({
            "evaluation_source_commit": commit,
            "source_sha256": {"src/a.rs": m.sha(self.frozen)},
        }))

        receipt = self.run_materialize()
        self.assertEqual((self.tmp / "out/src/a.rs").read_bytes(), self.frozen)
        self.assertEqual((self.tmp / "out/src/b.rs").read_bytes(),
                         b"historical dependency\n")
        self.assertFalse((self.tmp / "out/src/current_only.rs").exists())
        self.assertEqual(receipt["live"], [])
        self.assertEqual(receipt["replayed_from_snapshot"], {
            "src/a.rs": {"frozen_sha256": m.sha(self.frozen),
                         "live_sha256": m.sha(self.frozen)}})

        self.freeze.write_text(json.dumps({
            "source_base_main": commit,
            "source_sha256": {"src/a.rs": m.sha(self.frozen)},
        }))
        fallback = m.materialize(
            self.root, [self.freeze], self.tmp / "out-base", self.snapshots)
        self.assertEqual((self.tmp / "out-base/src/b.rs").read_bytes(),
                         b"historical dependency\n")
        self.assertEqual(fallback["replayed_from_snapshot"], receipt["replayed_from_snapshot"])


if __name__ == "__main__":
    unittest.main()

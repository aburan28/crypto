#!/usr/bin/env python3
"""The committed panel's source pins: live, snapshotted, or failed closed."""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import shutil
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("verify_committed", HERE / "verify_committed.py")
vc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vc)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class CheckSourcesTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = Path(tempfile.mkdtemp())
        self.root = self.tmp / "repo"
        (self.root / "src").mkdir(parents=True)
        self.snapshots = self.tmp / "store" / "historical_snapshots"
        self.snapshots.mkdir(parents=True)
        self.frozen = b"fn frozen() {}\n"
        self.pins = {"src/a.rs": sha(self.frozen)}
        self.write_manifest([])

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp)

    def write_manifest(self, entries: list[dict]) -> None:
        (self.snapshots.parent / "SNAPSHOTS.json").write_text(json.dumps({
            "schema": "compact-frozen-source-snapshots-v1", "snapshots": entries}))

    def snapshot(self, raw: bytes, path: str = "src/a.rs", stored: bytes | None = None) -> None:
        digest = sha(raw)
        packed = gzip.compress(raw if stored is None else stored, mtime=0)
        (self.snapshots / f"{digest}.gz").write_bytes(packed)
        self.write_manifest([{"path": path, "sha256": digest, "bytes": len(raw),
                              "gzip_sha256": sha(packed)}])

    def check(self) -> dict:
        return vc.check_sources(self.root, self.pins, self.snapshots)

    def test_live_file_matching_its_pin_passes(self) -> None:
        (self.root / "src/a.rs").write_bytes(self.frozen)
        self.assertEqual(self.check(), {"live": ["src/a.rs"], "snapshot": []})

    def test_drifted_file_with_its_snapshot_passes(self) -> None:
        (self.root / "src/a.rs").write_bytes(b"fn edited_later() {}\n")
        self.snapshot(self.frozen)
        self.assertEqual(self.check(), {"live": [], "snapshot": ["src/a.rs"]})

    def test_drifted_file_without_a_snapshot_fails_closed(self) -> None:
        (self.root / "src/a.rs").write_bytes(b"fn edited_later() {}\n")
        with self.assertRaises(AssertionError):
            self.check()

    def test_deleted_file_without_a_snapshot_fails_closed(self) -> None:
        with self.assertRaises(AssertionError):
            self.check()

    def test_snapshot_with_other_bytes_fails_closed(self) -> None:
        (self.root / "src/a.rs").write_bytes(b"fn edited_later() {}\n")
        self.snapshot(self.frozen, stored=b"fn tampered() {}\n")
        with self.assertRaises(AssertionError):
            self.check()

    def test_snapshot_filed_under_another_path_fails_closed(self) -> None:
        (self.root / "src/a.rs").write_bytes(b"fn edited_later() {}\n")
        self.snapshot(self.frozen, path="src/b.rs")
        with self.assertRaises(AssertionError):
            self.check()


class CommittedStoreTest(unittest.TestCase):
    def test_every_committed_snapshot_reproduces_its_pin(self) -> None:
        manifest = json.loads((vc.SNAPSHOT_DIR.parent / "SNAPSHOTS.json").read_text())
        loader = vc._snapshot_loader()
        index = loader.load_manifest(vc.SNAPSHOT_DIR)
        for entry in manifest["snapshots"]:
            raw = loader.snapshot_bytes(entry["path"], entry["sha256"], index, vc.SNAPSHOT_DIR)
            self.assertEqual(sha(raw), entry["sha256"])


if __name__ == "__main__":
    unittest.main()

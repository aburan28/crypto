#!/usr/bin/env python3
"""Fail-closed controls for the pinned-snapshot n53 rotation replay."""
from __future__ import annotations

import copy
import gzip
import shutil
import tempfile
import unittest
from pathlib import Path

import replay


class ReplayControls(unittest.TestCase):
    def setUp(self) -> None:
        self.frozen = replay.verify_freeze()
        self.check_protocol, _ = replay.original_modules()

    def test_each_snapshot_is_checked_against_its_own_frozen_hash(self) -> None:
        for key in replay.SNAPSHOTS:
            with self.subTest(key=key):
                frozen = copy.deepcopy(self.frozen)
                frozen["source_sha256"][key] = "0" * 64
                with self.assertRaisesRegex(AssertionError, key):
                    with replay.pinned_sources(frozen):
                        pass

    def test_corrupted_or_oversized_snapshot_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            for name in replay.SNAPSHOTS.values():
                shutil.copy(replay.SNAPSHOT_DIR / name, directory / name)
            (directory / replay.SNAPSHOTS["rust_producer"]).write_bytes(
                gzip.compress(b"fn main() {}\n", mtime=0))
            with self.assertRaisesRegex(AssertionError, "rust_producer"):
                with replay.pinned_sources(self.frozen, directory):
                    pass
            (directory / replay.SNAPSHOTS["rust_producer"]).write_bytes(
                gzip.compress(b"\0" * (replay.SNAPSHOT_MAX_BYTES + 1), mtime=0))
            with self.assertRaisesRegex(AssertionError, "rust_producer"):
                with replay.pinned_sources(self.frozen, directory):
                    pass

    def test_only_snapshot_keys_are_substituted_and_then_restored(self) -> None:
        live_function = self.check_protocol.sources
        live = live_function()
        with replay.pinned_sources(self.frozen) as pinned:
            substituted = self.check_protocol.sources()
            self.assertEqual(set(substituted), set(live))
            for key, path in substituted.items():
                if key in replay.SNAPSHOTS:
                    self.assertEqual(path, pinned[key])
                    self.assertNotEqual(path, live[key])
                else:
                    self.assertEqual(path, live[key])
        self.assertIs(self.check_protocol.sources, live_function)

    def test_original_checker_still_rejects_unpinned_producer_bytes(self) -> None:
        live_function = self.check_protocol.sources
        with tempfile.TemporaryDirectory() as temporary:
            drifted = Path(temporary) / "koblitz_s5_sat_instance.rs"
            drifted.write_bytes(b"// not the released producer\n")

            def sources():
                paths = live_function()
                paths["rust_producer"] = drifted
                return paths

            self.check_protocol.sources = sources
            try:
                with self.assertRaisesRegex(AssertionError, "source hash drift: rust_producer"):
                    self.check_protocol.preflight()
            finally:
                self.check_protocol.sources = live_function


if __name__ == "__main__":
    unittest.main()

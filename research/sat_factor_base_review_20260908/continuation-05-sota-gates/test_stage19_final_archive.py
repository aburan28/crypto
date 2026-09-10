#!/usr/bin/env python3
"""Regression checks for offline, relocated historical evidence replay."""
import json
from pathlib import Path, PurePosixPath
import unittest
from unittest.mock import patch

import verify_stage19_final_archive as archive


class HistoricalArchiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.context = archive.detached_archive(archive.REPO)
        cls.snapshot = cls.context.__enter__()

    @classmethod
    def tearDownClass(cls):
        cls.context.__exit__(None, None, None)

    def test_full_replay_without_original_host_ca(self):
        original = Path.read_bytes

        def no_live_ca(path):
            if str(path) == "/etc/ssl/cert.pem":
                raise AssertionError("historical replay read the live CA bundle")
            return original(path)

        with patch.object(Path, "read_bytes", no_live_ca):
            summary = archive.replay(self.snapshot)
        self.assertEqual("nine_clean_requests_complete", summary["artifact_status"])
        self.assertEqual(15, summary["verified_public_service_f4_terminal_coverage"])
        self.assertFalse(summary["full_solver_matrix_gate_passed"])

    def test_wrong_python_and_checkout_are_rejected(self):
        for field, replacement in (
            ("HISTORICAL_PYTHON", "/usr/bin/python3"),
            ("HISTORICAL_ROOT", PurePosixPath("/tmp/substituted-checkout")),
        ):
            with self.subTest(field=field), patch.object(archive, field, replacement):
                with self.assertRaises(archive.ArchiveError):
                    archive.verify_historical_commands(self.snapshot)

    def test_wrong_task_name_is_rejected(self):
        plan = json.loads((self.snapshot / archive.ARTIFACT / "plan.json").read_text())
        task = plan["tasks"][0]
        task["named_input"]["path"] = "inputs/../../substituted.magma"
        with self.assertRaises(archive.ArchiveError):
            archive.historical_command(task)

    def test_same_size_response_mutation_is_rejected(self):
        record = next((self.snapshot / archive.ARTIFACT / "attempts").glob("*/response.xml"))
        original = record.read_bytes()
        try:
            record.write_bytes(original.replace(b"SAT", b"BAD", 1))
            with self.assertRaises(archive.ArchiveError):
                archive.verify_current_archive(self.snapshot)
        finally:
            record.write_bytes(original)
        self.assertEqual(96, archive.verify_current_archive(self.snapshot)["files"])

    def test_unexpected_artifact_is_rejected(self):
        extra = self.snapshot / archive.ARTIFACT / "unexpected.txt"
        try:
            extra.write_text("extra")
            with self.assertRaises(archive.ArchiveError):
                archive.verify_current_archive(self.snapshot)
        finally:
            extra.unlink()

    def test_worktree_cleanup_on_failure(self):
        before = archive.git(archive.REPO, "worktree", "list", "--porcelain")
        with self.assertRaisesRegex(RuntimeError, "synthetic failure"):
            with archive.detached_archive(archive.REPO) as snapshot:
                self.assertTrue(snapshot.is_dir())
                raise RuntimeError("synthetic failure")
        self.assertFalse(snapshot.exists())
        self.assertEqual(before, archive.git(archive.REPO, "worktree", "list", "--porcelain"))


if __name__ == "__main__":
    unittest.main()

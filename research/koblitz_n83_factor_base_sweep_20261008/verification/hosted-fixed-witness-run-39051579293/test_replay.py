"""Tamper controls for the retained hosted fixed-witness artifact replay."""

import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import replay  # noqa: E402


class ReplayTamperTests(unittest.TestCase):
    def setUp(self):
        default_artifact = (
            "/private/tmp/n83-hosted-witness-run-39051579293/"
            "n83-hosted-witness-39051579293-1"
        )
        self.artifact = Path(os.environ.get("N83_WITNESS_ARTIFACT", default_artifact))
        self.s3_tar = Path(os.environ.get(
            "N83_WITNESS_S3_TAR", "/private/tmp/n83-hosted-witness-run-39051579293-s3.tar.gz"))
        if not self.artifact.is_dir() or not self.s3_tar.is_file():
            self.skipTest("downloaded run and S3 roundtrip required")

    def test_saved_replay_is_exact(self):
        result = replay.replay(self.artifact, self.s3_tar)
        self.assertEqual(result, json.loads((HERE / "audit.json").read_text()))
        self.assertEqual(result["status"], "PASS_FIXED_WITNESS_SAT_AND_S3_ROUNDTRIP")

    def test_false_solver_status_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            copied = Path(temporary) / self.artifact.name
            shutil.copytree(self.artifact, copied)
            result_path = copied / "solve.stdout"
            result = json.loads(result_path.read_text())
            result["status"] = "SAT_unverified"
            result_path.write_text(json.dumps(result))
            with self.assertRaisesRegex(ValueError, "raw output hash differs"):
                replay.replay(copied, self.s3_tar)

    def test_changed_s3_archive_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            copied = Path(temporary) / "changed.tar.gz"
            data = bytearray(self.s3_tar.read_bytes())
            data[-1] ^= 1
            copied.write_bytes(data)
            with self.assertRaisesRegex(ValueError, "S3 roundtrip receipt differs"):
                replay.replay(self.artifact, copied)


if __name__ == "__main__":
    unittest.main()

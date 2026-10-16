"""Fail-closed byte checks for the pending N83 v7 control inputs."""

import gzip
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import blake3
import preflight


class ControlInputPreflightTests(unittest.TestCase):
    def test_saved_receipt_keeps_run_unexecuted(self):
        row = json.loads(preflight.OUTPUT.read_text())
        self.assertEqual(row["status"], "PASS_local_byte_inputs_and_prior_evidence")
        self.assertEqual(row["prior_conservative_charge_seconds"], 19800)
        self.assertEqual(row["proposed_additional_cap_seconds"], 1200)
        self.assertFalse(row["control_run_executed_by_preflight"])
        self.assertFalse(row["s3_head_is_content_hash_proof"])
        self.assertEqual(row["small"]["compressed_sha256"],
                         preflight.gate.prior.SMALL_SHA)
        self.assertEqual(row["envelope"]["compressed_sha256"],
                         preflight.gate.prior.ENVELOPE_SHA)

    def test_corrupt_compressed_or_plain_archive_fails(self):
        plain = b"pinned control archive fixture\n"
        compressed = gzip.compress(plain, mtime=0)
        digest = hashlib.sha256(compressed).hexdigest()
        manifest = {"completed_base_count": 1, "bases": [{
            "bytes": len(compressed),
            "columns": 2,
            "points": 332,
            "s3_uri": "s3://example/pinned",
            "compressed_blake3": blake3.blake3(compressed).hexdigest(),
            "plain_blake3": blake3.blake3(plain).hexdigest(),
        }]}
        with TemporaryDirectory() as directory:
            path = Path(directory) / "archive.gz"
            path.write_bytes(compressed)
            with patch.dict(preflight.HEAD_SNAPSHOT,
                            {"K1182": {"content_length": len(compressed),
                                       "etag": "synthetic"}}):
                result = preflight.archive(path, manifest, "K1182", digest)
                self.assertEqual(result["plain_bytes"], len(plain))
                path.write_bytes(compressed[:-1] + bytes([compressed[-1] ^ 1]))
                with self.assertRaisesRegex(ValueError,
                                            "compressed archive changed"):
                    preflight.archive(path, manifest, "K1182", digest)
                path.write_bytes(compressed)
                manifest["bases"][0]["plain_blake3"] = "0" * 64
                with self.assertRaisesRegex(ValueError, "plain archive changed"):
                    preflight.archive(path, manifest, "K1182", digest)


if __name__ == "__main__":
    unittest.main()

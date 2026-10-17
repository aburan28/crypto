from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from run_remaining_census import (
    SegmentPaths,
    require_exact_file_set,
    segment_paths,
    sha256_file,
    validate_verification,
)


class RemainingCensusTests(unittest.TestCase):
    def test_segment_paths_are_fixed_width_and_disjoint(self) -> None:
        root = Path("/evidence")
        paths = segment_paths(root, "campaign-v1", 16_704, 33_088)
        self.assertEqual(paths.segment.name, "campaign-v1-segment-000000016704-000000033088.jsonl.gz")
        self.assertEqual(paths.checkpoint.name, "campaign-v1-checkpoint-000000033088.json")
        self.assertEqual(paths.verification_receipt.parent, root / "verification")

    def test_partial_generation_set_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            present = root / "segment"
            present.write_bytes(b"x")
            with self.assertRaisesRegex(RuntimeError, "partial generation evidence"):
                require_exact_file_set([present, root / "checkpoint", root / "receipt"], "generation")

    def test_verification_binds_range_settings_and_digests(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            segment = root / "segment.gz"
            checkpoint = root / "checkpoint.json"
            receipt = root / "verify.json"
            segment.write_bytes(b"segment")
            checkpoint.write_text("{}\n", encoding="utf-8")
            payload = {
                "status": "pass",
                "segment_digest": {"sha256": sha256_file(segment)},
                "checkpoint_digest": {"sha256": sha256_file(checkpoint)},
                "receipt": {
                    "audit_points": 2,
                    "audit_seed_x": 7,
                    "summary": {"x_start": 10, "x_end": 12, "node_records": 36},
                },
            }
            receipt.write_text(json.dumps(payload), encoding="utf-8")
            paths = SegmentPaths(segment, checkpoint, root / "generate.json", receipt, root / "manifest.json")
            self.assertEqual(validate_verification(paths, 10, 12)["status"], "pass")
            payload["receipt"]["audit_seed_x"] = 8
            receipt.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "audit settings mismatch"):
                validate_verification(paths, 10, 12)


if __name__ == "__main__":
    unittest.main()

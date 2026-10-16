"""Mutation checks for the planted N83 chained-S3 preflight."""

import json
import unittest

import replay


class ChainPreflightTests(unittest.TestCase):
    def test_saved_receipt_and_all_affine_chain(self):
        row = replay.replay()
        self.assertEqual(replay.canonical(row), replay.RECEIPT.read_bytes())
        self.assertEqual(row["chain"]["identity_mask"], 0)
        self.assertEqual(row["chain"]["s3_residues"], ["0"] * 4)
        self.assertEqual(row["ordered_global_point_indices"][:2], [166, 0])
        self.assertFalse(row["sat_search_executed"])

    def test_wrong_s3_output_is_rejected(self):
        fixture = json.loads(replay.FIXTURE.read_text())
        points = [replay.parse_point(term["point"]) for term in fixture["terms"]]
        points[:2] = points[1::-1]
        target = replay.parse_point(fixture["equation"]["target"])
        self.assertNotEqual(replay.semaev_s3(points[0][0], points[1][0],
                                             replay.point_add(points[0], points[1])[0] ^ 1), 0)
        with self.assertRaisesRegex(ValueError, "first pair violates lex order"):
            replay.check_chain(points[1::-1] + points[2:], target)

    def test_changed_fixture_bytes_are_rejected(self):
        from tempfile import TemporaryDirectory
        from pathlib import Path
        with TemporaryDirectory() as directory:
            changed = Path(directory) / "fixture.json"
            changed.write_bytes(replay.FIXTURE.read_bytes() + b" ")
            with self.assertRaisesRegex(ValueError, "fixture bytes changed"):
                replay.replay(changed)


if __name__ == "__main__":
    unittest.main()

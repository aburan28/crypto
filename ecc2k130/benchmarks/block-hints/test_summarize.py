import importlib.util
import pathlib
import tempfile
import unittest


HERE = pathlib.Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("block_summary", HERE / "summarize.py")
SUMMARY = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(SUMMARY)


class BlockHintParserTests(unittest.TestCase):
    def verify_text(self, block):
        return f"""packed table split forward: 1
packed table batch hints: 1
packed table block hints: {block}, queue 512
packed cycle fast2: 1
backend cuda-packed131: 96256 threads x 16 slots x 1 lanes = 1540096 walks, dp weight 48, 96 steps per launch
  finished: 1000.0 M it/s, 300 distinguished points (300 verified against the reference, 0 dropped)
"""

    def test_identity_and_geometry_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "verify.log"
            path.write_text(self.verify_text(1))
            parsed = SUMMARY.parse_verify(path, 1)
            self.assertTrue(parsed["valid"])
            self.assertEqual(parsed["blockHints"], 1)
            self.assertEqual(parsed["hintQueue"], 512)
            self.assertEqual(parsed["liveSlots"], 1540096)
            self.assertFalse(SUMMARY.parse_verify(path, 0)["valid"])
            path.write_text(self.verify_text(1).replace("queue 512", "queue 256"))
            self.assertFalse(SUMMARY.parse_verify(path, 1)["valid"])

    def test_timing_markers_follow_variant(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            control = root / "warmup-0-1-control.log"
            candidate = root / "warmup-0-2-candidate.log"
            control.write_text(self.verify_text(0))
            candidate.write_text(self.verify_text(1))
            rows = [
                {"phase": "warmup", "pair": 0, "order": 1, "variant": "control"},
                {"phase": "warmup", "pair": 0, "order": 2, "variant": "candidate"},
            ]
            self.assertEqual(SUMMARY.validate_block_timing_markers(root, rows), [])
            candidate.write_text(self.verify_text(0))
            self.assertTrue(SUMMARY.validate_block_timing_markers(root, rows))
            candidate.write_text(self.verify_text(1).replace("split forward: 1", "split forward: 0"))
            self.assertTrue(SUMMARY.validate_block_timing_markers(root, rows))
            candidate.write_text(self.verify_text(1).replace(
                "96256 threads x 16 slots", "48128 threads x 32 slots"))
            self.assertTrue(SUMMARY.validate_block_timing_markers(root, rows))


if __name__ == "__main__":
    unittest.main()

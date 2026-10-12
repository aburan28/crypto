import importlib.util
import pathlib
import struct
import tempfile
import unittest


HERE = pathlib.Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("fast2_summary", HERE / "summarize.py")
SUMMARY = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(SUMMARY)


class Fast2ParserTests(unittest.TestCase):
    def verify_text(self, fast2):
        return f"""packed table split forward: 1
packed table batch hints: 1
packed cycle fast2: {fast2}
backend cuda-packed131: 96256 threads x 16 slots x 1 lanes = 1540096 walks, dp weight 48, 96 steps per launch
  finished: 1000.0 M it/s, 300 distinguished points (300 verified against the reference, 0 dropped)
"""

    def test_fast2_identity_and_geometry_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "verify.log"
            path.write_text(self.verify_text(1))
            parsed = SUMMARY.parse_verify(path, 1)
            self.assertTrue(parsed["valid"])
            self.assertEqual(parsed["fast2"], 1)
            self.assertEqual(parsed["liveSlots"], 1540096)
            self.assertFalse(SUMMARY.parse_verify(path, 0)["valid"])
            path.write_text(self.verify_text(1).replace("x 16 slots", "x 32 slots"))
            self.assertFalse(SUMMARY.parse_verify(path, 1)["valid"])

    def test_v3_corpus_identity_ignores_record_order(self):
        a = bytes(range(32))
        b = bytes(reversed(range(32)))
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            one, two = root / "one.bin", root / "two.bin"
            header = struct.pack("<8sII", SUMMARY.BASE.MAGIC, 3, 32)
            one.write_bytes(header + a + b)
            two.write_bytes(header + b + a)
            x = SUMMARY.BASE.read_corpus(one)
            y = SUMMARY.BASE.read_corpus(two)
            self.assertEqual(x["records"], 2)
            self.assertEqual(x["sortedSha256"], y["sortedSha256"])

    def test_bounded_screen_gate_is_inherited(self):
        def row(phase, variant, rate, order):
            return {
                "phase": phase,
                "variant": variant,
                "rateMps": float(rate),
                "pair": 0,
                "order": order,
                "logSha256": "fixture",
                "gpuState": "fixture",
            }

        rows = [
            row("warmup", "control", 100, 1),
            row("warmup", "candidate", 100, 2),
            row("screen", "control", 100, 1),
            row("screen", "candidate", 100.5, 2),
            row("screen", "control", 99, 3),
        ]
        result = SUMMARY.BASE.decide_samples(rows)
        self.assertTrue(result["valid"])
        self.assertTrue(result["qualified"])
        rows[3]["rateMps"] = 100.49
        self.assertFalse(SUMMARY.BASE.decide_samples(rows)["qualified"])

    def test_preflight_requires_exact_reference_replayed_corpora(self):
        build = """ptxas info : Function properties for _ZN12eccPacked1314walkE10WalkParamsIjEPj
    400 bytes stack frame, 0 bytes spill stores, 0 bytes spill loads
ptxas info : Used 128 registers, used 1 barriers, 400 bytes cumulative stack size
ptxas info : Compile time = 1.0 ms
"""
        record = bytes(range(32))
        header = struct.pack("<8sII", SUMMARY.BASE.MAGIC, 3, 32)
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            for name, fast2 in (("control", 0), ("candidate", 1)):
                (root / f"build-{name}.log").write_text(build)
                (root / f"verify-{name}.log").write_text(self.verify_text(fast2))
                (root / f"dp-{name}.bin").write_bytes(header + record)
            for name in ("host.txt", "source-files.sha256", "binary-sha256.txt", "dp-identity.txt"):
                (root / name).write_text("fixture\n")
            result = SUMMARY.summarize(root, include_timing=False)
            self.assertTrue(result["valid"])
            self.assertTrue(result["corpusIdentity"])
            self.assertEqual(result["decision"], "preflight_pass")
            (root / "dp-candidate.bin").write_bytes(header + bytes(reversed(record)))
            mismatch = SUMMARY.summarize(root, include_timing=False)
            self.assertFalse(mismatch["valid"])
            self.assertIn("header-aware table-v3 corpus identities differ", mismatch["errors"])


if __name__ == "__main__":
    unittest.main()

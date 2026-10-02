import pathlib
import struct
import tempfile
import unittest

import summarize


class CorpusTests(unittest.TestCase):
    def write_corpus(self, path, records, magic=summarize.MAGIC, version=3, stride=32):
        path.write_bytes(struct.pack("<8sII", magic, version, stride) + b"".join(records))

    def test_header_is_removed_and_record_order_is_canonical(self):
        a = bytes(range(32))
        b = bytes(reversed(range(32)))
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            one, two = root / "one.bin", root / "two.bin"
            self.write_corpus(one, [a, b])
            self.write_corpus(two, [b, a])
            x, y = summarize.read_corpus(one), summarize.read_corpus(two)
            self.assertEqual(x["records"], 2)
            self.assertEqual(x["payloadBytes"], 64)
            self.assertEqual(x["sortedSha256"], y["sortedSha256"])

    def test_bad_header_and_partial_record_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "bad.bin"
            for data in (
                b"short",
                struct.pack("<8sII", b"WRONGHDR", 3, 32),
                struct.pack("<8sII", summarize.MAGIC, 2, 32),
                struct.pack("<8sII", summarize.MAGIC, 3, 72),
                struct.pack("<8sII", summarize.MAGIC, 3, 32) + b"x",
            ):
                path.write_bytes(data)
                with self.subTest(size=len(data)), self.assertRaises(ValueError):
                    summarize.read_corpus(path)


class QualificationTests(unittest.TestCase):
    @staticmethod
    def row(phase, variant, rate, pair=0, order=0):
        return {
            "phase": phase,
            "variant": variant,
            "rateMps": float(rate),
            "pair": pair,
            "order": order,
            "logSha256": "fixture",
            "gpuState": "fixture",
        }

    def base_rows(self, candidate):
        return [
            self.row("warmup", "control", 100, order=1),
            self.row("warmup", "candidate", 100, order=2),
            self.row("screen", "control", 100, order=1),
            self.row("screen", "candidate", candidate, order=2),
            self.row("screen", "control", 99, order=3),
        ]

    def test_half_percent_gate_and_unqualified_stop(self):
        qualified = summarize.decide_samples(self.base_rows(100.5))
        self.assertTrue(qualified["valid"])
        self.assertTrue(qualified["qualified"])
        self.assertEqual(qualified["decision"], "qualified_pending_confirmation")
        unqualified = summarize.decide_samples(self.base_rows(100.49))
        self.assertTrue(unqualified["valid"])
        self.assertFalse(unqualified["qualified"])
        self.assertEqual(unqualified["decision"], "unqualified_stop")

    def test_confirmation_order_and_ratios(self):
        rows = self.base_rows(101)
        rows += [
            self.row("confirm", "control", 100, 1, 1),
            self.row("confirm", "candidate", 102, 1, 2),
            self.row("confirm", "candidate", 101, 2, 1),
            self.row("confirm", "control", 100, 2, 2),
            self.row("confirm", "control", 100, 3, 1),
            self.row("confirm", "candidate", 103, 3, 2),
        ]
        result = summarize.decide_samples(rows)
        self.assertTrue(result["valid"])
        self.assertEqual(result["decision"], "confirmation_complete")
        self.assertEqual(result["pairedConfirmationRatios"], [1.02, 1.01, 1.03])

    def test_unqualified_confirmation_is_invalid(self):
        rows = self.base_rows(99)
        rows += [
            self.row("confirm", "control", 100, 1, 1),
            self.row("confirm", "candidate", 101, 1, 2),
        ]
        result = summarize.decide_samples(rows)
        self.assertFalse(result["valid"])
        self.assertEqual(result["decision"], "invalid")


class EvidenceParserTests(unittest.TestCase):
    def test_build_and_replay_identity_are_exact(self):
        build = """ptxas info : Function properties for _ZN12eccPacked1314walkE10WalkParamsIjEPj
    400 bytes stack frame, 0 bytes spill stores, 0 bytes spill loads
ptxas info : Used 128 registers, used 1 barriers, 400 bytes cumulative stack size
ptxas info : Compile time = 1.0 ms
"""
        verify = """packed table split forward: 1
  finished: 1000.0 M it/s, 300 distinguished points (300 verified against the reference, 0 dropped)
"""
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            (root / "build.log").write_text(build)
            (root / "verify.log").write_text(verify)
            resources = summarize.parse_build(root / "build.log")
            self.assertEqual(resources["registers"], 128)
            self.assertEqual(resources["stackFrameBytes"], 400)
            self.assertEqual(resources["spillStoreBytes"], 0)
            replay = summarize.parse_verify(root / "verify.log", 1)
            self.assertTrue(replay["valid"])
            self.assertEqual(replay["verified"], 300)
            self.assertEqual(replay["dropped"], 0)
            self.assertFalse(summarize.parse_verify(root / "verify.log", 0)["valid"])


if __name__ == "__main__":
    unittest.main()

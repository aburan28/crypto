import pathlib
import tempfile
import unittest

import summarize


def row(phase, arm, rate, pair=0, order=0):
    return {"phase": phase, "variant": arm, "rateMps": float(rate), "pair": pair, "order": order,
            "logSha256": "fixture", "gpuState": "fixture"}


class DecisionTests(unittest.TestCase):
    def screen(self, b16=100.0, b32=102.0, b64=101.0):
        return [
            row("warmup", "b16", b16), row("warmup", "b32", b32), row("warmup", "b64", b64),
            row("screen", "b16", b16, pair=1, order=1), row("screen", "b32", b32, pair=1, order=2),
            row("screen", "b64", b64, pair=1, order=3), row("screen", "b64", b64, pair=2, order=1),
            row("screen", "b32", b32, pair=2, order=2), row("screen", "b16", b16, pair=2, order=3),
        ]

    def test_top_candidate_and_gate(self):
        result = summarize.decide_samples(self.screen())
        self.assertTrue(result["valid"])
        self.assertEqual(result["topCandidate"], "b32")
        self.assertTrue(result["qualified"])
        self.assertEqual(result["decision"], "qualified_pending_confirmation")
        retained = summarize.decide_samples(self.screen(b32=100.4, b64=99.0))
        self.assertFalse(retained["qualified"])
        self.assertEqual(retained["decision"], "reference_retained")

    def test_dynamic_confirmation_order(self):
        rows = self.screen()
        rows += [
            row("confirm", "b16", 100, 1, 1), row("confirm", "b32", 102, 1, 2),
            row("confirm", "b32", 101, 2, 1), row("confirm", "b16", 100, 2, 2),
            row("confirm", "b16", 100, 3, 1), row("confirm", "b32", 103, 3, 2),
        ]
        result = summarize.decide_samples(rows)
        self.assertTrue(result["valid"])
        self.assertEqual(result["decision"], "confirmation_complete")
        self.assertEqual(result["pairedConfirmationRatios"], [1.02, 1.01, 1.03])

    def test_wrong_order_fails(self):
        rows = self.screen()
        rows += [row("confirm", "b32", 102, 1, 1)]
        self.assertFalse(summarize.decide_samples(rows)["valid"])


class GeometryParserTests(unittest.TestCase):
    def test_verify_geometry_and_markers(self):
        text = """packed table split forward: 1
packed table batch hints: 1
backend cuda-packed131: 48128 threads x 32 slots x 1 lanes = 1540096 walks, dp weight 48, 96 steps per launch
  finished: 1.0 M it/s, 300 distinguished points (300 verified against the reference, 0 dropped)
"""
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "verify.log"
            path.write_text(text)
            parsed = summarize.parse_verify(path, "b32")
            self.assertTrue(parsed["valid"])
            self.assertEqual(parsed["liveSlots"], 1540096)
            self.assertFalse(summarize.parse_verify(path, "b16")["valid"])


if __name__ == "__main__":
    unittest.main()

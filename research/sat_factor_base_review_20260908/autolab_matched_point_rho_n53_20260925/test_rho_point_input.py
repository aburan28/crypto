#!/usr/bin/env python3
"""Check point-only parity with a tiny fixture and exact-n53 rejection controls."""

import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
EXE = Path(os.environ.get("RHO_TEST_EXE",
                          REPO / "target/release/examples/koblitz_rho_batch_ks"))


class PointInputTest(unittest.TestCase):
    def invoke(self, n, point_path=None):
        env = {key: value for key, value in os.environ.items() if not key.startswith("KIC_")}
        env.update({"KIC_RHO_DP_BITS": "4", "KIC_RHO_PRECOMPUTE_WALKS": "0"})
        if point_path is not None:
            env["KIC_RHO_TARGET_POINTS_JSONL"] = str(point_path)
        return subprocess.run([str(EXE), str(n), "0", "signed_frobenius", "1", "42"],
                              cwd=REPO, env=env, capture_output=True, text=True)

    def test_point_only_matches_derived_tiny_fixture(self):
        reference = self.invoke(13)
        self.assertEqual(reference.returncode, 0, reference.stderr)
        expected = json.loads(reference.stdout.splitlines()[0])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "point.jsonl"
            path.write_text(json.dumps(expected["published_q"]) + "\n")
            observed = self.invoke(13, path)
        self.assertEqual(observed.returncode, 0, observed.stderr)
        row = json.loads(observed.stdout.splitlines()[0])
        self.assertEqual(row["published_q"], expected["published_q"])
        self.assertIsNone(row["published_fixture_scalar"])
        self.assertEqual(row["recovered_fixture_scalar"],
                         expected["published_fixture_scalar"])

    def test_n53_rejects_invalid_public_points(self):
        cases = (
            ([0, 0], "rho point target must be on the curve"),
            ([0, 1], "rho point target must belong to the prime-order subgroup"),
            ([1 << 53, 1], "rho point coordinates must be field elements"),
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "point.jsonl"
            for point, error in cases:
                with self.subTest(point=point):
                    path.write_text(json.dumps(point) + "\n")
                    observed = self.invoke(53, path)
                    self.assertNotEqual(observed.returncode, 0)
                    self.assertIn(error, observed.stderr)


if __name__ == "__main__":
    unittest.main()

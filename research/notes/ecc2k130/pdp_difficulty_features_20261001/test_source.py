"""Small closed-form controls; no archived target feature outcomes are read."""
from __future__ import annotations

import tempfile
from pathlib import Path
import unittest

import produce
import verify


class SourceControls(unittest.TestCase):
    def test_independent_square_laws_and_frobenius_invariance(self):
        # GF(2^3) / (z^3+z+1), an irreducible toy unrelated to held-out Q.
        for x in range(8):
            for y in range(8):
                a = produce.features(x, y, 3, [0, 1])
                b = verify.point_features(x, y, 3, [0, 1])
                self.assertEqual(a, b)
                self.assertEqual(a["signed_y_weight"],
                                 produce.features(x, x ^ y, 3, [0, 1])["signed_y_weight"])
                squared = produce.square(x, 3, [0, 1])
                self.assertEqual(a["frobenius_x_min_weight"],
                                 produce.features(squared, y, 3, [0, 1])["frobenius_x_min_weight"])

    def test_average_ties_and_rank_directions(self):
        self.assertEqual(produce.tied_ranks([1, 1, 3, 4]), [1.5, 1.5, 3.0, 4.0])
        self.assertEqual(verify.average_ranks([1, 1, 3, 4]), [1.5, 1.5, 3.0, 4.0])
        self.assertAlmostEqual(produce.spearman([1, 2, 3, 4], [10, 20, 30, 40]), 1.0)
        self.assertAlmostEqual(verify.spearman([1, 2, 3, 4], [40, 30, 20, 10]), -1.0)
        self.assertEqual(produce.spearman([1, 1, 1, 1], [1, 2, 3, 4]), 0.0)

    def test_independent_statistical_replay_on_toy_rows(self):
        rows = [{"x": i, "y": 0, "probes": (i * 7) % 13 + 1,
                 "features": {"x_weight": i % 4}} for i in range(12)]
        choice = {"feature": "x_weight", "sign": 1}
        spec = {"permutation_seed": 11}
        saved, scores = produce.evaluate(rows, choice, spec, 20)
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            produce.write_json(root / "toy_permutations.json", scores)
            verify.check_evaluation(rows, choice, spec, {"permutations": 20},
                                    root, "toy", saved)

    def test_permutation_seed_and_quartile_tie(self):
        rows = [
            {"x": x, "y": 0, "probes": p, "features": {"x_weight": f}}
            for x, f, p in ((4, 0, 1), (1, 0, 3), (2, 1, 5), (3, 2, 7))
        ]
        choice = {"feature": "x_weight", "sign": 1}
        spec = {"permutation_seed": 2026100141}
        summary, scores = produce.evaluate(rows, choice, spec, 20)
        self.assertEqual(summary["selected"], 1)
        self.assertEqual(summary["selected_median_probes"], 3)  # (x,y) tie-breaks the two zeroes.
        self.assertEqual(summary["all_median_probes"], 4.0)
        again, again_scores = produce.evaluate(rows, choice, spec, 20)
        self.assertEqual(summary, again)
        self.assertEqual(scores, again_scores)
        self.assertEqual(len(scores), 20)
        self.assertGreaterEqual(summary["permutation_p"], 1 / 21)


if __name__ == "__main__":
    unittest.main()

import copy
import random
import unittest

from boolean_closure import compute, equal_spans, normalize, solutions, verify


class ClosureTests(unittest.TestCase):
    def check_exact(self, n, generators):
        candidate = compute(n, generators)
        reference = compute(n, generators, strategy="exhaustive")
        self.assertEqual(verify(n, generators, candidate)["verified"], "exact_ideal")
        self.assertEqual(verify(n, generators, reference)["verified"], "exact_ideal")
        self.assertTrue(equal_spans(candidate, reference))
        points = solutions(n, generators)
        self.assertEqual(candidate["stats"]["rank"], (1 << n) - len(points))
        self.assertEqual(solutions(n, candidate["certificate"]["rows"]), points)
        return candidate

    def test_cancellation_and_idempotence(self):
        self.assertEqual(normalize(2, [[1, 1, 2]]), [[2]])
        result = self.check_exact(2, [[1]])
        self.assertEqual(result["stats"]["rank"], 2)

    def test_zero_coordinates_survive(self):
        self.check_exact(2, [[3]])  # x*y = 0
        self.assertEqual(solutions(2, [[3]]), [0, 1, 2])

    def test_high_degree_is_retained(self):
        result = self.check_exact(3, [[1, 0]])  # x + 1
        self.assertEqual(result["stats"]["max_scheduled_degree"], 3)
        self.assertTrue(any(7 in row for row in result["certificate"]["rows"]))

    def test_missing_closure_is_rejected(self):
        result = {"complete": True, "certificate": {
            "n": 2, "rows": [[1]], "origins": [["input", 0]]}}
        with self.assertRaisesRegex(ValueError, "variable closure"):
            verify(2, [[1]], result)

    def test_forged_provenance_is_rejected(self):
        result = copy.deepcopy(compute(2, [[1]]))
        result["certificate"]["rows"][0] = [0]
        with self.assertRaisesRegex(ValueError, "provenance"):
            verify(2, [[1]], result)

    def test_omitted_input_is_rejected(self):
        result = compute(2, [[1]])
        with self.assertRaisesRegex(ValueError, "generator containment"):
            verify(2, [[1], [0]], result)

    def test_budget_means_incomplete(self):
        for limits in ({"max_submissions": 1}, {"max_pending_terms": 1},
                       {"max_row_terms": 1}, {"max_total_terms": 1}):
            with self.subTest(limits=limits):
                result = compute(3, [[1, 2], [4]], **limits)
                self.assertFalse(result["complete"])
                with self.assertRaisesRegex(ValueError, "incomplete"):
                    verify(3, [[1, 2], [4]], result)
                verify(3, [[1, 2], [4]], result, require_complete=False)

    def test_edge_cases(self):
        for generators in ([], [[]], [[0]], [[1], [1, 0]], [[1], [1]]):
            with self.subTest(generators=generators):
                self.check_exact(3, generators)

    def test_seeded_controls(self):
        rng = random.Random(20260928)
        for n in range(1, 7):
            for sample in range(10):
                generators = [rng.sample(range(1 << n), min(5, 1 << n))
                              for _ in range(rng.randrange(1, 5))]
                with self.subTest(n=n, sample=sample):
                    self.check_exact(n, generators)


if __name__ == "__main__":
    unittest.main(verbosity=2)

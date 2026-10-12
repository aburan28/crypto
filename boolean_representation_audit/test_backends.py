import copy
import json
from pathlib import Path
import random
import unittest

from backends import compute, semantic_result, solutions, structure, verify


class BackendTests(unittest.TestCase):
    def check_system(self, n, generators):
        zeros = solutions(n, generators)
        for strategy in ('frontier', 'exhaustive'):
            sparse = compute(n, generators, strategy=strategy)
            packed = compute(n, generators, backend='packed', strategy=strategy)
            self.assertEqual(semantic_result(sparse), semantic_result(packed))
            self.assertEqual(verify(n, generators, packed)['verified'], 'exact_ideal')
            self.assertEqual(packed['stats']['rank'], (1 << n) - len(zeros))
            self.assertEqual(solutions(n, packed['certificate']['rows']), zeros)

    def test_frozen_cases(self):
        for case in json.loads(Path(__file__).with_name('inputs.json').read_text()):
            with self.subTest(case=case['name']):
                self.check_system(case['n'], case['generators'])

    def test_random_controls(self):
        rng = random.Random(91220260929)
        for n in range(1, 7):
            for sample in range(10):
                rows = [rng.sample(range(1 << n), min(5, 1 << n))
                        for _ in range(rng.randrange(1, 5))]
                with self.subTest(n=n, sample=sample):
                    self.check_system(n, rows)

    def test_edge_cases_and_scope(self):
        for rows in ([], [[]], [[0]], [[1, 1]], [[3]], [[1], [1, 0]]):
            self.check_system(3, rows)
        for n in (0, 11, True):
            with self.assertRaises(ValueError):
                compute(n, [[1]], backend='packed')
        with self.assertRaises(ValueError):
            compute(3, [[8]], backend='packed')

    def test_budgets_preserve_partial_semantics(self):
        for strategy in ('frontier', 'exhaustive'):
            for limits in ({'max_submissions': 1}, {'max_row_terms': 1},
                           {'max_total_terms': 1}, {'max_pending_terms': 1}):
                a = compute(3, [[1, 2], [4]], strategy=strategy, **limits)
                b = compute(3, [[1, 2], [4]], backend='packed', strategy=strategy, **limits)
                self.assertEqual(semantic_result(a), semantic_result(b))
                if not b['complete']:
                    with self.assertRaisesRegex(ValueError, 'incomplete'):
                        verify(3, [[1, 2], [4]], b)
                    verify(3, [[1, 2], [4]], b, require_complete=False)

    def test_verifier_rejects_mutation(self):
        result = copy.deepcopy(compute(3, [[1]], backend='packed'))
        result['certificate']['rows'][0] = [0]
        with self.assertRaisesRegex(ValueError, 'provenance'):
            verify(3, [[1]], result)

    def test_structural_controls(self):
        self.assertEqual(structure(4, [[1, 2], [2, 4], [4, 8]])['min_fill_width_upper_bound'], 1)
        self.assertEqual(structure(4, [[1, 2], [2, 4], [4, 8], [8, 1]])['min_fill_width_upper_bound'], 2)
        self.assertEqual(structure(4, [[1, 2, 4, 8]])['min_fill_width_upper_bound'], 3)
        self.assertEqual(structure(4, [[1], [2]])['min_fill_width_upper_bound'], 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)

"""Independent low-width controls for the production-path mathematical auditor."""
import itertools
import unittest

from f5_boolean_control import evaluate
from f5_production_control import decisive_reference, substitute


class ProductionControlTests(unittest.TestCase):
    def test_specialization_matches_truth_table_with_cancellation(self):
        rows = [[0, 1, 2, 3, 7], [2, 3], [0, 4, 5], []]
        for variable, value in itertools.product(range(3), (False, True)):
            specialized = substitute(rows, variable, value)
            for assignment in range(8):
                pinned = ((assignment | (1 << variable)) if value
                          else (assignment & ~(1 << variable)))
                self.assertEqual([evaluate(row, pinned) for row in rows],
                                 [evaluate(row, assignment) for row in specialized])

    def test_forces_require_linear_combinations(self):
        _, _, _, contradiction, forced, belongs = decisive_reference([[1, 2], [0, 2]], 3, 3)
        self.assertFalse(contradiction)
        self.assertEqual(sorted(forced), [[0, 1], [0, 2]])
        self.assertTrue(belongs([0, 1]))
        self.assertFalse(belongs([4]))

    def test_nonconstant_consequence_does_not_become_a_force(self):
        _, _, _, contradiction, forced, _ = decisive_reference([[1, 2]], 3, 3)
        self.assertFalse(contradiction)
        self.assertEqual(forced, [])

    def test_constant_refutation_is_detected(self):
        _, _, _, contradiction, _, _ = decisive_reference([[0]], 3, 3)
        self.assertTrue(contradiction)

    def test_nonlinear_decisive_rows_are_sound_on_every_exact_model(self):
        for rows in ([[3, 2], [0, 1]], [[3, 1, 2]], [[7, 1], [0, 2]]):
            models = [a for a in range(8) if all(evaluate(row, a) == 0 for row in rows)]
            self.assertTrue(models)
            _, _, _, contradiction, forced, _ = decisive_reference(rows, 3, 3)
            self.assertFalse(contradiction)
            self.assertTrue(all(evaluate(row, a) == 0 for row in forced for a in models))


if __name__ == '__main__':
    unittest.main()

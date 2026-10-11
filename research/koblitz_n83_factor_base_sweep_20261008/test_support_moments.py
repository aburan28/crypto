"""Small-group controls for the uniform-target support-count statement."""

import itertools
import unittest
from collections import Counter
from fractions import Fraction

from support_moments import support_count, vacuous_ceiling_column_threshold


class SupportMomentTests(unittest.TestCase):
    def test_multisets_and_nonidentity_target_bound(self):
        # Signed base in Z/7Z; the pair 1,-1 creates identity-sum supports.
        order = 7
        base = (1, 2, 6)
        for summands in (2, 3, 4):
            supports = list(itertools.combinations_with_replacement(base, summands))
            self.assertEqual(len(supports), support_count(len(base), summands))
            image = Counter(sum(support) % order for support in supports)
            identity = image[0]
            expected_nonidentity = Fraction(sum(image[q] for q in range(1, order)), order - 1)
            self.assertEqual(expected_nonidentity, Fraction(len(supports) - identity, order - 1))
            actual_hit_fraction = Fraction(sum(image[q] > 0 for q in range(1, order)), order - 1)
            ceiling = min(Fraction(1), Fraction(len(supports), order - 1))
            self.assertLessEqual(actual_hit_fraction, ceiling)
            self.assertEqual(Fraction(sum(image.values()), order), Fraction(len(supports), order))

    def test_first_vacuous_column_is_minimal(self):
        for order in (7, 8569786107849059, 2417851639230796216685689):
            for summands in range(2, 7):
                columns = vacuous_ceiling_column_threshold(order, summands)
                self.assertGreaterEqual(support_count(166 * columns, summands), order - 1)
                if columns > 1:
                    self.assertLess(support_count(166 * (columns - 1), summands), order - 1)


if __name__ == "__main__":
    unittest.main()

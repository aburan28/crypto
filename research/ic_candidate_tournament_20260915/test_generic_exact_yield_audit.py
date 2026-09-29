"""Pre-analysis controls for the exact three-summand group oracle."""
import unittest

from run_generic_exact_yield_audit import (
    PANEL, PANEL_SHA256, digest_bytes, exact_three_sum, pair_index,
    validate_registration,
)
from tournament import read


class AdditiveGroup:
    """Tiny group with the oracle's None-as-identity point convention."""
    @staticmethod
    def add(a, b):
        value = ((a or 0) + (b or 0)) % 11
        return value or None

    @staticmethod
    def neg(a):
        return None if a is None else (-a) % 11 or None


class ExactYieldAuditTests(unittest.TestCase):
    def test_registered_archive_and_scope_are_fixed(self):
        self.assertEqual(digest_bytes(PANEL.read_bytes()), PANEL_SHA256)
        panel = read(PANEL)
        validate_registration(panel)
        self.assertEqual(panel['summands'], 3)
        self.assertEqual(len(panel['cells']) * len(panel['solvers']), 20)

    def test_pair_index_covers_repeated_and_distinct_triples(self):
        group = AdditiveGroup()
        base = (1, 2)
        pairs = pair_index(group, base)
        for target in (3, 4, 5, 6):
            witness = exact_three_sum(group, base, pairs, target)
            self.assertIsNotNone(witness)
            self.assertEqual(len(witness), 3)
            self.assertEqual(sum(base[i] for i in witness) % 11, target)
        for target in (1, 2, 7, 8, 9, 10, None):
            self.assertIsNone(exact_three_sum(group, base, pairs, target))


if __name__ == '__main__':
    unittest.main()

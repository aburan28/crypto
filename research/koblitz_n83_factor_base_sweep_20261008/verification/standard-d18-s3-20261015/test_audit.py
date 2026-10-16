"""Small independent checks for the d18 archive audit's sign semantics."""

import unittest

import audit


class StandardSubspaceAuditTests(unittest.TestCase):
    def test_sign_pair_uses_binary_curve_negation(self):
        audit.check_sign_pair([(7, 11, False), (7, 11 ^ 7, True)], 0)

    def test_broken_labels_and_coordinates_are_rejected(self):
        for rows in (
            [(7, 11, False)],
            [(7, 11, False), (7, 11 ^ 7, False)],
            [(7, 11, False), (8, 11 ^ 7, True)],
            [(7, 11, False), (7, 11 ^ 6, True)],
        ):
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                audit.check_sign_pair(rows, 3)

    def test_canonical_field_coordinates(self):
        self.assertTrue(audit.canonical_decimal("0"))
        self.assertTrue(audit.canonical_decimal(str((1 << 83) - 1)))
        for value in ("00", "01", "-1", str(1 << 83), 1, "1.0"):
            with self.subTest(value=value):
                self.assertFalse(audit.canonical_decimal(value))


if __name__ == "__main__":
    unittest.main()

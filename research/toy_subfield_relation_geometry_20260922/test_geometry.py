import random
import unittest

from run import ToyCurve, ToyField, basis, inspect_base, poly_mod, span


def independent_mul(a, b, modulus):
    # Form the unreduced polynomial product, then divide by the modulus.
    out = 0
    for i in range(b.bit_length()):
        if b & (1 << i):
            out ^= a << i
    return poly_mod(out, modulus)


class GeometryTests(unittest.TestCase):
    def test_field_arithmetic_against_polynomial_product(self):
        rng = random.Random(981)
        for n in (6, 9, 10, 12):
            f = ToyField(n)
            for _ in range(256):
                a, b = rng.randrange(f.size), rng.randrange(f.size)
                self.assertEqual(f.mul(a, b), independent_mul(a, b, f.modulus))
            for a in range(1, f.size):
                self.assertEqual(f.mul(a, f.inv(a)), 1)

    def test_point_enumeration_against_all_coordinate_pairs(self):
        f = ToyField(6)
        for b in (1, 2, 13):
            c = ToyCurve(f, b)
            exact = {(x, y) for x in range(f.size) for y in range(f.size)
                     if independent_mul(y, y, f.modulus) ^ independent_mul(x, y, f.modulus)
                     == independent_mul(independent_mul(x, x, f.modulus), x, f.modulus) ^ b}
            self.assertEqual(c.affine, exact)

    def test_group_identities_and_associativity(self):
        rng = random.Random(438)
        for n in (6, 9, 10, 12):
            c = ToyCurve(ToyField(n), 2)
            points = [None] + sorted(c.affine)
            for _ in range(128):
                p, q, r = [rng.choice(points) for _ in range(3)]
                self.assertEqual(c.add(p, c.neg(p)), None)
                self.assertEqual(c.add(p, None), p)
                self.assertEqual(c.add(p, q), c.add(q, p))
                self.assertEqual(c.add(c.add(p, q), r), c.add(p, c.add(q, r)))

    def test_subfield_counts_and_product_span(self):
        for n, k in ((6, 2), (9, 3), (10, 2), (12, 4)):
            f = ToyField(n)
            fixed = [x for x in range(f.size) if f.frobenius(x, k) == x]
            u = basis(fixed)
            self.assertEqual(len(u), k)
            self.assertEqual(span(u), fixed)
            self.assertEqual(len(basis(f.mul(x, y) for x in u for y in u)), k)

    def test_subfield_curve_closure_and_empty_pair_workload(self):
        f = ToyField(6)
        c = ToyCurve(f, 1)
        u = basis(x for x in range(f.size) if f.frobenius(x, 2) == x)
        r = inspect_base(c, u, 2)
        self.assertTrue(r['all_base_points_in_subfield'])
        self.assertTrue(r['all_pair_sums_in_subfield'])
        empty = inspect_base(c, [], 2)
        self.assertEqual(empty['unordered_distinct_pairs'], 0)
        self.assertEqual(empty['distinct_nonidentity_targets'], 0)
        self.assertIsNone(empty['target_count_over_upper_bound'])

    def test_deduplicates_targets_and_excludes_infinity(self):
        c = ToyCurve(ToyField(6), 2)
        r = inspect_base(c, [1 << i for i in range(6)], 2)
        self.assertGreater(r['nonidentity_pairs_verified'], r['distinct_nonidentity_targets'])
        self.assertGreater(r['infinity_pairs_excluded'], 0)
        self.assertLessEqual(r['distinct_nonidentity_targets'], c.order - 1)

    def test_rejects_other_field_sizes(self):
        for n in (-1, 0, 5, 13, 131, 155, 185):
            with self.assertRaises(ValueError):
                ToyField(n)


if __name__ == '__main__':
    unittest.main()

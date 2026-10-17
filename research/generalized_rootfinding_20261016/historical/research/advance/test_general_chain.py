"""Independent full-coordinate tiny-field oracles for the general S3 join."""

from functools import lru_cache
import itertools
import unittest
from unittest.mock import patch

from rootfinder import polynomial as poly
from rootfinder.binaryfield import BinaryField
from rootfinder.fields import PrimeField, ExtensionField, make_field
from research.advance.general_chain import GeneralS3Fiber, solve_general_chain, _Unresolved


def independent_s3(field, x, y, z, curve):
    """Expanded symmetric polynomial, independent of the fiber formulas."""
    q = (x, y, z)
    square_products = [field.mul(field.square(q[i]), field.square(q[j]))
                       for i in range(3) for j in range(i)]
    total = 0
    for term in square_products:
        total = field.add(total, term)
    xyz = field.mul(field.mul(x, y), z)
    if field.characteristic == 2:
        return field.add(field.add(total, xyz), curve.get('b', 1))
    a2, a4, a6 = (curve.get(name, default) for name, default in
                  [('a2', 0), ('a4', 0), ('a6', 1)])
    sum_xyz = field.add(field.add(x, y), z)
    pair_products = 0
    for i in range(3):
        for j in range(i):
            pair_products = field.add(pair_products, field.mul(q[i], q[j]))
    two, four = 2 % field.characteristic, 4 % field.characteristic
    total = field.sub(total, field.mul(two, field.mul(xyz, sum_xyz)))
    total = field.sub(total, field.mul(two, field.mul(a4, pair_products)))
    total = field.sub(total, field.mul(four, field.mul(a2, xyz)))
    total = field.sub(total, field.mul(four, field.mul(a6, field.add(sum_xyz, a2))))
    return field.add(total, field.square(a4))


def zero_triples(field, curve):
    return {triple for triple in itertools.product(range(field.order), repeat=3)
            if independent_s3(field, *triple, curve) == 0}


def full_oracle(field, domain, target, curve):
    """Explicitly enumerate D^5 times K^3, without root solving or joining."""
    zeros = zero_triples(field, curve)
    answers = set()
    for x1, x2, x3, x4, x5 in itertools.product(domain, repeat=5):
        for u1, u2, u3 in itertools.product(range(field.order), repeat=3):
            if ((x1, x2, u1) in zeros and (u1, x3, u2) in zeros
                    and (u2, x4, u3) in zeros and (u3, x5, target) in zeros):
                answers.add((x1, x2, x3, x4, x5, u1, u2, u3))
    return answers


class GeneralChainTests(unittest.TestCase):
    def check_oracle(self, field, domain, target, curve, **options):
        expected = full_oracle(field, domain, target, curve)
        result = solve_general_chain(field, domain, target, curve=curve, **options)
        self.assertTrue(result['search_complete'], result)
        self.assertEqual(result['solution_count_exact'], len(expected))
        actual = {tuple(row) for row in result['solutions']}
        self.assertEqual(len(actual), len(result['solutions']))
        self.assertLessEqual(actual, expected)
        if result['output_complete']:
            self.assertEqual(actual, expected)
        planned, stats = result['stats']['planned'], result['stats']
        for name in ('root_calls', 'right_records', 'source_assignments'):
            self.assertLessEqual(stats[name], planned[name])
        return result

    def test_full_chains_prime_fields_including_characteristic_three(self):
        for field in (PrimeField(3), PrimeField(5), PrimeField(7)):
            curve = {'a2': 0, 'a4': 1, 'a6': 1}
            for target in range(min(field.order, 3)):
                with self.subTest(field=repr(field), target=target):
                    self.check_oracle(field, (0, 1), target, curve)

    def test_full_chains_odd_extensions_and_arbitrary_coefficients(self):
        field = make_field({'kind': 'extension', 'characteristic': 3, 'degree': 2,
                            'modulus': [1, 0, 1]})
        curve = {'a2': 3, 'a4': 1, 'a6': 4}
        for target in (0, 5):
            self.check_oracle(field, (1, 3), target, curve)
        field = ExtensionField(5, 2, (2, 0, 1))
        self.check_oracle(field, (6,), 7, {'a2': 6, 'a4': 1, 'a6': 1})

    def test_binary_fields_and_generic_characteristic_two_backends(self):
        for field in (BinaryField(1), BinaryField(2), PrimeField(2),
                      ExtensionField(2, 2, (1, 1, 1))):
            for target in range(min(2, field.order)):
                curve = {'b': 1}
                self.check_oracle(field, (0, 1), target, curve)

    def test_all_quadratic_fibers_match_direct_field_enumeration(self):
        for field, curve in [(PrimeField(3), {'a2': 1, 'a4': 0, 'a6': 1}),
                             (PrimeField(5), {'a2': 2, 'a4': 1, 'a6': 1}),
                             (ExtensionField(3, 2, (1, 0, 1)), {'a2': 3, 'a4': 1, 'a6': 4}),
                             (BinaryField(3), {'b': 3})]:
            fiber = GeneralS3Fiber(field, curve)
            for x, y in itertools.product(range(field.order), repeat=2):
                expected = tuple(z for z in range(field.order)
                                 if independent_s3(field, x, y, z, curve) == 0)
                with self.subTest(field=repr(field), x=x, y=y):
                    self.assertEqual(fiber.roots(x, y), expected)
                    c, b, a = fiber.coefficients(x, y)
                    for z in range(field.order):
                        specialized = field.add(field.add(c, field.mul(b, z)),
                                                field.mul(a, field.square(z)))
                        self.assertEqual(specialized, independent_s3(field, x, y, z, curve))

    def test_full_width_binary_coefficients_have_frobenius_count_certificate(self):
        field = BinaryField(83, 0xd84fae489bcdf275f7b61)
        curve = {'b': (1 << 79) + 1, 'a2': 1 << 71}
        fiber = GeneralS3Fiber(field, curve)
        for x, y in [(0, 1), ((1 << 65) + 7, (1 << 77) + 11),
                     ((1 << 70) + 1, (1 << 70) + 1)]:
            roots = fiber.roots(x, y)
            polynomial = poly.trim(fiber.coefficients(x, y), field)
            frobenius = poly.powmod_poly((0, 1), field.order, polynomial, field)
            restricted = poly.gcd_poly(polynomial, poly.sub(frobenius, (0, 1), field), field)
            self.assertEqual(len(roots), len(restricted) - 1)
            self.assertTrue(all(independent_s3(field, x, y, root, curve) == 0 for root in roots))

    def test_equal_coordinates_cubic_root_have_nonzero_constant_fiber(self):
        # C(X)=X^3+X has the root0, derivative1 there, hence S3(0,0,Z)=1.
        for field in (PrimeField(3), PrimeField(5), ExtensionField(3, 2, (1, 0, 1))):
            fiber = GeneralS3Fiber(field, {'a2': 0, 'a4': 1, 'a6': 0})
            self.assertEqual(fiber.coefficients(0, 0), (1, 0, 0))
            self.assertEqual(fiber.roots(0, 0), ())

    def test_singular_curves_rejected_even_when_no_sources(self):
        for field, curve in [(PrimeField(3), {'a2': 0, 'a4': 0, 'a6': 1}),
                             (PrimeField(5), {'a2': 0, 'a4': 0, 'a6': 0}),
                             (ExtensionField(3, 2, (1, 0, 1)), {'a2': 0, 'a4': 0, 'a6': 4}),
                             (BinaryField(3), {'b': 0})]:
            with self.subTest(field=repr(field)), self.assertRaises(ValueError):
                solve_general_chain(field, (), 0, curve=curve)

    def test_output_caps_preserve_all_duplicate_bucket_multiplicities(self):
        field, curve = BinaryField(2), {'b': 1}
        full = self.check_oracle(field, (0, 1), 0, curve)
        self.assertGreater(full['solution_count_exact'], 1)
        for cap in (0, 1, 3):
            result = self.check_oracle(field, (0, 1), 0, curve, max_output=cap)
            self.assertEqual(result['solution_count_exact'], full['solution_count_exact'])
            self.assertEqual(len(result['solutions']), min(cap, full['solution_count_exact']))

    def test_preflight_refusal_has_no_root_calls_or_exact_count(self):
        field, curve = PrimeField(5), {'a2': 0, 'a4': 1, 'a6': 1}
        for options in ({'max_domain': 1}, {'max_source_assignments': 0},
                        {'max_root_calls': 0}, {'max_records': 0}):
            with self.subTest(options=options):
                result = solve_general_chain(field, (0, 1), 0, curve=curve, **options)
                self.assertEqual(result['status'], 'unresolved')
                self.assertEqual(result['stats']['root_calls'], 0)
                self.assertIsNone(result['solution_count_exact'])

    def test_exhausted_quadratic_split_is_unresolved_not_empty(self):
        field, curve = PrimeField(5), {'a2': 0, 'a4': 1, 'a6': 1}
        fiber = GeneralS3Fiber(field, curve, max_split_trials=0)
        pair = next((x, y) for x, y in itertools.product(range(5), repeat=2)
                    if len([z for z in range(1, 5)
                            if independent_s3(field, x, y, z, curve) == 0]) == 2)
        with self.assertRaises(_Unresolved):
            fiber.roots(*pair)
        # Make the unresolved pair the first backwards fiber (x5,target).
        result = solve_general_chain(field, (pair[0],), pair[1], curve=curve,
                                     max_split_trials=0)
        self.assertEqual(result['status'], 'unresolved')
        self.assertIsNone(result['solution_count_exact'])
        self.assertIn('splitting', result['unresolved_reason'])
        self.assertEqual(result['solutions'], [])

    def test_zero_time_budget_is_explicitly_unresolved(self):
        result = solve_general_chain(BinaryField(2), (0, 1), 0, max_seconds=0)
        self.assertEqual(result['status'], 'unresolved')
        self.assertEqual(result['stats']['root_calls'], 0)
        self.assertIsNone(result['solution_count_exact'])

    def test_empty_domain_is_exactly_empty(self):
        result = solve_general_chain(PrimeField(5), (), 0,
                                     curve={'a2': 0, 'a4': 1, 'a6': 1},
                                     max_domain=0, max_source_assignments=0,
                                     max_root_calls=0, max_records=0)
        self.assertEqual(result['status'], 'complete')
        self.assertEqual(result['solution_count_exact'], 0)
        self.assertEqual(result['stats']['root_calls'], 0)

    def test_validation_and_unknown_curve_parameters(self):
        field, curve = PrimeField(5), {'a2': 0, 'a4': 1, 'a6': 1}
        for domain, target, options in [((0, 0), 0, {}), ((5,), 0, {}),
                                         ((0,), 5, {}), ((0,), 0, {'max_output': -1}),
                                         ((0,), 0, {'max_seconds': float('nan')})]:
            with self.subTest(domain=domain, target=target, options=options), self.assertRaises(ValueError):
                solve_general_chain(field, domain, target, curve=curve, **options)
        with self.assertRaises(ValueError):
            GeneralS3Fiber(field, {'b': 1})
        with self.assertRaises(ValueError):
            GeneralS3Fiber(BinaryField(2), {'a4': 1})


if __name__ == '__main__':
    unittest.main()

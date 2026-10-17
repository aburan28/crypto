"""Small independent norm-division and exhaustive parameter-family oracles."""

import itertools
import unittest
from unittest.mock import patch

from rootfinder import polynomial as poly
from rootfinder.binaryfield import BinaryField
from rootfinder.fields import PrimeField, ExtensionField
from lab.research.odd_family import (prepare_model, source_sums, domain_contains,
    coefficient_maps, discriminant, parameter_domain_remainder, conditional_family,
    cubic_from_parameter, supported_cubic)


def domain_values(model):
    f, descriptor = model['field'], model['domain']
    if descriptor['kind'] == 'interval':
        return tuple(f.add(model['offset'], i) for i in range(descriptor['size']))
    result = []
    for coordinate in itertools.product(range(f.characteristic), repeat=len(model['basis'])):
        value = model['offset']
        for scalar, vector in zip(coordinate, model['basis']):
            value = f.add(value, f.mul(scalar, vector))
        result.append(value)
    return tuple(result)


def literal_norm_cubic(model, s, b):
    f = model['field']
    z, y = model['fixed_point']
    c = f.mul(f.sub(f.sub(f.square(b), z), s), f.inv(2 % f.characteristic))
    d = f.neg(f.add(f.square(z), f.add(f.mul(c, z), f.mul(b, y))))
    ap = (d, c, 1)
    norm = poly.mul(ap, ap, f)
    curve = (model['curve']['a6'], model['curve']['a4'], model['curve']['a2'], 1)
    norm = poly.sub(norm, tuple(f.mul(f.square(b), value) for value in curve), f)
    cubic, remainder = poly.divmod_poly(norm, (f.neg(z), 1), f)
    assert remainder == () and cubic[-1] == 1
    return cubic


def oracle_supported(cubic, domain, field):
    residue = cubic
    for root in domain:
        while len(residue) > 1 and poly.evaluate_poly(residue, root, field) == 0:
            residue, remainder = poly.divmod_poly(residue, (field.neg(root), 1), field)
            assert not remainder
    return residue == (1,)


def exhaustive_parameters(model, s):
    f, domain = model['field'], domain_values(model)
    return {b for b in range(f.order) if oracle_supported(literal_norm_cubic(model, s, b), domain, f)}


class OddFamilyTests(unittest.TestCase):
    def models(self):
        for f in (PrimeField(3), PrimeField(5), PrimeField(7)):
            yield prepare_model(f, {'a2': 0, 'a4': 1, 'a6': 1}, (0, 1),
                                interval_length=2)
        f = ExtensionField(3, 2, (1, 0, 1))
        yield prepare_model(f, {'a2': 0, 'a4': 1, 'a6': 1}, (0, 1), offset=3, basis=(1,))
        f = ExtensionField(3, 3)
        yield prepare_model(f, {'a2': 1, 'a4': 0, 'a6': 1}, (0, 1), basis=(1, 3))
        f = ExtensionField(5, 2, (2, 0, 1))
        yield prepare_model(f, {'a2': 0, 'a4': 1, 'a6': 1}, (0, 1), basis=(1,))

    def test_coefficients_match_literal_norm_division(self):
        for model in self.models():
            f = model['field']
            for s in source_sums(model):
                maps = coefficient_maps(model, s)
                for b in range(f.order):
                    wanted = literal_norm_cubic(model, s, b)
                    got = cubic_from_parameter(model, s, b)
                    self.assertEqual(got, wanted)
                    self.assertLessEqual(len(maps['e2'])-1, 4)
                    self.assertLessEqual(len(maps['e3'])-1, 4)

    def test_all_parameter_families_match_independent_complete_oracle(self):
        for model in self.models():
            for s in source_sums(model):
                with self.subTest(field=repr(model['field']), s=s):
                    result = conditional_family(model, s)
                    self.assertEqual(result['status'], 'complete', result)
                    expected = exhaustive_parameters(model, s)
                    self.assertEqual(set(result['parameters']), expected)
                    self.assertEqual(result['solution_count_exact'], len(expected))
                    self.assertTrue(result['search_complete'])
                    self.assertFalse(result['stats']['domain_materialized'])
                    self.assertFalse(result['stats']['field_materialized'])

    def test_compact_remainder_specializes_to_direct_domain_remainder(self):
        for model in self.models():
            f = model['field']
            domain_polynomial = (1,)
            for x in domain_values(model):
                domain_polynomial = poly.mul(domain_polynomial, (f.neg(x), 1), f)
            for s in source_sums(model):
                formal = parameter_domain_remainder(model, s)
                self.assertTrue(all(d <= 2*model['domain']['size']
                                    for d in formal['stats']['actual_parameter_degrees']))
                for b in range(f.order):
                    cubic = literal_norm_cubic(model, s, b)
                    expected = poly.divmod_poly(domain_polynomial, cubic, f)[1]
                    actual = poly.trim(tuple(poly.evaluate_poly(row, b, f)
                                             for row in formal['coefficients']))
                    self.assertEqual(actual, expected)

    def test_repeated_root_branch_is_required_for_singleton_domain(self):
        model = prepare_model(PrimeField(5), {'a2': 0, 'a4': 1, 'a6': 1}, (0, 1))
        result = conditional_family(model, 0)
        self.assertEqual(result['status'], 'complete')
        self.assertIn(0, result['parameters'])
        cubic = literal_norm_cubic(model, 0, 0)
        self.assertEqual(cubic, (0, 0, 0, 1))
        self.assertEqual(result['stats']['simple_parameter_polynomial_degree'], 0)
        self.assertEqual(result['stats']['repeated_parameter_polynomial_degree'], 12)
        self.assertTrue(supported_cubic(model, cubic))

    def test_two_torsion_fixed_point_is_supported(self):
        model = prepare_model(PrimeField(7), {'a2': 0, 'a4': 1, 'a6': 0}, (0, 0),
                              interval_length=3)
        for s in source_sums(model):
            result = conditional_family(model, s)
            self.assertEqual(result['status'], 'complete')
            self.assertEqual(set(result['parameters']), exhaustive_parameters(model, s))

    def test_discriminant_degree_is_exact_in_characteristic_three(self):
        for model in self.models():
            f = model['field']
            for s in source_sums(model):
                delta = discriminant(model, s)
                self.assertEqual(len(delta)-1, 12)
                self.assertEqual(delta[-1], f.neg(f.inv(16 % f.characteristic)))

    def test_exact_domain_membership_and_sum_cover_with_modular_wrap(self):
        for model in self.models():
            f, domain = model['field'], domain_values(model)
            self.assertEqual(len(set(domain)), model['domain']['size'])
            self.assertEqual({x for x in range(f.order) if domain_contains(model, x)}, set(domain))
            expected = {f.add(f.add(x, y), z) for x, y, z in itertools.product(domain, repeat=3)}
            self.assertEqual(set(source_sums(model)), expected)

    def test_output_limit_preserves_exact_family_count(self):
        model = prepare_model(PrimeField(5), {'a2': 0, 'a4': 1, 'a6': 1}, (0, 1))
        result = conditional_family(model, 0, max_parameters=0)
        self.assertEqual(result['status'], 'output_truncated')
        self.assertGreaterEqual(result['solution_count_exact'], 1)
        self.assertEqual(result['parameters'], [])
        self.assertTrue(result['search_complete'])

    def test_degree_and_construction_caps_do_not_claim_empty_family(self):
        model = prepare_model(ExtensionField(3, 3), {'a2': 1, 'a4': 0, 'a6': 1}, (0, 1), basis=(1, 3))
        for budget in ({'max_parameter_degree': 0}, {'max_scalar_products': 0}):
            result = conditional_family(model, 0, **budget)
            self.assertEqual(result['status'], 'unresolved')
            self.assertIsNone(result['solution_count_exact'])
            self.assertFalse(result['search_complete'])

    def test_splitting_cap_does_not_promote_partial_roots_to_exact_count(self):
        model = prepare_model(PrimeField(5), {'a2': 0, 'a4': 1, 'a6': 1}, (0, 1))
        partial = {'status': 'unresolved', 'roots': [], 'restricted_degree': 2, 'split_trials': 0}
        with patch('lab.research.odd_family.poly.field_roots', return_value=partial):
            result = conditional_family(model, 0)
        self.assertEqual(result['status'], 'unresolved')
        self.assertIsNone(result['solution_count_exact'])
        self.assertFalse(result['search_complete'])

    def test_model_validation(self):
        f = PrimeField(5)
        curve = {'a2': 0, 'a4': 1, 'a6': 1}
        for options in ({'basis': (1, 1)}, {'basis': (1,), 'interval_length': 2},
                        {'interval_length': 6}, {'interval_length': 0}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                prepare_model(f, curve, (0, 1), **options)
        with self.assertRaises(ValueError):
            prepare_model(f, curve, (0, 0))
        with self.assertRaises(ValueError):
            prepare_model(BinaryField(2), {'b': 1}, (0, 1))


if __name__ == '__main__':
    unittest.main()

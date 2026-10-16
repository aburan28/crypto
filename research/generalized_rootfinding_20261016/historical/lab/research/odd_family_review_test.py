"""Independent review oracles for the odd compact-domain parameter component."""

from itertools import combinations_with_replacement
import unittest

from rootfinder import polynomial as P
from rootfinder.fields import PrimeField, ExtensionField
from rootfinder.binaryfield import BinaryField
from lab.research import odd_family as component
from lab.research.odd_support_math_check import domain_poly, remainder_power_recurrence


def explicit_domain(model):
    f = model['field']
    return tuple(value for value in range(f.order) if component.domain_contains(model, value))


def oracle_parameters(model, source_sum):
    f = model['field']
    supported = set()
    for roots in combinations_with_replacement(explicit_domain(model), 3):
        h = domain_poly(roots, f)
        if f.neg(h[2]) == source_sum:
            supported.add(h)
    return {parameter for parameter in range(f.order)
            if component.cubic_from_parameter(model, source_sum, parameter) in supported}


class OddFamilyIndependentReview(unittest.TestCase):
    def test_compact_frobenius_matches_ordinary_power_recurrence(self):
        for field in (PrimeField(3), PrimeField(5), ExtensionField(3, 2), ExtensionField(5, 2)):
            for basis in ((), (1,)) + (((1, field.characteristic),) if field.degree == 2 else ()):
                offset = 0 if field.degree == 1 else field.characteristic
                model = component.prepare_model(field, {'a2': 0, 'a4': 1, 'a6': 1},
                                                (0, 1), offset=offset, basis=basis)
                domain = explicit_domain(model)
                self.assertEqual(len(domain), model['domain']['size'])
                polynomial = domain_poly(domain, field)
                for source_sum in list(component.source_sums(model))[:3]:
                    construction = component.parameter_domain_remainder(model, source_sum)
                    maps = component.coefficient_maps(model, source_sum)
                    expected = remainder_power_recurrence(polynomial, maps['e2'], maps['e3'],
                                                          source_sum, field)
                    self.assertEqual(construction['coefficients'], expected,
                                     (field.order, basis, source_sum))
                    for coordinate, remainder in enumerate(expected):
                        if remainder:
                            self.assertLessEqual(len(remainder)-1, 2*(len(domain)-coordinate))

    def test_families_match_all_domain_multisets_in_four_fields(self):
        for field in (PrimeField(3), PrimeField(5), ExtensionField(3, 2), ExtensionField(5, 2)):
            offset = 0 if field.degree == 1 else field.characteristic
            model = component.prepare_model(field, {'a2': 0, 'a4': 1, 'a6': 1},
                                            (0, 1), offset=offset, basis=(1,))
            for source_sum in component.source_sums(model):
                report = component.conditional_family(model, source_sum)
                self.assertEqual(report['status'], 'complete', report)
                self.assertEqual(set(report['parameters']), oracle_parameters(model, source_sum))

    def test_full_extension_field_and_two_torsion_target(self):
        models = [component.prepare_model(ExtensionField(3, 2),
                                           {'a2': 0, 'a4': 1, 'a6': 1}, (0, 1), basis=(1, 3)),
                  component.prepare_model(PrimeField(5),
                                           {'a2': 0, 'a4': 4, 'a6': 0}, (0, 0), basis=(1,))]
        nonempty = 0
        for model in models:
            for source_sum in component.source_sums(model):
                report = component.conditional_family(model, source_sum)
                self.assertEqual(report['status'], 'complete', report)
                expected = oracle_parameters(model, source_sum)
                self.assertEqual(set(report['parameters']), expected)
                nonempty += len(expected)
        self.assertGreater(nonempty, 0)

    def test_interval_domain_offset_and_wraparound(self):
        for field in (PrimeField(5), ExtensionField(5, 2)):
            offset = 4 if field.degree == 1 else 5
            for length in (1, 2, 3, 5):
                model = component.prepare_model(field, {'a2': 0, 'a4': 1, 'a6': 1},
                                                (0, 1), offset=offset, interval_length=length)
                domain = explicit_domain(model)
                self.assertEqual(len(domain), length)
                actual_sums = {field.add(field.add(x, y), z)
                               for x in domain for y in domain for z in domain}
                self.assertEqual(actual_sums, set(component.source_sums(model)))
                for source_sum in component.source_sums(model):
                    report = component.conditional_family(model, source_sum)
                    self.assertEqual(report['status'], 'complete', report)
                    self.assertEqual(set(report['parameters']), oracle_parameters(model, source_sum))

    def test_budget_exhaustion_and_output_caps_remain_explicit(self):
        field = PrimeField(5)
        model = component.prepare_model(field, {'a2': 0, 'a4': 1, 'a6': 1}, (0, 1), basis=(1,))
        low_degree = component.conditional_family(model, 0, max_parameter_degree=1)
        self.assertEqual(low_degree['status'], 'unresolved')
        self.assertFalse(low_degree['search_complete'])
        exhausted = component.conditional_family(model, 0, max_scalar_products=0)
        self.assertEqual(exhausted['status'], 'unresolved')
        self.assertIsNone(exhausted['solution_count_exact'])
        for source_sum in component.source_sums(model):
            exact = oracle_parameters(model, source_sum)
            if exact:
                clipped = component.conditional_family(model, source_sum, max_parameters=0)
                self.assertEqual(clipped['status'], 'output_truncated')
                self.assertEqual(clipped['solution_count_exact'], len(exact))
                self.assertTrue(clipped['search_complete'])
                self.assertFalse(clipped['output_complete'])
                self.assertEqual(clipped['parameters'], [])
                break
        else:
            self.fail('nonempty output-cap fixture required')

    def test_rejected_assumptions(self):
        with self.assertRaises(ValueError):
            component.prepare_model(BinaryField(3), {'b': 1}, (0, 1))
        with self.assertRaises(ValueError):
            component.prepare_model(PrimeField(5), {'a2': 0, 'a4': 1, 'a6': 1},
                                    (0, 1), basis=(1, 2))
        with self.assertRaises(ValueError):
            component.prepare_model(PrimeField(3), {'a2': 0, 'a4': 0, 'a6': 1}, (0, 1))
        with self.assertRaises(ValueError):
            component.prepare_model(PrimeField(5), {'a2': 0, 'a4': 1, 'a6': 1},
                                    (0, 1), interval_length=6)


if __name__ == '__main__':
    unittest.main()

"""Independent S4-resultant all-x oracles, including nonrational twist cases."""

from itertools import product
import unittest

from rootfinder.fields import PrimeField, ExtensionField
from lab.research.odd_solver import prepare_x_model, solve_all_x


def quadratic(field, curve, x, y):
    # Numeric original S3 coefficients, independent of the norm construction.
    two, four = 2 % field.characteristic, 4 % field.characteristic
    xy, total = field.mul(x, y), field.add(x, y)
    a = field.square(field.sub(x, y))
    b = field.neg(field.mul(two, field.add(field.mul(total, field.add(xy, curve['a4'])),
             field.add(field.mul(field.mul(two, curve['a2']), xy), field.mul(two, curve['a6'])))))
    c = field.sub(field.square(field.sub(curve['a4'], xy)),
                  field.mul(field.mul(four, curve['a6']), field.add(total, curve['a2'])))
    return a, b, c


def resultant_s4(field, curve, x, y, w, z):
    # Exact homogeneous resultant of the two degree-two S3 polynomials.
    a, b, c = quadratic(field, curve, x, y)
    d, e, f = quadratic(field, curve, w, z)
    af_cd = field.sub(field.mul(a, f), field.mul(c, d))
    ae_bd = field.sub(field.mul(a, e), field.mul(b, d))
    bf_ce = field.sub(field.mul(b, f), field.mul(c, e))
    return field.sub(field.square(af_cd), field.mul(ae_bd, bf_ce))


def oracle(field, curve, domain, z):
    return {row for row in product(domain, repeat=3) if resultant_s4(field, curve, *row, z) == 0}


class OddAllXTests(unittest.TestCase):
    def check(self, field, curve, domain, z, **model_options):
        model = prepare_x_model(field, curve, z, **model_options)
        result = solve_all_x(model)
        self.assertEqual(result['status'], 'complete', result)
        expected = oracle(field, curve, domain, z)
        self.assertEqual({tuple(row) for row in result['solutions']}, expected)
        self.assertEqual(result['solution_count_exact'], len(expected))
        self.assertEqual(result['target_x'], z)
        return result

    def test_all_tiny_targets_prime_and_odd_extension(self):
        curve = {'a2': 0, 'a4': 1, 'a6': 1}
        for field, length in [(PrimeField(3), 2), (PrimeField(5), 3), (PrimeField(7), 2),
                              (ExtensionField(3, 2, (1, 0, 1)), 2),
                              (ExtensionField(5, 2, (2, 0, 1)), 2)]:
            for z in range(min(field.order, 9)):
                with self.subTest(field=repr(field), z=z):
                    self.check(field, curve, range(length), z, interval_length=length)

    def test_affine_extension_offsets_and_nonbase_target(self):
        field, curve = ExtensionField(3, 2, (1, 0, 1)), {'a2': 0, 'a4': 1, 'a6': 1}
        for z in (0, 3, 5, 7):
            self.check(field, curve, (3, 4, 5), z, offset=3, basis=(1,))

    def test_nonrational_pair_closure_supplies_omitted_roots(self):
        field, curve = PrimeField(5), {'a2': 0, 'a4': 1, 'a6': 1}
        result = self.check(field, curve, (0, 1), 0, interval_length=2)
        # C(1)=3 is nonsquare. The two nonrational conjugate points cancel,
        # while the third source has target abscissa zero.
        self.assertIn([1, 1, 0], result['solutions'])
        self.assertGreater(result['stats']['mixed_domain_elements_scanned'], 0)

    def test_all_nonrational_triples_at_two_torsion_target(self):
        field, curve = PrimeField(7), {'a2': 0, 'a4': 1, 'a6': 0}
        result = self.check(field, curve, (2, 4, 6), 0, offset=2,
                            interval_length=3, interval_step=2)
        # Every abscissa in this source domain has nonsquare C(x).
        self.assertGreater(result['solution_count_exact'], 0)
        self.assertGreater(result['stats']['point_triples_replayed'], 0)

    def test_nonrational_target_is_twisted_with_exact_output_mapping(self):
        field, curve = PrimeField(5), {'a2': 0, 'a4': 1, 'a6': 1}
        model = prepare_x_model(field, curve, 1, interval_length=2)
        self.assertNotEqual(model['scale_to_original'], 1)
        result = solve_all_x(model)
        self.assertEqual({tuple(row) for row in result['solutions']}, oracle(field, curve, (0, 1), 1))
        self.assertEqual(result['target_x'], 1)

    def test_output_cap_keeps_exact_unique_x_count(self):
        model = prepare_x_model(PrimeField(5), {'a2': 0, 'a4': 1, 'a6': 1}, 0, interval_length=2)
        full = solve_all_x(model)
        clipped = solve_all_x(model, max_output=0)
        self.assertGreater(full['solution_count_exact'], 0)
        self.assertEqual(clipped['solution_count_exact'], full['solution_count_exact'])
        self.assertEqual(clipped['solutions'], [])
        self.assertEqual(clipped['status'], 'output_truncated')
        self.assertTrue(clipped['search_complete'])

    def test_global_limits_do_not_claim_complete_empty_results(self):
        model = prepare_x_model(PrimeField(5), {'a2': 0, 'a4': 1, 'a6': 1}, 0, interval_length=2)
        for limits in ({'max_sums': 0}, {'max_domain_scan': 0}, {'max_nonresidue_trials': 0},
                       {'max_unique_solutions': 0}, {'max_parameter_degree': 0}):
            with self.subTest(limits=limits):
                result = solve_all_x(model, **limits)
                self.assertEqual(result['status'], 'unresolved')
                self.assertIsNone(result['solution_count_exact'])
                self.assertFalse(result['search_complete'])


if __name__ == '__main__':
    unittest.main()

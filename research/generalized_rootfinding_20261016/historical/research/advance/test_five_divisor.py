"""Exact five-divisor algebra/profile checks; no performance experiments."""

import itertools
import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from rootfinder.binaryfield import BinaryField
from rootfinder import polynomial as poly
from research.advance.divisor_triple import enumerate_family, negate_point
from research.advance.five_divisor import (conditional_five_family, five_coefficients_formula,
                          five_coefficients_norm, prepare_five_model,
                          profile_instance, support_conditions)
from research.advance.test_divisor_triple import oracle_add, point_list


def interpolate_function(field, points):
    """Independent field Gaussian interpolation of x³+r*x²+s*x+t+(h*x+k)y."""
    rows = [[field.mul(x, y), y, field.mul(x, x), x, 1,
             field.mul(field.mul(x, x), x)] for x, y in points]
    pivots, row_index = [], 0
    for column in range(5):
        chosen = next((i for i in range(row_index, len(rows)) if rows[i][column]), None)
        if chosen is None:
            continue
        rows[row_index], rows[chosen] = rows[chosen], rows[row_index]
        inverse = field.inv(rows[row_index][column])
        rows[row_index] = [field.mul(value, inverse) for value in rows[row_index]]
        for index in range(len(rows)):
            if index != row_index and rows[index][column]:
                factor = rows[index][column]
                rows[index] = [a ^ field.mul(factor, b)
                               for a, b in zip(rows[index], rows[row_index])]
        pivots.append(column)
        row_index += 1
    if any(not any(row[:5]) and row[5] for row in rows) or len(pivots) != 5:
        return None
    answer = [0] * 5
    for row, column in zip(rows, pivots):
        answer[column] = row[5]
    return answer


class FiveDivisorTests(unittest.TestCase):
    def test_formula_vs_literal_norm_every_parameter_small_field(self):
        f = BinaryField(2)
        for a, b in ((0, 1), (1, 2)):
            for target in point_list(f, a, b)[::2]:
                model = prepare_five_model(f, a, b, target, 1, (1,))
                for h, k, r, s in itertools.product(range(f.order), repeat=4):
                    self.assertEqual(five_coefficients_norm(model, h, k, r, s),
                                     five_coefficients_formula(model, h, k, r, s))

    def test_reconstructed_five_point_divisors_against_independent_group_oracle(self):
        f, rng = BinaryField(4), random.Random(616103)
        verified = 0
        for a, b in ((0, 1), (1, 7), (3, 2)):
            points = point_list(f, a, b)
            local = 0
            for _ in range(200):
                sources = tuple(rng.sample(points, 5))
                total = None
                for source in sources:
                    total = oracle_add(f, a, total, source)
                if total is None or negate_point(total) in sources:
                    continue
                tpoint = negate_point(total)
                coefficients = interpolate_function(f, sources + (tpoint,))
                if coefficients is None:
                    continue
                h, k, r, s, t = coefficients
                model = prepare_five_model(f, a, b, total, 0, (1, 2, 4, 8))
                actual = five_coefficients_norm(model, h, k, r, s)
                expected = (1,)
                for x, _ in sources:
                    expected = poly.mul(expected, (x, 1), f)
                self.assertEqual(actual, expected)
                z, yt = tpoint
                self.assertEqual(t, f.mul(f.square(z), z) ^ f.mul(r, f.square(z))
                                 ^ f.mul(s, z) ^ f.mul(f.mul(h, z) ^ k, yt))
                self.assertTrue(all(support_conditions(model, actual)))
                verified += 1
                local += 1
                if local == 12:
                    break
            self.assertEqual(local, 12)
        self.assertEqual(verified, 36)

    def test_affine_conditional_families_are_exact_necessary_constraints(self):
        f = BinaryField(3)
        model = prepare_five_model(f, 1, 1, (0, 1), 3, (1, 4))
        for h, k in itertools.product(range(f.order), repeat=2):
            family = conditional_five_family(model, h, k)
            expected = {r | (s << f.degree) for r, s in itertools.product(range(f.order), repeat=2)
                        if all(support_conditions(model, five_coefficients_norm(model, h, k, r, s)))}
            self.assertEqual(set(enumerate_family(family)), expected)
            self.assertEqual(family["cardinality"], len(expected))

    def test_bilinear_products_prevent_joint_linearity(self):
        f = BinaryField(5)
        model = prepare_five_model(f, 0, 1, (0, 1), 0, (1,))
        base = five_coefficients_formula(model, 1, 0, 0, 0)
        k_only = five_coefficients_formula(model, 1, 1, 0, 0)
        r_only = five_coefficients_formula(model, 1, 0, 1, 0)
        joint = five_coefficients_formula(model, 1, 1, 1, 0)
        self.assertNotEqual(tuple(a ^ b ^ c for a, b, c in zip(base, k_only, r_only)), joint)

    def test_actual_n83_target_support_ranks_and_exponential_obstruction(self):
        path = Path(__file__).parent / "instances/n83-d18-m5-T001-offset0.json"
        receipt = profile_instance(path, outer_samples=1)
        self.assertEqual(receipt["support_dimensions"], [18, 35, 52, 69, 83])
        self.assertEqual(receipt["conditional_constraint_rank_ceiling"], 93)
        self.assertEqual(receipt["conditional_nullity_floor"], 73)
        self.assertEqual(int(receipt["negative_target"][1], 16),
                         int(receipt["target"][0], 16) ^ int(receipt["target"][1], 16))
        self.assertFalse(receipt["solver_complete"])
        for row in receipt["profiles"]:
            self.assertLessEqual(row["rank"], 93)
            if row["consistent"]:
                self.assertGreaterEqual(row["nullity"], 73)
                self.assertGreaterEqual(row["relaxed_cardinality_exact"], 1 << 73)


if __name__ == "__main__":
    unittest.main()

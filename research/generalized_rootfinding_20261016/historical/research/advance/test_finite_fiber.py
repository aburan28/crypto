"""Independent exhaustive small-field oracles for finite-fiber recovery."""

import itertools
import random
import unittest

from rootfinder.binaryfield import BinaryField
from research.advance.finite_fiber import (
    BinaryQuadraticFiber, evaluate_s3, five_source_intermediate_fiber,
)


class FiniteFiberTests(unittest.TestCase):
    def test_quadratics_against_every_field_element(self):
        for degree in (1, 2, 3):
            field = BinaryField(degree)
            solver = BinaryQuadraticFiber(field)
            for a, b, c in itertools.product(range(field.order), repeat=3):
                if a == b == c == 0:
                    with self.assertRaises(ValueError):
                        solver.roots(a, b, c)
                    continue
                expected = tuple(x for x in range(field.order)
                                 if field.mul(a, field.square(x))
                                 ^ field.mul(b, x) ^ c == 0)
                self.assertEqual(solver.roots(a, b, c), expected)

    def test_full_chains_against_independent_triple_enumeration(self):
        rng = random.Random(20261012)
        for degree in (2, 3):
            field = BinaryField(degree)
            solver = BinaryQuadraticFiber(field)
            for _ in range(20):
                source = tuple(rng.randrange(field.order) for _ in range(5))
                target = rng.randrange(field.order)
                curve_b = rng.randrange(1, field.order)
                expected = []
                for u1, u2, u3 in itertools.product(range(field.order), repeat=3):
                    triples = ((source[0], source[1], u1), (u1, source[2], u2),
                               (u2, source[3], u3), (u3, source[4], target))
                    if all(evaluate_s3(field, *triple, curve_b=curve_b) == 0
                           for triple in triples):
                        expected.append((u1, u2, u3))
                result = five_source_intermediate_fiber(
                    field, source, target, curve_b=curve_b, quadratic_solver=solver)
                self.assertEqual(result["status"], "complete")
                self.assertEqual(result["intermediate_chains"], expected)
                self.assertEqual(result["solution_count_exact"], len(expected))

    def test_zero_curve_constant_does_not_claim_bounded_fibers(self):
        field = BinaryField(3)
        with self.assertRaises(ValueError):
            five_source_intermediate_fiber(field, (0,) * 5, 0, curve_b=0)

    def test_even_degree_artin_schreier_map(self):
        field = BinaryField(4)
        solver = BinaryQuadraticFiber(field)
        for rhs in range(field.order):
            expected = tuple(x for x in range(field.order) if field.square(x) ^ x == rhs)
            self.assertEqual(solver.artin_schreier(rhs), expected)


if __name__ == "__main__":
    unittest.main()

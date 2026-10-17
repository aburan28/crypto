"""Independent small-field point-triple oracles; these are not benchmarks."""

import itertools
import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from rootfinder.binaryfield import BinaryField
from rootfinder.fields import PrimeField
from rootfinder import polynomial as poly
from research.advance.divisor_triple import (add_points, conditional_family, cubic_from_parameters,
                            cubic_in_domain, cubic_domain_roots, domain_contains, enumerate_family,
                            negate_point, on_curve, prepare_model, solve_triples)


def point_list(field, curve_a, curve_b, domain=None):
    """Brute evaluate every ordinate, without using root isolation."""
    if domain is None:
        domain = range(field.order)
    points = []
    for x in domain:
        for y in range(field.order):
            x2 = field.mul(x, x)
            if field.mul(y, y) ^ field.mul(x, y) == (
                    field.mul(x2, x) ^ field.mul(curve_a, x2) ^ curve_b):
                points.append((x, y))
    return points


def domain_list(offset, basis):
    result = []
    for bits in itertools.product((0, 1), repeat=len(basis)):
        value = offset
        for bit, vector in zip(bits, basis):
            if bit:
                value ^= vector
        result.append(value)
    return result


def oracle_add(field, curve_a, left, right):
    """Derive the third intersection from a line, independently of backend."""
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and (y1 != y2 or x1 == 0):
        return None
    # On the tangent, implicit differentiation gives x*y' = y+x².
    numerator, denominator = (y1 ^ field.mul(x1, x1), x1) if left == right else (
        y1 ^ y2, x1 ^ x2)
    m = field.mul(numerator, field.inv(denominator))
    intercept = y1 ^ field.mul(m, x1)
    # Substitution y=m*x+n in the curve gives monic intersection cubic
    # x³+(A+m²+m)*x²+n*x+(B+n²).
    xthird = curve_a ^ field.mul(m, m) ^ m ^ x1 ^ x2
    ythird = field.mul(m, xthird) ^ intercept
    return xthird, ythird ^ xthird


def brute_triples(field, curve_a, curve_b, target, offset, basis):
    allowed = point_list(field, curve_a, curve_b, domain_list(offset, basis))
    return {triple for triple in itertools.product(allowed, repeat=3)
            if oracle_add(field, curve_a,
                          oracle_add(field, curve_a, triple[0], triple[1]),
                          triple[2]) == target}


class DivisorTripleTests(unittest.TestCase):
    def check_oracle(self, field, curve_a, curve_b, target, offset, basis,
                     **kwargs):
        expected = brute_triples(field, curve_a, curve_b, target, offset, basis)
        result = solve_triples(field, curve_a, curve_b, target, offset, basis,
                               **kwargs)
        self.assertTrue(result["complete"], result)
        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["solution_count_exact"], len(expected))
        self.assertEqual(result["solution_count_lower_bound"], len(expected))
        actual = set(result["solutions"])
        self.assertEqual(len(actual), len(result["solutions"]))
        self.assertLessEqual(actual, expected)
        if result["output_complete"]:
            self.assertEqual(actual, expected)
        self.assertEqual(result["stats"]["candidates_visited"],
                         result["stats"]["relaxed_candidates_planned"])
        self.assertEqual(result["stats"]["outer_visited"],
                         result["stats"]["outer_cardinality_exact"])
        return result, expected

    def test_all_targets_small_fields_and_nonzero_curve_a(self):
        for degree in (1, 2, 3):
            f = BinaryField(degree)
            domains = [(0, ()), (1, (1,)),
                       (0, tuple(1 << i for i in range(degree)))]
            if degree > 1:
                domains.append((2, (1,)))
            for a in (0, 1):
                for b in sorted(set((1, f.order - 1))):
                    for target in point_list(f, a, b):
                        for offset, basis in domains:
                            with self.subTest(degree=degree, a=a, b=b,
                                              target=target, domain=(offset, basis)):
                                self.check_oracle(f, a, b, target, offset, basis)

    def test_affine_scaled_nonprefix_subspaces(self):
        f = BinaryField(4)
        domains = ((3, (2, 9)), (7, (3, 8)), (15, (1, 2, 4)))
        for a, b in ((0, 1), (1, 2), (9, 7)):
            points = point_list(f, a, b)
            for target in points[::max(1, len(points) // 4)]:
                for offset, basis in domains:
                    self.check_oracle(f, a, b, target, offset, basis)

    def test_repeated_point_and_torsion_target_are_covered(self):
        f = BinaryField(3)
        torsion = (0, 1)
        result, expected = self.check_oracle(f, 0, 1, torsion, 0, (1, 2, 4))
        self.assertIn((torsion, torsion, torsion), expected)
        self.assertIn((torsion, torsion, torsion), result["solutions"])
        self.assertTrue(any(len(set(point[0] for point in triple)) < 3
                            for triple in expected))

    def test_repeated_cubic_support_uses_cube_of_vanishing_polynomial(self):
        f = BinaryField(4)
        model = prepare_model(f, 0, 1, (0, 1), 0, (1, 2))
        cubic = poly.mul(poly.mul((1, 1), (1, 1), f), (1, 1), f)
        self.assertTrue(cubic_in_domain(model, cubic))
        # L_D itself cannot have the triple root because L_D is squarefree.
        ld = (0, 1)
        for x in (1, 2, 3):
            ld = poly.mul(ld, (x, 1), f)
        self.assertTrue(poly.divmod_poly(ld, cubic, f)[1])
        outside = poly.mul(cubic[:], (8, 1), f)
        outside = poly.divmod_poly(outside, (1, 1), f)[0]
        self.assertFalse(cubic_in_domain(model, outside))

    def test_quotient_coefficients_match_literal_norm_for_every_parameter(self):
        f = BinaryField(3)
        for a, b in ((0, 1), (1, 7), (3, 2)):
            target = point_list(f, a, b)[-1]
            model = prepare_model(f, a, b, target, 3, (1, 4))
            for parameter_b in range(f.order):
                for c in range(f.order):
                    h, d = cubic_from_parameters(model, parameter_b, c)
                    z, yt = model["negative_target"]
                    self.assertEqual(d ^ f.mul(c, z) ^ f.mul(z, z)
                                     ^ f.mul(parameter_b, yt), 0)
                    # Reconstruct the norm independently by evaluation at
                    # every field x; no conditional coefficient formulas.
                    for x in range(f.order):
                        ap = f.mul(x, x) ^ f.mul(c, x) ^ d
                        rhs = f.mul(f.mul(x, x), x) ^ f.mul(a, f.mul(x, x)) ^ b
                        norm = f.mul(ap, ap) ^ f.mul(parameter_b, f.mul(x, ap))
                        norm ^= f.mul(f.mul(parameter_b, parameter_b), rhs)
                        self.assertEqual(f.mul(x ^ z, poly.evaluate_poly(h, x, f)), norm)

    def test_conditional_families_match_brute_product_span_constraints(self):
        f = BinaryField(4)
        model = prepare_model(f, 7, 3, point_list(f, 7, 3)[-1], 9, (2, 7))
        w2 = set(domain_list(0, model["w2_basis"]))
        w3 = set(domain_list(0, model["w3_basis"]))
        outer = set(enumerate_family(model["outer_family"]))
        brute_outer = set()
        for b in range(f.order):
            family = conditional_family(model, b)
            expected = set()
            for c in range(f.order):
                e3, e2, e1, _ = cubic_from_parameters(model, b, c)[0]
                a = model["offset"]
                a2 = f.mul(a, a)
                if domain_contains(model, e1):
                    brute_outer.add(b)
                    if e2 ^ a2 in w2 and (
                        e3 ^ f.mul(a2, e1) ^ f.mul(a, e2) ^ f.mul(a2, a)) in w3:
                        expected.add(c)
            self.assertEqual(set(enumerate_family(family)), expected)
            self.assertEqual(family["cardinality"], len(expected))
        self.assertEqual(outer, brute_outer)

    def test_output_cap_preserves_exact_count(self):
        f = BinaryField(3)
        result, expected = self.check_oracle(f, 0, 1, (0, 1), 0, (1, 2, 4),
                                            max_output=1)
        self.assertGreater(len(expected), 1)
        self.assertFalse(result["output_complete"])
        self.assertEqual(len(result["solutions"]), 1)
        zero, expected = self.check_oracle(f, 0, 1, (0, 1), 0, (1, 2, 4),
                                          max_output=0)
        self.assertFalse(zero["solutions"])
        self.assertEqual(zero["solution_count_exact"], len(expected))

    def test_outer_and_candidate_caps_apply_before_family_enumeration(self):
        f = BinaryField(3)
        outer = solve_triples(f, 0, 1, (0, 1), 0, (1, 2, 4), max_outer=0)
        self.assertFalse(outer["complete"])
        self.assertEqual(outer["stats"]["outer_visited"], 0)
        inner = solve_triples(f, 0, 1, (0, 1), 0, (1, 2, 4), max_candidates=0)
        self.assertFalse(inner["complete"])
        self.assertEqual(inner["stats"]["candidates_visited"], 0)
        self.assertIsNone(inner["solution_count_exact"])

    def test_memory_cap_returns_verified_lower_bound(self):
        f = BinaryField(3)
        result = solve_triples(f, 0, 1, (0, 1), 0, (1, 2, 4), max_unique=1)
        self.assertFalse(result["complete"])
        self.assertIsNone(result["solution_count_exact"])
        self.assertEqual(result["solution_count_lower_bound"], 1)
        self.assertEqual(len(result["solutions"]), 1)
        expected = brute_triples(f, 0, 1, (0, 1), 0, (1, 2, 4))
        self.assertLessEqual(set(result["solutions"]), expected)

    def test_compact_subspace_splitting_needs_no_randomized_trials(self):
        f = BinaryField(4)
        result, expected = self.check_oracle(f, 0, 1, (0, 1), 0, (1, 2, 4, 8),
                                             max_split_trials=0)
        self.assertTrue(result["complete"])
        self.assertEqual(result["solution_count_exact"], len(expected))
        self.assertEqual(result["stats"]["split_trials"], 0)

    def test_compact_domain_cubic_roots_against_direct_field_evaluation(self):
        f = BinaryField(5)
        model = prepare_model(f, 0, 1, (0, 1), 19, (3, 8, 16))
        rng = random.Random(7283)
        domain = domain_list(19, (3, 8, 16))
        for _ in range(128):
            chosen = tuple(rng.choice(domain) for _ in range(3))
            cubic = (1,)
            for x in chosen:
                cubic = poly.mul(cubic, (x, 1), f)
            result = cubic_domain_roots(model, cubic)
            expected = {x for x in range(f.order) if poly.evaluate_poly(cubic, x, f) == 0}
            self.assertEqual(set(result["roots"]), expected)
            self.assertEqual(result["distinct_root_count_exact"], len(set(chosen)))
        unsupported = poly.mul(poly.mul((domain[0], 1), (domain[1], 1), f),
                               (next(x for x in range(f.order) if x not in domain), 1), f)
        self.assertIsNone(cubic_domain_roots(model, unsupported))

    def test_validation(self):
        f = BinaryField(3)
        for kwargs in ({"curve_b": 0}, {"target": None},
                       {"target": (0, 0)}, {"basis": (1, 1)},
                       {"basis": (0,)}, {"offset": 8},
                       {"max_candidates": True}, {"max_unique": -1}):
            arguments = dict(field=f, curve_a=0, curve_b=1, target=(0, 1),
                             offset=0, basis=(1,))
            arguments.update(kwargs)
            with self.assertRaises((TypeError, ValueError)):
                solve_triples(**arguments)
        with self.assertRaises(TypeError):
            prepare_model(PrimeField(2), 0, 1, (0, 1), 0, ())

    def test_group_law_against_independent_line_oracle(self):
        f = BinaryField(4)
        for a, b in ((0, 1), (1, 7), (3, 2)):
            points = [None] + point_list(f, a, b)
            for left, right in itertools.product(points, repeat=2):
                self.assertEqual(add_points(f, a, left, right),
                                 oracle_add(f, a, left, right))
                total = add_points(f, a, left, right)
                self.assertTrue(on_curve(f, a, b, total))
            for point in points:
                self.assertIsNone(add_points(f, a, point, negate_point(point)))

    def test_full_width_support_spans_without_field_enumeration(self):
        f = BinaryField(83, 0xd84fae489bcdf275f7b61)
        # Torsion target is rational without any full-field ordinate search.
        model = prepare_model(f, 0, 1, (0, 1), 0,
                              tuple(1 << index for index in range(18)))
        self.assertEqual(len(model["w2_basis"]), 35)
        self.assertEqual(len(model["w3_basis"]), 52)
        self.assertLessEqual(model["outer_family"]["cardinality"], 1 << 19)
        family = conditional_family(model, 0)
        self.assertTrue(family["consistent"])
        self.assertGreaterEqual(family["nullity"], 4)
        self.assertEqual(family["cardinality"], 1 << family["nullity"])

    def test_full_width_affine_cubic_reconstruction_and_deep_coset_splits(self):
        f = BinaryField(83, 0xd84fae489bcdf275f7b61)
        scale = 0x3211b4df3bc4ff9f0808b
        offset = 0x574cdce6672c80407c4ec
        basis = tuple(f.mul(scale, 1 << index) for index in range(18))
        model = prepare_model(f, 1, 1, (0, 1), offset, basis)
        # These roots differ only at the bottom of the 18-dimensional
        # nested-subspace flag, exercising long same-coset paths.
        for chosen in ((offset, offset ^ basis[0], offset ^ basis[1]),
                       (offset ^ basis[0], offset ^ basis[0], offset ^ basis[1])):
            self.assertTrue(all(root.bit_length() > 70 for root in chosen))
            cubic = (1,)
            for root in chosen:
                cubic = poly.mul(cubic, (root, 1), f)
            result = cubic_domain_roots(model, cubic)
            self.assertEqual(set(result["roots"]), set(chosen))
            self.assertGreaterEqual(result["gcd_split_steps"], 17)


if __name__ == "__main__":
    unittest.main()

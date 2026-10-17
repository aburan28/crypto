"""Independent finite-coordinate oracles for the separator-fiber join."""

from functools import lru_cache
import itertools
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from rootfinder.binaryfield import BinaryField
from research.advance.chain_join import solve_chain_join


def instance(degree, dimension, target, basis=None):
    field = BinaryField(degree)
    return {"extension_field": {"characteristic": 2, "degree": degree,
                                "modulus_hex": hex(field.modulus)},
            "source_dimension": dimension, "summands": 5,
            "n_vars": 5 * dimension + 3 * degree,
            "source_basis_hex": [hex(value) for value in
                                 (basis if basis is not None else [1 << i for i in range(dimension)])],
            "target": {"x_hex": hex(target)}, "curve": {"b_hex": "0x1"}}


def _s3(field, x, y, z, curve_b):
    # Independently expanded original polynomial, not the circuit identity.
    product = field.mul(x, y)
    return (field.mul(field.square(x ^ y), field.square(z))
            ^ field.mul(field.mul(x, y), z) ^ field.square(product) ^ curve_b)


def full_boolean_oracle(problem):
    field = BinaryField(problem["extension_field"]["degree"],
                        int(problem["extension_field"]["modulus_hex"], 16))
    d, n = problem["source_dimension"], field.degree
    basis = [int(value, 16) for value in problem["source_basis_hex"]]
    target = int(problem["target"]["x_hex"], 16)
    curve_b = int(problem["curve"]["b_hex"], 16)
    domain = [0] * (1 << d)
    for mask in range(1 << d):
        for i, vector in enumerate(basis):
            if mask & (1 << i):
                domain[mask] ^= vector
    result = set()
    for mask in range(1 << problem["n_vars"]):
        source = [domain[(mask >> (i * d)) & ((1 << d) - 1)] for i in range(5)]
        u = [(mask >> (5 * d + i * n)) & (field.order - 1) for i in range(3)]
        links = ((source[0], source[1], u[0]), (u[0], source[2], u[1]),
                 (u[1], source[3], u[2]), (u[2], source[4], target))
        if all(_s3(field, *link, curve_b) == 0 for link in links):
            result.add(mask)
    return result


def source_sweep_halftrace_oracle(problem):
    """Independent source sweep with half-trace, no cached matrix or join."""
    field = BinaryField(problem["extension_field"]["degree"],
                        int(problem["extension_field"]["modulus_hex"], 16))
    if field.degree % 2 != 1:
        raise ValueError("this independent oracle uses odd-degree half-trace")
    dimension = problem["source_dimension"]
    basis = [int(value, 16) for value in problem["source_basis_hex"]]
    target = int(problem["target"]["x_hex"], 16)
    curve_b = int(problem["curve"]["b_hex"], 16)
    domain = []
    for mask in range(1 << dimension):
        value = 0
        for bit, vector in enumerate(basis):
            if mask & (1 << bit):
                value ^= vector
        domain.append(value)

    @lru_cache(maxsize=None)
    def roots(left, right):
        a, b = field.square(left ^ right), field.mul(left, right)
        c = field.square(b) ^ curve_b
        if not a:
            result = () if not b else (field.mul(c, field.inv(b)),)
        elif not b:
            result = (field.pow(field.mul(c, field.inv(a)), 1 << (field.degree - 1)),)
        else:
            rhs = field.mul(field.mul(a, c), field.inv(field.square(b)))
            trace, power = rhs, rhs
            for _ in range(1, field.degree):
                power = field.square(power)
                trace ^= power
            if trace:
                result = ()
            else:
                halftrace, power = rhs, rhs
                for _ in range(1, (field.degree + 1) // 2):
                    power = field.square(field.square(power))
                    halftrace ^= power
                scale = field.mul(b, field.inv(a))
                result = tuple(sorted((field.mul(scale, halftrace),
                                       field.mul(scale, halftrace ^ 1))))
        if any(_s3(field, left, right, root, curve_b) for root in result):
            raise ArithmeticError("independent half-trace oracle produced an invalid root")
        return result

    result = set()
    for source_bits in itertools.product(range(len(domain)), repeat=5):
        x1, x2, x3, x4, x5 = (domain[bits] for bits in source_bits)
        for u1 in roots(x1, x2):
            for u2 in roots(u1, x3):
                for u3 in roots(u2, x4):
                    if _s3(field, u3, x5, target, curve_b):
                        continue
                    mask = sum(bits << (index * dimension)
                               for index, bits in enumerate(source_bits))
                    mask |= u1 << (5 * dimension)
                    mask |= u2 << (5 * dimension + field.degree)
                    mask |= u3 << (5 * dimension + 2 * field.degree)
                    result.add(mask)
    return result


class ChainJoinTests(unittest.TestCase):
    def test_every_tiny_target_matches_full_boolean_enumeration(self):
        for degree in (1, 2, 3):
            for target in range(1 << degree):
                problem = instance(degree, 1, target)
                expected = full_boolean_oracle(problem)
                result = solve_chain_join(problem, max_output=100_000)
                self.assertTrue(result["output_complete"], result)
                self.assertEqual(result["solution_count_exact"], len(expected))
                self.assertEqual({int(value, 16) for value in result["assignments_hex"]}, expected)
                self.assertEqual(len(result["assignments_hex"]), len(expected))

    def test_nonprefix_basis_matches_independent_source_sweep(self):
        problem = instance(3, 2, 5, basis=[7, 3])
        expected = source_sweep_halftrace_oracle(problem)
        result = solve_chain_join(problem, max_output=100_000)
        self.assertTrue(result["output_complete"])
        self.assertEqual({int(value, 16) for value in result["assignments_hex"]}, expected)

    def test_primary_count_can_be_derived_when_omitted(self):
        problem = instance(3, 1, 1)
        expected = solve_chain_join(problem)
        del problem["n_vars"]
        actual = solve_chain_join(problem)
        self.assertEqual(actual["n_vars"], 14)
        self.assertEqual(actual["assignments_hex"], expected["assignments_hex"])

    def test_preflight_rejects_before_materializing_domain(self):
        problem = instance(3, 2, 1)
        for budget in ("max_source_assignments", "max_root_calls", "max_records"):
            with patch("research.advance.chain_join._decode_domain_value",
                       side_effect=AssertionError("domain must not be materialized")):
                result = solve_chain_join(problem, **{budget: 0})
            self.assertFalse(result["search_complete"])
            self.assertIsNone(result["solution_count_exact"])
            self.assertEqual(result["stats"]["root_calls"], 0)
            self.assertFalse(result["stats"]["domain_materialized"])

    def test_output_cap_preserves_exact_count_without_full_output_claim(self):
        problem = instance(3, 1, 1)
        expected = full_boolean_oracle(problem)
        self.assertGreater(len(expected), 0)
        result = solve_chain_join(problem, max_output=0)
        self.assertTrue(result["search_complete"])
        self.assertFalse(result["output_complete"])
        self.assertEqual(result["status"], "output_truncated")
        self.assertEqual(result["solution_count_exact"], len(expected))
        self.assertEqual(result["assignments_hex"], [])

    def test_zero_time_budget_has_no_exact_count(self):
        result = solve_chain_join(instance(3, 1, 1), max_seconds=0)
        self.assertFalse(result["search_complete"])
        self.assertIsNone(result["solution_count_exact"])
        self.assertFalse(result["stats"]["domain_materialized"])

    def test_dependent_basis_rejected(self):
        with self.assertRaises(ValueError):
            solve_chain_join(instance(3, 2, 1, basis=[3, 3]))

    def test_nonzero_curve_constant_required(self):
        problem = instance(3, 1, 1)
        problem["curve"]["b_hex"] = "0x0"
        with self.assertRaises(ValueError):
            solve_chain_join(problem)

    def test_wrong_layout_receipts_and_equation_hash_rejected(self):
        problem = instance(3, 1, 1)
        problem["chain_links"] = []
        with self.assertRaises(ValueError):
            solve_chain_join(problem)
        problem = instance(3, 1, 1)
        problem["native_coordinate_evaluation_checks"] = [
            {"assignment_hex": "0x0", "coordinate_output_hex": "0x0"}]
        with self.assertRaises(ValueError):
            solve_chain_join(problem)
        problem = instance(3, 1, 1)
        problem["equations"] = [[] for _ in range(12)]
        problem["equations_sha256"] = "bad"
        with self.assertRaises(ValueError):
            solve_chain_join(problem)

    def test_n11_export_matches_source_sweep_and_original_equations(self):
        path = Path(__file__).parent / "instances/n11-d3-m5-target17.json"
        if not path.exists():
            self.skipTest("optional native export not installed")
        problem = json.loads(path.read_text())
        expected = source_sweep_halftrace_oracle(problem)
        result = solve_chain_join(problem, max_output=100_000)
        self.assertTrue(result["output_complete"], result)
        self.assertEqual({int(value, 16) for value in result["assignments_hex"]}, expected)
        self.assertEqual(result["solution_count_exact"], len(expected))

    def test_optional_projected_sat_matches_tiny_join(self):
        from research.advance.projected_sat import load_backend, solve_projected_chain
        try:
            load_backend()
        except RuntimeError:
            self.skipTest("optional projected SAT backend unavailable")
        problem = instance(3, 1, 1)
        expected = solve_projected_chain(problem, max_seconds=10)
        actual = solve_chain_join(problem)
        self.assertTrue(expected["output_complete"])
        self.assertEqual(actual["assignments_hex"], expected["assignments_hex"])


if __name__ == "__main__":
    unittest.main()

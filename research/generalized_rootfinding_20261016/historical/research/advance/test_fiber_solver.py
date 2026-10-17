"""Independent exhaustive checks for the residual-span DAG research reference."""

import itertools
import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from rootfinder.binaryfield import BinaryField
from rootfinder.fields import PrimeField, ExtensionField
from fiber_solver import solve_fiber_dag


def oracle(equations, n, domain, field):
    result = set()
    for assignment in itertools.product(domain, repeat=n):
        zero = True
        for equation in equations:
            total = 0
            for coefficient, powers in equation:
                term = coefficient
                for value, exponent in zip(assignment, powers):
                    for _ in range(exponent):
                        term = field.mul(term, value)
                total = field.add(total, term)
            if total:
                zero = False
                break
        if zero:
            result.add(assignment)
    return result


def powers(n, *indices):
    result = [0] * n
    for index in indices:
        result[index] += 1
    return tuple(result)


class FiberDAGTests(unittest.TestCase):
    def check(self, equations, n, domain, field, **options):
        expected = oracle(equations, n, domain, field)
        result = solve_fiber_dag(equations, n, domain, field, **options)
        self.assertTrue(result['count_complete'], result)
        self.assertEqual(result['solution_count_exact'], len(expected))
        actual = {tuple(row) for row in result['solutions']}
        self.assertEqual(len(actual), len(result['solutions']))
        self.assertLessEqual(actual, expected)
        if result['output_complete']:
            self.assertEqual(actual, expected)
        return result

    def test_random_cubic_prime_binary_and_odd_extension(self):
        fields = [PrimeField(2), PrimeField(3), PrimeField(5), BinaryField(3, 0b1011),
                  ExtensionField(3, 2, (1, 0, 1)), ExtensionField(5, 2, (2, 0, 1))]
        n = 4
        monomials = [exponents for exponents in itertools.product(range(4), repeat=n)
                     if sum(exponents) <= 3]
        for field_index, field in enumerate(fields):
            randomizer = random.Random(51187 + field_index)
            for trial in range(12):
                equations = [[(randomizer.randrange(field.order), randomizer.choice(monomials))
                              for _ in range(12)] for _ in range(1 + trial % 4)]
                with self.subTest(field=repr(field), trial=trial):
                    self.check(equations, n, tuple(range(min(field.order, 3))), field,
                               variable_order=(3, 1, 0, 2) if trial % 2 else None)

    def test_gf_2_83_arbitrary_coefficients_and_non_boolean_domain(self):
        field = BinaryField(83, 0xd84fae489bcdf275f7b61)
        n = 3
        equations = [[(1 << 76, powers(n, 0, 1, 2)),
                      (1 << 65, powers(n, 0, 0, 0, 0)),
                      (17, powers(n, 1, 1)), (17, powers(n, 1, 1)),
                      (9, powers(n, 2)), (9, powers(n))]]
        self.check(equations, n, (0, 1, 1 << 77), field)

    def test_domain_normal_form_collapses_high_degree(self):
        # On GF(5), x^5=x; all variables therefore become free after cancellation.
        field = PrimeField(5)
        n = 4
        equation = [(1, powers(n, 0, 0, 0, 0, 0)), (4, powers(n, 0))]
        result = self.check([equation], n, range(5), field, max_solutions=2)
        self.assertEqual(result['solution_count_exact'], 625)
        self.assertEqual(result['metrics']['states'], 1)

    def test_non_subfield_domain_remainder_is_exact(self):
        field = PrimeField(7)
        self.check([[(1, (8, 0, 0)), (3, (0, 5, 1)), (1, (0, 0, 0))]],
                   3, (1, 3, 6), field)

    def test_unary_forcing_and_linear_span_prune(self):
        field = PrimeField(5)
        equations = [[(1, (1, 0, 0)), (4, (0, 0, 0))],
                     [(1, (0, 1, 0)), (3, (0, 0, 0))],
                     [(1, (0, 0, 1)), (2, (0, 0, 0))]]
        result = self.check(equations, 3, range(5), field)
        self.assertEqual(result['solutions'], [[1, 2, 3]])
        self.assertEqual(result['metrics']['transitions'], 2)
        self.assertEqual(result['metrics']['unary_restrictions'], 2)
        inconsistent = equations + [[(1, (1, 0, 0)), (3, (0, 0, 0))]]
        result = self.check(inconsistent, 3, range(5), field)
        self.assertEqual(result['metrics']['states'], 1)
        self.assertEqual(result['metrics']['constant_prunes'], 1)

    def test_row_span_equivalence_merges_distinct_specializations(self):
        # For x0=1 and 2, the remaining equation differs only by nonzero scalar.
        field = PrimeField(3)
        equations = [[(1, (1, 1, 0)), (1, (1, 0, 1))]]
        result = self.check(equations, 3, (1, 2), field)
        self.assertGreaterEqual(result['metrics']['cache_hits'], 1)
        self.assertGreaterEqual(result['metrics']['domain_root_cache_hits'], 1)

    def test_connected_cubic_cumulative_sum_merges_states(self):
        field = PrimeField(3)
        n = 6
        # Dense polynomial in S=sum(x): S^3=S in characteristic three, so
        # use S^2-1, plus cubic interaction S*x0*x1. All variables are coupled.
        equation = [(2, powers(n))]
        for i in range(n):
            equation.append((1, powers(n, i, i)))
            for j in range(i):
                equation.append((2, powers(n, i, j)))
        cubic = [(1, powers(n, i, 0, 1)) for i in range(n)]
        result = self.check([equation, cubic], n, range(3), field)
        self.assertLess(result['metrics']['states'], result['exhaustive_assignments'])
        self.assertGreater(result['metrics']['cache_hits'], 0)

    def test_free_variables_and_output_caps_preserve_exact_count(self):
        for cap in [0, 1, 3, 64]:
            result = self.check([[(1, (1, 0, 0, 0)), (2, (0, 0, 0, 0))]],
                                4, range(3), PrimeField(3), max_solutions=cap)
            self.assertEqual(result['solution_count_exact'], 27)
            self.assertEqual(len(result['solutions']), min(cap, 27))
            self.assertEqual(result['metrics']['free_variable_collapses'], 3)

    def test_count_only_does_not_materialize_terminal_domain_roots(self):
        field = PrimeField(7)
        result = self.check([[(1, (2,)), (6, (0,))]], 1, range(7), field, max_solutions=0)
        self.assertEqual(result['solution_count_exact'], 2)
        self.assertEqual(result['metrics']['domain_root_tests'], 0)

    def test_singleton_and_empty_domains(self):
        self.check([[(1, (3, 1)), (4, (0, 0))]], 2, (1,), PrimeField(5))
        result = solve_fiber_dag([], 3, (), PrimeField(3), max_states=0, max_arithmetic=0)
        self.assertEqual(result['solution_count_exact'], 0)

    def test_budgets_never_promote_incomplete_count(self):
        field = PrimeField(3)
        equation = [(1, (2, 0, 1)), (2, (0, 1, 0))]
        for options in [{'max_states': 0}, {'max_transitions': 0}, {'max_arithmetic': 0},
                        {'max_terms': 0}, {'max_cached_terms': 0}]:
            with self.subTest(options=options):
                result = solve_fiber_dag([equation], 3, range(3), field, **options)
                self.assertEqual(result['status'], 'budget_exceeded')
                self.assertIsNone(result['solution_count_exact'])
                self.assertFalse(result['count_complete'])

    def test_arithmetic_budget_during_extraction_preserves_completed_count(self):
        field = PrimeField(7)
        equation = [(1, (2,)), (6, (0,))]
        count_result = solve_fiber_dag([equation], 1, range(7), field, max_solutions=0)
        budget = count_result['metrics']['arithmetic_primitives']
        result = solve_fiber_dag([equation], 1, range(7), field, max_arithmetic=budget)
        self.assertEqual(result['status'], 'output_truncated')
        self.assertEqual(result['solution_count_exact'], 2)
        self.assertTrue(result['count_complete'])
        self.assertFalse(result['output_complete'])

    def test_validation(self):
        field = PrimeField(3)
        bad = [dict(domain=(0, 0)), dict(variable_order=(0, 0)), dict(max_states=-1)]
        for options in bad:
            params = {'domain': (0, 1)} | options
            with self.subTest(options=options), self.assertRaises(ValueError):
                solve_fiber_dag([[(1, (1, 0))]], 2, field=field, **params)
        with self.assertRaises(ValueError):
            solve_fiber_dag([[(1, (1,))]], 2, (0, 1), field)
        with self.assertRaises(ValueError):
            solve_fiber_dag([[(3, (1, 0))]], 2, (0, 1), field)


if __name__ == '__main__':
    unittest.main()

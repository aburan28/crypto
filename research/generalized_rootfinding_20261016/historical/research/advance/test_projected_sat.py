import unittest

from rootfinder.binaryfield import BinaryField
from research.advance.finite_fiber import evaluate_s3
from research.advance.projected_sat import solve_projected_chain, load_backend


def instance(degree, dimension, target):
    field = BinaryField(degree)
    return {"extension_field": {"characteristic": 2, "degree": degree,
                                "modulus_hex": hex(field.modulus)},
            "source_dimension": dimension, "summands": 5,
            "n_vars": 5 * dimension + 3 * degree,
            "source_basis_hex": [hex(1 << i) for i in range(dimension)],
            "target": {"x_hex": hex(target)}, "curve": {"b_hex": "0x1"}}


def oracle(problem):
    f = BinaryField(problem["extension_field"]["degree"], int(problem["extension_field"]["modulus_hex"], 16))
    d, n = problem["source_dimension"], f.degree
    target = int(problem["target"]["x_hex"], 16)
    result = set()
    for mask in range(1 << problem["n_vars"]):
        sources = [(mask >> (i * d)) & ((1 << d) - 1) for i in range(5)]
        u = [(mask >> (5 * d + i * n)) & ((1 << n) - 1) for i in range(3)]
        links = [(sources[0], sources[1], u[0]), (u[0], sources[2], u[1]),
                 (u[1], sources[3], u[2]), (u[2], sources[4], target)]
        if all(evaluate_s3(f, *link) == 0 for link in links):
            result.add(mask)
    return result


class ProjectedSATTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        load_backend()

    def test_every_tiny_target_matches_independent_full_boolean_oracle(self):
        for degree in (1, 2, 3):
            for target in range(1 << degree):
                problem = instance(degree, 1, target)
                expected = oracle(problem)
                for projection in ("sources", "all_variables"):
                    result = solve_projected_chain(problem, projection=projection, max_seconds=10)
                    self.assertTrue(result["output_complete"], result)
                    self.assertEqual(result["solution_count_exact"], len(expected))
                    self.assertEqual({int(v, 16) for v in result["assignments_hex"]}, expected)

    def test_source_projection_uses_one_query_per_fiber(self):
        problem = instance(3, 1, 1)
        projected = solve_projected_chain(problem, projection="sources", max_seconds=10)
        full = solve_projected_chain(problem, projection="all_variables", max_seconds=10)
        self.assertEqual(projected["assignments_hex"], full["assignments_hex"])
        self.assertEqual(projected["stats"]["sat_calls"], projected["stats"]["source_tuples"] + 1)
        self.assertEqual(full["stats"]["sat_calls"], full["solution_count_exact"] + 1)
        self.assertLessEqual(projected["stats"]["sat_calls"], full["stats"]["sat_calls"])

    def test_output_cap_retains_exact_count_but_not_full_output_claim(self):
        problem = instance(3, 1, 1)
        result = solve_projected_chain(problem, max_seconds=10, max_output=0)
        self.assertTrue(result["search_complete"])
        self.assertGreater(result["solution_count_exact"], 0)
        self.assertFalse(result["output_complete"])
        self.assertEqual(result["assignments_hex"], [])

    def test_zero_query_budget_cannot_claim_unsat_or_exact_count(self):
        result = solve_projected_chain(instance(3, 1, 1), max_seconds=10, max_sat_calls=0)
        self.assertFalse(result["search_complete"])
        self.assertIsNone(result["solution_count_exact"])
        self.assertEqual(result["stats"]["sat_calls"], 0)

    def test_native_exported_ladder_receipts_before_budget_refusal(self):
        import json
        from pathlib import Path
        path = Path(__file__).parent / "instances/n11-d3-m5-target17.json"
        if not path.exists():
            self.skipTest("optional upstream export fixture not installed")
        result = solve_projected_chain(json.loads(path.read_text()), max_seconds=20,
                                       max_sat_calls=0, max_output=10_000)
        self.assertFalse(result["search_complete"])
        self.assertIsNone(result["solution_count_exact"])
        self.assertEqual(result["n_vars"], 48)
        self.assertEqual(result["stats"]["sat_calls"], 0)
        self.assertGreater(result["stats"]["xor_rows"], 0)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Pre-panel independent arithmetic and rank checks; no frozen targets."""
import json
from pathlib import Path
import unittest

import verify


HERE = Path(__file__).resolve().parent
PILOT_RESULT = (verify.ROOT /
    "research/ecc2k130_factor_base_pilot_20260924/results_final.json")


class ReferenceArithmeticTests(unittest.TestCase):
    def test_unmetered_degree_seven_geometry(self):
        verify.restore_bare_curve()
        field = verify.pilot.FastGF2m(21, verify.pilot.IRR)
        source = verify.pilot.Koblitz(field, 0, 1)
        twist = verify.pilot.Koblitz(field, 1, 1)
        lines, selected, _ = verify.pilot.order_seven_lines(field, twist)
        complement = next(lines[key]["generator"] for key in sorted(lines)
                          if key != selected)
        generator, _ = verify.pilot.choose_generator(source)
        transport = verify.pilot.BinaryVeluMap.from_generator(
            source, twist, lines[selected]["generator"], 7)
        dual = verify.DualTransport(source, twist, lines[selected]["generator"],
                                    complement, 7, generator)
        self.assertEqual(transport(generator), dual.forward(generator))
        self.assertEqual(dual.compose(generator), source.mul(generator, 7))

    def test_bit_polynomial_field_matches_existing_field(self):
        reference = verify.ReferenceField()
        fast = verify.pilot.FastGF2m(21, verify.pilot.IRR)
        values = [1, 2, 3, 0x12345, 0x1fffff]
        state = 0x12ab34
        for _ in range(16):
            state = (state * 1103515245 + 12345) & reference.mask
            values.append(state)
        for a in values:
            self.assertEqual(reference.sqr(a), fast.sqr(a))
            if a:
                self.assertEqual(reference.inv(a), fast.inv(a))
            for b in values:
                self.assertEqual(reference.mul(a, b), fast.mul(a, b))

    def test_source_and_leaf_point_laws(self):
        saved = json.loads(PILOT_RESULT.read_text())
        field = verify.pilot.FastGF2m(21, verify.pilot.IRR)
        source = verify.pilot.Koblitz(field, 0, 1)
        leaf_b = saved["geometry"]["selected_codomain_b"]
        leaf = verify.pilot.Koblitz(field, 0, leaf_b)
        source_generator = tuple(saved["challenge"]["generator"])
        leaf_generator = None
        for x in range(1, 128):
            lifts = leaf.points_over(x)
            if lifts:
                candidate = leaf.mul(lifts[0], verify.pilot.COFACTOR)
                if candidate is not None:
                    leaf_generator = candidate
                    break
        self.assertIsNotNone(leaf_generator)
        for curve, generator in ((source, source_generator),
                                 (leaf, leaf_generator)):
            reference = verify.ReferenceCurve(curve.b)
            self.assertTrue(reference.on_curve(generator))
            self.assertIsNone(reference.mul(generator, 421))
            for a in (0, 1, 2, 3, 7, 21, 100):
                P = curve.mul(generator, a)
                self.assertEqual(reference.mul(generator, a), P)
                self.assertTrue(reference.on_curve(P))
                for b in (0, 1, 2, 13):
                    Q = curve.mul(generator, b)
                    self.assertEqual(reference.add(P, Q), curve.add(P, Q))

    def test_rank_solver_checks_full_solution_and_dependencies(self):
        expected = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        rows = []
        for i in range(9):
            row = [0] * 9
            row[i] = 1
            row[(i + 1) % 9] = 2
            rows.append((row, sum(a * b for a, b in zip(row, expected)) % 421))
        rank, solution = verify.independent_linear_system(rows, 9)
        self.assertEqual(rank, 9)
        self.assertEqual(solution, expected)
        self.assertEqual(verify.independent_linear_system(rows + [rows[0]], 9),
                         (9, expected))


if __name__ == "__main__":
    unittest.main()

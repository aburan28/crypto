"""Independent finite-space and scalar-elimination checks."""
import time
import unittest

import experiment
import e2e
import pullback


class ExactAlgebraTests(unittest.TestCase):
    def test_restricted_linear_map(self):
        oracle = experiment.previous.PairOracle(5, 3)
        f = oracle.f
        target = experiment.nagaocompare.affinePoints(f, oracle.curve)[1]
        search = pullback.Search(f, oracle.curve, 3, target, time.perf_counter()+60)
        for basis in search.imageSquares.values():
            space = [0]
            for value, _ in basis:
                space += [f.add(v, value) for v in space]
            for c2bits in [1, 3, 17, 31]:
                for c1bits in [0, 1, 7, 29]:
                    c2, c1 = f.fromCoords(c2bits), f.fromCoords(c1bits)
                    fibers = {}
                    for v in space:
                        mapped = f.add(f.mul(c2, f.sqr(v)), f.mul(c1, v))
                        fibers.setdefault(mapped, set()).add(v)
                    for bits in range(32):
                        rhs = f.fromCoords(bits)
                        self.assertEqual(set(search.imagePreimages(basis, c2, c1, rhs)),
                                         fibers.get(rhs, set()))

    def test_modular_rows_and_dependent_relations(self):
        ops = e2e.ScalarOps(127)
        rows = e2e.Rows(3, ops)
        logs = [13, 42, 97]
        for coefficients in [[0, 1, 2], [0, 2, 4], [3, 1, 0], [1, 0, 2]]:
            rhs = sum(a*b for a, b in zip(coefficients, logs)) % 127
            rows.add(coefficients+[rhs])
        self.assertEqual(rows.solve(), logs)
        with self.assertRaises(ArithmeticError):
            rows.add([0, 0, 0, 1])
        with self.assertRaises(ValueError):
            e2e.Rows(3, ops).solve()


if __name__ == '__main__':
    unittest.main()

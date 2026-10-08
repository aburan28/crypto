import random
import unittest

import xl


class TestXL(unittest.TestCase):
    def test_planted_satisfies_every_equation(self):
        F, rng = xl.lr.Field(17), random.Random(1)
        l, b = 5, 2
        for _ in range(20):
            V, W, w0, t, sol = xl.planted_instance(F, l, b, rng)
            c, d = sol[:l], sol[l:]
            cm = sum(1 << i for i in range(l) if c[i])
            dm = sum(1 << j for j in range(b) if d[j])
            for p in xl.equation_polys(F, V, W, w0, t, l, b):
                val = sum(1 for (pc, pd) in p if (pc & cm) == pc and (pd & dm) == pd) & 1
                self.assertEqual(val, 0)

    def test_extracted_linear_polys_vanish_on_planted(self):
        F, rng = xl.lr.Field(17), random.Random(2)
        l, b = 6, 3
        for D1, D2 in ((2, 1), (1, 2)):
            V, W, w0, t, sol = xl.planted_instance(F, l, b, rng)
            rows, order, nn = xl.macaulay(xl.equation_polys(F, V, W, w0, t, l, b), l, b, D1, D2)
            for r in xl.extract_linear(xl.reduce_rows(rows), nn, l, b):
                val = (r >> (l + b)) & 1
                for v in range(l + b):
                    if (r >> v) & 1:
                        val ^= sol[v]
                self.assertEqual(val, 0)

    def test_linearization_cliff(self):
        # (1,1) at n=17, l=5: resolves b=1, fails b=3 (count 18 > 17)
        rng = random.Random(3)
        self.assertGreaterEqual(xl.xl_cell(17, 5, 1, 1, 1, 10, rng)["resolved"], 9)
        self.assertEqual(xl.xl_cell(17, 5, 3, 1, 1, 10, rng)["resolved"], 0)


if __name__ == "__main__":
    unittest.main()

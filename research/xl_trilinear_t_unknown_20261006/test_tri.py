import random
import unittest

import tri


def evaluate(eq, x):
    mask = sum(1 << i for i, v in enumerate(x) if v)
    return sum(1 for m in eq if (m & mask) == m) & 1


class TestTri(unittest.TestCase):
    def test_planted_satisfies_both_blocks(self):
        F, rng = tri.lr.Field(11), random.Random(1)
        for a1, l, b in ((1, 3, 1), (2, 4, 2), (0, 3, 1)):
            inst = tri.planted(F, a1, l, b, rng)
            _, e1, e2 = tri.build(F, inst, a1, l, b)
            for eq in e1 + e2:
                self.assertEqual(evaluate(eq, inst["sol"]), 0)

    def test_x1_lies_in_factor_base(self):
        F, rng = tri.lr.Field(11), random.Random(2)
        inst = tri.planted(F, 2, 4, 1, rng)
        span = {0}
        for v in inst["V"]:
            span |= {s ^ v for s in span}
        self.assertTrue(all(u in span for u in inst["U"]))
        self.assertIn(inst["w1"], span)

    def test_extracted_linear_polys_vanish_on_planted(self):
        F, rng = tri.lr.Field(11), random.Random(3)
        a1, l, b, caps = 1, 3, 1, (1, 1, 1, 2)
        inst = tri.planted(F, a1, l, b, rng)
        Bk, e1, e2 = tri.build(F, inst, a1, l, b)
        rows, order, nn = tri.macaulay(Bk, e1, e2, caps)
        piv = tri.reduce_rows(rows)
        nv = Bk.nvar
        for c, r in piv.items():
            if c >= nn:
                lin = r >> nn
                val = (lin >> nv) & 1
                for v in range(nv):
                    if (lin >> v) & 1:
                        val ^= inst["sol"][v]
                self.assertEqual(val, 0)


if __name__ == "__main__":
    unittest.main()

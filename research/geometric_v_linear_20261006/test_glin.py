import random
import unittest

import glin


class TestGlin(unittest.TestCase):
    def test_identity(self):
        # S3(t, X2, X3) == P^2 + t^2 A^2 + t P + 1 for random field elements
        F, rng = glin.lr.Field(19), random.Random(1)
        for _ in range(200):
            t, x2, x3 = (rng.getrandbits(19) for _ in range(3))
            A, P = x2 ^ x3, F.mul(x2, x3)
            self.assertEqual(glin.lr.S3(F, t, x2, x3),
                             F.sq(P) ^ F.mul(F.sq(t), F.sq(A)) ^ F.mul(t, P) ^ 1)

    def test_geometric_product_span(self):
        F, rng = glin.lr.Field(23), random.Random(2)
        base = glin.Base(F, 6, "geometric", rng)
        self.assertLessEqual(len(base.P), 2 * 6 - 1)
        prods = [F.mul(a, b) for a in base.V for b in base.V]
        self.assertEqual(glin.rank_of(base.P + prods), len(base.P))

    def test_planted_recovered(self):
        F, rng = glin.lr.Field(19), random.Random(3)
        base = glin.Base(F, 5, "geometric", rng)
        for _ in range(10):
            T, pair = glin.planted_target(F, base, rng)
            _, cand = glin.solve(base, T[0])
            good, _ = glin.decompositions(base, T, cand)
            self.assertTrue(any(set(g) == pair for g in good))


if __name__ == "__main__":
    unittest.main()

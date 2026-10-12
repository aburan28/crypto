import random
import unittest

import lr


class TestLR(unittest.TestCase):
    def test_s3_vanishes_on_sums(self):
        F, rng = lr.Field(13), random.Random(1)
        n = 0
        while n < 200:
            P, Q = lr.lift(F, rng.getrandbits(13)), lr.lift(F, rng.getrandbits(13))
            if not P or not Q or P[0] == Q[0]:
                continue
            S = lr.add(F, P, Q)
            if S is None:
                continue
            self.assertEqual(lr.S3(F, P[0], Q[0], S[0]), 0)
            n += 1

    def test_points_on_curve(self):
        F, rng = lr.Field(17), random.Random(2)
        for _ in range(200):
            P = lr.lift(F, rng.getrandbits(17))
            if P:
                x, y = P
                self.assertEqual(F.sq(y) ^ F.mul(x, y), F.mul(F.sq(x), x) ^ 1)

    def test_linearized_columns_match_evaluation(self):
        # sum of selected columns + constant == S3(X2, X3, t) for every assignment
        F, rng = lr.Field(13), random.Random(3)
        l, b = 4, 2
        V = [rng.getrandbits(13) for _ in range(l)]
        W = V[:b]
        w0, t = lr.rand_in(V, rng), rng.getrandbits(13)
        cols, rhs = lr.columns(F, V, W, w0, t)
        for _ in range(100):
            c = [rng.getrandbits(1) for _ in range(l)]
            d = [rng.getrandbits(1) for _ in range(b)]
            mono = [c[i] & d[j] for i in range(l) for j in range(b)] + c + d
            acc = rhs
            for m, cv in zip(mono, cols):
                if m:
                    acc ^= cv
            X2 = 0
            for i in range(l):
                if c[i]:
                    X2 ^= V[i]
            X3 = w0
            for j in range(b):
                if d[j]:
                    X3 ^= W[j]
            self.assertEqual(acc, lr.S3(F, X2, X3, t))

    def test_symmetry_defect(self):
        # W < V: the linear columns of a shared direction coincide when w0 = 0
        F, rng = lr.Field(13), random.Random(4)
        V = [rng.getrandbits(13) for _ in range(4)]
        cols, _ = lr.columns(F, V, V[:1], 0, rng.getrandbits(13))
        self.assertEqual(cols[4], cols[8])  # c_0 column == d_0 column


if __name__ == "__main__":
    unittest.main()

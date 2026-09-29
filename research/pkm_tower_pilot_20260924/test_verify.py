"""Checks of verify.py's exhaustive counts on a small prime.

verify.py decides, without F4, whether a tuple of V^m solves the summation
system. These tests check its three counters against planted decompositions
on curves over F_1009. For m = 4 they also check it against the other order
of eliminating the chain's two free unknowns.

    python3 -m unittest discover -s research/pkm_tower_pilot_20260924 -p 'test_*.py'
"""

import random
import unittest

from verify import pmul, psub, res2, res_zero, s3, s3_coeffs, s4_in_u

P = 1009


def add(pt, qt, a, p):
    if pt is None:
        return qt
    if qt is None:
        return pt
    (x1, y1), (x2, y2) = pt, qt
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if pt == qt:
        slope = (3 * x1 * x1 + a) * pow(2 * y1, p - 2, p) % p
    else:
        slope = (y2 - y1) * pow(x2 - x1, p - 2, p) % p
    x3 = (slope * slope - x1 - x2) % p
    return (x3, (slope * (x1 - x3) - y1) % p)


def curves(rng, count):
    """Random curves over F_P with their points."""
    out = []
    while len(out) < count:
        a, b = rng.randrange(P), rng.randrange(1, P)
        if (4 * a**3 + 27 * b * b) % P == 0:
            continue
        pts = [(x, y) for x in range(P) for y in range(P) if (y * y - x**3 - a * x - b) % P == 0]
        out.append((a, b, pts))
    return out


def planted(rng, a, pts, m):
    """m points with distinct abscissae and the abscissa of their sum, or None."""
    while True:
        chosen = rng.sample(pts, m)
        if len({q[0] for q in chosen}) == m:
            break
    total = None
    for q in chosen:
        total = add(total, q, a, P)
    if total is None:
        return None
    return [q[0] for q in chosen], total[0]


def s4_other_order(x3, x4, xr, a, b, p):
    """Res_w(S3(U, x3, w), S3(w, x4, xR)) as a polynomial in U: the chain's
    other end eliminated first."""
    A, B, C = s3_coeffs(x4, xr, a, b, p)
    a2 = [x3 * x3 % p, (p - 2) * x3 % p, 1]
    b2 = [(p - 2) * ((a * x3 + 2 * b) % p) % p, (p - 2) * ((x3 * x3 + a) % p) % p, (p - 2) * x3 % p]
    c2 = [(a * a - 4 * b * x3) % p, (p - 2) * ((a * x3 + 2 * b) % p) % p, x3 * x3 % p]
    u = psub([A * c % p for c in c2], [C * c % p for c in a2], p)
    v = psub([A * c % p for c in b2], [B * c % p for c in a2], p)
    w = psub([B * c % p for c in c2], [C * c % p for c in b2], p)
    return psub(pmul(u, u, p), pmul(v, w, p), p)


class PlantedDecompositions(unittest.TestCase):
    def setUp(self):
        self.rng = random.Random(20260924)
        self.curves = curves(self.rng, 4)

    def test_m2_planted_pair_is_a_zero_of_s3(self):
        for a, b, pts in self.curves:
            got = planted(self.rng, a, pts, 2)
            if got is None:
                continue
            (x1, x2), xr = got
            self.assertEqual(s3(x1, x2, xr, a, b, P), 0)

    def test_m3_planted_triple_is_a_zero_of_s4(self):
        for a, b, pts in self.curves:
            got = planted(self.rng, a, pts, 3)
            if got is None:
                continue
            (x1, x2, x3), xr = got
            self.assertEqual(res2(s3_coeffs(x1, x2, a, b, P), s3_coeffs(x3, xr, a, b, P), P), 0)

    def test_m4_planted_quadruple_is_a_zero_of_s5(self):
        for a, b, pts in self.curves:
            got = planted(self.rng, a, pts, 4)
            if got is None:
                continue
            (x1, x2, x3, x4), xr = got
            A, B, C = s3_coeffs(x4, xr, a, b, P)
            self.assertTrue(res_zero(s4_in_u(x1, x2, x3, a, b, P), [C, B, A], P))

    def test_m4_both_elimination_orders_agree(self):
        zeros = 0
        for a, b, pts in self.curves[:2]:
            v = sorted({q[0] for q in self.rng.sample(pts, 7)})
            xr = self.rng.choice(pts)[0]
            for x1 in v:
                for x2 in v:
                    f = s3_coeffs(x1, x2, a, b, P)
                    for x3 in v:
                        s4 = s4_in_u(x1, x2, x3, a, b, P)
                        for x4 in v:
                            A, B, C = s3_coeffs(x4, xr, a, b, P)
                            one = res_zero(s4, [C, B, A], P)
                            other = res_zero([f[2], f[1], f[0]], s4_other_order(x3, x4, xr, a, b, P), P)
                            self.assertEqual(one, other, (x1, x2, x3, x4))
                            zeros += one
        # The grid is small against P, but not so small that no tuple solves it.
        self.assertGreater(zeros, 0)


if __name__ == "__main__":
    unittest.main()

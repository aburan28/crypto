import itertools, random, unittest

import pycryptosat

from encode import Curve, Lin, Model, _xor_clauses, build, decode, verify
from run import factor_space, targets


class T(unittest.TestCase):
    def test_s3_formula(self):
        rng = random.Random(1)
        for n in (13, 19):
            E, f = Curve(n), Curve(n).f
            for _ in range(50):
                P, Q = E.random_point(rng), E.random_point(rng)
                S = E.add(P, Q)
                if S is None:
                    continue
                x1, x2, x3 = P[0], Q[0], S[0]
                e = f.sq(f.mul(x1, x2) ^ f.mul(x1, x3) ^ f.mul(x2, x3)) \
                    ^ f.mul(f.mul(x1, x2), x3) ^ 1
                self.assertEqual(e, 0)

    def test_xor_expansion(self):
        for k in range(1, 7):
            for r in (0, 1):
                cl = _xor_clauses(list(range(1, k + 1)), r)
                for bits in itertools.product((0, 1), repeat=k):
                    ok = all(any((b if l > 0 else 1 - b) for l in c
                                 for b in [bits[abs(l) - 1]]) for c in cl)
                    self.assertEqual(ok, sum(bits) % 2 == r)

    def test_field_product_circuit(self):
        rng = random.Random(2)
        n = 13
        M = Model(n)
        x, y = M.elem(), M.elem()
        z = M.mul(x, y)
        zv = [M.as_var(a) for a in z]
        for _ in range(5):
            a, b = rng.getrandbits(n), rng.getrandbits(n)
            s = pycryptosat.Solver()
            M.load(s)
            assum = [(xx.v and next(iter(xx.v))) * (1 if (a >> k) & 1 else -1)
                     for k, xx in enumerate(x)]
            assum += [next(iter(yy.v)) * (1 if (b >> k) & 1 else -1)
                      for k, yy in enumerate(y)]
            sat, sol = s.solve(assum)
            self.assertTrue(sat)
            got = sum(sol[v] << k for k, v in enumerate(zv))
            self.assertEqual(got, M.f.mul(a, b))

    def test_planted_witness_satisfies_all_formulations(self):
        from run import FORMULATIONS
        for n, s in ((13, 4), (19, 6)):
            rng = random.Random(3 + n)
            E = Curve(n)
            basis, valid = factor_space(E, s, rng)
            for R, w in targets(E, valid, "planted", 2, rng):
                for form, (final, sub, nx) in FORMULATIONS.items():
                    M, cvars = build(n, basis, R, final, sub, nx)
                    # coordinates of each witness x in the basis
                    assum = []
                    for x, ci in zip(w, cvars):
                        coords = next(c for c in itertools.product((0, 1), repeat=s)
                                      if _span(basis, c) == x)
                        assum += [v if b else -v for v, b in zip(ci, coords)]
                    solver = pycryptosat.Solver()
                    M.load(solver)
                    sat, sol = solver.solve(assum)
                    self.assertTrue(sat, (n, form))
                    self.assertIsNotNone(verify(E, decode(sol, cvars, basis), R))


def _span(basis, coords):
    x = 0
    for b, c in zip(basis, coords):
        if c:
            x ^= b
    return x


if __name__ == "__main__":
    unittest.main()

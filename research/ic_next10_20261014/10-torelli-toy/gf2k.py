"""Small binary fields F_{2^k}, k <= 8, as integer tables (numpy)."""
import itertools
import numpy as np

MODULI = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101,
          6: 0b1000011, 7: 0b10000011, 8: 0b100011011}


def _mul(a, b, k, mod):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> k & 1:
            a ^= mod
    return r


class Field:
    def __init__(self, k):
        self.k, self.q = k, 1 << k
        q, mod = self.q, MODULI[k]
        if k == 1:
            self.mul = np.array([[0, 0], [0, 1]], dtype=np.int64)
        else:
            self.mul = np.array([[_mul(a, b, k, mod) for b in range(q)] for a in range(q)],
                                dtype=np.int64)
        self.inv = np.zeros(q, dtype=np.int64)
        for a in range(1, q):
            self.inv[a] = int(np.nonzero(self.mul[a] == 1)[0][0])
        # absolute trace to F_2
        tr = np.zeros(q, dtype=np.int64)
        for a in range(q):
            s, x = 0, a
            for _ in range(k):
                s ^= x
                x = int(self.mul[x, x])
            assert s in (0, 1)
            tr[a] = s
        self.trace = tr
        # powers table pw[a, e] for e <= 16
        self.pw = np.zeros((q, 25), dtype=np.int64)
        self.pw[:, 0] = 1
        for e in range(1, 25):
            self.pw[:, e] = self.mul[self.pw[:, e - 1], np.arange(q)]

    def proj_points(self, n):
        """Normalised points of P^n(F_q): first nonzero coordinate equal to 1."""
        pts = []
        q = self.q
        for lead in range(n + 1):
            for tail in itertools.product(range(q), repeat=n - lead):
                pts.append((0,) * lead + (1,) + tail)
        return np.array(pts, dtype=np.int64)

    def monomial_values(self, pts, exps):
        """V[P, m] = prod_i P_i^{e_mi} in F_q."""
        V = np.ones((len(pts), len(exps)), dtype=np.int64)
        for m, e in enumerate(exps):
            v = np.ones(len(pts), dtype=np.int64)
            for i, ei in enumerate(e):
                if ei:
                    v = self.mul[v, self.pw[pts[:, i], ei]]
            V[:, m] = v
        return V


def monomials(nvars, d):
    out = []
    def rec(prefix, left, slots):
        if slots == 1:
            out.append(tuple(prefix + [left]))
            return
        for a in range(left, -1, -1):
            rec(prefix + [a], left - a, slots - 1)
    rec([], d, nvars)
    return out

"""All smooth hyperelliptic curves y^2 + h y = f over F_2 of genus g
(deg h <= g+1, deg f <= 2g+2, weighted projective model), reduced modulo
y -> y + e(x): each representative stands for 2^{g+1} models.

Smoothness: exact gcd criterion at finite points and at infinity.
Output: dict {(N_1..N_g): number of models}.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from gf2k import Field


def pdeg(a):
    return a.bit_length() - 1


def pmulmod2(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def pgcd(a, b):
    while b:
        while a and pdeg(a) >= pdeg(b):
            a ^= b << (pdeg(a) - pdeg(b))
        a, b = b, a
    return a


def pderiv(a):
    # d/dx sum a_i x^i = sum_{i odd} a_i x^{i-1}
    r, i = 0, 1
    while (a >> i):
        if (a >> i) & 1:
            r |= 1 << (i - 1)
        i += 2
    return r


def reverse(a, n):
    r = 0
    for i in range(n + 1):
        if (a >> i) & 1:
            r |= 1 << (n - i)
    return r


def smooth(h, f, g):
    if h == 0:
        return False
    hp, fp = pderiv(h), pderiv(f)
    c = pmulmod2(fp, fp) ^ pmulmod2(pmulmod2(hp, hp), f)
    if pgcd(h, c) != 1:
        return False
    H, F = reverse(h, g + 1), reverse(f, 2 * g + 2)
    if H & 1:
        return True  # no point at infinity over the root of H at 0
    Hp, Fp = pderiv(H), pderiv(F)
    C = pmulmod2(Fp, Fp) ^ pmulmod2(pmulmod2(Hp, Hp), F)
    return (C & 1) == 1


def complement_reps(h, g):
    """Bit positions spanning a complement of V_h = {e^2 + h e} in deg<=2g+2."""
    nf = 2 * g + 3
    rows = []
    for i in range(g + 2):
        rows.append((1 << (2 * i)) ^ (h << i))
    # row-reduce over F2, pivot = highest bit
    piv = {}
    for r in rows:
        for p in sorted(piv, reverse=True):
            if (r >> p) & 1:
                r ^= piv[p]
        if r:
            p = pdeg(r)
            # reduce existing
            for q in list(piv):
                if (piv[q] >> p) & 1:
                    piv[q] ^= r
            piv[p] = r
    assert len(piv) == g + 1, (h, len(piv))
    return [i for i in range(nf) if i not in piv]


def enumerate_hyper(g, fields=None):
    fields = fields or [Field(k) for k in range(1, g + 1)]
    out = {}
    nf = 2 * g + 3
    total_smooth = 0
    # per field: powers of x for x in F_q, i = 0..2g+2
    xp = [F.pw[:, :nf] for F in fields]
    for h in range(1, 1 << (g + 2)):
        free = complement_reps(h, g)
        reps = []
        for m in range(1 << len(free)):
            f = 0
            for j, b in enumerate(free):
                if (m >> j) & 1:
                    f |= 1 << b
            reps.append(f)
        sm = [smooth(h, f, g) for f in reps]
        reps = [f for f, s in zip(reps, sm) if s]
        if not reps:
            continue
        fb = np.array([[(f >> i) & 1 for i in range(nf)] for f in reps], dtype=np.float32)
        hb = [(h >> i) & 1 for i in range(g + 2)]
        counts = np.zeros((len(reps), g), dtype=np.int64)
        for ki, F in enumerate(fields):
            q = F.q
            P = xp[ki]
            hx = np.zeros(q, dtype=np.int64)
            for i, b in enumerate(hb):
                if b:
                    hx ^= P[:, i]
            nz = hx != 0
            w = F.mul[F.inv[hx], F.inv[hx]]  # 1/h(x)^2 (garbage where h=0)
            Tm = F.trace[F.mul[P, w[:, None]]].astype(np.float32)  # q x nf
            par = (fb @ Tm.T) % 2  # reps x q
            cnt = (2 - 2 * par[:, nz]).sum(1) + (~nz).sum()
            # infinity: H(0) = coeff x^{g+1} of h, F(0) = coeff x^{2g+2} of f
            if (h >> (g + 1)) & 1:
                Finf = fb[:, 2 * g + 2]
                trv = (F.k * Finf) % 2
                cnt = cnt + 2 - 2 * trv
            else:
                cnt = cnt + 1
            counts[:, ki] = cnt.astype(np.int64)
        total_smooth += len(reps)
        for row in counts:
            key = tuple(int(x) for x in row)
            out[key] = out.get(key, 0) + (1 << (g + 1))
    return out, total_smooth * (1 << (g + 1))


if __name__ == "__main__":
    import sys
    g = int(sys.argv[1])
    out, tot = enumerate_hyper(g)
    print(g, tot, len(out))

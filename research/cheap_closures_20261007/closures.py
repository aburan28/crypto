"""Two cheap closures from the geometric-V idea screen (../geometric_v_linear_20261006).

Part I (identities, exhaustive over all points for small n) on y^2 + xy = x^3 + 1,
with tau the 2-Frobenius (x, y) -> (x^2, y^2):
  I1  x((1 + tau) P) = x + 1/x
  I2  x(P + T2) = 1/x for the 2-torsion point T2 = (0, 1)
  I3  (1 - tau)(P + T) = (1 - tau) P for every T in E(F_2)
  I4  lambda^2 + lambda = x(2P), lambda = x + y/x (Knudsen lambda-coordinate)
If they hold, the 2-torsion, E(F_2)-translation and lambda invariants of a factor
base are pullbacks by the separable endomorphisms 1 + tau, 1 - tau and [2]: the
same decomposition problem up to their kernels (at most 4 points), so at most
log2 4 = 2 bits per summand.

Part II (additive energy).  For a factor base V (geometric or random, dim l), let
S be the points whose abscissa lies in V (both signs).  The additive energy
E(S) = #{(a, b, c, d) in S^4 : a + b = c + d}.  Compared against a random point set
of the same size (closed under negation).  Trivial solutions number exactly
3|S|^2 - 3|S| for a negation-closed S; the nontrivial excess is compared with the
random expectation |S|^4 / #E.
"""
import json
import os
import random
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "linearization_reach_20260930"))
import lr  # noqa: E402

lr.IRR.setdefault(11, (1 << 11) | 0b101)
lr.IRR.setdefault(7, (1 << 7) | 0b11)


def all_points(F):
    pts = [None]
    for x in range(1 << F.n):
        if x == 0:
            pts.append((0, 1))      # y^2 = 1 -> y = 1 (the 2-torsion point)
            continue
        P = lr.lift(F, x)
        if P:
            pts.extend([P, lr.neg(P)])
    return pts


def frob(F, P):
    return None if P is None else (F.sq(P[0]), F.sq(P[1]))


def sub(F, P, Q):
    return lr.add(F, P, lr.neg(Q))


def identities(n):
    F = lr.Field(n)
    pts = all_points(F)
    EF2 = [None, (0, 1), (1, 0), (1, 1)]      # E(F_2) for y^2 + xy = x^3 + 1
    for T in EF2[1:]:
        assert F.sq(T[1]) ^ F.mul(T[0], T[1]) == F.mul(F.sq(T[0]), T[0]) ^ 1
    c = Counter()
    for P in pts:
        if P is None or P[0] == 0:
            continue
        x, y = P
        Q = lr.add(F, P, frob(F, P))
        c["I1_total"] += 1
        c["I1_ok"] += (Q is not None and Q[0] == x ^ F.inv(x))
        S = lr.add(F, P, (0, 1))
        c["I2_total"] += 1
        c["I2_ok"] += (S is not None and S[0] == F.inv(x))
        base = sub(F, P, frob(F, P))
        for T in EF2:
            PT = lr.add(F, P, T)
            c["I3_total"] += 1
            c["I3_ok"] += (sub(F, PT, frob(F, PT)) == base)
        D = lr.add(F, P, P)
        lam = x ^ F.mul(y, F.inv(x))
        c["I4_total"] += 1
        c["I4_ok"] += (D is not None and F.sq(lam) ^ lam == D[0])
    return dict(n=n, points=len(pts), **c)


def point_set_from_V(F, V):
    span = [0]
    for v in V:
        span = span + [s ^ v for s in span]
    S = []
    for x in span:
        if x == 0:
            continue
        P = lr.lift(F, x)
        if P:
            S.extend([P, lr.neg(P)])
    return S


def energy(F, S):
    r = Counter()
    m = len(S)
    for i in range(m):
        a = S[i]
        r[a and lr.add(F, a, a)] += 1                     # ordered (a, a)
        for j in range(i + 1, m):
            r[lr.add(F, a, S[j])] += 2                    # (a, b) and (b, a)
    return sum(v * v for v in r.values())


def base(F, l, kind, rng):
    n = F.n
    if kind == "geometric":
        theta, g = rng.getrandbits(n) or 1, rng.getrandbits(n) or 2
        V, p = [], theta
        for _ in range(l):
            V.append(p)
            p = F.mul(p, g)
        return V
    return [rng.getrandbits(n) for _ in range(l)]


def energy_cell(n, l, seeds, rng):
    """Per seed: a geometric set, a random-V set, and a random point set of exactly
    the geometric set's size (closed under negation).  Nontrivial energy is
    E - (3m^2 - 3m): solutions other than (a,b)=(c,d), (a,b)=(d,c) and a+b=O=c+d."""
    F = lr.Field(n)
    pts = [P for P in all_points(F) if P is not None]
    order = len(pts) + 1
    xs_all = sorted({P[0] for P in pts if P[0] != 0})
    out = []
    for s in range(seeds):
        Sg = point_set_from_V(F, base(F, l, "geometric", rng))
        Sr = point_set_from_V(F, base(F, l, "randomV", rng))
        Sx = []
        for x in rng.sample(xs_all, len(Sg) // 2):
            P = lr.lift(F, x)
            Sx.extend([P, lr.neg(P)])
        for kind, S in (("geometric", Sg), ("randomV", Sr), ("randomset", Sx)):
            E = energy(F, S)
            m = len(S)
            nontriv = E - (3 * m * m - 3 * m)
            exp_nt = m ** 4 / order
            out.append(dict(n=n, l=l, kind=kind, seed=s, size=m, energy=E, nontrivial=nontriv,
                            expected_nontrivial=round(exp_nt, 1),
                            ratio_nontrivial=round(nontriv / exp_nt, 4) if exp_nt else None))
            print(json.dumps(out[-1]), flush=True)
    return out


def main(out_dir, seed):
    rng = random.Random(seed)
    ids = [identities(n) for n in (7, 11, 13)]
    for r in ids:
        print(json.dumps(r), flush=True)
    cells = []
    for n, ls in ((13, (6,)), (17, (7, 8)), (19, (8, 9)), (23, (10, 11))):
        for l in ls:
            cells += energy_cell(n, l, 3, rng)
    with open(f"{out_dir}/closures.json", "w") as fh:
        json.dump(dict(seed=seed, identities=ids, energy=cells), fh, indent=1)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))

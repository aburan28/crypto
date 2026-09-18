"""The torsion-point route: an elliptic curve E'/F_p as the auxiliary group.

Kim–Cheon (Math. Comp. 2014; survey §5) replace the multiplicative group
F_p^* of Cheon's p − 1 case by an auxiliary elliptic curve E' over the
*scalar* field F_p, with #E'(F_p) = d·δ.  A point P with x(P) = rα is
mapped by the multiplication-by-d map into the subgroup H of order δ,
using the division polynomials:

    x([d]P) = φ_d(rα) / ψ_d²(rα),       deg φ_d = d²,  deg ψ_d² = d² − 1.

So `d²` auxiliary inputs G_i = [α^i]G let one compute both
[φ_d(rα)]G and [ψ_d²(rα)]G.  A birthday collision between
{ x(P̃) : P̃ ∈ H } and { φ_d(r_j α)/ψ_d²(r_j α) } would take m ≈ √δ
samples on each side, i.e. √(p/d): the DLPwAI floor again, but now for
*any* d dividing some #E'(F_p) in the Hasse interval, not only divisors
of p ± 1.  That is the promise.

The obstacle (survey §5: "computing list L2 seems hard"): the quotient
cannot be formed in the exponent.  One can only test the cross-multiplied
equality  [x(P̃_i)]·[ψ_d²(r_j α)]G  =  [φ_d(r_j α)]G,  and that is one
scalar multiplication *per pair* (i, j): m² of them, i.e. δ = p/d, not
√(p/d).  This script measures exactly that, honestly, and also runs an
"oracle" version that is handed the quotient for free, to show that the
√δ collision count is real and the exponent division is the only wall.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
from typing import Dict, List, Optional, Set, Tuple

import ec
from ec import Counter, Curve, Point

HERE = os.path.dirname(os.path.abspath(__file__))


# -- polynomial helpers over F_p ---------------------------------------------


def pmul(f: List[int], g: List[int], p: int) -> List[int]:
    out = [0] * (len(f) + len(g) - 1)
    for i, fi in enumerate(f):
        if fi:
            for j, gj in enumerate(g):
                out[i + j] = (out[i + j] + fi * gj) % p
    return out


def padd(f: List[int], g: List[int], p: int) -> List[int]:
    m = max(len(f), len(g))
    return [((f[i] if i < len(f) else 0) + (g[i] if i < len(g) else 0)) % p for i in range(m)]


def psub(f: List[int], g: List[int], p: int) -> List[int]:
    return padd(f, [(-x) % p for x in g], p)


def psc(f: List[int], c: int, p: int) -> List[int]:
    return [c * x % p for x in f]


def pdiv_exact(f: List[int], g: List[int], p: int) -> List[int]:
    """f / g, asserting zero remainder."""
    f = list(f)
    while f and f[-1] == 0:
        f.pop()
    while g and g[-1] == 0:
        g = g[:-1]
    if not f:
        return [0]
    inv = pow(g[-1], -1, p)
    q = [0] * (len(f) - len(g) + 1)
    for i in range(len(f) - len(g), -1, -1):
        c = f[i + len(g) - 1] * inv % p
        q[i] = c
        if c:
            for j, gj in enumerate(g):
                f[i + j] = (f[i + j] - c * gj) % p
    assert all(x == 0 for x in f), "inexact division"
    return q


def evalpoly(f: List[int], x: int, p: int) -> int:
    r = 0
    for c in reversed(f):
        r = (r * x + c) % p
    return r


# -- division polynomials -----------------------------------------------------
# ψ_n is stored as (f, e): the element f(x)·y^e of F_p[x, y]/(y² − F),
# e ∈ {0, 1}, F = x³ + a x + b.  Odd n have e = 0, even n have e = 1.


def divpolys(p: int, a: int, b: int, n: int):
    F = [b % p, a % p, 0, 1]

    def mul(u, v):
        f, e = pmul(u[0], v[0], p), u[1] + v[1]
        if e >= 2:
            f, e = pmul(f, F, p), e - 2
        return (f, e)

    def sub(u, v):
        assert u[1] == v[1]
        return (psub(u[0], v[0], p), u[1])

    psi = [None] * (n + 3)
    psi[0] = ([0], 0)
    psi[1] = ([1], 0)
    psi[2] = ([2], 1)
    psi[3] = ([(-a * a) % p, 12 * b % p, 6 * a % p, 0, 3], 0)
    psi[4] = (psc([(-8 * b * b - a**3) % p, (-4 * a * b) % p, (-5 * a * a) % p, 20 * b % p,
                   5 * a % p, 0, 1], 4, p), 1)
    for k in range(5, n + 3):
        m = k // 2
        if k % 2 == 1:
            t1 = mul(psi[m + 2], mul(psi[m], mul(psi[m], psi[m])))
            t2 = mul(psi[m - 1], mul(psi[m + 1], mul(psi[m + 1], psi[m + 1])))
            psi[k] = sub(t1, t2)
            assert psi[k][1] == 0
        else:
            t1 = mul(psi[m + 2], mul(psi[m - 1], psi[m - 1]))
            t2 = mul(psi[m - 2], mul(psi[m + 1], psi[m + 1]))
            num = mul(psi[m], sub(t1, t2))
            # divide by 2y: y^e / y = y^{e-1}; e = 0 means the polynomial
            # carries an F = y² factor to borrow from.
            f, e = num
            if e == 1:
                psi[k] = (psc(f, pow(2, -1, p), p), 0)
            else:
                psi[k] = (psc(pdiv_exact(f, F, p), pow(2, -1, p), p), 1)
            assert psi[k][1] == 1
    return psi, F


def x_of_dP_polys(p: int, a: int, b: int, d: int) -> Tuple[List[int], List[int]]:
    """φ_d, ψ_d² ∈ F_p[x] with x([d]P) = φ_d(x)/ψ_d²(x)."""
    psi, F = divpolys(p, a, b, d + 1)

    def sq(u):  # (f y^e)² = f² F^e
        f = pmul(u[0], u[0], p)
        return pmul(f, F, p) if u[1] else f

    def mul0(u, v):  # product that must land in F_p[x]
        f, e = pmul(u[0], v[0], p), u[1] + v[1]
        assert e % 2 == 0
        return pmul(f, F, p) if e == 2 else f

    psi2 = sq(psi[d])
    phi = psub(pmul([0, 1], psi2, p), mul0(psi[d - 1], psi[d + 1]), p)
    return phi, psi2


# -- the experiment ----------------------------------------------------------


def _multi_scalar(E: Curve, coeffs: List[int], pts: List[Point]) -> Point:
    R: Point = None
    for c, P in zip(coeffs, pts):
        if c:
            R = E.add(R, E.mul(c, P))
    return R


def find_aux_curve(p: int, d_target: int, rng: random.Random):
    """E'/F_p with #E' = d·δ, d near d_target, gcd(d, δ) = 1, and H = [d]E'
    cyclic of order δ (a generator is returned): (E', #E', d, Hgen)."""
    for _ in range(100000):
        a, b = rng.randrange(p), rng.randrange(p)
        if (4 * a**3 + 27 * b**2) % p == 0:
            continue
        E = Curve(p, a, b)
        P = E.random_point(rng)
        ks = ec.point_order_candidates(E, P)
        if len(ks) != 1:
            continue
        n = ks[0]
        fac = ec.factorise(n)
        d = ec.best_divisor(fac, d_target)
        if not (d_target // 2 <= d <= 2 * d_target and d >= 3 and math.gcd(d, n // d) == 1):
            continue
        delta = n // d
        delta_fac = ec.factorise(delta)
        for _ in range(8):
            Hgen = E.mul(d, E.random_point(rng))
            if Hgen is not None and all(E.mul(delta // q, Hgen) is not None for q in delta_fac):
                E.ctr = Counter()
                return E, n, d, Hgen
    raise RuntimeError("no aux curve")


def d_torsion_points(Eaux: Curve, n_aux: int, d: int, rng: random.Random) -> Set[Point]:
    """E'(F_p)[d] by closure of a few [δ]R samples (small groups only)."""
    delta = n_aux // d
    gens = set()
    for _ in range(40):
        T = Eaux.mul(delta, Eaux.random_point(rng))
        if T is not None:
            gens.add(T)
    group: Set[Point] = {None}
    frontier = [None]
    while frontier:
        nxt = []
        for X in frontier:
            for T in gens:
                Y = Eaux.add(X, T)
                if Y not in group:
                    group.add(Y)
                    nxt.append(Y)
        frontier = nxt
    return group


def run(bits: int, seed: int, d_exp: float = 0.25) -> Dict:
    rng = random.Random(seed)
    inst = ec.find_curve(bits, seed=5000 + bits, want="p-1", exponent=0.5, tolerance=10)
    p, G = inst.p, inst.G
    alpha = rng.randrange(2, p - 1)
    Eaux, n_aux, d, Hgen = find_aux_curve(p, max(3, int(round(p**d_exp))), rng)
    delta = n_aux // d
    phi, psi2 = x_of_dP_polys(p, Eaux.a, Eaux.b, d)
    deg = max(len(phi), len(psi2)) - 1
    # sanity: x([d]P) = φ/ψ² on a random point
    P = Eaux.random_point(rng)
    dP = Eaux.mul(d, P)
    assert dP is not None
    assert evalpoly(phi, P[0], p) * pow(evalpoly(psi2, P[0], p), -1, p) % p == dP[0]

    # auxiliary inputs G_i = [α^i]G, i ≤ deg (given by the problem, not charged)
    E0 = inst.curve()
    G_pows = [G]
    for _ in range(deg):
        G_pows.append(E0.mul(alpha, G_pows[-1]))

    m = math.isqrt(delta) + 1
    out = {"bits": bits, "p": p, "aux_curve": {"a": Eaux.a, "b": Eaux.b, "order": n_aux},
           "d": d, "delta": delta, "aux_inputs_needed": deg + 1, "m": m, "alpha": alpha,
           "sqrt_p": math.sqrt(p), "floor_dlpwai": math.sqrt(p / d), "rows": []}

    # ---- honest version: cross-multiplied pairwise test ------------------
    ctr = Counter()
    E = inst.curve(ctr)
    Eaux.ctr = ctr  # aux-curve work is charged too
    ctr.begin("torsion:L1")
    L1: List[Tuple[int, int]] = []          # (x(P̃), k) with P̃ = [k]Hgen
    # L1 is the cheap side (one aux-curve multiplication per entry); take
    # 3m entries so the expected number of matching pairs is ≈ 3.
    for _ in range(3 * m):
        k = rng.randrange(1, delta)
        L1.append((Eaux.mul(k, Hgen)[0], k))
    ctr.end()
    ctr.begin("torsion:L2")
    L2 = []
    for _ in range(m):
        r = rng.randrange(1, p)
        cphi = [c * pow(r, i, p) % p for i, c in enumerate(phi)]
        cpsi = [c * pow(r, i, p) % p for i, c in enumerate(psi2)]
        L2.append((r, _multi_scalar(E, cphi, G_pows), _multi_scalar(E, cpsi, G_pows)))
    ctr.end()
    ctr.begin("torsion:pairwise")
    found = None
    pairs = 0
    for i, (xi, k) in enumerate(L1):
        for j, (r, Gphi, Gpsi) in enumerate(L2):
            pairs += 1
            if E.mul(xi, Gpsi) == Gphi:
                found = (i, j)
                break
        if found:
            break
    ctr.end()
    alpha_rec = None
    if found:
        i, j = found
        ctr.begin("torsion:recover")
        xi, k = L1[i]
        r = L2[j][0]
        # [d]P = ±P̃ with P̃ = [k]Hgen; one preimage is [k d^{-1} mod δ]Hgen,
        # the rest differ by E'(F_p)[d].
        P0 = Eaux.mul(k * pow(d, -1, delta) % delta, Hgen)
        rinv = pow(r, -1, p)
        for T in d_torsion_points(Eaux, n_aux, d, rng):
            for Q in (Eaux.add(P0, T), Eaux.add(Eaux.neg(P0), T)):
                if Q is None:
                    continue
                cand = Q[0] * rinv % p
                if E.mul(cand, G) == G_pows[1]:
                    alpha_rec = cand
                    break
            if alpha_rec is not None:
                break
        ctr.end()
    snap = ctr.snapshot()
    out["rows"].append({"variant": "torsion_embedding_honest", "correct": alpha_rec == alpha,
                        "ops": snap["ops"], "phases": snap["phases"],
                        "S": snap["ops"] / math.sqrt(p),
                        "ops_over_floor": snap["ops"] / math.sqrt(p / d),
                        "pairs_tested": pairs, "m": m, "L1": 3 * m, "L2": m})

    # ---- oracle version: quotient handed over for free -------------------
    ctr = Counter()
    E = inst.curve(ctr)
    Eaux.ctr = ctr
    ctr.begin("oracle:L1")
    T1: Dict[Point, int] = {}
    for _ in range(m):
        k = rng.randrange(1, delta)
        x = Eaux.mul(k, Hgen)[0]
        T1[E.mul(x, G)] = x
    ctr.end()
    ctr.begin("oracle:L2")
    hit = None
    tries = 0
    for _ in range(6 * m):
        tries += 1
        r = rng.randrange(1, p)
        x = r * alpha % p                     # the oracle: uses α to form the quotient
        den = evalpoly(psi2, x, p)
        if den == 0:
            continue
        q = evalpoly(phi, x, p) * pow(den, -1, p) % p
        if E.mul(q, G) in T1:
            hit = r
            break
    ctr.end()
    snap = ctr.snapshot()
    out["rows"].append({"variant": "torsion_embedding_oracle_quotient", "collision": hit is not None,
                        "samples_to_collision": tries, "ops": snap["ops"], "phases": snap["phases"],
                        "S": snap["ops"] / math.sqrt(p),
                        "ops_over_floor": snap["ops"] / math.sqrt(p / d)})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bits", type=int, nargs="+", default=[16, 18, 20])
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for bits in args.bits:
        res = []
        for s in range(args.seeds):
            r = run(bits, s)
            res.append(r)
            for row in r["rows"]:
                print(f"  {bits}b seed {s} d={r['d']} delta={r['delta']} inputs={r['aux_inputs_needed']} "
                      f"m={r['m']} {row['variant']:34s} ops={row['ops']:>10,d} S={row['S']:9.2f} "
                      f"floor*={row['ops_over_floor']:8.1f} "
                      f"{'OK' if row.get('correct') else ('collision' if row.get('collision') else 'FAIL')}",
                      flush=True)
        path = os.path.join(args.out, f"torsion_{bits:02d}.json")
        with open(path, "w") as fh:
            json.dump(res, fh, indent=1)
        print("wrote", path)


if __name__ == "__main__":
    main()

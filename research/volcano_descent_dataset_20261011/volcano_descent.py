#!/usr/bin/env python3
"""Isogeny-volcano descent over prime fields, built on PARI (cypari).

Algorithm (Kohel 1996; Fouquet-Morain 2002; Sutherland 2013, "Isogeny volcanoes")
-----------------------------------------------------------------------------
For an ordinary E / F_p with Frobenius trace t, write t^2 - 4p = D_K * f^2 with
D_K a fundamental discriminant.  For a prime ell != p the rational ell-isogeny
graph component of E is a volcano of depth k = v_ell(f):

* level i (distance from the crater) holds curves with End(E) of conductor
  c * ell^i, where c = f / ell^k;
* a curve at level i with 0 < i <= k has ell + 1 rational ell-isogenies,
  one ascending and ell descending;
* a floor curve (i = k) has exactly one rational ell-isogeny (the ascending one);
* a crater curve (i = 0) has ell + 1 - (D_c/ell) descending and 1 + (D_c/ell)
  horizontal isogenies (D_c = D_K c^2), all rational when k >= 1.

Neighbours of j are the roots in F_p of the classical modular polynomial
Phi_ell(X, j).  Two procedures are implemented:

``explore_from_crater``  exact: multi-source BFS from the known crater set
    (roots of the Hilbert class polynomial of D_c mod p).  BFS distance to the
    crater set equals the level, so every vertex is labelled with ground truth.

``descend_to_floor``     blind: needs only a curve E, its trace, and ell.  From
    the start vertex it classifies neighbours by one non-backtracking walk each
    (a descending walk reaches a degree-1 vertex in exactly k - i steps, any
    other walk takes at least k - i + 2), takes the shortest, then keeps
    stepping to any neighbour other than the one it came from until the vertex
    has no other neighbour.  The level of the start vertex is k minus the
    length of the shortest walk, so the End conductor is read off exactly.

Everything is exact integer arithmetic; PARI supplies polynomial factoring,
class polynomials and point counting only.
"""
from __future__ import annotations

import math
import random
from collections import deque
from dataclasses import dataclass, field

from cypari import pari

pari.allocatemem(2 * 10**9)

# ---------------------------------------------------------------- PARI helpers
pari(
    r"""
nbrs(PHI, j, p) = {
  my(P = subst(PHI, 'y, j), F = factormod(P, p), R = List(), tot = 0);
  for(i = 1, #F~,
    my(g = lift(F[i,1]), e = F[i,2], d = poldegree(g));
    if(d == 1,
      my(r = (-polcoef(g,0) * lift(1/Mod(polcoef(g,1), p))) % p);
      listput(R, [r, e]);
      tot += e));
  [Vec(R), tot];
}
"""
)
_nbrs = pari("nbrs")
_phi_cache: dict[int, object] = {}


def phi(ell: int):
    """Classical modular polynomial Phi_ell(x, y) over Z (cached)."""
    if ell not in _phi_cache:
        _phi_cache[ell] = pari.polmodular(ell)
    return _phi_cache[ell]


def neighbors(ell: int, j: int, p: int):
    """Rational ell-neighbours of j: ({root: multiplicity}, rational multiplicity)."""
    res, tot = _nbrs(phi(ell), j, p)
    return {int(r): int(e) for r, e in res}, int(tot)


# ---------------------------------------------------------- CM parameter algebra
def omega_data(D: int):
    """omega = (D + sqrt D)/2:  Tr = D,  Nm = (D^2 - D)/4."""
    return D, (D * D - D) // 4


def frobenius_from(D: int, f: int, t: int):
    """Return (p, x) for pi = x + f*omega of trace t, given parity condition."""
    assert (t - f * D) % 2 == 0
    x = (t - f * D) // 2
    p = x * x + x * f * D + f * f * ((D * D - D) // 4)
    assert 4 * p == t * t - D * f * f
    return p, x


def class_number(D: int) -> int:
    return int(pari.qfbclassno(D))


def v_ell(n: int, ell: int) -> int:
    v = 0
    while n % ell == 0:
        n //= ell
        v += 1
    return v


def decompose(p: int, t: int):
    """t^2 - 4p = D_K f^2  ->  (D_K, f)."""
    delta = t * t - 4 * p
    assert delta < 0
    fac = pari.coredisc(delta, 1)  # [D_K, f]
    return int(fac[0]), int(fac[1])


def predicted_group(p: int, t: int, D: int, d: int):
    """Structure of E(F_p) when End(E) has conductor d: Z/n1 x Z/(N/n1).

    E(F_p) = End(E)/(pi - 1) (Lenstra); n1 = gcd(A, f/d) with pi - 1 = A + f*omega.
    """
    Dk, f = decompose(p, t)
    assert Dk == D and f % d == 0
    x = (t - f * D) // 2
    A = x - 1
    n1 = math.gcd(A, f // d)
    N = p + 1 - t
    assert N % n1 == 0 and (N // n1) % n1 == 0
    return n1, N // n1


@dataclass
class Instance:
    """A tall-volcano instance: curves with Frobenius pi = x + f*omega over F_p."""
    ell: int
    D: int            # fundamental discriminant D_K
    c: int            # conductor of the crater order (gcd(c, ell) = 1)
    k: int            # depth = v_ell(f)
    m: int            # v_ell(A) clipped to k: how far the rational ell^m torsion extends
    p: int
    t: int
    f: int = field(init=False)
    N: int = field(init=False)

    def __post_init__(self):
        self.f = self.c * self.ell ** self.k
        self.N = self.p + 1 - self.t


def find_instance(ell, D, c, k, m, bits, rng, max_tries=200000, big_prime_slack=None):
    """Random (p, t) with t^2 - 4p = D (c ell^k)^2, v_ell(A) >= m, p ~ 2^bits."""
    f = c * ell ** k
    assert math.gcd(c, ell) == 1 and 0 <= m <= k
    if 2 ** (bits + 1) <= -D * f * f:
        raise ValueError("bits too small for this depth")
    t_hi = math.isqrt(4 * 2 ** bits + D * f * f)
    t_lo = math.isqrt(max(4 * 2 ** (bits - 1) + D * f * f, 1)) + 1
    # A = (t - f D)/2 - 1 == 0 mod ell^m  <=>  t == f D + 2 mod 2 ell^m
    mod = 2 * ell ** m
    r0 = (f * D + 2) % mod
    for _ in range(max_tries):
        t = rng.randrange(t_lo, t_hi)
        t += (r0 - t) % mod
        if t >= t_hi:
            continue
        if (t - f * D) % 2:
            continue
        p = (t * t - D * f * f) // 4
        if p < 5 or p % ell == 0 or t % p == 0:
            continue
        if not pari.isprime(p):
            continue
        N = p + 1 - t
        if big_prime_slack is not None:
            # require N = (smooth part < 2^slack) * (one large prime)
            n = N
            for q in pari.primes(200):
                q = int(q)
                while n % q == 0:
                    n //= q
            if not (n == 1 or pari.isprime(n)):
                continue
            if N // max(n, 1) >= 2 ** big_prime_slack:
                continue
        return Instance(ell, D, c, k, m, p, t)
    raise RuntimeError("no instance found")


# --------------------------------------------------------------- curve models
def curve_for_j(j: int, p: int, t: int, rng=None):
    """Short-Weierstrass (a, b) over F_p with j-invariant j and trace t."""
    target = p + 1 - t
    rng = rng or random.Random(0)
    if j % p not in (0, 1728 % p):
        a = (3 * j * (1728 - j)) % p
        b = (2 * j * (1728 - j) ** 2) % p
        cands = [(a, b)]
        while True:
            d = rng.randrange(2, p)
            if pow(d, (p - 1) // 2, p) == p - 1:
                break
        cands.append(((a * d * d) % p, (b * d ** 3) % p))
    else:
        cands = []
        for s in range(1, 4000):
            cands.append((0, s) if j % p == 0 else (s, 0))
    for a, b in cands:
        if (4 * a ** 3 + 27 * b * b) % p == 0:
            continue
        E = pari.ellinit([a, b], p)
        if int(pari.ellcard(E)) == target:
            return a, b
    raise RuntimeError("no model with the requested trace")


def j_of(a: int, b: int, p: int) -> int:
    return int(pari.ellj(pari.ellinit([a, b], p)).lift())


# ------------------------------------------------------- exact BFS from crater
@dataclass
class Volcano:
    inst: Instance
    level: dict           # j -> level (distance to crater set)
    edges: list           # (src_j, dst_j, multiplicity, rational_total_at_src)
    comp: dict            # j -> component id (ell-isogeny component)
    totals: dict          # j -> total rational multiplicity of Phi_ell(X, j)


def crater_js(inst: Instance):
    Dc = inst.D * inst.c ** 2
    roots = pari.polrootsmod(pari.polclass(Dc), inst.p)
    js = sorted(int(r.lift()) for r in roots)
    assert len(js) == class_number(Dc), (len(js), class_number(Dc))
    return js


def explore_from_crater(inst: Instance, max_vertices=2_000_000):
    p, ell = inst.p, inst.ell
    level, totals, edges = {}, {}, []
    dq = deque()
    for j in crater_js(inst):
        level[j] = 0
        dq.append(j)
    uf = {j: j for j in level}

    def find(x):
        while uf[x] != x:
            uf[x] = uf[uf[x]]
            x = uf[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            uf[ra] = rb

    while dq:
        j = dq.popleft()
        nb, tot = neighbors(ell, j, p)
        totals[j] = tot
        for r, e in nb.items():
            edges.append((j, r, e, tot))
            if r not in level:
                level[r] = level[j] + 1
                uf[r] = r
                dq.append(r)
                if len(level) > max_vertices:
                    raise RuntimeError("vertex cap exceeded")
            union(j, r)
    comp = {j: find(j) for j in level}
    return Volcano(inst, level, edges, comp, totals)


def level_counts(v: Volcano):
    out = {}
    for lv in v.level.values():
        out[lv] = out.get(lv, 0) + 1
    return dict(sorted(out.items()))


def expected_level_counts(inst: Instance):
    return {i: class_number(inst.D * (inst.c * inst.ell ** i) ** 2) for i in range(inst.k + 1)}


# --------------------------------------------------------------- blind descent
def _avail(nb, prev):
    """Kernels available for a non-backtracking step: the root `prev` loses one
    multiplicity (the dual of the arriving isogeny); the rest are usable.  At
    j = 0, 1728 several kernels share one j-invariant, so multiplicities matter."""
    return {r: (e - 1 if r == prev else e) for r, e in nb.items() if (e - 1 if r == prev else e) > 0}


def _walk_len(ell, p, start, prev, rng, cap):
    """Length of one random non-backtracking walk from `start` (reached from
    `prev`) until a vertex with no other kernel; None if longer than `cap`."""
    cur, pv = start, prev
    for n in range(cap + 1):
        nb, _ = neighbors(ell, cur, p)
        av = _avail(nb, pv)
        if not av:
            return n
        roots = list(av)
        weights = [av[r] for r in roots]
        pv, cur = cur, rng.choices(roots, weights)[0]
    return None


def descend_to_floor(j0: int, p: int, t: int, ell: int, rng=None, walks=3):
    """Blind descent.  Returns dict(path=[j0, ...], floor_j, start_level, k)."""
    rng = rng or random.Random(1)
    D, f = decompose(p, t)
    k = v_ell(f, ell)
    if k == 0:
        return dict(path=[j0], floor_j=j0, start_level=0, k=0)
    nb, tot = neighbors(ell, j0, p)
    path = [j0]
    if tot == 1:
        return dict(path=path, floor_j=j0, start_level=k, k=k)     # degree 1: the floor
    best, best_len = None, None
    for r, e in nb.items():
        if r == j0:      # horizontal loop(s)
            continue
        ls = [_walk_len(ell, p, r, j0, rng, 2 * k + 3) for _ in range(walks)]
        ls = [x for x in ls if x is not None]
        if not ls:
            continue
        m = min(ls) + 1
        if best_len is None or m < best_len:
            best, best_len = r, m
    assert best is not None, "no descending neighbour found"
    start_level = k - best_len
    prev, cur = j0, best
    path.append(cur)
    while True:
        nb, _ = neighbors(ell, cur, p)
        av = _avail(nb, prev)
        if not av:
            break
        prev, cur = cur, next(iter(av))
        path.append(cur)
    return dict(path=path, floor_j=cur, start_level=start_level, k=k)


# ------------------------------------------------- fast twist selection (no SEA)
def curve_for_j_fast(j: int, p: int, t: int, rng):
    """As curve_for_j but picks the twist with random-point order tests (a few
    scalar multiplications instead of a point count).  Falls back to ellcard."""
    if j % p in (0, 1728 % p):
        return curve_for_j(j, p, t, rng)
    N, N2 = p + 1 - t, p + 1 + t
    a = (3 * j * (1728 - j)) % p
    b = (2 * j * (1728 - j) ** 2) % p
    E = pari.ellinit([a, b], p)
    ok1 = ok2 = True
    for _ in range(4):
        P = list(rand_point(E, p, rng))
        if pari.ellmul(E, P, N) != [0]:
            ok1 = False
        if pari.ellmul(E, P, N2) != [0]:
            ok2 = False
    if ok1 and not ok2:
        return a, b
    if ok2 and not ok1:
        while True:
            d = rng.randrange(2, p)
            if pow(d, (p - 1) // 2, p) == p - 1:
                break
        return (a * d * d) % p, (b * d ** 3) % p
    return curve_for_j(j, p, t, rng)


# ------------------------------------------- deterministic points / group structure
def rand_point(E, p: int, rng):
    """A random affine point of E drawn from `rng` (no PARI RNG involved)."""
    while True:
        x = rng.randrange(p)
        ys = pari.ellordinate(E, x)
        if len(ys):
            y = int(ys[rng.randrange(len(ys))].lift()) if hasattr(ys[0], "lift") else int(ys[0])
            return x, y


def group_structure(E, p: int, N: int, Nfac, rng, npoints=24):
    """(n1, n2) with E(F_p) = Z/n1 x Z/n2, n1 | n2, from the exponent of 24 random points."""
    e = 1
    for _ in range(npoints):
        x, y = rand_point(E, p, rng)
        o = int(pari.ellorder(E, [x, y], Nfac))
        e = e * o // math.gcd(e, o)
    n2 = e
    assert N % n2 == 0
    n1 = N // n2
    assert n2 % n1 == 0, "exponent too small or group not 2-generated"
    return n1, n2

"""Cheon's algorithms for the discrete logarithm with auxiliary inputs,
with every group operation counted.

Problem (DLPwAI).  Given G of prime order p and the auxiliary inputs
G_i = [α^i]G, recover α.  Cheon (EUROCRYPT 2006, J. Cryptology 2010):

* p − 1 case: from G, G_1, G_d with d | p − 1, in O(√((p−1)/d) + √d)
  *exponentiations*, by embedding α^d into the subgroup of F_p^* of order
  (p − 1)/d.
* p + 1 case: from G, G_1, …, G_{2d} with d | p + 1, in O(√((p+1)/d) + d)
  exponentiations, by embedding ((α − θ)/(α + θ))^d, θ² = a non-residue,
  into the subgroup of order (p + 1)/d of the norm-one torus of F_{p²}^*.

Both are matched by the generic lower bound Ω(√(p/d)) for the DLPwAI.

Every "exponentiation" is a scalar multiplication of a *fixed* base (G,
G_d or a derived point) by a scalar known in F_p, so it can use a
fixed-base comb table; that table is charged.  Three variants of the
p − 1 case are exposed so that the price of the exponentiations is visible:

  bsgs_naive   every step is a full double-and-add scalar multiplication
  bsgs_comb    fixed-base comb tables (Kozaki–Kutsuma–Matsuo style)
  kangaroo     Cheon's memoryless variant, comb tables, distinguished points

and the p + 1 case with comb tables.  A ``plain_bsgs`` row reproduces what
the repository's ``cheon_attack.rs`` actually computes (baby steps of
length d on G, giant steps of length (p−1)/d), which ignores the
auxiliary input and is not Cheon's algorithm.
"""

from __future__ import annotations

import math
import random
from typing import Dict, List, Optional, Tuple

from ec import Counter, Curve, FixedBase, Point, factorise, primitive_root, sqrt_mod


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _comb_width(nbits: int, steps: int) -> int:
    """Pick the comb width minimising table + per-step cost for `steps`
    multiplications of `nbits`-bit scalars.  Returns w in 2..12."""
    best = None
    for w in range(2, 13):
        chunks = -(-nbits // w)
        cost = chunks * ((1 << w) - 1) + steps * (chunks - 1)
        if best is None or cost < best[0]:
            best = (cost, w)
    return best[1]


class Scalar:
    """Fixed-base multiplier: either naive double-and-add or a comb."""

    def __init__(self, E: Curve, P: Point, nbits: int, steps: int, comb: bool):
        self.E, self.P, self.comb = E, P, None
        if comb:
            self.comb = FixedBase(E, P, nbits, _comb_width(nbits, steps))

    def mul(self, k: int) -> Point:
        if self.comb is not None:
            return self.comb.mul(k)
        return self.E.mul(k, self.P)

    @property
    def table_size(self) -> int:
        return self.comb.size if self.comb else 0


# ---------------------------------------------------------------------------
# what cheon_attack.rs does: BSGS with baby length d, ignoring G_d
# ---------------------------------------------------------------------------


def plain_bsgs(E: Curve, G: Point, G1: Point, p: int, d: int) -> Dict:
    """Baby steps [i]G for i < d, giant stride [d]G, up to (p−1)/d giants.
    Cost d + (p−1)/d group operations.  This is the algorithm in the
    repository's cheon_attack.rs (the `d2_g` argument is unused there)."""
    ctr = E.ctr
    start = ctr.ops
    ctr.begin("plain:baby")
    table: Dict[Point, int] = {}
    R: Point = None
    for i in range(d):
        table[R] = i
        if R == G1:
            ctr.end()
            return {"alpha": i, "ops": ctr.ops - start}
        R = E.add(R, G)
    ctr.end()
    ctr.begin("plain:giant")
    stride = E.neg(E.mul(d, G))
    cur = G1
    for j in range((p - 1) // d + 2):
        if cur in table:
            ctr.end()
            return {"alpha": (table[cur] + j * d) % p, "ops": ctr.ops - start}
        cur = E.add(cur, stride)
    ctr.end()
    return {"alpha": None, "ops": ctr.ops - start}


# ---------------------------------------------------------------------------
# p − 1 case
# ---------------------------------------------------------------------------


def _bsgs_in_exponent(E: Curve, base: Scalar, target: Scalar, zeta_hat: int, N: int,
                      p: int, phase: str) -> Optional[int]:
    """Find k in [0, N) with target.P = [ζ̂^k] base.P, where ζ̂ has order N
    in F_p^*.  Baby steps [ζ̂^i] base.P (i < m), giant steps
    [ζ̂^{−m j}] target.P.  Every step is one fixed-base multiplication."""
    m = math.isqrt(N) + 1
    ctr = E.ctr
    ctr.begin(phase + ":baby")
    table: Dict[Point, int] = {}
    z = 1
    for i in range(m):
        table.setdefault(base.mul(z), i)
        z = z * zeta_hat % p
    ctr.end()
    ctr.begin(phase + ":giant")
    step = pow(zeta_hat, -m, p)
    z = 1
    for j in range(m + 1):
        R = target.mul(z)
        if R in table:
            ctr.end()
            return (table[R] + m * j) % N
        z = z * step % p
    ctr.end()
    return None


def cheon_p_minus_1_bsgs(E: Curve, G: Point, G1: Point, Gd: Point, p: int, d: int,
                         comb: bool = True) -> Dict:
    """Cheon, Theorem 1.  Inputs G, G_1 = [α]G, G_d = [α^d]G, d | p − 1."""
    assert (p - 1) % d == 0
    ctr = E.ctr
    start = ctr.ops
    fac = factorise(p - 1)
    zeta = primitive_root(p, fac)
    N = (p - 1) // d
    zeta_hat = pow(zeta, d, p)          # order N
    nbits = p.bit_length()

    # Step 1: α^d = ζ̂^{k0}.
    ctr.begin("p-1:precompute")
    m1 = math.isqrt(N) + 1
    base = Scalar(E, G, nbits, m1, comb)
    targ = Scalar(E, Gd, nbits, m1, comb)
    ctr.end()
    k0 = _bsgs_in_exponent(E, base, targ, zeta_hat, N, p, "p-1:step1")
    if k0 is None:
        return {"alpha": None, "ops": ctr.ops - start}

    # Step 2: α = ζ^{k0 + k1 N}, k1 < d.  Let ζ̌ = ζ^N (order d),
    # α̌ = α ζ^{−k0} = ζ̌^{k1}, so [α̌]G = [ζ^{−k0}] G_1.
    ctr.begin("p-1:step2:precompute")
    m2 = math.isqrt(d) + 1
    G1_shift = E.mul(pow(zeta, -k0, p), G1)
    zeta_check = pow(zeta, N, p)
    base2 = Scalar(E, G, nbits, m2, comb)
    targ2 = Scalar(E, G1_shift, nbits, m2, comb)
    ctr.end()
    k1 = _bsgs_in_exponent(E, base2, targ2, zeta_check, d, p, "p-1:step2")
    if k1 is None:
        return {"alpha": None, "ops": ctr.ops - start}
    alpha = pow(zeta, k0 + k1 * N, p)
    ctr.begin("p-1:verify")
    ok = E.mul(alpha, G) == G1
    ctr.end()
    return {"alpha": alpha if ok else None, "ops": ctr.ops - start, "verified": ok,
            "k0": k0, "k1": k1, "N": N, "table_points": base.table_size + targ.table_size
            + base2.table_size + targ2.table_size, "comb_w": base.comb.w if comb else None}


def cheon_p_minus_1_kangaroo(E: Curve, G: Point, G1: Point, Gd: Point, p: int, d: int,
                             seed: int, theta_bits: Optional[int] = None) -> Dict:
    """Cheon §3.1: two kangaroos in the exponent domain of ζ̂.  The tame
    kangaroo sits at [ζ̂^{w}]G with w known; the wild one at [ζ̂^{w}]G_d =
    [ζ̂^{w + k0}]G.  Jump sizes depend only on the current point, so a
    collision of points is a collision of exponents.  Every jump is one
    fixed-base multiplication (comb tables charged)."""
    assert (p - 1) % d == 0
    rng = random.Random(seed)
    ctr = E.ctr
    start = ctr.ops
    fac = factorise(p - 1)
    zeta = primitive_root(p, fac)
    N = (p - 1) // d
    nbits = p.bit_length()

    def solve(base_pt: Point, target_pt: Point, gen: int, order: int, phase: str) -> Optional[int]:
        """Return k in [0, order) with target_pt = [gen^k] base_pt."""
        if order <= 64:  # tiny: just try them all (counted)
            ctr.begin(phase + ":exhaust")
            z = 1
            for k in range(order):
                if E.mul(z, base_pt) == target_pt:
                    ctr.end()
                    return k
                z = z * gen % p
            ctr.end()
            return None
        ctr.begin(phase + ":precompute")
        expected = 2 * math.isqrt(order)
        base = Scalar(E, base_pt, nbits, expected, True)
        targ = Scalar(E, target_pt, nbits, expected, True)
        ctr.end()
        # jump set: powers of two with mean ≈ √order / 2 (Pollard)
        mean = max(2, math.isqrt(order) // 2)
        r = 16
        jumps = [1 << (i * max(1, mean.bit_length() - 1) // r) for i in range(r)]
        # rescale so the mean is right
        scale = max(1, mean * r // max(1, sum(jumps)))
        jumps = [j * scale for j in jumps]
        tb = theta_bits if theta_bits is not None else max(1, order.bit_length() // 4)
        mask = (1 << tb) - 1
        ctr.begin(phase + ":walk")
        # The exponent space is the cyclic group Z/order, so a kangaroo can
        # wrap around and fall into its own cycle.  A distinguished point
        # already in a kangaroo's *own* table means exactly that: restart it
        # from a fresh random position (van Oorschot–Wiener style).
        tame: Dict[Point, int] = {}
        wild: Dict[Point, int] = {}
        wt = rng.randrange(order)   # tame exponent (absolute)
        ww = rng.randrange(order)   # wild offset (relative to the unknown k0)
        T = base.mul(pow(gen, wt, p))
        W = targ.mul(pow(gen, ww, p))
        restarts = 0
        for _ in range(64 * order):
            # tame
            if T[0] & mask == 0:
                if T in wild:
                    ctr.end()
                    return (wt - wild[T]) % order
                if T in tame:
                    restarts += 1
                    wt = rng.randrange(order)
                    T = base.mul(pow(gen, wt, p))
                    continue
                tame[T] = wt
            wt = (wt + jumps[T[0] % r]) % order
            T = base.mul(pow(gen, wt, p))
            # wild
            if W[0] & mask == 0:
                if W in tame:
                    ctr.end()
                    return (tame[W] - ww) % order
                if W in wild:
                    restarts += 1
                    ww = rng.randrange(order)
                    W = targ.mul(pow(gen, ww, p))
                    continue
                wild[W] = ww
            ww = (ww + jumps[W[0] % r]) % order
            W = targ.mul(pow(gen, ww, p))
        ctr.end()
        return None

    zeta_hat = pow(zeta, d, p)
    k0 = solve(G, Gd, zeta_hat, N, "p-1:kang1")
    if k0 is None:
        return {"alpha": None, "ops": ctr.ops - start}
    ctr.begin("p-1:kang2:shift")
    G1_shift = E.mul(pow(zeta, -k0, p), G1)
    ctr.end()
    k1 = solve(G, G1_shift, pow(zeta, N, p), d, "p-1:kang2")
    if k1 is None:
        return {"alpha": None, "ops": ctr.ops - start}
    alpha = pow(zeta, k0 + k1 * N, p)
    ctr.begin("p-1:verify")
    ok = E.mul(alpha, G) == G1
    ctr.end()
    return {"alpha": alpha if ok else None, "ops": ctr.ops - start, "verified": ok,
            "k0": k0, "k1": k1, "N": N}


# ---------------------------------------------------------------------------
# p + 1 case
# ---------------------------------------------------------------------------


class Fp2:
    """F_{p²} = F_p[θ]/(θ² − a), a a non-residue.  Elements are (s, t) = s + tθ."""

    def __init__(self, p: int, a: int):
        self.p, self.a = p, a

    def mul(self, x, y):
        p, a = self.p, self.a
        return ((x[0] * y[0] + a * x[1] * y[1]) % p, (x[0] * y[1] + x[1] * y[0]) % p)

    def pow(self, x, e: int):
        e %= (self.p * self.p - 1)
        r = (1, 0)
        while e:
            if e & 1:
                r = self.mul(r, x)
            x = self.mul(x, x)
            e >>= 1
        return r

    def inv(self, x):
        p, a = self.p, self.a
        n = (x[0] * x[0] - a * x[1] * x[1]) % p
        ni = pow(n, -1, p)
        return (x[0] * ni % p, (-x[1]) * ni % p)


def _poly_mul(f: List[int], g: List[int], p: int) -> List[int]:
    out = [0] * (len(f) + len(g) - 1)
    for i, fi in enumerate(f):
        if fi:
            for j, gj in enumerate(g):
                out[i + j] = (out[i + j] + fi * gj) % p
    return out


def _poly_pow(f: List[int], e: int, p: int) -> List[int]:
    r = [1]
    while e:
        if e & 1:
            r = _poly_mul(r, f, p)
        f = _poly_mul(f, f, p)
        e >>= 1
    return r


def _fp2_poly_pow(A: List[int], B: List[int], e: int, p: int, a: int):
    """(A(x) + B(x)θ)^e as a pair of polynomials."""
    RA, RB = [1], [0]
    while e:
        if e & 1:
            RA, RB = (_add(_poly_mul(RA, A, p), _scale(_poly_mul(RB, B, p), a, p), p),
                      _add(_poly_mul(RA, B, p), _poly_mul(RB, A, p), p))
        A, B = (_add(_poly_mul(A, A, p), _scale(_poly_mul(B, B, p), a, p), p),
                _scale(_poly_mul(A, B, p), 2, p))
        e >>= 1
    return RA, RB


def _add(f, g, p):
    n = max(len(f), len(g))
    return [((f[i] if i < len(f) else 0) + (g[i] if i < len(g) else 0)) % p for i in range(n)]


def _scale(f, c, p):
    return [c * x % p for x in f]


def _multi_scalar(E: Curve, coeffs: List[int], pts: List[Point]) -> Point:
    """Σ [c_i] P_i by plain double-and-add on each term (counted)."""
    R: Point = None
    for c, P in zip(coeffs, pts):
        if c:
            R = E.add(R, E.mul(c, P))
    return R


def cheon_p_plus_1_bsgs(E: Curve, G_pows: List[Point], p: int, d: int) -> Dict:
    """Cheon, Theorem 2.  Inputs G_i = [α^i]G for 0 ≤ i ≤ 2d, d | p + 1.

    β := ((α − θ)/(α + θ))^d = (A(α) + B(α)θ) / C(α) with
    A + Bθ = (x − θ)^{2d} and C = (x² − a)^d, all of degree ≤ 2d, so
    [A(α)]G, [B(α)]G, [C(α)]G are multi-scalar sums over the auxiliary
    inputs.  β lies in the subgroup of order N = (p+1)/d of the norm-one
    torus H ⊂ F_{p²}^*; BSGS on β = ζ̂^{k0} compares
        (A + Bθ)(s' + t'θ) = C (s_u + t_u θ)
    coordinate-wise in the exponent: 2 fixed-base multiplications per baby
    step, 4 per giant step."""
    assert (p + 1) % d == 0 and len(G_pows) >= 2 * d + 1
    ctr = E.ctr
    start = ctr.ops
    G, G1 = G_pows[0], G_pows[1]
    nbits = p.bit_length()
    # non-residue a
    a = 2
    while sqrt_mod(a, p) is not None:
        a += 1
    F = Fp2(p, a)
    # generator ζ of the torus H (order p+1): take random x, ζ = x^{p−1}
    fac = factorise(p + 1)
    rng = random.Random(p)
    while True:
        x = (rng.randrange(p), rng.randrange(1, p))
        z = F.pow(x, p - 1)
        if all(F.pow(z, (p + 1) // q) != (1, 0) for q in fac):
            zeta = z
            break
    N = (p + 1) // d

    # Polynomials A, B, C in the indeterminate x (= α), by the binomial
    # theorem: (x − θ)^{2d} = Σ C(2d,i) x^i (−θ)^{2d−i} with θ^{2j} = a^j,
    # θ^{2j+1} = a^j θ;  (x² − a)^d = Σ C(d,j) x^{2j} (−a)^{d−j}.
    ctr.begin("p+1:poly")
    A = [0] * (2 * d + 1)
    B = [0] * (2 * d + 1)
    C = [0] * (2 * d + 1)
    for i in range(2 * d + 1):
        e = 2 * d - i                      # power of (−θ)
        c = math.comb(2 * d, i) * pow(-1, e) * pow(a, e // 2, p) % p
        if e % 2 == 0:
            A[i] = c
        else:
            B[i] = c
    for j in range(d + 1):
        C[2 * j] = math.comb(d, j) * pow(-a, d - j, p) % p
    GA = _multi_scalar(E, A, G_pows[: 2 * d + 1])
    GB = _multi_scalar(E, B, G_pows[: 2 * d + 1])
    GC = _multi_scalar(E, C, G_pows[: 2 * d + 1])
    ctr.end()

    def bsgs(zeta_hat, order: int, GA, GB, GC, phase: str) -> Optional[int]:
        """k with (A+Bθ)/C = ζ̂^k, using only [A]G, [B]G, [C]G."""
        m = math.isqrt(order) + 1
        ctr.begin(phase + ":precompute")
        mC = Scalar(E, GC, nbits, 2 * m, True)
        mA = Scalar(E, GA, nbits, 2 * m, True)
        mB = Scalar(E, GB, nbits, 2 * m, True)
        ctr.end()
        ctr.begin(phase + ":baby")
        table: Dict[Tuple[Point, Point], int] = {}
        z = (1, 0)
        for u in range(m):
            table.setdefault((mC.mul(z[0]), mC.mul(z[1])), u)
            z = F.mul(z, zeta_hat)
        ctr.end()
        ctr.begin(phase + ":giant")
        step = F.pow(zeta_hat, -m)
        z = (1, 0)
        for v in range(m + 1):
            s, t = z
            # (A + Bθ)(s + tθ) = (A s + a B t) + (A t + B s)θ
            L0 = E.add(mA.mul(s), mB.mul(a * t % p))
            L1 = E.add(mA.mul(t), mB.mul(s))
            k = table.get((L0, L1))
            if k is not None:
                ctr.end()
                return (k + m * v) % order
            z = F.mul(z, step)
        ctr.end()
        return None

    zeta_hat = F.pow(zeta, d)                 # order N
    k0 = bsgs(zeta_hat, N, GA, GB, GC, "p+1:step1")
    if k0 is None:
        return {"alpha": None, "ops": ctr.ops - start}
    # Step 2.  γ := (α − θ)/(α + θ) = ζ^{k}, k = k0 + k1 N, k1 < d (β = γ^d
    # fixed k mod N).  γ = (α − θ)²/(α² − a) = (A' + B'θ)/C' with
    # A' = x² + a, B' = −2x, C' = x² − a: degree 2, so G_0, G_1, G_2 suffice.
    # Find k1 from γ ζ^{−k0} = ζ̌^{k1}, ζ̌ = ζ^N of order d.
    ctr.begin("p+1:step2:shift")
    GA1 = _multi_scalar(E, [a, 0, 1], G_pows[:3])
    GB1 = _multi_scalar(E, [0, (-2) % p], G_pows[:2])
    GC1 = _multi_scalar(E, [(-a) % p, 0, 1], G_pows[:3])
    s, t = F.pow(zeta, -k0)
    GA2 = E.add(E.mul(s, GA1), E.mul(a * t % p, GB1))
    GB2 = E.add(E.mul(t, GA1), E.mul(s, GB1))
    ctr.end()
    k1 = bsgs(F.pow(zeta, N), d, GA2, GB2, GC1, "p+1:step2")
    if k1 is None:
        return {"alpha": None, "ops": ctr.ops - start}
    k = (k0 + k1 * N) % (p + 1)
    gamma = F.pow(zeta, k)
    # α = θ(1 + γ)/(1 − γ), which must land in F_p.
    num = F.mul((0, 1), ((1 + gamma[0]) % p, gamma[1]))
    den = ((1 - gamma[0]) % p, (-gamma[1]) % p)
    alpha = None
    if den != (0, 0):
        cand = F.mul(num, F.inv(den))
        if cand[1] == 0:
            alpha = cand[0]
    ctr.begin("p+1:verify")
    ok = alpha is not None and E.mul(alpha, G) == G1
    ctr.end()
    return {"alpha": alpha if ok else None, "ops": ctr.ops - start, "verified": ok,
            "k0": k0, "k1": k1, "N": N}

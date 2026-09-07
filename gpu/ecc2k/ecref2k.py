#!/usr/bin/env python3
"""Pure-Python reference for the ECC2K GPU kernels.

Target is ECC2K-95, the Certicom Koblitz challenge: the curve

    E_0 :  y^2 + x y = x^3 + 1     over  F_{2^97},  f(t) = t^97 + t^6 + 1

whose group has order 4r with r a 95-bit prime.  The challenge was solved in
1998 by Harley et al.; it is here as a realistic, retired target for the
kernels, not as an attack on anything live.

Koblitz curves are worth special treatment because the Frobenius map
tau(x, y) = (x^2, y^2) is a cheap endomorphism of the group.  On the prime
subgroup it acts as multiplication by an integer s with s^2 + s + 2 = 0
(mod r), and the 2m points {+-tau^i(P)} form an equivalence class, so a
rho search over classes is sqrt(2m) ~ 13.9 times shorter.

This file is the oracle for everything in gpu/ecc2k/: field and curve
arithmetic, the walk definition, and the generators for the C headers.

    python3 ecref2k.py info    --curve ecc2k95
    python3 ecref2k.py params  --curve ecc2k95 > curve_ecc2k95.h
    python3 ecref2k.py vectors --curve ecc2k95 > vec_ecc2k95.h
"""
import argparse
import random
import sys

WORDS = 4           # field elements are 4 x 32-bit words (m <= 128)
LIMB_BITS = 32
SC_WORDS = 4        # scalar ring mod r is also 4 x 32-bit words


# ---------------------------------------------------------------------------
# F_2[t] helpers.  A polynomial is a Python int; bit i is the coefficient of t^i.
# ---------------------------------------------------------------------------

def pdeg(a):
    return a.bit_length() - 1


def pmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
    return r


def pmod(a, f):
    df = pdeg(f)
    while a.bit_length() - 1 >= df and a:
        a ^= f << (a.bit_length() - 1 - df)
    return a


def pgcd(a, b):
    while b:
        a = pmod(a, b)
        a, b = b, a
    return a


def is_irreducible(f):
    m = pdeg(f)
    x = 2
    r = x
    for _ in range(m):
        r = pmod(pmul(r, r), f)
    if r != x:
        return False
    n = m
    primes = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            primes.append(d)
            n //= d
        d += 1
    if n > 1:
        primes.append(n)
    for p in set(primes):
        r = x
        for _ in range(m // p):
            r = pmod(pmul(r, r), f)
        if pgcd(r ^ x, f) != 1:
            return False
    return True


class F2m:
    """F_2[t]/(f), f an irreducible trinomial t^m + t^k + 1."""

    def __init__(self, m, k):
        self.m, self.k = m, k
        self.f = (1 << m) | (1 << k) | 1
        assert is_irreducible(self.f), "t^%d+t^%d+1 is reducible" % (m, k)

    def mul(self, a, b):
        return pmod(pmul(a, b), self.f)

    def sqr(self, a):
        return self.mul(a, a)

    def inv(self, a):
        assert a != 0
        return self.pow(a, (1 << self.m) - 2)

    def pow(self, a, e):
        r = 1
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.sqr(a)
            e >>= 1
        return r

    def frob(self, a, n=1):
        for _ in range(n % self.m):
            a = self.sqr(a)
        return a

    def trace(self, a):
        """Tr(a) = a + a^2 + ... + a^(2^(m-1)), an element of F_2."""
        t, cur = 0, a
        for _ in range(self.m):
            t ^= cur
            cur = self.sqr(cur)
        assert t in (0, 1)
        return t

    def solve_quad(self, c):
        """Solve z^2 + z = c.  Requires Tr(c) = 0; returns one root."""
        if c == 0:
            return 0
        if self.trace(c) != 0:
            return None
        if self.m % 2 == 1:
            # half-trace: z = sum c^(4^i), i = 0..(m-1)/2
            z, cur = 0, c
            for _ in range((self.m - 1) // 2):
                cur = self.sqr(self.sqr(cur))
                z ^= cur
            z ^= c
            # the half-trace formula gives z with z^2+z = c
            z = 0
            cur = c
            for i in range((self.m - 1) // 2 + 1):
                z ^= cur
                cur = self.sqr(self.sqr(cur))
            return z
        raise NotImplementedError("even m")

    def random(self, rng):
        return rng.getrandbits(self.m)


# ---------------------------------------------------------------------------
# Koblitz curve  y^2 + x y = x^3 + a x^2 + 1  over F_2^m
# ---------------------------------------------------------------------------

class Koblitz:
    def __init__(self, name, m, k, a):
        self.name, self.F, self.a = name, F2m(m, k), a
        self.m = m
        self.b = 1
        self.n = None      # full group order
        self.r = None      # prime subgroup order
        self.h = None      # cofactor
        self.G = None
        self.s = None      # tau acts as multiplication by s on the subgroup

    # -- group law ---------------------------------------------------------
    def neg(self, P):
        return None if P is None else (P[0], P[0] ^ P[1])

    def add(self, P, Q):
        F = self.F
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if y1 ^ y2 == x1:      # Q == -P
                return None
            return self.dbl(P)
        lam = F.mul(y1 ^ y2, F.inv(x1 ^ x2))
        x3 = F.sqr(lam) ^ lam ^ x1 ^ x2 ^ self.a
        y3 = F.mul(lam, x1 ^ x3) ^ x3 ^ y1
        return (x3, y3)

    def dbl(self, P):
        F = self.F
        if P is None or P[0] == 0:
            return None
        x1, y1 = P
        lam = x1 ^ F.mul(y1, F.inv(x1))
        x3 = F.sqr(lam) ^ lam ^ self.a
        y3 = F.sqr(x1) ^ F.mul(lam ^ 1, x3)
        return (x3, y3)

    def mul(self, kk, P):
        if kk < 0:
            return self.mul(-kk, self.neg(P))
        R = None
        while kk:
            if kk & 1:
                R = self.add(R, P)
            P = self.dbl(P)
            kk >>= 1
        return R

    def frob(self, P, n=1):
        """tau^n(P) = (x^(2^n), y^(2^n)) -- a group endomorphism."""
        if P is None:
            return None
        return (self.F.frob(P[0], n), self.F.frob(P[1], n))

    def on_curve(self, P):
        if P is None:
            return True
        x, y = P
        F = self.F
        lhs = F.sqr(y) ^ F.mul(x, y)
        rhs = F.mul(F.sqr(x), x) ^ F.mul(self.a, F.sqr(x)) ^ self.b
        return lhs == rhs

    def lift_x(self, x):
        """A point with this x, or None if x is not on the curve."""
        F = self.F
        if x == 0:
            return None
        # y^2 + xy = x^3 + a x^2 + b  =>  z^2 + z = (x^3+ax^2+b)/x^2, y = xz
        c = F.mul(F.mul(F.sqr(x), x) ^ F.mul(self.a, F.sqr(x)) ^ self.b,
                  F.inv(F.sqr(x)))
        z = F.solve_quad(c)
        if z is None:
            return None
        return (x, F.mul(x, z))

    def random_point(self, rng):
        while True:
            P = self.lift_x(self.F.random(rng))
            if P is not None and self.on_curve(P):
                return P

    # -- parameters --------------------------------------------------------
    def compute_order(self):
        mu = -1 if self.a == 0 else 1
        t0, t1 = 2, mu
        for _ in range(self.m - 1):
            t0, t1 = t1, mu * t1 - 2 * t0
        self.n = (1 << self.m) + 1 - t1
        return self.n


def is_probable_prime(n, rounds=32):
    if n < 2:
        return False
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def sqrt_mod(a, p):
    """Tonelli-Shanks."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p
            i += 1
        bb = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, bb * bb % p, t * bb * bb % p, r * bb % p
    return r


def setup_curve(E, seed=1):
    """Fill in n, r, h, G and s.  s is verified against tau numerically."""
    rng = random.Random(seed)
    E.compute_order()
    for h in (2, 4, 6, 8, 12, 16):
        if E.n % h == 0 and is_probable_prime(E.n // h):
            E.h, E.r = h, E.n // h
            break
    else:
        raise SystemExit("no prime subgroup for %s (n = %d)" % (E.name, E.n))

    # generator of the prime-order subgroup
    while True:
        P = E.random_point(rng)
        G = E.mul(E.h, P)
        if G is not None and E.mul(E.r, G) is None:
            E.G = G
            break

    # tau acts as multiplication by s: s^2 - mu*s + 2 = 0 (mod r)
    mu = -1 if E.a == 0 else 1
    disc = sqrt_mod((mu * mu - 8) % E.r, E.r)
    assert disc is not None, "tau eigenvalue does not exist mod r"
    inv2 = pow(2, -1, E.r)
    for cand in ((mu + disc) * inv2 % E.r, (mu - disc) * inv2 % E.r):
        if E.frob(E.G) == E.mul(cand, E.G):
            E.s = cand
            break
    else:
        raise SystemExit("could not identify the Frobenius eigenvalue")
    return E


# ---------------------------------------------------------------------------
# The Frobenius-invariant class function
# ---------------------------------------------------------------------------

def normal_masks(E, gamma=None):
    """Masks M_j with  bit k of M_j  =  Tr(t^k * gamma^(2^j)).

    Then parity(popcount(x & M_j)) = Tr(x * gamma^(2^j)), and

        g(x) = #{ j : Tr(x * gamma^(2^j)) = 1 }

    is invariant under x -> x^2, because squaring merely rotates j.  With
    gamma a normal element this is exactly the Hamming weight of x in the
    normal basis generated by gamma; correctness of the invariance does not
    depend on gamma being normal, only the quality of the distribution does.
    """
    F = E.F
    m = F.m
    if gamma is None:
        gamma = pick_normal_element(F)
    masks = []
    g = gamma
    for _ in range(m):
        mask = 0
        for kk in range(m):
            if F.trace(F.mul(1 << kk, g)):
                mask |= 1 << kk
        masks.append(mask)
        g = F.sqr(g)
    return gamma, masks


def pick_normal_element(F):
    """gamma whose conjugates form a basis: the m x m matrix of coordinates
    of gamma^(2^j) must be invertible over F_2."""
    for cand in range(2, 1 << min(F.m, 20)):
        rows = []
        g = cand
        for _ in range(F.m):
            rows.append(g)
            g = F.sqr(g)
        if gf2_rank(rows, F.m) == F.m:
            return cand
    raise SystemExit("no normal element found")


def gf2_rank(rows, n):
    rows = list(rows)
    rank = 0
    for col in range(n):
        piv = None
        for i in range(rank, len(rows)):
            if rows[i] >> col & 1:
                piv = i
                break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(len(rows)):
            if i != rank and (rows[i] >> col & 1):
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def hw_class(x, masks):
    """The Frobenius-invariant weight g(x) defined above."""
    return sum(1 for M in masks if bin(x & M).count("1") & 1)


# ---------------------------------------------------------------------------
# The walk
# ---------------------------------------------------------------------------

class FrobWalk:
    """Frobenius-class r-adding walk, the ECC2K-130 iteration shape:

        j(P) = (g(x_P) mod NJ) + JMIN
        P   -> P + tau^j(P)

    j depends only on the Frobenius-invariant g, so the map descends to the
    classes {+-tau^i(P)}: f(tau P) = tau f(P) and f(-P) = -f(P).  Walking on
    classes of size 2m shortens the search by sqrt(2m).

    A step multiplies the coefficients: P = aP0 + bQ0 becomes
    (1 + s^j)(aP0 + bQ0), so a <- a(1+s^j), b <- b(1+s^j) mod r.
    """

    NJ = 8
    JMIN = 3

    def __init__(self, E, masks):
        self.E, self.masks = E, masks

    def j_of(self, P):
        return (hw_class(P[0], self.masks) % self.NJ) + self.JMIN

    def step(self, P):
        j = self.j_of(P)
        return self.E.add(P, self.E.frob(P, j)), j

    def is_dp(self, P, dp_threshold):
        return hw_class(P[0], self.masks) <= dp_threshold

    def canonical(self, P):
        """Representative of the class: the Frobenius image with the
        smallest x.  Returns (x_min, rotation)."""
        best, beste = P[0], 0
        cur = P[0]
        for e in range(1, self.E.m):
            cur = self.E.F.sqr(cur)
            if cur < best:
                best, beste = cur, e
        return best, beste


# ---------------------------------------------------------------------------
# Curve catalogue
# ---------------------------------------------------------------------------

# (m, k, a) with t^m + t^k + 1 irreducible
CURVES = {
    "ecc2k95": (97, 6, 0),     # the Certicom challenge curve
    "k23": (23, 5, 0),         # toy: 21-bit subgroup, solvable end to end
    "k41": (41, 3, 0),         # toy: 40-bit subgroup
}

_CACHE = {}


def get_curve(name):
    if name not in _CACHE:
        if name not in CURVES:
            raise SystemExit("unknown curve %s (have: %s)" %
                             (name, ", ".join(sorted(CURVES))))
        m, k, a = CURVES[name]
        _CACHE[name] = setup_curve(Koblitz(name, m, k, a))
    return _CACHE[name]


# ---------------------------------------------------------------------------
# C header emission
# ---------------------------------------------------------------------------

def limbs(x, n):
    return [(x >> (LIMB_BITS * i)) & 0xFFFFFFFF for i in range(n)]


def c_limbs(x, n=WORDS):
    return "{" + ", ".join("0x%08xu" % v for v in limbs(x, n)) + "}"


def c_point(E, P):
    if P is None:
        return "{%s, %s, 1}" % (c_limbs(0), c_limbs(0))
    return "{%s, %s, 0}" % (c_limbs(P[0]), c_limbs(P[1]))


def emit_params(E, out):
    w = out.write
    F = E.F
    gamma, masks = normal_masks(E)
    R256 = 1 << (LIMB_BITS * SC_WORDS)
    r = E.r
    w("/* Generated by ecref2k.py -- %s parameters.  Do not edit. */\n" % E.name)
    w("#ifndef GPU_ECC2K_CURVE_PARAMS_H\n#define GPU_ECC2K_CURVE_PARAMS_H\n")
    w('#define CURVE2K_NAME "%s"\n' % E.name)
    w("#define F2M_M %d\n" % F.m)
    w("#define F2M_K %d\n" % F.k)
    w("#define F2M_WORDS %d\n" % WORDS)
    w("#define CURVE2K_A %d\n" % E.a)
    w("#define CURVE2K_COFACTOR %d\n" % E.h)
    w("/* full group order and prime subgroup order */\n")
    w("#define CURVE2K_N_LIMBS %s\n" % c_limbs(E.n))
    w("#define CURVE2K_R_LIMBS %s\n" % c_limbs(E.r))
    w("#define CURVE2K_R_BITS %d\n" % E.r.bit_length())
    w("/* tau acts on the subgroup as multiplication by s */\n")
    w("#define CURVE2K_S_LIMBS %s\n" % c_limbs(E.s))
    w("#define CURVE2K_GX_LIMBS %s\n" % c_limbs(E.G[0]))
    w("#define CURVE2K_GY_LIMBS %s\n" % c_limbs(E.G[1]))
    w("/* Montgomery constants for the scalar ring mod r */\n")
    w("#define SC_WORDS %d\n" % SC_WORDS)
    w("#define SC_NPRIME 0x%08xu\n" % ((-pow(r, -1, 1 << 32)) % (1 << 32)))
    w("#define SC_R1_LIMBS %s\n" % c_limbs(R256 % r))
    w("#define SC_R2_LIMBS %s\n" % c_limbs(R256 * R256 % r))
    w("#define SC_RM2_LIMBS %s\n" % c_limbs(r - 2))
    w("/* Frobenius-invariant class weight: g(x) = #{j : parity(x & M_j)},\n"
      "   with gamma = 0x%x a normal element, so g is the normal-basis\n"
      "   Hamming weight.  Invariance needs only the conjugate structure. */\n" % gamma)
    # windowed change-of-basis table: 4-bit windows over the m bits of x
    nwin = (F.m + 3) // 4
    w("#define F2M_CB_WINDOWS %d\n" % nwin)
    w("/* cb_table[w][v] = the class-weight bit-vector contributed by nibble v\n"
      "   at window w; g(x) = popcount(XOR of the selected entries). */\n")
    # Entries are padded to WORDS+1 so that, staged in shared memory, the
    # stride is odd and the 16 possible nibbles of a warp land in 16
    # different banks.  A power-of-two stride would cost a bank conflict on
    # the hottest read in the walk.
    stride = WORDS + 1
    w("#define F2M_CB_STRIDE %d\n" % stride)
    parts = []
    for wnd in range(nwin):
        rows = []
        for v in range(16):
            acc = 0
            for bit in range(4):
                kk = wnd * 4 + bit
                if kk < F.m and (v >> bit) & 1:
                    for j, M in enumerate(masks):
                        if (M >> kk) & 1:
                            acc ^= 1 << j
            rows.append(c_limbs(acc, stride))
        parts.append("{%s}" % ", ".join(rows))
    w("#define F2M_CB_TABLE {%s}\n" % ", ".join(parts))
    w("#endif\n")


def emit_vectors(E, out, count=64, seed=4242):
    rng = random.Random(seed)
    F = E.F
    gamma, masks = normal_masks(E)
    w = out.write
    w("/* Generated by ecref2k.py -- %s test vectors. */\n" % E.name)
    w("#ifndef GPU_ECC2K_TEST_VECTORS_H\n#define GPU_ECC2K_TEST_VECTORS_H\n")
    w('#define VEC2K_CURVE_NAME "%s"\n' % E.name)

    w("#define VEC2K_FIELD_COUNT %d\n" % count)
    w("/* a, b, a+b, a*b, a^2, a^-1, tau^3(a), and g(a) in the last slot */\n")
    w("static const uint32_t vec2k_field[VEC2K_FIELD_COUNT][7][%d] = {\n" % WORDS)
    gvals = []
    for i in range(count):
        if i == 0:
            a, b = 0, 0
        elif i == 1:
            a, b = 1, 1
        elif i == 2:
            a, b = (1 << F.m) - 1, (1 << F.m) - 1
        elif i == 3:
            a, b = 1 << (F.m - 1), 1 << (F.m - 1)
        else:
            a, b = F.random(rng), F.random(rng)
        inv = F.inv(a) if a else 0
        vals = [a, b, a ^ b, F.mul(a, b), F.sqr(a), inv, F.frob(a, 3)]
        gvals.append(hw_class(a, masks))
        w("  {" + ", ".join(c_limbs(v) for v in vals) + "},\n")
    w("};\n")
    w("static const uint32_t vec2k_field_g[VEC2K_FIELD_COUNT] = {%s};\n" %
      ", ".join(str(v) for v in gvals))

    w("#define VEC2K_POINT_COUNT %d\n" % count)
    w("typedef struct { uint32_t x[%d], y[%d]; uint32_t inf; } vec2k_pt_t;\n"
      % (WORDS, WORDS))
    w("static const struct { vec2k_pt_t P, Q, sum, dbl, kP, tauP; uint32_t k[%d]; }\n"
      "  vec2k_point[VEC2K_POINT_COUNT] = {\n" % SC_WORDS)
    for i in range(count):
        P = E.random_point(rng)
        if i == 0:
            Q = E.neg(P)
        elif i == 1:
            Q = None
        elif i == 2:
            Q = P
        else:
            Q = E.random_point(rng)
        if i == 3:
            kk = 0
        elif i == 4:
            kk = 1
        elif i == 5:
            kk = E.n
        elif i == 6:
            kk = E.n - 1
        else:
            kk = rng.getrandbits(LIMB_BITS * SC_WORDS)
        w("  {%s, %s, %s, %s, %s, %s, %s},\n" % (
            c_point(E, P), c_point(E, Q), c_point(E, E.add(P, Q)),
            c_point(E, E.dbl(P)), c_point(E, E.mul(kk, P)),
            c_point(E, E.frob(P)), c_limbs(kk, SC_WORDS)))
    w("};\n")

    # tau(G) == s*G, the identity the whole attack rests on
    w("/* tau(G) and s*G must agree -- the Frobenius eigenvalue check */\n")
    w("static const vec2k_pt_t vec2k_tau_G = %s;\n" % c_point(E, E.frob(E.G)))
    w("static const vec2k_pt_t vec2k_sG = %s;\n" % c_point(E, E.mul(E.s, E.G)))

    # walk traces
    walk = FrobWalk(E, masks)
    nsteps = 128
    starts = [E.mul(E.h, E.random_point(rng)) for _ in range(4)]
    starts = [P for P in starts if P is not None][:4]
    w("#define VEC2K_WALK_STEPS %d\n" % nsteps)
    w("#define VEC2K_WALK_COUNT %d\n" % len(starts))
    w("#define VEC2K_WALK_NJ %d\n" % FrobWalk.NJ)
    w("#define VEC2K_WALK_JMIN %d\n" % FrobWalk.JMIN)
    w("static const vec2k_pt_t vec2k_walk_start[%d] = {\n" % len(starts))
    for P in starts:
        w("  %s,\n" % c_point(E, P))
    w("};\n")
    ends, traces = [], []
    for idx, P in enumerate(starts):
        cur = P
        tr = []
        for _ in range(nsteps):
            cur, j = walk.step(cur)
            if cur is None:
                break
            tr.append(j)
        ends.append(cur)
        if idx == 0:
            traces = tr
    w("static const vec2k_pt_t vec2k_walk_end[%d] = {\n" % len(starts))
    for P in ends:
        w("  %s,\n" % c_point(E, P))
    w("};\n")
    w("static const uint8_t vec2k_walk_trace[%d] = {%s};\n" %
      (len(traces), ", ".join(str(t) for t in traces)))

    # canonical class representatives
    w("#define VEC2K_CANON_COUNT 16\n")
    w("static const uint32_t vec2k_canon_in[16][%d] = {\n" % WORDS)
    canon = []
    for _ in range(16):
        P = E.random_point(rng)
        canon.append(walk.canonical(P))
        w("  %s,\n" % c_limbs(P[0]))
    w("};\n")
    w("static const uint32_t vec2k_canon_out[16][%d] = {\n" % WORDS)
    for xmin, _ in canon:
        w("  %s,\n" % c_limbs(xmin))
    w("};\n")
    w("#endif\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["info", "params", "vectors"])
    ap.add_argument("--curve", default="ecc2k95")
    ap.add_argument("--count", type=int, default=64)
    args = ap.parse_args()
    E = get_curve(args.curve)
    if args.what == "info":
        print("curve  %s : y^2 + xy = x^3 + %d x^2 + 1  over F_2^%d" %
              (E.name, E.a, E.m))
        print("field  f(t) = t^%d + t^%d + 1" % (E.F.m, E.F.k))
        print("order  #E = %d  = %d * %d" % (E.n, E.h, E.r))
        print("       r is %d bits, prime" % E.r.bit_length())
        print("tau    s = %d" % E.s)
        print("       s^2 %s s + 2 = 0 mod r : %s" %
              ("+" if E.a == 0 else "-",
               (E.s * E.s + (1 if E.a == 0 else -1) * E.s + 2) % E.r == 0))
        print("       tau(G) == s*G : %s" % (E.frob(E.G) == E.mul(E.s, E.G)))
        print("G      x = 0x%x" % E.G[0])
        print("       y = 0x%x" % E.G[1])
        gamma, masks = normal_masks(E)
        print("gamma  0x%x (normal element)" % gamma)
        rng = random.Random(7)
        xs = [E.random_point(rng)[0] for _ in range(200)]
        gs = [hw_class(x, masks) for x in xs]
        print("g(x)   mean %.1f over %d samples (m/2 = %.1f)" %
              (sum(gs) / len(gs), len(gs), E.m / 2))
        inv = all(hw_class(x, masks) == hw_class(E.F.sqr(x), masks) for x in xs)
        print("       Frobenius invariant on all samples: %s" % inv)
        print("class  size 2m = %d, rho speedup sqrt(2m) = %.1fx"
              % (2 * E.m, (2 * E.m) ** 0.5))
    elif args.what == "params":
        emit_params(E, sys.stdout)
    else:
        emit_vectors(E, sys.stdout, count=args.count)


if __name__ == "__main__":
    main()

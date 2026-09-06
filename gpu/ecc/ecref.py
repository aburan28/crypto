#!/usr/bin/env python3
"""Pure-Python reference for the GPU elliptic-curve kernels.

This file is the oracle: everything the CUDA headers compute is checked
against it (on the host, via the CPU-emulation build, and on the device via
`bench`).  It deliberately uses only Python's built-in big integers so it has
no dependencies and is easy to audit.

Two things live here:

  1. A minimal short-Weierstrass curve model (affine arithmetic, scalar
     multiplication, Montgomery-form helpers) plus the secp256k1 parameters
     and a deterministic search for a small prime-order "toy" curve that the
     end-to-end Pollard-rho test can actually solve.

  2. Generators for the C headers consumed by the kernels and tests:

        python3 ecref.py params --curve secp256k1 > curve_secp256k1.h
        python3 ecref.py params --curve toy40     > curve_toy40.h
        python3 ecref.py vectors --curve secp256k1 --mode mont > vec_secp256k1_mont.h
        python3 ecref.py vectors --curve secp256k1 --mode fast > vec_secp256k1_fast.h
        python3 ecref.py vectors --curve toy40     --mode mont > vec_toy40_mont.h

The r-adding walk used by the rho kernel is defined here as well (see
`RhoWalk`); the C code must match it bit-for-bit, including the choice of
hashing the *internal* field representation rather than the canonical one
(see gpu/README.md, "Walk definition").
"""
import argparse
import hashlib
import random
import sys

LIMBS = 8          # 256-bit numbers as 8 x 32-bit little-endian limbs
LIMB_BITS = 32
R256 = 1 << 256    # Montgomery radix


# ----------------------------------------------------------------------------
# Number theory helpers
# ----------------------------------------------------------------------------

def is_probable_prime(n, rounds=32):
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for q in small:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    rng = random.Random(0xC0FFEE ^ n)
    for _ in range(rounds):
        a = rng.randrange(2, n - 1)
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


def isqrt(n):
    x = int(n ** 0.5)
    while x * x > n:
        x -= 1
    while (x + 1) * (x + 1) <= n:
        x += 1
    return x


# ----------------------------------------------------------------------------
# Curve model
# ----------------------------------------------------------------------------

class Curve:
    """y^2 = x^3 + a x + b over F_p.  Points are (x, y) tuples or None (O)."""

    def __init__(self, name, p, a, b, n=None, G=None):
        self.name, self.p, self.a, self.b, self.n, self.G = name, p, a, b, n, G

    # -- affine group law ---------------------------------------------------
    def neg(self, P):
        return None if P is None else (P[0], (-P[1]) % self.p)

    def add(self, P, Q):
        p = self.p
        if P is None:
            return Q
        if Q is None:
            return P
        if P[0] == Q[0]:
            if (P[1] + Q[1]) % p == 0:
                return None
            lam = (3 * P[0] * P[0] + self.a) * pow(2 * P[1], -1, p) % p
        else:
            lam = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p) % p
        x3 = (lam * lam - P[0] - Q[0]) % p
        y3 = (lam * (P[0] - x3) - P[1]) % p
        return (x3, y3)

    def dbl(self, P):
        return self.add(P, P)

    def mul(self, k, P):
        R = None
        while k > 0:
            if k & 1:
                R = self.add(R, P)
            P = self.dbl(P)
            k >>= 1
        return R

    def on_curve(self, P):
        if P is None:
            return True
        x, y = P
        return (y * y - (x * x * x + self.a * x + self.b)) % self.p == 0

    def random_point(self, rng):
        p = self.p
        while True:
            x = rng.randrange(p)
            rhs = (x * x * x + self.a * x + self.b) % p
            if p % 4 == 3:
                y = pow(rhs, (p + 1) // 4, p)
            else:
                y = tonelli_shanks(rhs, p)
                if y is None:
                    continue
            if y * y % p == rhs:
                return (x, y)

    # -- Montgomery helpers --------------------------------------------------
    @property
    def r_mod_p(self):
        return R256 % self.p

    def to_mont(self, x):
        return x * R256 % self.p

    def from_mont(self, x):
        return x * pow(R256, -1, self.p) % self.p


def tonelli_shanks(n, p):
    n %= p
    if n == 0:
        return 0
    if pow(n, (p - 1) // 2, p) != 1:
        return None
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


SECP256K1 = Curve(
    "secp256k1",
    p=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F,
    a=0,
    b=7,
    n=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141,
    G=(0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
       0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8),
)


def point_order_bsgs(E, P):
    """Order of P for a curve with p large enough that the Hasse interval
    (width 4 sqrt p) is shorter than p.  Returns the unique m in the Hasse
    interval with m*P = O, or None if there are several candidates."""
    p = E.p
    w = 2 * isqrt(p) + 1
    lo = p + 1 - w
    Q = E.mul(lo, P)                  # want t >= 0 with (lo + t) P = O
    msz = isqrt(2 * w) + 1
    baby = {}
    B = None                          # B = j*P
    for j in range(msz):
        key = B
        baby.setdefault(key, j)
        B = E.add(B, P)
    # giant: Q + i*msz*P == -(j P)  <=>  (lo + i*msz + j) P = O
    step = E.mul(msz, P)
    cands = []
    cur = Q
    for i in range(msz + 2):
        negcur = E.neg(cur)
        if negcur in baby:
            cands.append(lo + i * msz + baby[negcur])
        cur = E.add(cur, step)
    cands = sorted(set(c for c in cands if lo <= c <= p + 1 + w and E.mul(c, P) is None))
    return cands[0] if len(cands) == 1 else None


def make_toy_curve(bits=40, seed=1):
    """Deterministically find y^2 = x^3 + a x + b over a `bits`-bit prime with
    prime group order.  a != 0 so the generic doubling path is exercised."""
    rng = random.Random(seed)
    while True:
        p = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if not is_probable_prime(p):
            continue
        for _ in range(40):
            a, b = rng.randrange(1, p), rng.randrange(1, p)
            if (4 * a ** 3 + 27 * b ** 2) % p == 0:
                continue
            E = Curve("toy%d" % bits, p, a, b)
            P = E.random_point(rng)
            m = point_order_bsgs(E, P)
            if m is None or not is_probable_prime(m):
                continue
            E.n, E.G = m, P
            assert E.mul(m, P) is None
            return E


_CURVES = {}


def get_curve(name):
    if name not in _CURVES:
        if name == "secp256k1":
            _CURVES[name] = SECP256K1
        elif name.startswith("toy"):
            _CURVES[name] = make_toy_curve(int(name[3:]))
        else:
            raise SystemExit("unknown curve " + name)
    return _CURVES[name]


# ----------------------------------------------------------------------------
# Field representation ("mode") and the r-adding walk
# ----------------------------------------------------------------------------

class Repr:
    """How the kernels store field elements.  'fast' = canonical integers
    (secp256k1 special reduction); 'mont' = Montgomery form x*2^256 mod p."""

    def __init__(self, E, mode):
        self.E, self.mode = E, mode
        if mode == "fast" and E.name != "secp256k1":
            raise SystemExit("fast mode is secp256k1-only")

    def internal(self, x):
        return x if self.mode == "fast" else self.E.to_mont(x)


class RhoWalk:
    """r-adding walk with optional negation map, defined on the *internal*
    representation of the affine coordinates.

    partition(P) = limb0(internal(x)) & (R-1)
    DP test      = ((limb0(internal(x)) >> 8) & dp_mask) == 0
    negation map = replace (x, y) by (x, -y) if internal(y) > (p-1)/2
    step         = P <- P + M[partition(P)]   (M[j] = c_j P + d_j Q)
    """

    def __init__(self, E, repr_, r_bits, neg_map, table_seed=0):
        self.E, self.repr, self.R = E, repr_, 1 << r_bits
        self.neg_map = neg_map
        self.table_seed = table_seed

    def build_table(self, P, Q):
        rng = random.Random(self.table_seed)
        self.coef = [(rng.randrange(1, self.E.n), rng.randrange(1, self.E.n))
                     for _ in range(self.R)]
        self.M = [self.E.add(self.E.mul(c, P), self.E.mul(d, Q)) for c, d in self.coef]
        return self.M

    def partition(self, P):
        return self.repr.internal(P[0]) & (self.R - 1)

    def is_dp(self, P, dp_mask):
        return (((self.repr.internal(P[0]) & 0xFFFFFFFF) >> 8) & dp_mask) == 0

    def canonical(self, P):
        if P is None or not self.neg_map:
            return P, False
        if self.repr.internal(P[1]) > (self.E.p - 1) // 2:
            return self.E.neg(P), True
        return P, False

    def step(self, P):
        """One step; returns (new point, partition index, negated flag)."""
        j = self.partition(P)
        Pn = self.E.add(P, self.M[j])
        Pn, negated = self.canonical(Pn)
        return Pn, j, negated

    def walk(self, P, steps):
        for _ in range(steps):
            P, _, _ = self.step(P)
            if P is None:
                break
        return P


# ----------------------------------------------------------------------------
# C header emission
# ----------------------------------------------------------------------------

def limbs(x, n=LIMBS):
    return [(x >> (LIMB_BITS * i)) & 0xFFFFFFFF for i in range(n)]


def c_limbs(x, n=LIMBS):
    return "{" + ", ".join("0x%08xu" % v for v in limbs(x, n)) + "}"


def c_point(P):
    if P is None:
        return "{%s, %s, 1}" % (c_limbs(0), c_limbs(0))
    return "{%s, %s, 0}" % (c_limbs(P[0]), c_limbs(P[1]))


def emit_params(E, out):
    p = E.p
    nprime32 = (-pow(p, -1, 1 << 32)) % (1 << 32)
    r2 = R256 * R256 % p
    r1 = R256 % p
    w = out.write
    w("/* Generated by ecref.py -- curve parameters for %s.  Do not edit. */\n" % E.name)
    w("#ifndef GPU_ECC_CURVE_PARAMS_H\n#define GPU_ECC_CURVE_PARAMS_H\n")
    w('#define CURVE_NAME "%s"\n' % E.name)
    w("#define CURVE_BITS %d\n" % p.bit_length())
    w("#define CURVE_A_IS_ZERO %d\n" % (1 if E.a == 0 else 0))
    w("#define CURVE_IS_SECP256K1 %d\n" % (1 if E.name == "secp256k1" else 0))
    w("/* prime p, little-endian 32-bit limbs */\n")
    w("#define FP_P_LIMBS %s\n" % c_limbs(p))
    w("/* -p^{-1} mod 2^32 (Montgomery CIOS constant) */\n")
    w("#define FP_NPRIME 0x%08xu\n" % nprime32)
    w("/* R = 2^256 mod p and R^2 mod p */\n")
    w("#define FP_R1_LIMBS %s\n" % c_limbs(r1))
    w("#define FP_R2_LIMBS %s\n" % c_limbs(r2))
    w("/* (p-1)/2, for the negation-map canonical choice of y */\n")
    w("#define FP_HALF_LIMBS %s\n" % c_limbs((p - 1) // 2))
    w("/* p - 2, Fermat inversion exponent */\n")
    w("#define FP_PM2_LIMBS %s\n" % c_limbs(p - 2))
    w("/* curve coefficients (canonical integers) and group order */\n")
    w("#define CURVE_A_LIMBS %s\n" % c_limbs(E.a))
    w("#define CURVE_B_LIMBS %s\n" % c_limbs(E.b))
    w("#define CURVE_N_LIMBS %s\n" % c_limbs(E.n))
    w("#define CURVE_GX_LIMBS %s\n" % c_limbs(E.G[0]))
    w("#define CURVE_GY_LIMBS %s\n" % c_limbs(E.G[1]))
    n = E.n
    w("/* Montgomery constants for the scalar ring mod n (host-side DLP solve) */\n")
    w("#define N_BITS %d\n" % n.bit_length())
    w("#define N_NPRIME 0x%08xu\n" % ((-pow(n, -1, 1 << 32)) % (1 << 32)))
    w("#define N_R1_LIMBS %s\n" % c_limbs(R256 % n))
    w("#define N_R2_LIMBS %s\n" % c_limbs(R256 * R256 % n))
    w("#define N_HALF_LIMBS %s\n" % c_limbs((n - 1) // 2))
    w("#define N_PM2_LIMBS %s\n" % c_limbs(n - 2))
    w("#endif\n")


def emit_vectors(E, mode, out, count=64, seed=12345):
    rng = random.Random(seed)
    rp = Repr(E, mode)
    p = E.p
    w = out.write
    w("/* Generated by ecref.py -- test vectors for %s (%s mode). */\n" % (E.name, mode))
    w("#ifndef GPU_ECC_TEST_VECTORS_H\n#define GPU_ECC_TEST_VECTORS_H\n")
    w('#define VEC_CURVE_NAME "%s"\n#define VEC_MODE "%s"\n' % (E.name, mode))

    # ---- field vectors ----
    w("#define VEC_FIELD_COUNT %d\n" % count)
    w("static const uint32_t vec_field[VEC_FIELD_COUNT][7][8] = {\n")
    for i in range(count):
        if i == 0:
            a, b = 0, 0
        elif i == 1:
            a, b = p - 1, p - 1
        elif i == 2:
            a, b = 1, p - 1
        elif i == 3:
            a, b = p - 1, 1
        else:
            a, b = rng.randrange(p), rng.randrange(p)
        inv = pow(a, -1, p) if a else 0
        vals = [a, b, (a + b) % p, (a - b) % p, a * b % p, a * a % p, inv]
        w("  {" + ", ".join(c_limbs(v) for v in vals) + "},\n")
    w("};\n")

    # ---- point vectors: P, Q, P+Q, 2P, kP, k ----
    w("#define VEC_POINT_COUNT %d\n" % count)
    w("typedef struct { uint32_t x[8], y[8]; uint32_t inf; } vec_pt_t;\n")
    w("static const struct { vec_pt_t P, Q, sum, dbl, kP; uint32_t k[8]; } vec_point[VEC_POINT_COUNT] = {\n")
    for i in range(count):
        P = E.random_point(rng)
        if i == 0:
            Q = E.neg(P)              # P + (-P) = O
        elif i == 1:
            Q = None                  # P + O
        elif i == 2:
            Q = P                     # add(P,P) must equal dbl
        else:
            Q = E.random_point(rng)
        if i == 3:
            k = 0
        elif i == 4:
            k = 1
        elif i == 5:
            k = E.n                   # n P = O
        elif i == 6:
            k = E.n - 1               # -P
        elif i == 7:
            k = (1 << 256) - 1
        else:
            k = rng.getrandbits(256)
        w("  {%s, %s, %s, %s, %s, %s},\n" % (
            c_point(P), c_point(Q), c_point(E.add(P, Q)), c_point(E.dbl(P)),
            c_point(E.mul(k, P)), c_limbs(k)))
    w("};\n")

    # ---- rho walk vectors ----
    for neg in (0, 1):
        r_bits = 5
        walk = RhoWalk(E, rp, r_bits, neg_map=bool(neg), table_seed=777 + neg)
        P0 = E.G
        Q0 = E.mul(rng.randrange(2, E.n), E.G)
        M = walk.build_table(P0, Q0)
        tag = "neg" if neg else "plain"
        w("#define VEC_WALK_%s_RBITS %d\n" % (tag.upper(), r_bits))
        w("static const vec_pt_t vec_walk_%s_P = %s;\n" % (tag, c_point(P0)))
        w("static const vec_pt_t vec_walk_%s_Q = %s;\n" % (tag, c_point(Q0)))
        w("static const vec_pt_t vec_walk_%s_table[%d] = {\n" % (tag, 1 << r_bits))
        for T in M:
            w("  %s,\n" % c_point(T))
        w("};\n")
        w("static const uint32_t vec_walk_%s_coef[%d][2][8] = {\n" % (tag, 1 << r_bits))
        for c, d in walk.coef:
            w("  {%s, %s},\n" % (c_limbs(c), c_limbs(d)))
        w("};\n")
        nsteps = 200
        starts = [E.random_point(rng) for _ in range(4)]
        w("#define VEC_WALK_%s_STEPS %d\n" % (tag.upper(), nsteps))
        w("#define VEC_WALK_%s_COUNT %d\n" % (tag.upper(), len(starts)))
        w("static const vec_pt_t vec_walk_%s_start[%d] = {\n" % (tag, len(starts)))
        for S in starts:
            w("  %s,\n" % c_point(S))
        w("};\n")
        w("static const vec_pt_t vec_walk_%s_end[%d] = {\n" % (tag, len(starts)))
        for S in starts:
            w("  %s,\n" % c_point(walk.walk(S, nsteps)))
        w("};\n")
        # per-step trace of the first walk (partition index + negated flag)
        trace = []
        S = starts[0]
        for _ in range(nsteps):
            S, j, negated = walk.step(S)
            trace.append(j | (negated << 8))
        w("static const uint16_t vec_walk_%s_trace[%d] = {%s};\n" % (
            tag, nsteps, ", ".join(str(t) for t in trace)))
    w("#endif\n")


def emit_vhdl_vectors(E, out, count=256, seed=999):
    """Plain-text vectors for the GHDL testbenches in hdl/ecc/.

    Format, one record per line, all values 64 hex digits:
        MUL <a> <b> <a*b mod p>
        ADD <a> <b> <a+b mod p>
        SUB <a> <b> <a-b mod p>
        PADD <x1> <y1> <x2> <y2> <inv(x2-x1)> <x3> <y3>
    """
    rng = random.Random(seed)
    p = E.p
    h = lambda v: "%064x" % v
    w = out.write
    w("# generated by ecref.py vhdl --curve %s -- do not edit\n" % E.name)
    for i in range(count):
        if i == 0:
            a, b = 0, 0
        elif i == 1:
            a, b = p - 1, p - 1
        elif i == 2:
            a, b = 1, p - 1
        elif i == 3:
            a, b = p - 1, 2
        elif i == 4:
            a, b = (1 << 255), (1 << 255)
        elif i == 5:
            a, b = p - 1, p - 2
        else:
            a, b = rng.randrange(p), rng.randrange(p)
        w("MUL %s %s %s\n" % (h(a), h(b), h(a * b % p)))
        w("ADD %s %s %s\n" % (h(a), h(b), h((a + b) % p)))
        w("SUB %s %s %s\n" % (h(a), h(b), h((a - b) % p)))
    for i in range(count // 4):
        P = E.random_point(rng)
        Q = E.random_point(rng)
        if P[0] == Q[0]:
            continue
        inv = pow(Q[0] - P[0], -1, p)
        R = E.add(P, Q)
        w("PADD %s %s %s %s %s %s %s\n" %
          (h(P[0]), h(P[1]), h(Q[0]), h(Q[1]), h(inv), h(R[0]), h(R[1])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["params", "vectors", "info", "vhdl"])
    ap.add_argument("--curve", default="secp256k1")
    ap.add_argument("--mode", default="mont", choices=["mont", "fast"])
    ap.add_argument("--count", type=int, default=64)
    args = ap.parse_args()
    E = get_curve(args.curve)
    if args.what == "info":
        print("curve", E.name, "p =", hex(E.p), "a =", hex(E.a), "b =", hex(E.b))
        print("n =", hex(E.n), "prime" if is_probable_prime(E.n) else "COMPOSITE")
        print("G =", tuple(hex(c) for c in E.G))
    elif args.what == "vhdl":
        emit_vhdl_vectors(E, sys.stdout)
    elif args.what == "params":
        emit_params(E, sys.stdout)
    else:
        emit_vectors(E, args.mode, sys.stdout, count=args.count)


if __name__ == "__main__":
    main()

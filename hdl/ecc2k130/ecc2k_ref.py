#!/usr/bin/env python3
# Reference model and test-vector generator for hdl/ecc2k130.
#
#   python3 ecc2k_ref.py check              # verify the hardware algorithms
#   python3 ecc2k_ref.py vectors > vectors_ecc2k130.txt
#   python3 ecc2k_ref.py cost               # gate counts of the constant matrices
#
# The field is GF(2^131) in the permuted type-II optimal normal basis used by
# the ECC2K-130 client (ecc2k130/codegen/field.py).  This script imports that
# model -- the same one every generated CUDA routine is verified against -- so
# the VHDL is checked against the oracle the GPU code is checked against.
#
# Coordinates: an element is a 131-bit integer whose bit (i-1) is the
# coefficient of gamma_i = zeta^i + zeta^-i, i = 1..131, zeta a primitive 263rd
# root of unity.  Hex fields in the vector file are 33 digits (132 bits, top
# bit always zero), most significant first, matching gf_t(130 downto 0).
#
# The hardware multiplier does not convolve in the normal basis.  It converts
# both operands to the optimal polynomial basis {c, c^2, ..., c^131}, c =
# gamma_1 = zeta + zeta^-1 (Bernstein-Lange), multiplies them as ordinary
# polynomials over GF(2), and converts the 261-coefficient product back.  Both
# conversions are fixed GF(2)-linear maps:
#
#   gamma_i = T_i(c), T_0 = 0, T_1 = c, T_i = c T_{i-1} + T_{i-2}   (Dickson)
#   c^k     = sum_j C(k, j) zeta^(k-2j),  C(k, j) odd  iff  j & (k-j) = 0
#
# and zeta^e + zeta^-e = gamma_{fold(e)}, fold(e) = min(e mod 263, 263 - e mod 263).
# `check` proves both maps against the field model on random elements before
# any vector is written.
#
# No type hints, camelCase identifiers, no itertools (project convention).

import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'ecc2k130', 'codegen'))

import curves   # noqa: E402
import field    # noqa: E402

M = 131
N = 2 * M + 1
ELL = 680564733841876926932320129493409985129
DP_WEIGHT = 34            # the challenge cutoff (ecc2k130/generated/eccF131.h)
WALK_DP_WEIGHT = 56       # loose cutoff so simulation walks actually report

onb = field.Onb(M)
curve = curves.Curve(onb)


# --------------------------------------------------------------------------
# coordinate helpers
# --------------------------------------------------------------------------
def fold(e):
    e %= N
    return e if e <= M else N - e


def coords(u):
    return onb.toCoords(u)


def elem(a):
    return onb.fromCoords(a)


def hexs(a):
    assert 0 <= a < (1 << M)
    return '%033x' % a


def weight(a):
    return bin(a).count('1')


def frobCoords(a, k):
    """sigma^k on coordinates: bit (i-1) -> bit (fold(i * 2^k) - 1)."""
    e = pow(2, k, N)
    r = 0
    for i in range(1, M + 1):
        if (a >> (i - 1)) & 1:
            r |= 1 << (fold(i * e) - 1)
    return r


def sigmaJ(a, s):
    """sigma^(3+s) built the way the hardware does it: sigma^3 then three
    conditional squarings by 1, 2 and 4."""
    r = frobCoords(a, 3)
    if s & 1:
        r = frobCoords(r, 1)
    if s & 2:
        r = frobCoords(r, 2)
    if s & 4:
        r = frobCoords(r, 4)
    return r


def jOf(hw):
    return 3 + ((hw >> 1) & 7)


# --------------------------------------------------------------------------
# the hardware multiplier, modelled step by step
# --------------------------------------------------------------------------
def dicksonRows():
    """T_i(c) for i = 1..M as bit-ints, bit k = coefficient of c^k."""
    rows = [0] * (M + 1)
    tm2, tm1 = 0, 2          # T_0 = 0, T_1 = c
    rows[1] = tm1
    for i in range(2, M + 1):
        rows[i] = (tm1 << 1) ^ tm2
        tm2, tm1 = tm1, rows[i]
    return rows


DICKSON = dicksonRows()


def prep(a):
    """ONB coordinates -> polynomial a'(c) with A(c) = c * a'(c); bit t of the
    result is the coefficient of c^(t+1)."""
    p = 0
    for i in range(1, M + 1):
        if (a >> (i - 1)) & 1:
            p ^= DICKSON[i]
    assert p & 1 == 0
    return p >> 1


def polyMul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def toOnbRows():
    """For t = 0..2M-2 the ONB coordinates of c^(t+2)."""
    rows = []
    for t in range(2 * M - 1):
        k = t + 2
        r = 0
        for j in range((k - 1) // 2 + 1):
            if j & (k - j) == 0:
                r ^= 1 << (fold(k - 2 * j) - 1)
        # C(k, k/2) is even for every k >= 2, so c^k never has a constant term
        rows.append(r)
    return rows


TOONB = toOnbRows()


def toOnb(h):
    r = 0
    t = 0
    while h:
        if h & 1:
            r ^= TOONB[t]
        h >>= 1
        t += 1
    return r


def hwMul(a, b):
    return toOnb(polyMul(prep(a), prep(b)))


def refMul(a, b):
    return coords(onb.mul(elem(a), elem(b)))


def symMul(a, b):
    """The direct normal-basis product (gf_mul_ref in the VHDL package):
    gamma_i gamma_j = gamma_fold(i+j) + gamma_|i-j|, gamma_0 = 0."""
    r = 0
    for i in range(1, M + 1):
        if not (a >> (i - 1)) & 1:
            continue
        for j in range(1, M + 1):
            if (b >> (j - 1)) & 1:
                r ^= 1 << (fold(i + j) - 1)
                if i != j:
                    r ^= 1 << (abs(i - j) - 1)
    return r


def hwInv(d):
    """Itoh-Tsujii with the chain 1,2,4,8,16,32,64,128,130: eight multiplies."""
    b2 = hwMul(d, frobCoords(d, 1))
    acc = b2
    for k in (2, 4, 8, 16, 32, 64):
        acc = hwMul(acc, frobCoords(acc, k))
    b130 = hwMul(b2, frobCoords(acc, 2))
    return frobCoords(b130, 1)


def hwStep(x, y):
    """One iteration R -> R + sigma^j(R), j = 3 + ((HW(x) >> 1) & 7), in the
    exact order of operations the step pipe uses."""
    hw = weight(x)
    s = (hw >> 1) & 7
    x2 = sigmaJ(x, s)
    y2 = sigmaJ(y, s)
    d = x ^ x2
    inv = hwInv(d)
    lam = hwMul(y ^ y2, inv)
    x3 = frobCoords(lam, 1) ^ lam ^ d
    y3 = hwMul(lam, x ^ x3) ^ x3 ^ y
    return x3, y3


def refStep(x, y):
    p = (elem(x), elem(y))
    j = jOf(weight(x))
    r = curve.add(p, curve.frob(p, j))
    return coords(r[0]), coords(r[1])


def randomSubgroupPoint(rng):
    while True:
        x = elem(rng.getrandbits(M))
        p = curve.pointFromX(x)
        if p is None:
            continue
        p = curve.mul(p, 4)
        if p is None:
            continue
        return coords(p[0]), coords(p[1])


# --------------------------------------------------------------------------
# checks
# --------------------------------------------------------------------------
def check(rng, rounds=48):
    one = (1 << M) - 1
    assert coords(onb.one()) == one
    assert hwMul(one, one) == one
    assert hwMul(1, 1) == 2, 'gamma_1^2 must be gamma_2'
    for _ in range(rounds):
        a = rng.getrandbits(M)
        b = rng.getrandbits(M)
        r = refMul(a, b)
        assert hwMul(a, b) == r
        assert symMul(a, b) == r
        assert hwMul(a, one) == a
        for k in (1, 2, 3, 4, 8, 64):
            assert frobCoords(a, k) == coords(onb.frob(elem(a), k))
        for s in range(8):
            assert sigmaJ(a, s) == coords(onb.frob(elem(a), 3 + s))
        if a:
            assert hwInv(a) == coords(onb.inv(elem(a)))
    for _ in range(4):
        x, y = randomSubgroupPoint(rng)
        assert hwStep(x, y) == refStep(x, y)
    return True


# --------------------------------------------------------------------------
# vectors
# --------------------------------------------------------------------------
def emitVectors(rng, out):
    one = (1 << M) - 1
    out.write('# generated by hdl/ecc2k130/ecc2k_ref.py vectors -- do not edit\n')
    specials = [0, one, 1, 2, 1 << (M - 1), (1 << M) - 2]
    pairs = []
    for a in specials:
        for b in specials:
            pairs.append((a, b))
    for i in range(1, 12):
        pairs.append((1 << (i - 1), 1 << (i - 1)))            # gamma_i^2 = gamma_2i
    while len(pairs) < 256:
        pairs.append((rng.getrandbits(M), rng.getrandbits(M)))
    for a, b in pairs:
        out.write('MUL %s %s %s\n' % (hexs(a), hexs(b), hexs(refMul(a, b))))

    # 200 = 12 full batches of 16 plus a partial one, to exercise the flush
    for _ in range(200):
        x, y = randomSubgroupPoint(rng)
        x3, y3 = refStep(x, y)
        out.write('STEP %s %s %d %s %s %d\n'
                  % (hexs(x), hexs(y), weight(x), hexs(x3), hexs(y3), weight(x3)))

    for _ in range(64):
        x, y = randomSubgroupPoint(rng)
        x0, y0 = x, y
        k = 0
        while True:
            x, y = refStep(x, y)
            k += 1
            if weight(x) <= WALK_DP_WEIGHT:
                break
        out.write('WALK %d %s %s %d %s %s\n'
                  % (WALK_DP_WEIGHT, hexs(x0), hexs(y0), k, hexs(x), hexs(y)))


def cost():
    """Two-input gate counts of the multiplier, read off the constant
    matrices.  A k-input XOR is k-1 XOR2 gates."""
    prepIn = [sum((DICKSON[i] >> k) & 1 for i in range(1, M + 1)) for k in range(1, M + 1)]
    toOnbIn = [sum((r >> k) & 1 for r in TOONB) for k in range(M)]
    prepXor = sum(k - 1 for k in prepIn if k)
    toOnbXor = sum(k - 1 for k in toOnbIn if k)
    print('prep  (gamma -> c-powers): %5d XOR2, widest output %d inputs, x2 operands'
          % (prepXor, max(prepIn)))
    print('toOnb (c-powers -> gamma): %5d XOR2, widest output %d inputs'
          % (toOnbXor, max(toOnbIn)))
    print('Karatsuba levels (gf2_kmul), n x n product over GF(2):')
    for levels in range(5):
        a, x, leaf = kmulCost(M, levels)
        print('  %d level(s): leaf %3d bits, %5d AND, %5d XOR2; per multiply %5d AND, %5d XOR2,'
              ' latency %d clk' % (levels, leaf, a, x, a, 2 * prepXor + x + toOnbXor,
                                   2 * levels + 4))


def kmulCost(n, levels):
    """(AND, XOR2, leaf width) of gf2_kmul: schoolbook leaf, unequal halves
    zero-extended, post-combine counted bit by bit."""
    if levels == 0:
        return n * n, (n - 1) * (n - 1), n
    h = (n + 1) // 2
    a, x, leaf = kmulCost(h, levels - 1)
    pre = 2 * h                                  # a0 + a1, b0 + b1
    mid = 2 * (2 * h - 1)                        # p0 + p1 + p2
    # overlaps of p0, mid << h, p2 << 2h within the 2n-1 output bits
    terms = [0] * (2 * n - 1)
    for i in range(2 * h - 1):
        for lo in (0, h, 2 * h):
            if i + lo < 2 * n - 1:
                terms[i + lo] += 1
    post = sum(t - 1 for t in terms if t)
    return 3 * a, 3 * x + pre + mid + post, leaf


def main():
    rng = random.Random(0x2C130)
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'check'
    if cmd == 'check':
        check(rng)
        print('ecc2k_ref: hardware model agrees with ecc2k130/codegen field model')
    elif cmd == 'vectors':
        check(rng, rounds=8)
        emitVectors(rng, sys.stdout)
    elif cmd == 'cost':
        cost()
    else:
        sys.stderr.write('usage: ecc2k_ref.py check|vectors|cost\n')
        sys.exit(2)


if __name__ == '__main__':
    main()

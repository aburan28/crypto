"""lambda-projective coordinates for the ECC2K-130 walk: formulas, invariance, price.

Companion to ../LAMBDA-PROJECTIVE.md.  Three things, each checked against the
generator's ONB model of GF(2^m) on y^2 + xy = x^3 + 1 (codegen/curves.py):

  1. the lambda-projective mixed (8M + 2S) and full (11M + 2S) addition
     formulas of Oliveira, Lopez, Aranha and Rodriguez-Henriquez, the
     lambda-affine addition reduced to one inversion (6M + 2S + I) and the
     Lopez-Dahab mixed addition (8M + 5S) agree with affine
     addition on random points, at m = 23, 41 and 131, and the operation
     counts are what the note prices;
  2. a walk driven from a projective representative is a function of the
     point iff its branch selector and distinguished-point predicate are
     invariant under scaling of (X, L, Z); with the weight of X in place of
     the weight of X/Z two representatives of the same point part within a
     step or two, and the predicate fires spuriously at its own rate;
  3. the price of every variant in the repository's unit (carry-less and
     ALU lane-instructions per scalar update, from THROUGHPUT-20B.md), as
     the single table the note carries.

Usage: python3 lambda_projective.py [--trials N] [--seed S] [--table]

No type hints, camelCase identifiers, no itertools (project convention).
"""
import argparse
import random
import sys

import curves
import field


# ---------------------------------------------------------------------------
# a counting wrapper around the field so the formulas price themselves

class Counted:
    """Field wrapper that counts multiplications, squarings and inversions."""

    def __init__(self, f):
        self.f = f
        self.reset()

    def reset(self):
        self.muls = 0
        self.sqrs = 0
        self.invs = 0

    def add(self, a, b):
        return self.f.add(a, b)

    def mul(self, a, b):
        self.muls += 1
        return self.f.mul(a, b)

    def sqr(self, a):
        self.sqrs += 1
        return self.f.sqr(a)

    def inv(self, a):
        self.invs += 1
        return self.f.inv(a)

    def frob(self, a, k):
        return self.f.frob(a, k)

    def one(self):
        return self.f.one()

    def counts(self):
        return (self.muls, self.sqrs, self.invs)


# ---------------------------------------------------------------------------
# lambda coordinates

def toLambdaAffine(f, p):
    """(x, y) -> (x, lambda) with lambda = x + y/x.  x = 0 is the 2-torsion point."""
    x, y = p
    if x == 0:
        raise ValueError('lambda coordinates are undefined at x = 0')
    return (x, f.add(x, f.mul(y, f.inv(x))))


def fromLambdaAffine(f, p):
    x, lam = p
    return (x, f.mul(x, f.add(lam, x)))


def toLambdaProjective(f, p, z):
    """A representative (X, L, Z) = (z x, z lambda, z) of the lambda-affine point."""
    x, lam = p
    return (f.mul(x, z), f.mul(lam, z), z)


def fromLambdaProjective(f, r):
    X, L, Z = r
    zi = f.inv(Z)
    return (f.mul(X, zi), f.mul(L, zi))


def lambdaMixedAdd(f, r, q):
    """(X_P, L_P, Z_P) + (x_Q, lambda_Q), 8M + 2S.

    A = L_P + lambda_Q Z_P,  B = (X_P + x_Q Z_P)^2,
    X = (A X_P)(A x_Q Z_P),  Z = A B Z_P,  L = (A x_Q Z_P + B)^2 + A B (L_P + Z_P).
    Requires x_P != x_Q, which the walk guarantees except on a 2^-131 event.
    """
    X1, L1, Z1 = r
    x2, l2 = q
    t1 = f.mul(l2, Z1)                       # M
    a = f.add(L1, t1)
    t2 = f.mul(x2, Z1)                       # M
    b = f.sqr(f.add(X1, t2))                 # S
    t3 = f.mul(a, t2)                        # M   A x_Q Z_P
    c = f.sqr(f.add(t3, b))                  # S
    ab = f.mul(a, b)                         # M
    z3 = f.mul(ab, Z1)                       # M
    l3 = f.add(c, f.mul(ab, f.add(L1, Z1)))  # M
    x3 = f.mul(f.mul(a, X1), t3)             # M M
    return (x3, l3, z3)


def lambdaFullAdd(f, r, s):
    """(X_P, L_P, Z_P) + (X_Q, L_Q, Z_Q), 11M + 2S, both projective.

    A = L_P Z_Q + L_Q Z_P,  B = (X_P Z_Q + X_Q Z_P)^2,
    X = (A X_P Z_Q)(A X_Q Z_P),  Z = A B Z_P Z_Q,
    L = (A X_Q Z_P + B)^2 + A B Z_Q (L_P + Z_P).
    """
    X1, L1, Z1 = r
    X2, L2, Z2 = s
    a = f.add(f.mul(L1, Z2), f.mul(L2, Z1))      # M M
    u = f.mul(X1, Z2)                            # M
    v = f.mul(X2, Z1)                            # M
    b = f.sqr(f.add(u, v))                       # S
    t3 = f.mul(a, v)                             # M   A X_Q Z_P
    c = f.sqr(f.add(t3, b))                      # S
    abz2 = f.mul(f.mul(a, b), Z2)                # M M   A B Z_Q, shared by Z and L
    z3 = f.mul(abz2, Z1)                         # M
    l3 = f.add(c, f.mul(abz2, f.add(L1, Z1)))    # M
    x3 = f.mul(f.mul(a, u), t3)                  # M M
    return (x3, l3, z3)


def lambdaAffineAdd(f, p, q):
    """(x_P, lambda_P) + (x_Q, lambda_Q) in lambda-affine form, one inversion.

    x = x_P x_Q A / B,  lambda = (x_Q A + B)^2 / (A B) + lambda_P + 1,
    A = lambda_P + lambda_Q,  B = (x_P + x_Q)^2; invert A B once.  6M + 2S + 1I.
    """
    x1, l1 = p
    x2, l2 = q
    a = f.add(l1, l2)
    b = f.sqr(f.add(x1, x2))                 # S
    inv = f.inv(f.mul(a, b))                 # M I   1 / (A B)
    invB = f.mul(a, inv)                     # M     1 / B
    t = f.mul(x2, a)                         # M     x_Q A
    x3 = f.mul(f.mul(x1, t), invB)           # M M
    l3 = f.add(f.add(f.mul(f.sqr(f.add(t, b)), inv), l1), f.one())  # S M
    return (x3, l3)


def ldMixedAdd(f, r, q):
    """Lopez-Dahab (X, Y, Z), x = X/Z, y = Y/Z^2, plus affine (x_Q, y_Q); a = 0.

    Lopez and Dahab 1999 as written.  8M + 5S on this curve: the a Z_1^2 term is
    a multiplication by the constant a (0 here, 1 on the other Koblitz curve),
    never a field product, so the count is the published best.
    """
    X1, Y1, Z1 = r
    x2, y2 = q
    z1sq = f.sqr(Z1)                         # S
    a = f.add(f.mul(y2, z1sq), Y1)           # M
    b = f.add(f.mul(x2, Z1), X1)             # M
    c = f.mul(Z1, b)                         # M
    d = f.mul(f.sqr(b), c)                   # S M   (a = 0 drops the a Z1^2 term)
    z3 = f.sqr(c)                            # S
    e = f.mul(a, c)                          # M
    x3 = f.add(f.add(f.sqr(a), d), e)        # S
    ff = f.add(x3, f.mul(x2, z3))            # M
    g = f.mul(f.add(x2, y2), f.sqr(z3))      # S M
    y3 = f.add(f.mul(f.add(e, z3), ff), g)   # M
    return (x3, y3, z3)


def fromLd(f, r):
    X, Y, Z = r
    zi = f.inv(Z)
    return (f.mul(X, zi), f.mul(Y, f.sqr(zi)))


def sigmaProjective(f, r, j):
    """The Frobenius acts coordinate-wise: (X, L, Z) -> (X^(2^j), L^(2^j), Z^(2^j))."""
    X, L, Z = r
    return (f.frob(X, j), f.frob(L, j), f.frob(Z, j))


def negProjective(f, r):
    """-(x, lambda) = (x, lambda + 1), i.e. (X, L + Z, Z)."""
    X, L, Z = r
    return (X, f.add(L, Z), Z)


def sameProjectivePoint(f, r, s):
    """[X:L:Z] = [X':L':Z'] iff the 2x2 minors vanish; no inversion needed to test."""
    X1, L1, Z1 = r
    X2, L2, Z2 = s
    return (f.mul(X1, Z2) == f.mul(X2, Z1)) and (f.mul(L1, Z2) == f.mul(L2, Z1))


# ---------------------------------------------------------------------------
# 1. formulas against affine addition

def randomAffinePoint(curve, rng):
    while True:
        p = curve.pointFromX(curve.f.randomElement(rng))
        if p is not None and p[0] != 0:
            return p


def checkFormulas(m, trials, rng):
    onb = field.Onb(m)
    curve = curves.Curve(onb)
    f = Counted(onb)
    mixedCounts = None
    fullCounts = None
    affineCounts = None
    ldCounts = None
    for _ in range(trials):
        p = randomAffinePoint(curve, rng)
        q = randomAffinePoint(curve, rng)
        if p[0] == q[0]:
            continue
        expected = curve.add(p, q)
        if expected is None or expected[0] == 0:
            continue
        lp = toLambdaAffine(onb, p)
        lq = toLambdaAffine(onb, q)

        f.reset()
        got = lambdaAffineAdd(f, lp, lq)
        c = f.counts()
        assert affineCounts in (None, c), (affineCounts, c)
        affineCounts = c
        assert fromLambdaAffine(onb, got) == expected, ('lambda-affine', m)

        f.reset()
        zl = onb.randomElement(rng) or onb.one()
        ld = (onb.mul(p[0], zl), onb.mul(p[1], onb.sqr(zl)), zl)
        got = ldMixedAdd(f, ld, q)
        c = f.counts()
        assert ldCounts in (None, c), (ldCounts, c)
        ldCounts = c
        assert fromLd(onb, got) == expected, ('lopez-dahab', m)
        zp = onb.randomElement(rng) or onb.one()
        zq = onb.randomElement(rng) or onb.one()
        rp = toLambdaProjective(onb, lp, zp)
        rq = toLambdaProjective(onb, lq, zq)
        want = toLambdaProjective(onb, toLambdaAffine(onb, expected), onb.one())

        f.reset()
        got = lambdaMixedAdd(f, rp, lq)
        c = f.counts()
        assert mixedCounts in (None, c), (mixedCounts, c)
        mixedCounts = c
        assert sameProjectivePoint(onb, got, want), ('mixed', m)
        assert fromLambdaAffine(onb, fromLambdaProjective(onb, got)) == expected

        f.reset()
        got = lambdaFullAdd(f, rp, rq)
        c = f.counts()
        assert fullCounts in (None, c), (fullCounts, c)
        fullCounts = c
        assert sameProjectivePoint(onb, got, want), ('full', m)
        assert fromLambdaAffine(onb, fromLambdaProjective(onb, got)) == expected

        # the shipping walk's addend: sigma^j(R) is projective too, and sigma
        # acts coordinate-wise on any representative
        j = rng.randrange(3, 11)
        sj = sigmaProjective(onb, rp, j)
        want2 = curve.add(p, curve.frob(p, j))
        got2 = lambdaFullAdd(onb, rp, sj)
        assert fromLambdaAffine(onb, fromLambdaProjective(onb, got2)) == want2

        # negation on a representative
        assert fromLambdaAffine(onb, fromLambdaProjective(onb, negProjective(onb, rp))) == curve.neg(p)
    return mixedCounts, fullCounts, affineCounts, ldCounts


# ---------------------------------------------------------------------------
# 2. the walk on representatives

def weightOf(onb, a):
    return onb.hammingWeight(a)


def walkStep(onb, r, selector):
    """R <- R + sigma^j(R), j = 3 + ((HW/2) mod 8), HW from `selector`."""
    hw = selector(r)
    j = 3 + ((hw // 2) % 8)
    return lambdaFullAdd(onb, r, sigmaProjective(onb, r, j))


def invariantSelector(onb):
    """HW(X/Z): the affine x recovered by an inversion, hence scale-invariant."""
    def sel(r):
        X, L, Z = r
        return weightOf(onb, onb.mul(X, onb.inv(Z)))
    return sel


def representativeSelector(onb):
    """HW(X) of whatever representative is in hand: cheap and not a class function."""
    def sel(r):
        return weightOf(onb, r[0])
    return sel


def splitStatistics(m, trials, steps, rng):
    """Walk the same point from two representatives; report when the points part."""
    onb = field.Onb(m)
    curve = curves.Curve(onb)
    invSel = invariantSelector(onb)
    repSel = representativeSelector(onb)
    invariantSplits = 0
    firstSplit = []
    # a predicate of the same flavour as HW(x) <= 34 on m = 131, with its
    # threshold set so that it fires on roughly a tenth of the points and the
    # spurious rate is visible in a few hundred trials
    dpWeight = int(m / 2 - 0.64 * m ** 0.5)
    dpAgree = 0
    dpTrue = 0
    dpSpurious = 0
    dpMissed = 0
    for _ in range(trials):
        p = randomAffinePoint(curve, rng)
        lp = toLambdaAffine(onb, p)
        c = onb.randomElement(rng) or onb.one()
        r1 = toLambdaProjective(onb, lp, onb.one())
        r2 = toLambdaProjective(onb, lp, c)

        # invariant selector: the two trails must stay on the same points
        a, b = r1, r2
        for _ in range(steps):
            a = walkStep(onb, a, invSel)
            b = walkStep(onb, b, invSel)
            if not sameProjectivePoint(onb, a, b):
                invariantSplits += 1
                break

        # representative selector: count steps until the points differ
        a, b = r1, r2
        split = None
        for t in range(1, steps + 1):
            a = walkStep(onb, a, repSel)
            b = walkStep(onb, b, repSel)
            if not sameProjectivePoint(onb, a, b):
                split = t
                break
        firstSplit.append(split if split is not None else steps + 1)

        # the predicate: HW(X) <= w on the representative against HW(x) <= w
        truth = weightOf(onb, lp[0]) <= dpWeight
        claim = weightOf(onb, r2[0]) <= dpWeight
        dpTrue += truth
        dpAgree += (truth == claim)
        dpSpurious += (claim and not truth)
        dpMissed += (truth and not claim)
    return {
        'm': m,
        'trials': trials,
        'steps': steps,
        'invariantSplits': invariantSplits,
        'meanFirstSplit': sum(firstSplit) / len(firstSplit),
        'splitWithin1': sum(1 for s in firstSplit if s == 1) / len(firstSplit),
        'splitWithin2': sum(1 for s in firstSplit if s <= 2) / len(firstSplit),
        'neverSplit': sum(1 for s in firstSplit if s > steps) / len(firstSplit),
        'dpWeight': dpWeight,
        'dpRate': dpTrue / trials,
        'dpAgree': dpAgree / trials,
        'dpSpurious': dpSpurious / trials,
        'dpMissed': dpMissed / trials,
    }


# ---------------------------------------------------------------------------
# 3. the price, in the repository's unit

# Per-routine prices from THROUGHPUT-20B.md section 3 (static SASS on the
# shipping preset, sm_120), the same source as the ITERATION-FUNCTION table:
#   product   mulPolynomial131        161 ALU slots, 6 clmad
#   squaring  spread32p + reduce       83 ALU slots, 6.25 clmad
#   inversion 8 products + 5 Frobenius networks: 8 * 360 + 1,312 ALU, 48 clmad
# and the pipe rates of ITERATION-FUNCTION section 1 on the RTX PRO 6000.
PRODUCT_ALU = 161.0
PRODUCT_CLMAD = 6.0
SQUARE_ALU = 83.0
SQUARE_CLMAD = 6.25
INVERSION_ALU = 8 * 360.0 + 1312.0
INVERSION_CLMAD = 48.0
SMS = 188
CLMAD_RATE = 1.65           # lanes per SM-clock
ALU_RATE = 62.1
CLOCKS_GHZ = (2.30, 2.43)


def price(muls, sqrs, invs):
    alu = muls * PRODUCT_ALU + sqrs * SQUARE_ALU + invs * INVERSION_ALU
    clmad = muls * PRODUCT_CLMAD + sqrs * SQUARE_CLMAD + invs * INVERSION_CLMAD
    return alu, clmad


def ceilingBs(clmad, alu):
    """Updates per second at 100% of whichever pipe binds, over the clock range."""
    out = []
    for ghz in CLOCKS_GHZ:
        byClmad = SMS * ghz * CLMAD_RATE / clmad if clmad else float('inf')
        byAlu = SMS * ghz * ALU_RATE / alu if alu else float('inf')
        out.append(min(byClmad, byAlu))
    return out


def affineBatched(batch, mixed=True):
    """Affine addition with Montgomery's trick over `batch` slots.

    lambda = e/d (1M), lambda (x + x') (1M), lambda^2 (1S), the trick's
    3 (B - 1) / B products and one inversion per batch.  The same for the
    mixed (table) and full (sigma^j R) addends: sigma is free in the normal
    basis and the addend is affine either way.
    """
    muls = 2.0 + 3.0 * (batch - 1) / batch
    return muls, 1.0, 1.0 / batch


def priceTable(mixedCounts, fullCounts, lambdaAffineCounts=(6, 2, 1), ldCounts=(8, 5, 0)):
    rows = []
    for batch in (1, 2, 3, 4, 8, 16, 32):
        mu, sq, iv = affineBatched(batch)
        alu, cl = price(mu, sq, iv)
        rows.append(('affine, Montgomery batch %d' % batch, mu, sq, iv, alu, cl))
    mu, sq, iv = affineBatched(10 ** 9)
    alu, cl = price(mu, sq, iv)
    rows.append(('affine, batch -> infinity', mu, sq, iv, alu, cl))
    am, as_, ai = lambdaAffineCounts
    mu = am + 3.0 * 15 / 16
    alu, cl = price(mu, as_, ai / 16.0)
    rows.append(('lambda-affine, Montgomery batch 16', mu, as_, ai / 16.0, alu, cl))
    lm, ls, li = ldCounts
    alu, cl = price(lm, ls, li)
    rows.append(('Lopez-Dahab mixed, hash free', lm, ls, li, alu, cl))
    mm, ms, mi = mixedCounts
    alu, cl = price(mm, ms, mi)
    rows.append(('lambda-projective mixed (table walk), hash free', mm, ms, mi, alu, cl))
    fm, fs, fi = fullCounts
    alu, cl = price(fm, fs, fi)
    rows.append(('lambda-projective full (shipping walk), hash free', fm, fs, fi, alu, cl))
    alu, cl = price(mm, ms, mi + 1)
    rows.append(('lambda-projective mixed + invariant hash by one inversion', mm, ms, mi + 1, alu, cl))
    mu = mm + 3.0 * 15 / 16
    alu, cl = price(mu, ms, mi + 1.0 / 16)
    rows.append(('lambda-projective mixed + invariant hash by batched inversion (16)', mu, ms, mi + 1.0 / 16, alu, cl))
    return rows


def formatTable(rows):
    out = []
    out.append('| variant | M | S | I | ALU / update | clmad / update | ceiling B/s (100%, 2.30 - 2.43 GHz) | binding pipe |')
    out.append('|---|---:|---:|---:|---:|---:|---:|---|')
    for name, mu, sq, iv, alu, cl in rows:
        lo, hi = ceilingBs(cl, alu)
        byClmad = SMS * CLOCKS_GHZ[0] * CLMAD_RATE / cl
        byAlu = SMS * CLOCKS_GHZ[0] * ALU_RATE / alu
        binding = 'clmad' if byClmad < byAlu else 'ALU'
        out.append('| %s | %.2f | %.2f | %.3f | %.0f | %.1f | %.1f - %.1f | %s |' % (
            name, mu, sq, iv, alu, cl, lo, hi, binding))
    return '\n'.join(out)


# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--trials', type=int, default=200)
    ap.add_argument('--seed', type=int, default=20260919)
    ap.add_argument('--steps', type=int, default=64)
    ap.add_argument('--table', action='store_true', help='print the price table only')
    args = ap.parse_args(argv)
    rng = random.Random(args.seed)

    mixedCounts = None
    fullCounts = None
    affineCounts = None
    ldCounts = None
    if not args.table:
        for m in (23, 41, 131):
            trials = args.trials if m < 131 else max(20, args.trials // 10)
            mixed, full, affine, ld = checkFormulas(m, trials, rng)
            print('m = %3d: %d random additions agree with affine; lambda-projective mixed %dM+%dS, '
                  'full %dM+%dS; lambda-affine %dM+%dS+%dI; Lopez-Dahab mixed %dM+%dS' % (
                      m, trials, mixed[0], mixed[1], full[0], full[1],
                      affine[0], affine[1], affine[2], ld[0], ld[1]))
            mixedCounts = mixedCounts or mixed
            fullCounts = fullCounts or full
            affineCounts = affineCounts or affine
            ldCounts = ldCounts or ld
            assert (mixed, full, affine, ld) == (mixedCounts, fullCounts, affineCounts, ldCounts)
        for m in (23, 41):
            s = splitStatistics(m, args.trials, args.steps, rng)
            print('m = %3d: invariant selector, trails on the same point after %d steps: %d of %d split' % (
                m, s['steps'], s['invariantSplits'], s['trials']))
            print('         representative selector HW(X): mean first split at step %.2f; '
                  'split at step 1: %.1f%%, within 2: %.1f%%, never in %d: %.1f%%' % (
                      s['meanFirstSplit'], 100 * s['splitWithin1'], 100 * s['splitWithin2'],
                      s['steps'], 100 * s['neverSplit']))
            print('         predicate HW(.) <= %d: true rate %.3f, representative agrees %.3f, '
                  'spurious %.3f, missed %.3f' % (
                      s['dpWeight'], s['dpRate'], s['dpAgree'], s['dpSpurious'], s['dpMissed']))
    if mixedCounts is None:
        mixedCounts, fullCounts, affineCounts, ldCounts = (8, 2, 0), (11, 2, 0), (6, 2, 1), (8, 5, 0)
    print()
    print(formatTable(priceTable(mixedCounts, fullCounts, affineCounts, ldCounts)))
    return 0


if __name__ == '__main__':
    sys.exit(main())

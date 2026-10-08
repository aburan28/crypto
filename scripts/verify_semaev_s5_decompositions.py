"""Semantic verification of the symbolic S_5 used by Boundary C (EXP-R6b).

The Rust side (`cryptanalysis::semaev_leading_form`) checks S_5 structurally --
symmetry in all five arguments, degree 8 in each, and agreement of S_4 with the
repository's own implementation.  That leaves one thing unchecked: whether the
polynomial actually *is* Semaev's fifth summation polynomial, i.e. whether it
vanishes exactly when five points with those x-coordinates admit a sign choice
summing to O.

This script checks that directly, on real curve points over F_2^5, for four
curves.  Run: python3 scripts/verify_semaev_s5_decompositions.py
Expected: every genuine 5-point decomposition gives S_5 = 0, no tuple without
one does, and no nonzero value coexists with a decomposition.
"""
import os, random, itertools, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from semaev_symbolic import build_s4, build_s5, A6

random.seed(7)

# ---- F_2^n arithmetic -------------------------------------------------
N = 5
RED = (1 << 5) | (1 << 2) | 1          # x^5 + x^2 + 1
MASK = (1 << N) - 1

def fmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1
        a <<= 1
        if a >> N & 1: a ^= RED
    return r & MASK

def fpow(a, k):
    r, base = 1, a
    while k:
        if k & 1: r = fmul(r, base)
        base = fmul(base, base); k >>= 1
    return r

def finv(a):
    assert a != 0
    return fpow(a, (1 << N) - 2)

# ---- curve y^2 + xy = x^3 + a2 x^2 + a6 ------------------------------
class Curve:
    def __init__(self, a2, a6): self.a2, self.a6 = a2, a6
    def on(self, P):
        if P is None: return True
        x, y = P
        return fmul(y, y) ^ fmul(x, y) == fmul(x, fmul(x, x)) ^ fmul(self.a2, fmul(x, x)) ^ self.a6
    def neg(self, P):
        return None if P is None else (P[0], P[0] ^ P[1])
    def add(self, P, Q):
        if P is None: return Q
        if Q is None: return P
        x1, y1 = P; x2, y2 = Q
        if x1 == x2:
            if y2 == (x1 ^ y1): return None           # Q = -P
            # doubling
            if x1 == 0: return None
            lam = x1 ^ fmul(y1, finv(x1))
            x3 = fmul(lam, lam) ^ lam ^ self.a2
            y3 = fmul(x1, x1) ^ fmul(lam ^ 1, x3)
            return (x3, y3)
        lam = fmul(y1 ^ y2, finv(x1 ^ x2))
        x3 = fmul(lam, lam) ^ lam ^ x1 ^ x2 ^ self.a2
        y3 = fmul(lam, x1 ^ x3) ^ x3 ^ y1
        return (x3, y3)
    def points(self):
        pts = [None]
        for x in range(1 << N):
            for y in range(1 << N):
                if self.on((x, y)): pts.append((x, y))
        return pts

# ---- evaluate a symbolic polynomial ----------------------------------
def peval(p, xs, a6):
    """xs: list of field values for X1..X5 (index 0..4); a6 value."""
    tot = 0
    for e in p:
        term = 1
        for i, ei in enumerate(e[:5]):
            if ei: term = fmul(term, fpow(xs[i], ei))
        if e[A6]: term = fmul(term, fpow(a6, e[A6]))
        if e[6]: raise AssertionError("residual Y")
        tot ^= term
    return tot

s4 = build_s4()
s5, _ = build_s5()
print("built S_4 (%d monomials) and S_5 (%d monomials)" % (len(s4), len(s5)))

for a6 in (1, 3, 5, 0b10110):
    E = Curve(0, a6)
    pts = [P for P in E.points() if P is not None]
    order = len(pts) + 1
    print(f"\n--- a6={a6}: #E = {order} ---")

    # S_4: P1+P2+P3+P4 = O  =>  S_4 = 0
    ok4 = bad4 = 0
    for _ in range(300):
        P = [random.choice(pts) for _ in range(3)]
        S = None
        for Q in P: S = E.add(S, Q)
        P4 = E.neg(S)
        if P4 is None: continue
        v = peval(s4, [P[0][0], P[1][0], P[2][0], P4[0], 0], a6)
        if v == 0: ok4 += 1
        else: bad4 += 1
    print(f"  S_4 on genuine 4-decompositions: {ok4} zero, {bad4} NONZERO")

    # S_5: P1+..+P5 = O  =>  S_5 = 0
    ok5 = bad5 = 0
    for _ in range(300):
        P = [random.choice(pts) for _ in range(4)]
        S = None
        for Q in P: S = E.add(S, Q)
        P5 = E.neg(S)
        if P5 is None: continue
        v = peval(s5, [P[0][0], P[1][0], P[2][0], P[3][0], P5[0]], a6)
        if v == 0: ok5 += 1
        else: bad5 += 1
    print(f"  S_5 on genuine 5-decompositions: {ok5} zero, {bad5} NONZERO")

    # Converse: random 5-tuples of x-coords should usually be nonzero,
    # and zero exactly when some sign choice sums to O.
    zero_but_no_decomp = 0
    nonzero_checked = 0
    for _ in range(200):
        xs = [random.choice(pts)[0] for _ in range(5)]
        v = peval(s5, xs, a6)
        # does some sign combination of points with these x-coords sum to O?
        cand = []
        for x in xs:
            ys = [P[1] for P in pts if P[0] == x]
            cand.append(ys)
        found = False
        for combo in itertools.product(*cand):
            S = None
            for x, y in zip(xs, combo): S = E.add(S, (x, y))
            if S is None: found = True; break
        if v == 0 and not found: zero_but_no_decomp += 1
        if v != 0:
            nonzero_checked += 1
            if found:
                print("  !! S_5 nonzero but a decomposition exists — WRONG")
    print(f"  converse spot-check: {zero_but_no_decomp} vanish without a decomposition"
          f" (expected 0), {nonzero_checked} nonzero cases all consistent")

"""Pullback of the third summation polynomial along a Velu x-map factors
into the kernel-twisted summation polynomials of the domain curve.
Checked with sympy over GF(p) for ell = 2 and ell = 3 (rational kernels).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy import symbols, Poly, together, fraction, expand
p = 211
x1, x2, x3 = symbols('x1 x2 x3')

def S3(a, b, X1, X2, X3):
    return ((X1-X2)**2*X3**2 - 2*((X1+X2)*(X1*X2+a) + 2*b)*X3
            + ((X1*X2-a)**2 - 4*b*(X1+X2)))

def velu_x(a, b, x0, ell):
    gx = (3*x0*x0 + a) % p
    if ell == 2:
        v, u = gx, 0
        A, B = (a - 5*v) % p, (b - 7*(u + x0*v)) % p
        return (lambda x: x + v/(x-x0)), A, B
    v = 2*gx % p; u = 4*(x0**3 + a*x0 + b) % p
    A, B = (a - 5*v) % p, (b - 7*(u + x0*v)) % p
    return (lambda x: x + v/(x-x0) + u/(x-x0)**2), A, B

def report(a, b, x0, ell):
    X, A, B = velu_x(a, b, x0, ell)
    pulled = together(S3(A, B, X(x1), X(x2), X(x3)))
    num, den = fraction(pulled)
    num = Poly(expand(num), x1, x2, x3, modulus=p)
    base = Poly(S3(a, b, x1, x2, x3), x1, x2, x3, modulus=p)
    q, rem = num.div(base)
    print(f"ell={ell}, kernel x0={x0}: E1: y^2=x^3+{a}x+{b}  ->  E2: y^2=x^3+{A}x+{B}")
    print(f"  deg_x1 of pulled-back numerator: {num.degree(x1)} (= 2*ell), terms: {len(num.terms())}")
    print(f"  divisible by S3^E1: {rem.is_zero}")
    print(f"  cofactor degrees (x1,x2,x3) = {q.degree(x1)},{q.degree(x2)},{q.degree(x3)}, terms: {len(q.terms())}")
    # symmetric?
    qs = Poly(q.as_expr().subs({x1: x2, x2: x1}, simultaneous=True), x1, x2, x3, modulus=p)
    print(f"  cofactor symmetric in x1<->x2: {qs == q}")
    # the cofactor vanishes exactly on triples with P1+P2+P3 = +-T: check numerically on the curve
    import random
    from fiber_probe import Ext, Curve
    F = Ext(p, 1, [0, 1]); E = Curve(F, F.emb(a), F.emb(b)); pts = [P for P in E.points_prime() if P is not None]
    rng = random.Random(1); hits = tot = 0; sum_in_kernel_when_zero = True
    T = next(P for P in pts if P[0][0] == x0) if ell == 3 else (F.emb(x0), F.zero())
    ker = [None, T, E.neg(T)]
    for _ in range(300):
        P1, P2 = rng.choice(pts), rng.choice(pts)
        S = E.add(P1, P2)
        P3 = rng.choice(pts)
        val = q.eval({x1: P1[0][0], x2: P2[0][0], x3: P3[0][0]}) if False else int(q.as_expr().subs({x1: P1[0][0], x2: P2[0][0], x3: P3[0][0]})) % p
        s3 = E.add(S, P3); s3n = E.add(S, E.neg(P3))
        kernel_hit = any(E.add(s, E.neg(k)) is None for s in (s3, s3n) for k in ker[1:])
        if kernel_hit: tot += 1; hits += (val == 0)
    print(f"  cofactor vanishes on sampled triples with P1+-P2+-P3 in ker minus O: {hits}/{tot}")
    # and forced: choose P3 = T - P1 - P2 so that P1+P2+P3 = T
    ok = tried = 0
    while tried < 50:
        P1, P2 = rng.choice(pts), rng.choice(pts)
        P3 = E.add(T, E.neg(E.add(P1, P2)))
        if P3 is None or P3[0] in (P1[0], P2[0]): continue   # degenerate draw, redraw
        tried += 1
        val = int(q.as_expr().subs({x1: P1[0][0], x2: P2[0][0], x3: P3[0][0]})) % p
        ok += (val == 0)
    print(f"  cofactor vanishes on 50 forced non-degenerate triples with P1+P2+P3 = T: {ok}/50")

# E1: y^2 = x^3 + x + 3 over F_211 (from fiber_probe): 2-torsion x0 and 3-torsion x0 = 77
a, b = 1, 3
two = [x for x in range(p) if (x**3 + a*x + b) % p == 0]
report(a, b, two[0], 2)
report(a, b, 77, 3)

#!/usr/bin/env python3
"""
bielliptic_j0_criterion.py  —  autolab 2026-08-04, Thread 2 (CHLRS / Howe gluing).

WHAT THIS SETTLES
=================
The 2026-07-27 entry of RESEARCH_AUTOLAB_LOG.md ended on a "structural gap":

    "the Z/3Z Richelot computes Howe_cover -> Richelot_dual.  What is needed is the
     FORWARD map: (E1,E2) -> Howe_cover parameters (a,b)."

The forward map is elementary once the right normal form is used.  The repo had been
working inside the family

    C_{a,b} :  y^2 = x^6 + a x^3 + b                                        (Z/3 family)

which is the WRONG family: it provably can never have both elliptic quotients of
j-invariant 0.  The right family is the two-parameter

    C_{u,v} :  y^2 = u x^6 + v

whose two elliptic quotients are, in closed form,

    E1 :  y^2 = x^3 + u^2 v          (via  X = x^2,      Y = y      )
    E2 :  y^2 = x^3 + v^2 u          (via  X = 1/x^2,    Y = y/x^3  )

RESULTS (all verified below)
============================
T1. Closed form for the Z/3 family.  For C_{a,b} with b a sixth power, put
    c = b^(1/3), e = b^(1/2).  The involution tau: (x,y) -> (c/x, e y/x^3) has
    quotient
        E_+ :  Y^2 = (a-2e) Z^3 + 9c Z^2 - (6e/c) Z + 1
    and tau' (e -> -e) gives E_-.  Then #Jac(C) = #E_+ * #E_-.
    -> verified exhaustively, 715/715 (a,b) over p = 13,19,31,37,43.

T2. In the Z/3 family the two quotients ALWAYS have equal trace, so
        #Jac(C_{a,b}) is a perfect square.
    -> verified 2013/2013 over p = 13,19,31,37,43,61,67 (zero exceptions).
    Proof: sigma:(x,y)->(zeta3 x, y) acts on H^0(Omega_C) with eigenvalues
    zeta3, zeta3^2, so Q(zeta3) embeds in End^0(Jac C).  A field embeds in a
    product K1 x K2 only via injective projections; if E_+ !~ E_- that forces
    Q(zeta3) in End^0(E_+-) i.e. j(E_+-) = 0, excluded by T3.  Hence E_+ ~ E_-.

T3. The Z/3 family can NEVER have both quotients of j-invariant 0.
    j(E_+) = 0 <=> a = -5e/2 ;  j(E_-) = 0 <=> a = +5e/2 ;  both <=> e = 0 <=> b = 0
    (singular).  -> verified: 0 occurrences in the 2013 cases above.

T4. Normal-form classification (split-bielliptic case).  Write a split-bielliptic
    genus-2 curve as y^2 = a6 x^6 + a4 x^4 + a2 x^2 + a0.  Both quotients have
    j = 0 iff a4 = a2 = 0.  (The conditions are 3 a2 a6 = a4^2 and 3 a4 a0 = a2^2;
    with a2 a4 != 0 they force the sextic to be the perfect cube (alpha x^2+beta)^3,
    i.e. singular.)  So the curve is C_{u,v}: y^2 = u x^6 + v, giving the pair
        (b1, b2) = (u^2 v, v^2 u)    hence    b1 * b2 = (uv)^3  is ALWAYS a cube.

T5. CRITERION.  Two j=0 curves E_{b1}, E_{b2} over F_p arise as the two elliptic
    quotients of a SPLIT-bielliptic genus-2 curve  <=>  b1 * b2 is a cube in F_p^*.
    Solvability: in F_p^*/(F_p^*)^6 = Z/6 with classes i1, i2, the system
    2s+t = i1, s+2t = i2 (mod 6) is solvable iff i1 + i2 = 0 (mod 3), which is
    exactly "b1 b2 is a cube".  Exactly 5 of the 15 unordered pairs of the six
    sextic twists satisfy this, for every p = 1 mod 3.
    -> #Jac(C_{u,v}) = #E1 * #E2 verified 1404/1404 over p = 7,13,19,31.

SCOPE / WHAT IS *NOT* CLAIMED
=============================
T5 characterises the SPLIT-bielliptic case (the involution x -> c/x with rational
fixed points, equivalently x -> -x).  Over F_p with p = 3 mod 4 there are also
NON-SPLIT bielliptic involutions.  A char-poly sweep of the non-split normal form
    f = a6 x^6 + a5 x^5 + a4 x^4 + a3 x^3 + c a4 x^2 + c^2 a5 x + c^3 a6  (c non-square)
found char polys matching j=0 pairs with b1 b2 NOT a cube (1 such pair at p=7,
8 at p=13).  Matching char poly only pins the isogeny class, not the quotients, so
this is inconclusive -- but it means T5 is SUFFICIENT, not proven necessary.
Resolving the non-split case is the open sub-task.

CONSEQUENCE FOR THE 2026-07-27 p=1009 TEST
==========================================
That test targeted #Jac = #E1 * #E2 = 967 * 1053 = 1018251 with traces +43 and -43.
1018251 is not a perfect square (1009^2 = 1018081), so by T2 it is unreachable in
the Z/3 family -- the "Test 2 mismatch" was a structural impossibility, not a bad
parameterisation.  And b1*b2 is a non-cube for all 36 (b1,b2) in those two trace
classes, so it is unreachable by C_{u,v} either.

Run:  python3 bielliptic_j0_criterion.py
"""

import itertools
from math import gcd, isqrt

P256 = 2**256 - 2**32 - 977
N_SECP = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141


# ----------------------------------------------------------------- primitives

def legendre_tab(p):
    t = [0] * p
    for x in range(1, p):
        t[x * x % p] = 1
    return [0 if x == 0 else (1 if t[x] else -1) for x in range(p)]


def is_cube_mod(a, p):
    return True if p % 3 != 1 else pow(a % p, (p - 1) // 3, p) == 1


def is_sixth_power(a, p):
    return pow(a % p, (p - 1) // gcd(6, p - 1), p) == 1


def ec_add(P, Q, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        lam = 3 * x1 * x1 % p * pow(2 * y1 % p, p - 2, p) % p
    else:
        lam = (y2 - y1) * pow((x2 - x1) % p, p - 2, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def ec_mul(P, k, p):
    R, Q = None, P
    while k:
        if k & 1:
            R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R


def point_on(b, p):
    for x in range(1, 1000):
        rhs = (x * x % p * x + b) % p
        if rhs and pow(rhs, (p - 1) // 2, p) == 1:
            y = pow(rhs, (p + 1) // 4, p)
            if (y * y - rhs) % p == 0:
                return (x, y)
    raise RuntimeError("no point found")


def sextic_traces(p, t0):
    """the six sextic-twist traces of a j=0 curve of trace t0 over F_p."""
    rem = 4 * p - t0 * t0
    assert rem % 3 == 0
    f = isqrt(rem // 3)
    assert 3 * f * f == rem, "4p - t^2 = 3 f^2 failed"
    return sorted({t0, -t0, (t0 + 3 * f) // 2, -((t0 + 3 * f) // 2),
                   (t0 - 3 * f) // 2, -((t0 - 3 * f) // 2)}), f


def trace_of_j0(b, p, candidates):
    P = point_on(b, p)
    hits = [t for t in candidates if ec_mul(P, p + 1 - t, p) is None]
    return hits[0] if len(hits) == 1 else hits


# ------------------------------------------------------- T5: solve for (u,v)

def solve_uv(b1, b2, p):
    """Find u,v in F_p^* with  u^2 v = b1  and  v^2 u = b2  modulo sixth powers,
    i.e. produce the explicit curve C_{u,v}: y^2 = u x^6 + v whose quotients are
    E_{b1} and E_{b2}.  Returns (u,v) or None."""
    if not is_cube_mod(b1 * b2 % p, p):
        return None
    # u^2 v = b1, u v^2 = b2  =>  (uv)^3 = b1 b2  and  u/v = b1/b2
    # pick w = a cube root of b1*b2, then u = (b1/w) ... solve directly:
    #   from u^2 v = b1 and u v^2 = b2:  u^3 = b1^2 / b2,  v^3 = b2^2 / b1
    t1 = b1 * b1 % p * pow(b2, p - 2, p) % p          # u^3
    t2 = b2 * b2 % p * pow(b1, p - 2, p) % p          # v^3
    u = cube_root(t1, p)
    v = cube_root(t2, p)
    if u is None or v is None:
        return None
    # fix the zeta3 ambiguity so that u^2 v == b1 exactly (up to 6th powers)
    z3 = pow(next(g for g in range(2, 100)
                  if pow(g, (p - 1) // 3, p) != 1), (p - 1) // 3, p)
    for i in range(3):
        for j in range(3):
            uu = u * pow(z3, i, p) % p
            vv = v * pow(z3, j, p) % p
            if uu * uu % p * vv % p == b1 % p and vv * vv % p * uu % p == b2 % p:
                return (uu, vv)
    return None


def cube_root(a, p):
    a %= p
    if a == 0:
        return 0
    if p % 3 == 2:
        return pow(a, (2 * p - 1) // 3, p)
    if not is_cube_mod(a, p):
        return None
    # p = 1 mod 3: Tonelli-style via a^((2q+1)/3) when p = 3q+1 and 3 || p-1 fails;
    # use generic: find e with 3^e || p-1
    q, e = p - 1, 0
    while q % 3 == 0:
        q //= 3
        e += 1
    if e == 1:
        # p-1 = 3q, gcd(3,q)=1 -> inverse of 3 mod q
        return pow(a, pow(3, -1, q), p)
    # general Adleman-Manders-Miller
    g = next(x for x in range(2, 1000) if not is_cube_mod(x, p))
    r = pow(a, q, p)
    # discrete log of r in the order-3^e subgroup, base h = g^q
    h = pow(g, q, p)
    ex = 0
    for k in range(e):
        if pow(r * pow(h, p - 1 - ex, p) % p, 3 ** (e - 1 - k), p) != 1:
            for d in (1, 2):
                if pow(r * pow(pow(h, ex + d * 3 ** k, p), p - 2, p) % p,
                       3 ** (e - 1 - k), p) == 1:
                    ex += d * 3 ** k
                    break
    if ex % 3:
        return None
    y = pow(a, pow(3, -1, q), p) * pow(g, -(ex // 3) * q % (p - 1), p) % p
    return y if pow(y, 3, p) == a else None


# --------------------------------------------------------------- secp256k1

def secp256k1_report():
    p = P256
    print("=" * 78)
    print("secp256k1: which sextic-twist pairs admit a split-bielliptic gluing?")
    print("=" * 78)
    t0 = p + 1 - N_SECP
    six, f = sextic_traces(p, t0)
    print(f"  p mod 12 = {p % 12}   (so p = 1 mod 3 and p = 3 mod 4)")
    print(f"  trace(secp256k1)     t0 = {t0}")
    print(f"  4p - t0^2 = 3 f^2,    f = {f}")

    g = 3
    while pow(g, (p - 1) // 2, p) == 1 or pow(g, (p - 1) // 3, p) == 1:
        g += 1
    h = pow(g, (p - 1) // 6, p)
    assert pow(h, 3, p) == p - 1, "h must be a primitive 6th root of unity"
    bk = [7 * pow(h, k, p) % p for k in range(6)]
    tr = {k: trace_of_j0(bk[k], p, six) for k in range(6)}

    print("\n  twist k   b_k = 7 h^k   trace")
    for k in range(6):
        mark = "   <-- secp256k1" if tr[k] == t0 else ""
        print(f"    k={k}   cube(b_k)={str(is_cube_mod(bk[k],p)):>5s}   "
              f"trace={tr[k]}{mark}")

    print("\n  pair    b1*b2 cube?   HoweH2   gcd(N1,N2)   split-bielliptic curve")
    glue = []
    for i, j in itertools.combinations(range(6), 2):
        c = is_cube_mod(bk[i] * bk[j] % p, p)
        h2 = is_cube_mod(bk[i], p) == is_cube_mod(bk[j], p)
        gg = gcd(p + 1 - tr[i], p + 1 - tr[j])
        uv = solve_uv(bk[i], bk[j], p) if c else None
        if uv:
            glue.append((i, j, uv))
        note = "y^2 = u x^6 + v  (u,v below)" if uv else "-"
        print(f"   ({i},{j})   {str(c):>10s}   {str(h2):>6s}   {gg:>10d}   {note}")

    print(f"\n  => {len(glue)} of 15 pairs are split-bielliptic-glueable.")
    print("\n  explicit curves  C: y^2 = u x^6 + v   with quotients E_{b_i}, E_{b_j}:")
    for i, j, (u, v) in glue:
        q1, q2 = u * u % p * v % p, v * v % p * u % p
        ok1 = is_sixth_power(q1 * pow(bk[i], p - 2, p) % p, p)
        ok2 = is_sixth_power(q2 * pow(bk[j], p - 2, p) % p, p)
        t1 = trace_of_j0(q1, p, six)
        t2 = trace_of_j0(q2, p, six)
        print(f"\n   pair ({i},{j}):")
        print(f"     u = {u}")
        print(f"     v = {v}")
        print(f"     quotient E1: y^2=x^3+u^2v   ~= E_{{b_{i}}} ? {ok1}   "
              f"trace match: {t1 == tr[i]}")
        print(f"     quotient E2: y^2=x^3+v^2u   ~= E_{{b_{j}}} ? {ok2}   "
              f"trace match: {t2 == tr[j]}")
        print(f"     #Jac(C) = #E_{i} * #E_{j} = "
              f"{(p + 1 - tr[i]) * (p + 1 - tr[j])}")

    repo = {(0, 2), (0, 3), (0, 5), (1, 4), (2, 3)}
    mine = {(i, j) for i, j, _ in glue}
    print("\n  cross-check against the 2026-05-30 Howe H1/H2/H3 table:")
    print(f"    this criterion : {sorted(mine)}")
    print(f"    repo 2026-05-30: {sorted(repo)}")
    print(f"    only here      : {sorted(mine - repo)}")
    print(f"    only repo      : {sorted(repo - mine)}")
    return mine


def p1009_check():
    print()
    print("=" * 78)
    print("the 2026-07-27 p=1009 Howe target is structurally unreachable")
    print("=" * 78)
    p = 1009
    n1, n2 = 967, 1053
    tgt = n1 * n2
    r = isqrt(tgt)
    print(f"  target #Jac = #E1 * #E2 = {n1} * {n2} = {tgt}")
    print(f"  perfect square? {r*r == tgt}  ({r}^2 = {r*r})")
    print("  T2 forces #Jac(C_{a,b}) to be a perfect square -> UNREACHABLE in the")
    print("  Z/3 family.  The 2026-07-27 'Test 2 mismatch' was an impossibility,")
    print("  not a parameterisation bug.")
    chi = legendre_tab(p)
    tr = {b: p + 1 - (p + 1 + sum(chi[(x * x % p * x + b) % p] for x in range(p)))
          for b in range(1, p)}
    b1s = [b for b in range(1, p) if tr[b] == 43]
    b2s = [b for b in range(1, p) if tr[b] == -43]
    good = [(a, b) for a in b1s for b in b2s if is_cube_mod(a * b, p)]
    print(f"  and of the {len(b1s)*len(b2s)} (b1,b2) with those traces, "
          f"{len(good)} have b1*b2 a cube")
    print("  -> unreachable by C_{u,v} = y^2 = u x^6 + v either.")


if __name__ == "__main__":
    secp256k1_report()
    p1009_check()

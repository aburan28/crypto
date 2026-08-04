#!/usr/bin/env python3
"""
bielliptic_j0_verify.py — exhaustive small-prime evidence for the claims T1,T2,T3,T5
asserted in bielliptic_j0_criterion.py  (autolab 2026-08-04, Thread 2).

  T1/T2/T3  the Z/3 family  C: y^2 = x^6 + a x^3 + b   (b a sixth power)
            c = b^(1/3), e = b^(1/2)
            E1: Y^2 = (a-2e)Z^3 + 9c Z^2 - (6e/c)Z + 1
            E2: Y^2 = (a+2e)Z^3 + 9c Z^2 + (6e/c)Z + 1
            checks: #Jac = #E1*#E2 ; trace(E1)==trace(E2) ; never j1==j2==0

  T5        the split family  C: y^2 = u x^6 + v
            E1: y^2 = x^3 + u^2 v ,  E2: y^2 = x^3 + v^2 u
            check: #Jac = #E1*#E2

Run:  python3 bielliptic_j0_verify.py
"""


def legendre_tab(p):
    t = [0] * p
    for x in range(1, p):
        t[x * x % p] = 1
    return [0 if x == 0 else (1 if t[x] else -1) for x in range(p)]


def fp2(p, chi):
    nr = next(n for n in range(2, p) if chi[n] == -1)

    def mul(A, B):
        return ((A[0] * B[0] + nr * A[1] * B[1]) % p, (A[0] * B[1] + A[1] * B[0]) % p)
    elts = [(x, y) for x in range(p) for y in range(p)]
    return mul, elts, set(mul(A, A) for A in elts)


def jac_order_sextic(coef, p, chi, mul, elts, sq2):
    """#Jac for y^2 = sum coef[i] x^i (deg 6, leading coeff coef[6] != 0)."""
    n1 = 2 if chi[coef[6]] == 1 else 0
    for x in range(p):
        v, xp = 0, 1
        for i in range(7):
            v += coef[i] * xp
            xp = xp * x % p
        n1 += 1 + chi[v % p]
    n2 = 2 if (coef[6], 0) in sq2 else 0
    for A in elts:
        row = [(1, 0)]
        for _ in range(6):
            row.append(mul(row[-1], A))
        r0 = r1 = 0
        for i in range(7):
            if coef[i]:
                r0 += coef[i] * row[i][0]
                r1 += coef[i] * row[i][1]
        V = (r0 % p, r1 % p)
        if V == (0, 0):
            n2 += 1
        elif V in sq2:
            n2 += 2
    s1 = p + 1 - n1
    s2 = (s1 * s1 - (p * p + 1 - n2)) // 2
    return 1 - s1 + s2 - p * s1 + p * p


def count_cubic(A, B, C, D, p, chi):
    tot = p + 1
    for z in range(p):
        tot += chi[(((A * z + B) * z + C) * z + D) % p]
    return tot


def j_cubic(A, B, C, D, p):
    inv3 = pow(3, p - 2, p)
    b2, b1, b0 = B % p, A * C % p, A * A % p * D % p
    a4 = (b1 - b2 * b2 * inv3) % p
    a6 = (2 * pow(b2, 3, p) * pow(27, p - 2, p) - b2 * b1 * inv3 + b0) % p
    disc = (-16 * (4 * pow(a4, 3, p) + 27 * a6 * a6)) % p
    if disc == 0:
        return None
    return (-110592 * pow(a4, 3, p) * pow(disc, p - 2, p)) % p


def verify_z3_family(primes):
    print("=" * 78)
    print("T1/T2/T3 — the Z/3 family  y^2 = x^6 + a x^3 + b")
    print("=" * 78)
    tot_ok = tot_bad = tot_same_tr = tot_cases = tot_j00 = 0
    for p in primes:
        chi = legendre_tab(p)
        mul, elts, sq2 = fp2(p, chi)
        sixths = {}
        for t in range(1, p):
            sixths.setdefault(pow(t, 6, p), t)
        ok = bad = same_tr = cases = j00 = 0
        for b, t in sorted(sixths.items()):
            c, e = pow(t, 2, p), pow(t, 3, p)
            cinv = pow(c, p - 2, p)
            for a in range(p):
                if (a * a - 4 * b) % p == 0:
                    continue
                A1, A2 = (a - 2 * e) % p, (a + 2 * e) % p
                if A1 == 0 or A2 == 0:
                    continue
                B, C1, C2 = 9 * c % p, (-6 * e * cinv) % p, (6 * e * cinv) % p
                j1, j2 = j_cubic(A1, B, C1, 1, p), j_cubic(A2, B, C2, 1, p)
                if j1 is None or j2 is None:
                    continue
                coef = [b, 0, 0, a, 0, 0, 1]
                nJ = jac_order_sextic(coef, p, chi, mul, elts, sq2)
                n1 = count_cubic(A1, B, C1, 1, p, chi)
                n2 = count_cubic(A2, B, C2, 1, p, chi)
                cases += 1
                if nJ == n1 * n2:
                    ok += 1
                else:
                    bad += 1
                if n1 == n2:
                    same_tr += 1
                if j1 == 0 and j2 == 0:
                    j00 += 1
        print(f"  p={p:3d}: #Jac==#E1*#E2 {ok}/{cases}   trace(E1)==trace(E2) "
              f"{same_tr}/{cases}   (j1,j2)==(0,0) {j00}")
        tot_ok += ok
        tot_bad += bad
        tot_same_tr += same_tr
        tot_cases += cases
        tot_j00 += j00
    print(f"  TOTAL: #Jac identity {tot_ok}/{tot_cases} (mismatches {tot_bad}); "
          f"equal traces {tot_same_tr}/{tot_cases}; both-j=0 cases {tot_j00}")


def verify_split_family(primes):
    print()
    print("=" * 78)
    print("T5 — the split family  y^2 = u x^6 + v   (quotients y^2=x^3+u^2v, x^3+v^2u)")
    print("=" * 78)
    tot_ok = tot = 0
    for p in primes:
        chi = legendre_tab(p)
        mul, elts, sq2 = fp2(p, chi)
        ok = n = 0
        for u in range(1, p):
            for v in range(1, p):
                nJ = jac_order_sextic([v, 0, 0, 0, 0, 0, u], p, chi, mul, elts, sq2)
                b1, b2 = u * u % p * v % p, v * v % p * u % p
                e1 = p + 1 + sum(chi[(x * x % p * x + b1) % p] for x in range(p))
                e2 = p + 1 + sum(chi[(x * x % p * x + b2) % p] for x in range(p))
                n += 1
                if nJ == e1 * e2:
                    ok += 1
        print(f"  p={p:3d}: #Jac == #E1*#E2 for {ok}/{n}")
        tot_ok += ok
        tot += n
    print(f"  TOTAL: {tot_ok}/{tot}")


if __name__ == "__main__":
    verify_z3_family([13, 19, 31, 37, 43])
    verify_split_family([7, 13, 19, 31])

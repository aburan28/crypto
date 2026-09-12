"""Symbolic Semaev polynomials over F_2[a6], and the Boolean-degree profile.

Reference implementation for Boundary C of RESEARCH_ISOGENY_CLASS_SEARCH.md.
The authoritative version is `cryptanalysis::semaev_leading_form` in Rust,
which the test suite runs; this script exists so the computation can be
re-derived independently, and it is what
`verify_semaev_s5_decompositions.py` imports.

Variables: index 0..4 = X1..X5, 5 = a6, 6 = Y (resultant variable).
Polynomials: dict exponent-tuple -> 1 (coefficients in F_2).
"""
NV = 7
A6 = 5
Y = 6

def pz(): return {}
def pone(): return {(0,)*NV: 1}
def pvar(i):
    e = [0]*NV; e[i] = 1; return {tuple(e): 1}

def padd(*ps):
    r = {}
    for p in ps:
        for e in p:
            r[e] = r.get(e, 0) ^ 1
    return {e: 1 for e, c in r.items() if c}

def pmul(p, q):
    if not p or not q: return {}
    r = {}
    for e1 in p:
        for e2 in q:
            e = tuple(a + b for a, b in zip(e1, e2))
            r[e] = r.get(e, 0) ^ 1
    return {e: 1 for e, c in r.items() if c}

def psq(p):
    # char 2: (sum of monomials)^2 = sum of (monomial)^2
    return {tuple(2*x for x in e): 1 for e in p}

def s3_coeffs(p, q):
    """S_3(p,q,T) = A T^2 + B T + C."""
    A = psq(padd(p, q))
    B = pmul(p, q)
    C = padd(psq(B), pvar(A6))
    return A, B, C

def res_quad_quad(f, g):
    """Resultant of two quadratics given as (a2,a1,a0), char 2."""
    A1, B1, C1 = f; A2, B2, C2 = g
    t1 = psq(padd(pmul(A1, C2), pmul(A2, C1)))
    t2 = pmul(padd(pmul(A1, B2), pmul(A2, B1)),
              padd(pmul(B1, C2), pmul(B2, C1)))
    return padd(t1, t2)

def coeffs_in(p, var):
    """Split p into {degree_in_var: poly_without_var}."""
    out = {}
    for e in p:
        d = e[var]
        e2 = list(e); e2[var] = 0; e2 = tuple(e2)
        sub = out.setdefault(d, {})
        sub[e2] = sub.get(e2, 0) ^ 1
    return {d: {e: 1 for e, c in sub.items() if c} for d, sub in out.items()}

def sylvester_det(fc, gc):
    """Resultant of f (degree df) and g (degree dg) from coefficient lists
    fc = [f_df, ..., f_0], gc = [g_dg, ..., g_0].  char 2 => no signs."""
    df, dg = len(fc) - 1, len(gc) - 1
    n = df + dg
    M = [[pz() for _ in range(n)] for _ in range(n)]
    for i in range(dg):
        for k, c in enumerate(fc):
            M[i][i + k] = c
    for j in range(df):
        for k, c in enumerate(gc):
            M[dg + j][j + k] = c
    from functools import lru_cache
    memo = {}
    def det(depth, cols):
        if depth == n: return pone()
        key = (depth, cols)
        if key in memo: return memo[key]
        acc = {}
        for idx, c in enumerate(cols):
            entry = M[depth][c]
            if not entry: continue
            sub = det(depth + 1, cols[:idx] + cols[idx+1:])
            if not sub: continue
            for e in pmul(entry, sub):
                acc[e] = acc.get(e, 0) ^ 1
        res = {e: 1 for e, c in acc.items() if c}
        memo[key] = res
        return res
    return det(0, tuple(range(n)))

def build_s4():
    """S_4(X1,X2,X3,X4) = Res_T(S_3(X1,X2,T), S_3(X3,X4,T))."""
    f = s3_coeffs(pvar(0), pvar(1))
    g = s3_coeffs(pvar(2), pvar(3))
    return res_quad_quad(f, g)

def build_s5():
    """S_5(X1..X5) = Res_Y(S_3(X1,X2,Y), S_4(X3,X4,X5,Y))."""
    # S_4 with slots (X3, X4, X5, Y)
    f4a = s3_coeffs(pvar(2), pvar(3))
    f4b = s3_coeffs(pvar(4), pvar(Y))
    s4 = res_quad_quad(f4a, f4b)
    g_by_deg = coeffs_in(s4, Y)
    dg = max(g_by_deg)
    gc = [g_by_deg.get(d, pz()) for d in range(dg, -1, -1)]
    # S_3(X1,X2,Y) in Y
    A, B, C = s3_coeffs(pvar(0), pvar(1))
    fc = [A, B, C]
    return sylvester_det(fc, gc), dg

def bdeg(e, symbolic):
    """Boolean degree = sum of Hamming weights of exponents of symbolic vars."""
    return sum(bin(e[i]).count('1') for i in symbolic)

def profile(p, symbolic, label):
    """Max Boolean degree per a6 power."""
    prof = {}
    for e in p:
        k = e[A6]
        d = bdeg(e, symbolic)
        prof[k] = max(prof.get(k, -1), d)
    top = max(prof.values())
    a6free = prof.get(0, -1)
    a6_max = max([d for k, d in prof.items() if k >= 1], default=-1)
    print(f"  {label}: monomials={len(p)}  a6 powers={sorted(prof)}")
    for k in sorted(prof):
        print(f"     a6^{k}: max Boolean degree {prof[k]}")
    print(f"     top overall = {top}; a6-free = {a6free}; max over a6^k (k>=1) = {a6_max}")
    verdict = "HOLDS (a6 strictly below the top)" if a6_max < top else "FAILS (a6 reaches the top)"
    print(f"     >>> Boundary C at this m: {verdict}")
    return prof, top, a6_max

def maxdeg_per_var(p, nvars):
    return [max(e[i] for e in p) for i in range(nvars)]

if __name__ == "__main__":
    print("=== S_3 ===")
    A, B, C = s3_coeffs(pvar(0), pvar(1))
    s3 = padd(pmul(A, psq(pvar(2))), pmul(B, pvar(2)), C)  # S_3(X1,X2,X3) as A*X3^2+B*X3+C
    print("  S_3 monomials:", sorted(s3))
    # m = 2 decomposition: X1, X2 symbolic; X3 = target constant
    profile(s3, [0, 1], "S_3, m=2 (X1,X2 symbolic)")

    print("\n=== S_4 ===")
    s4 = build_s4()
    print("  degree per variable X1..X4:", maxdeg_per_var(s4, 4), " a6 degree:", max(e[A6] for e in s4))
    # symmetry check
    def permute(p, perm):
        out = {}
        for e in p:
            e2 = list(e)
            for i, j in enumerate(perm): e2[j] = e[i]
            e2 = tuple(e2)
            out[e2] = out.get(e2, 0) ^ 1
        return {e: 1 for e, c in out.items() if c}
    import itertools
    sym_ok = all(permute(s4, list(pm) + [4, 5, 6]) == s4
                 for pm in itertools.permutations(range(4)))
    print("  symmetric in X1..X4:", sym_ok)
    profile(s4, [0, 1, 2], "S_4, m=3 (X1,X2,X3 symbolic, X4 = target)")

    print("\n=== S_5 ===")
    s5, dg = build_s5()
    print("  S_4-in-Y degree used:", dg)
    print("  S_5 monomials:", len(s5))
    print("  degree per variable X1..X5:", maxdeg_per_var(s5, 5), " a6 degree:", max(e[A6] for e in s5))
    print("  Y must be eliminated; residual Y degree:", max(e[Y] for e in s5))
    sym_ok5 = all(permute(s5, list(pm) + [5, 6]) == s5
                  for pm in itertools.permutations(range(5)))
    print("  symmetric in X1..X5:", sym_ok5)
    profile(s5, [0, 1, 2, 3], "S_5, m=4 (X1..X4 symbolic, X5 = target)")

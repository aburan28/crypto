"""Is there a Frobenius-invariant F_2-subspace of x-values lying inside 4E?

For K_0: y^2 + x y = x^3 + 1 over GF(2^n), n odd, 2-part of E(F) cyclic
of order 4.  Closed-form test (derived in section 14.2 of the note):
    P = (x, y) in 2E  <=>  Tr(x) = 0
    P in 4E           <=>  Tr(x) = 0  and  Tr(y) + Tr(x * H(x)) = 0
with H the half-trace.  Step 1 validates this against scalar
multiplication at n = 19.  Step 2 enumerates every invariant subspace of
GF(2^31) inside the trace-zero hyperplane spanned by one or two of the six
5-dimensional cyclotomic components, and reports the fraction of its
points lying in 4E (1/2 for a "random" subspace; 1 would be a linear
index-4 factor base usable by the algebraic oracles).
"""
import json, random, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from image_vs_fiber_census import K0, gmul, half_trace, factor

def tr(E, a):
    t, c = 0, a
    for _ in range(E.n):
        t ^= c
        c = gmul(c, c, E.f, E.n)
    return t & 1

def in_4E_closed_form(E, P):
    x, y = P
    if tr(E, x) != 0:
        return False
    return (tr(E, y) ^ tr(E, gmul(x, half_trace(x, E.f, E.n), E.f, E.n))) == 0

def in_4E_smul(E, P, odd):
    Q = E.smul(odd, P)          # 2-primary component, in a cyclic group of order 4
    return Q is None

def validate(n=19, samples=3000, seed=2):
    E = K0(n, 0)
    two = 1
    while E.N % (2 * two) == 0:
        two *= 2
    assert two == 4
    odd = E.N // 4
    rng = random.Random(seed)
    agree = tot = in4 = 0
    while tot < samples:
        L = E.lift(rng.randrange(1, 1 << n))
        if not L:
            continue
        for P in L:
            a, b = in_4E_closed_form(E, P), in_4E_smul(E, P, odd)
            agree += (a == b); tot += 1; in4 += b
    return dict(n=n, N=E.N, samples=tot, agree=agree, fraction_in_4E=in4 / tot)

# ---------------------------------------------------------------- GF(2) linear algebra
def kernel_of_linearised(E, coeffs):
    """Kernel of X -> sum_k coeffs[k] X^(2^k) on GF(2^n), as a list of basis ints."""
    n = E.n
    cols = []
    for i in range(n):
        x = 1 << i
        acc, cur = 0, x
        for k, c in enumerate(coeffs):
            if c:
                acc ^= cur
            cur = gmul(cur, cur, E.f, n)
        cols.append(acc)
    # rows of the matrix M (n x n) with M[:, i] = cols[i]; solve M v = 0
    # Gaussian elimination on the augmented identity
    rows = [(cols[i], 1 << i) for i in range(n)]   # (image, preimage-combination)
    basis = []
    pivots = {}
    for img, comb in rows:
        while img:
            p = img.bit_length() - 1
            if p in pivots:
                pimg, pcomb = pivots[p]
                img ^= pimg; comb ^= pcomb
            else:
                pivots[p] = (img, comb)
                break
        if img == 0:
            basis.append(comb)
    return basis

def irreducible_quintics():
    out = []
    for mid in range(16):
        f = (1 << 5) | (mid << 1) | 1
        # no roots in GF(2): f(0)=1 already, f(1)= parity of bits
        if bin(f).count("1") % 2 == 0:
            continue
        # not divisible by x^2 + x + 1 (0b111)
        a = f
        while a.bit_length() >= 3:
            a ^= 0b111 << (a.bit_length() - 3)
        if a == 0:
            continue
        out.append(f)
    return out

def span(basis):
    S = [0]
    for g in basis:
        S = S + [s ^ g for s in S]
    return S

def scan(n=31):
    E = K0(n, 0)
    fac = factor(E.N)
    quint = irreducible_quintics()
    assert len(quint) == 6, quint
    comps = []
    for f in quint:
        coeffs = [(f >> k) & 1 for k in range(6)]
        ker = kernel_of_linearised(E, coeffs)
        assert len(ker) == 5, (f, len(ker))
        comps.append((f, ker))
    results = []
    import itertools
    for size in (1, 2):
        for S in itertools.combinations(range(6), size):
            basis = [b for j in S for b in comps[j][1]]
            pts = tot4 = tot2 = 0
            for x in span(basis):
                if x == 0:
                    continue
                for P in E.lift(x):
                    pts += 1
                    tot2 += tr(E, P[0]) == 0
                    tot4 += in_4E_closed_form(E, P)
            results.append(dict(factors=[bin(comps[j][0]) for j in S], dim=5 * size,
                                points=pts, in_2E=tot2, in_4E=tot4,
                                fraction_in_4E=tot4 / max(pts, 1)))
    return dict(n=n, N=E.N, factorisation={str(k): v for k, v in fac.items()}, subspaces=results)

if __name__ == "__main__":
    v = validate()
    print(json.dumps(dict(validation=v), indent=1), flush=True)
    s = scan()
    print(json.dumps(s, indent=1))

# Measure first fall degree (operational), plain-XL solving degree, and the
# Huang–Kosters–Yeo last fall degree d_F of the Weil descent of Semaev S_3
# restricted to an F_2-subspace V of dim n' inside F_{2^n}.
# Works in the Boolean ring B = F_2[x]/(x_i^2 - x_i); by the reduction argument
# (pi(V_{F ∪ FE, c}) = Boolean closure) this equals HKY's d_F for F ∪ field eqs.
import sys, json, time
import random as pyrandom
from sage.all import *

def descend_S3(n, nprime, seed, family="random", sat=None):
    seed = int(seed); rnd = pyrandom.Random(seed)
    set_random_seed(seed)
    K = GF(2**n, 'z', modulus='minimal_weight')
    z = K.gen(); fpoly = K.modulus()  # over GF(2)
    fc = [int(c) for c in fpoly.list()]  # coefficients, length n+1
    B = BooleanPolynomialRing(2*nprime, 'x', order='deglex')
    X = B.gens()
    # subspace basis
    if family == "subfield":
        assert n % nprime == 0
        sub = K.subfield(nprime) if hasattr(K,'subfield') else None
        # elements of subfield: fixed by Frobenius^nprime
        elems = [e for e in K if e**(2**nprime) == e]
        Vb = []
        Vspan = set([K(0)])
        for e in elems:
            if e not in Vspan:
                Vb.append(e); Vspan = set(a+b for a in Vspan for b in [K(0),e])
        assert len(Vb) == nprime
    else:
        while True:
            Vb = [K.random_element() for _ in range(nprime)]
            M = matrix(GF(2), [v._vector_() for v in Vb])
            if M.rank() == nprime: break
    # coordinates of element with boolean-poly coefficients: list length n over B
    def coords(e):
        return [B(int(c)) for c in e._vector_()]
    def add(a,b): return [a[i]+b[i] for i in range(n)]
    def mul(a,b):
        prod = [B(0)]*(2*n-1)
        for i in range(n):
            if a[i]==0: continue
            for j in range(n):
                if b[j]==0: continue
                prod[i+j] += a[i]*b[j]
        # reduce mod fpoly (monic, degree n): z^n = sum_{k<n} fc[k] z^k
        for d in range(2*n-2, n-1, -1):
            c = prod[d]
            if c == 0: continue
            prod[d] = B(0)
            for k in range(n):
                if fc[k]: prod[d-n+k] += c
        return prod[:n]
    def lin(vars_):  # sum_j vars_[j] * Vb[j]
        out = [B(0)]*n
        for j,v in enumerate(Vb):
            vc = v._vector_()
            for i in range(n):
                if vc[i]: out[i] += vars_[j]
        return out
    x1 = lin(X[:nprime]); x2 = lin(X[nprime:])
    b = K.random_element()
    while b == 0: b = K.random_element()
    # target: decomposable (sat) or random
    if sat is True:
        # choose random x1,x2 in V and compute x3 from S_3? Not direct; instead pick
        # points: use curve y^2+xy=x^3+ax^2+b with a=0; need x3 = x(P1+P2). Simpler: pick
        # x3 as root of S_3(v1,v2,X)=0 for random v1,v2 in V.
        R = PolynomialRing(K,'T'); T = R.gen()
        while True:
            v1 = sum(rnd.choice([K(0),K(1)])*v for v in Vb); v2 = sum(rnd.choice([K(0),K(1)])*v for v in Vb)
            S = (v1*v2 + v1*T + v2*T)**2 + v1*v2*T + b
            rts = S.roots()
            if rts: x3 = rnd.choice(rts)[0]; break
    else:
        x3 = K.random_element()
    x3c = coords(x3)
    e = add(add(mul(x1,x2), mul(x1,x3c)), mul(x2,x3c))
    S3 = add(add(mul(e,e), mul(mul(x1,x2),x3c)), coords(b))
    F = [f for f in S3]  # n boolean polys
    global LAST
    LAST = dict(K=K, Vb=Vb, x3=x3, b=b, B=B)
    return B, F, dict(n=n, nprime=nprime, seed=seed, family=family, x3=str(x3), b=str(b))

def rref_polys(polys, B):
    """Return RREF basis (as polynomials) of span of polys."""
    polys = [p for p in polys if p != 0]
    if not polys: return []
    # collect monomials, sorted descending in deglex so RREF leading terms are max
    monset = set()
    for p in polys:
        for m in p.monomials(): monset.add(m)
    mons = sorted(monset, reverse=True)
    idx = {m:i for i,m in enumerate(mons)}
    M = matrix(GF(2), len(polys), len(mons), sparse=False)
    for i,p in enumerate(polys):
        for m in p.monomials(): M[i, idx[m]] = 1
    E = M.echelon_form()
    out = []
    for r in E.rows():
        nz = r.nonzero_positions()
        if not nz: break
        out.append(sum(mons[j] for j in nz))
    return out

def macaulay_rows(F, B, c):
    X = B.gens()
    rows = []
    for f in F:
        d = f.degree()
        if d > c: continue
        # multilinear monomials of degree <= c-d
        for k in range(0, c-d+1):
            for comb in Combinations(range(len(X)), k):
                m = B(1)
                for i in comb: m *= X[i]
                rows.append(m*f)
    return rows

def closure(F, B, c, mutants=True):
    """V_{F,c}: mutant closure. Returns RREF basis."""
    X = B.gens()
    basis = rref_polys(macaulay_rows(F, B, c), B)
    if not mutants: return basis
    while True:
        new = list(basis)
        for bpoly in basis:
            if bpoly.degree() < c:
                for x in X:
                    new.append(x*bpoly)
        nb = rref_polys(new, B)
        if len(nb) == len(basis):
            return nb
        basis = nb

def ideal_dims(F, B, cmax):
    I = B.ideal(F)
    G = I.groebner_basis()
    DG = max([g.degree() for g in G]) if G else 0
    lms = [B(g.lm()) for g in G]
    N = B.ngens()
    # count standard monomials of degree <= c
    dims = {}
    unsat = any(g == 1 for g in G)
    for c in range(cmax+1):
        tot = sum(binomial(N,k) for k in range(c+1))
        if unsat:
            dims[c] = tot; continue
        std = 0
        for k in range(c+1):
            for comb in Combinations(range(N), k):
                m = B(1)
                for i in comb: m *= B.gens()[i]
                if not any(m.reducible_by(l) for l in lms):
                    std += 1
        dims[c] = tot - std
    return G, DG, dims, unsat

def measure(n, nprime, seed, family="random", sat=None, cmax=8):
    t0 = time.time()
    B, F, meta = descend_S3(n, nprime, seed, family, sat)
    G, DG, idim, unsat = ideal_dims(F, B, cmax)
    meta.update(unsat=unsat, DG=DG, nsol=None)
    res = {}
    prev_plain = 0; ffd = None; xl_solve = None; e0 = None
    plain_dims = {}; clos_dims = {}
    for c in range(2, cmax+1):
        plain = closure(F, B, c, mutants=False)
        low = sum(1 for p in plain if p.degree() < c)
        if ffd is None and low > prev_plain: ffd = c
        prev_plain = len(plain)
        plain_dims[c] = len(plain)
        if xl_solve is None and len(plain) == idim[c] and c >= DG: xl_solve = c
        clo = closure(F, B, c, mutants=True)
        clos_dims[c] = len(clo)
        if e0 is None and len(clo) == idim[c] and c >= DG: e0 = c
        if e0 is not None: break
    dlast = None
    if e0 is not None:
        dlast = e0
        for c in range(e0-1, 1, -1):
            if clos_dims.get(c) == idim[c]: dlast = c
            else: break
    meta.update(ffd=ffd, xl_solving=xl_solve, d_last=dlast, e0=e0,
                plain_dims=plain_dims, closure_dims=clos_dims, ideal_dims={c: idim[c] for c in plain_dims},
                secs=round(time.time()-t0,1))
    return meta

if __name__ == "__main__" and not globals().get("LFD_LIBRARY"):
    n = int(sys.argv[1]); nprime = int(sys.argv[2]); seeds = int(sys.argv[3])
    family = sys.argv[4] if len(sys.argv) > 4 else "random"
    sat = None if len(sys.argv) <= 5 else (sys.argv[5] == "sat")
    cmax = int(sys.argv[6]) if len(sys.argv) > 6 else 8
    for s in range(seeds):
        r = measure(n, nprime, 1000*n + s, family, sat, cmax)
        print(json.dumps(r, default=lambda o: int(o) if isinstance(o, Integer) else float(o)), flush=True)

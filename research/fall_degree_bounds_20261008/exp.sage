# Experiments E1 (trace mechanism) and E2 (random-quadratic controls) for
# research/fall_degree_bounds_20261008/PREREGISTRATION.md.
#   sage exp.sage trace   n n' draws family [rand|sat] cmax
#   sage exp.sage control n n' draws kind  [rand|sat] cmax   kind in {plain, planted_linear}
import sys, json
import random as pyrandom
LFD_LIBRARY = True
load("lfd_fast.sage")

def last_fall(B, F, cmax):
    G, DG, idim, unsat = ideal_dims(F, B, cmax)
    ffd = None; prev = 0; clos = {}; e0 = None
    for c in range(2, cmax+1):
        plain = closure(F, B, c, mutants=False)
        if ffd is None and sum(1 for p in plain if p.degree() < c) > prev: ffd = c
        prev = len(plain)
        clos[c] = len(closure(F, B, c, mutants=True))
        if c >= DG and clos[c] == idim[c]:
            e0 = c; break
    dl = None
    if e0 is not None:
        dl = e0
        for c in range(e0-1, 1, -1):
            if clos.get(c) == idim[c]: dl = c
            else: break
    return dict(unsat=bool(unsat), DG=int(DG), ffd=ffd, d_last=dl, cap=cmax,
                closure_dims={int(k): int(v) for k, v in clos.items()},
                ideal_dims={int(c): int(idim[c]) for c in clos})

def trace_cell(n, nprime, seed, family, sat, cmax):
    B, F, meta = descend_S3(n, nprime, seed, family, sat)
    K, Vb, x3 = LAST["K"], LAST["Vb"], LAST["x3"]
    X = B.gens()
    trV = [int(v.trace()) for v in Vb]
    trx3 = int(x3.trace())
    # predicted Kosters–Yeo relation restricted to V
    L = B(trx3) + sum(B(trV[j])*(X[j] + X[nprime+j]) for j in range(nprime))
    span = rref_polys(F, B)
    low = [p for p in span if p.degree() <= 1]
    in_span = len(rref_polys(span + [L], B)) == len(span) if L != 0 else True
    r = last_fall(B, F, cmax)
    r.update(meta)
    r.update(experiment="E1-trace", trace_on_V_zero=all(t == 0 for t in trV),
             tr_x3=trx3, predicted_L=str(L), L_in_span=bool(in_span),
             deg_le1_dim=len(low), deg_le1=[str(p) for p in low])
    return r

def syz_cell(n, nprime, seed, family):
    """Degree-3 linear syzygies sum_i l_i f_i = 0 (l_i affine-linear) of the descended system."""
    B, F, meta = descend_S3(n, nprime, seed, family, None)
    K, Vb, x3, b = LAST["K"], LAST["Vb"], LAST["x3"], LAST["b"]
    X = B.gens(); N = len(X)
    trV = [int(v.trace()) for v in Vb]
    c = int((b/x3**2).trace()) if x3 != 0 else None
    L = B(c) + sum(B(trV[j])*(X[j] + X[nprime+j]) for j in range(nprime))
    mults = [B(1)] + list(X)
    prods = [m*f for f in F for m in mults]          # index = i*(N+1) + k
    monset = set()
    for p in prods:
        for mm in p.monomials(): monset.add(mm)
    mons = sorted(monset, reverse=True); idx = {mm: j for j, mm in enumerate(mons)}
    M = matrix(GF(2), len(prods), len(mons))
    for r, p in enumerate(prods):
        for mm in p.monomials(): M[r, idx[mm]] = 1
    ker = M.left_kernel().basis()
    out = []
    for v in ker:
        ells = [sum(B(int(v[i*(N+1)+k]))*mults[k] for k in range(N+1)) for i in range(len(F))]
        nz = [e for e in ells if e != 0]
        common = nz[0] if nz and all(e == nz[0] for e in nz) else None
        S = [i for i, e in enumerate(ells) if e != 0]
        FS = sum(F[i] for i in S) if S else B(0)
        out.append(dict(common_multiplier=str(common) if common is not None else None,
                        F_S=str(FS), F_S_equals_L=bool(FS == L), multiplier_equals_L_plus_1=bool(common == L + 1)))
    r = dict(experiment="E1b-syzygy", n=n, nprime=nprime, seed=int(seed), family=family,
             trace_on_V_zero=all(t == 0 for t in trV), predicted_L=str(L), syz_dim=len(ker), syzygies=out)
    return r

def random_quad(B, rnd):
    X = B.gens(); N = len(X); p = B(rnd.getrandbits(1))
    for i in range(N):
        if rnd.getrandbits(1): p += X[i]
        for j in range(i+1, N):
            if rnd.getrandbits(1): p += X[i]*X[j]
    return p

def control_cell(n, nvars, seed, kind, sat, cmax):
    seed = int(seed); rnd = pyrandom.Random(seed)
    B = BooleanPolynomialRing(nvars, 'x', order='deglex'); X = B.gens()
    F = [random_quad(B, rnd) for _ in range(n)]
    if kind == "planted_linear":
        # mimic the trace relation: one F_2-combination of the equations is affine-linear
        ell = B(rnd.getrandbits(1)) + sum(X[i] for i in range(nvars) if rnd.getrandbits(1))
        F[-1] = sum(F[:-1]) + ell
    if sat:
        pt = [rnd.getrandbits(1) for _ in range(nvars)]
        F = [f + B(f(*pt)) for f in F]   # shift constants so pt is a root
        if kind == "planted_linear":
            ell = F[-1] + sum(F[:-1])    # keep the linear relation consistent
    r = last_fall(B, F, cmax)
    r.update(experiment="E2-control", kind=kind, n=n, nprime=nvars//2, seed=seed,
             sat_requested=bool(sat))
    return r

if __name__ == "__main__":
    mode = sys.argv[1]; n = int(sys.argv[2]); np_ = int(sys.argv[3]); draws = int(sys.argv[4])
    arg = sys.argv[5]; sat = sys.argv[6] == "sat"; cmax = int(sys.argv[7])
    for s in range(draws):
        if mode == "syz":
            r = syz_cell(n, np_, 5000*n + s, arg)
        elif mode == "trace":
            r = trace_cell(n, np_, 7000*n + s, arg, True if sat else None, cmax)
        else:
            r = control_cell(n, 2*np_, 9000*n + s, arg, sat, cmax)
        print(json.dumps(r, default=lambda o: int(o) if isinstance(o, Integer) else str(o)), flush=True)

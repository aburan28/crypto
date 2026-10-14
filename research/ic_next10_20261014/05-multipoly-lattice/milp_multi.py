"""Exact finite-ℓ optimum of the multi-polynomial full-rank triangular family over ALL monomial
sets in a bounding box and EVERY admissible leading-vertex choice (MILP, HiGHS via scipy).

Family: shifts p^{ℓ−|k|} x^γ Π f_j^{k_j}, leading monomial β = γ + Σ k_j a_j (a_j a vertex of
Newt(f_j); the a_j must be simultaneous unique maxima of one weight w, which makes every product
triangular), at most one shift per leading monomial, support (generic Minkowski sum) inside M.
Success at p → ∞:  Σ_β |k(β)| > Σ_{β∈M} λ·β,  λ = δ on small unknowns, 1 on large ones.

MILP per δ:  max Σ |k| z_{β,k} − Σ_β λ(δ)·β m_β
             z_{β,k} ≤ m_{β'} for β' ∈ supp(x^γ f^k);  Σ_k z_{β,k} ≤ 1;  m, z binary.
Dinkelbach on δ; the returned δ is the exact rational ratio of the optimal M.

Usage: python3 -I milp_multi.py out.json [case ...]
"""
import itertools
import json
import sys
import time
from fractions import Fraction

import numpy as np
from scipy.optimize import LinearConstraint, linprog, milp, Bounds
from scipy.sparse import coo_matrix


def box_support(nv, var_deg):
    """dense box support: var_deg = {v: d}"""
    vs = sorted(var_deg)
    out = []
    for e in itertools.product(*[range(var_deg[v] + 1) for v in vs]):
        x = [0] * nv
        for v, ev in zip(vs, e):
            x[v] = ev
        out.append(tuple(x))
    return out


def curve_support(nv, xv, yv):
    """y² − x³ − a x − b, a, b ≠ 0"""
    out = []
    for ex, ey in ((0, 2), (3, 0), (1, 0), (0, 0)):
        x = [0] * nv
        x[xv] = ex
        x[yv] = ey
        out.append(tuple(x))
    return out


def vertices(supp):
    """vertices of the Newton polytope = support points that are unique maxima of some weight"""
    vs = []
    pts = np.array(supp, dtype=float)
    for i, a in enumerate(supp):
        others = [s for s in supp if s != a]
        if not others:
            vs.append(a)
            continue
        A = np.array([np.array(s) - np.array(a) for s in others], dtype=float)  # w·(s−a) ≤ −1
        res = linprog(np.zeros(len(a)), A_ub=A, b_ub=-np.ones(len(others)), bounds=[(None, None)] * len(a))
        if res.status == 0:
            vs.append(a)
    return vs


def common_weight(supps, leads):
    rows = []
    for supp, a in zip(supps, leads):
        for s in supp:
            if s != a:
                rows.append(np.array(s) - np.array(a))
    A = np.array(rows, dtype=float)
    res = linprog(np.zeros(A.shape[1]), A_ub=A, b_ub=-np.ones(len(rows)), bounds=[(None, None)] * A.shape[1])
    return res.status == 0


def minkowski(a, b):
    return {tuple(x + y for x, y in zip(p, q)) for p in a for q in b}


def solve_case(name, nv, large, supps, leads, ell, ext, verbose=False):
    r = len(supps)
    mons = list(itertools.product(*[range(e + 1) for e in ext]))
    idx = {m: i for i, m in enumerate(mons)}
    nm = len(mons)
    small_sum = np.array([sum(m[v] for v in range(nv) if not large[v]) for m in mons], dtype=float)
    large_sum = np.array([sum(m[v] for v in range(nv) if large[v]) for m in mons], dtype=float)
    # powers' supports
    kvecs = [k for k in itertools.product(range(ell + 1), repeat=r) if 1 <= sum(k) <= ell]
    zero = tuple([0] * nv)
    psupp = {}
    for k in kvecs:
        s = {zero}
        for j, kj in enumerate(k):
            for _ in range(kj):
                s = minkowski(s, supps[j])
        psupp[k] = s
    zs = []  # (beta index, |k|, list of support indices)
    for k in kvecs:
        lead = [sum(k[j] * leads[j][v] for j in range(r)) for v in range(nv)]
        sup = [tuple(np.array(s) - np.array(lead)) for s in psupp[k]]  # relative to β
        for bi, b in enumerate(mons):
            ok = True
            cols = []
            gamma = [b[v] - lead[v] for v in range(nv)]
            if min(gamma) < 0:
                continue
            for s in sup:
                t = tuple(b[v] + s[v] for v in range(nv))
                j = idx.get(t)
                if j is None:
                    ok = False
                    break
                cols.append(j)
            if ok:
                zs.append((bi, sum(k), cols, k))
    nz = len(zs)
    nvar = nm + nz
    # constraints: z - m_{β'} ≤ 0 ; Σ_k z_{β,k} ≤ 1
    rows, cols_, vals = [], [], []
    rc = 0
    for zi, (bi, ks, cols, _) in enumerate(zs):
        for j in cols:
            rows += [rc, rc]
            cols_ += [nm + zi, j]
            vals += [1.0, -1.0]
            rc += 1
    nimp = rc
    per_beta = {}
    for zi, (bi, _, _, _) in enumerate(zs):
        per_beta.setdefault(bi, []).append(zi)
    for bi, lst in per_beta.items():
        if len(lst) > 1:
            for zi in lst:
                rows.append(rc)
                cols_.append(nm + zi)
                vals.append(1.0)
            rc += 1
    A = coo_matrix((vals, (rows, cols_)), shape=(rc, nvar)).tocsr()
    ub = np.concatenate([np.zeros(nimp), np.ones(rc - nimp)])
    cons = LinearConstraint(A, -np.inf, ub)
    gain = np.array([z[1] for z in zs], dtype=float)
    delta = Fraction(0)
    best = None
    it = 0
    timed_out = False
    while True:
        it += 1
        c = np.concatenate([large_sum + float(delta) * small_sum, -gain])  # minimise −objective
        res = milp(c, constraints=cons, integrality=np.ones(nvar), bounds=Bounds(0, 1),
                   options={"disp": False, "time_limit": 400, "mip_rel_gap": 0})
        if res.status != 0:
            timed_out = True
            if res.x is None:
                break
        else:
            timed_out = False
        x = np.round(res.x).astype(int)
        msel = x[:nm]
        zsel = x[nm:]
        G = int(sum(zs[i][1] for i in range(nz) if zsel[i]))
        Lg = int(round(float(np.dot(msel, large_sum))))
        Sm = int(round(float(np.dot(msel, small_sum))))
        obj = Fraction(G - Lg) - delta * Sm
        if obj <= 0 or Sm == 0:
            break
        nd = Fraction(G - Lg, Sm)
        if best is not None and nd <= best[0]:
            break
        best = (nd, int(msel.sum()), G, [mons[i] for i in range(nm) if msel[i]])
        delta = nd
        if it > 30:
            break
    if best is None:
        return {"system": name, "ell": ell, "leads": [list(a) for a in leads], "ext": list(ext),
                "delta": 0.0, "delta_exact": "0", "dim": 0, "sum_k": 0, "binds": False, "nz": nz, "nm": nm}
    nd, dim, G, M = best
    maxes = [max(m[v] for m in M) for v in range(nv)]
    binds = any(maxes[v] == ext[v] for v in range(nv))
    return {"system": name, "ell": ell, "leads": [list(a) for a in leads], "ext": list(ext),
            "delta": float(nd), "delta_exact": f"{nd.numerator}/{nd.denominator}", "dim": dim, "sum_k": G,
            "maxes": maxes, "binds": binds, "milp_time_limit_hit": timed_out, "nz": nz, "nm": nm, "iters": it}


def run_system(name, nv, large, supps, ell, ext, lead_filter=None, out=None):
    vlists = [vertices(s) for s in supps]
    combos = []
    for leads in itertools.product(*vlists):
        if common_weight(supps, leads):
            if lead_filter is None or lead_filter(leads):
                combos.append(leads)
    best = None
    for leads in combos:
        t0 = time.time()
        r = solve_case(name, nv, large, supps, leads, ell, ext)
        r["seconds"] = round(time.time() - t0, 2)
        print(f"| {name} | {ell} | {r['leads']} | {r['ext']} | {r['delta_exact']} | {r['delta']:.6f} | {r['dim']} | "
              f"{r['sum_k']} | {r.get('maxes')} | {'BINDS' if r['binds'] else 'no'} | {r['nm']}+{r['nz']} | {r['seconds']} |",
              flush=True)
        out.append(r)
        json.dump(out, open(sys.argv[1], "w"), indent=1)
        if best is None or r["delta"] > best["delta"]:
            best = r
    print(f"| **{name} best over {len(combos)} leading-vertex choices** | {ell} | {best['leads']} | | "
          f"**{best['delta_exact']}** | **{best['delta']:.6f}** | | | | | | |", flush=True)
    out.append({"system": name, "ell": ell, "best_over_leads": best["delta"], "best_exact": best["delta_exact"],
                "n_lead_choices": len(combos)})
    return best


def main():
    out = []
    which = set(sys.argv[2:])
    print("| system | ℓ | leading vertices | bounding box | δ_ℓ exact | δ_ℓ | dim | Σk | max exps | box binds? | vars m+z | s |")
    print("|:--|--:|:--|:--|:--|--:|--:|--:|:--|:--|:--|--:|")
    # C1: S3 alone (x1, x2)
    if not which or "c1" in which:
        s3 = box_support(2, {0: 2, 1: 2})
        run_system("C1 S3 (all M)", 2, [False, False], [s3], 1, [6, 6], out=out)
    if "c1l2" in which:
        s3 = box_support(2, {0: 2, 1: 2})
        run_system("C1 S3 (all M)", 2, [False, False], [s3], 2, [9, 9], out=out)
    # V3: S3 pair sharing x1 (x1, x2, x3), all M — validates the down-closed restriction
    if not which or "v3" in which:
        f1 = box_support(3, {0: 2, 1: 2})
        f2 = box_support(3, {0: 2, 2: 2})
        run_system("V3 S3 pair (all M)", 3, [False] * 3, [f1, f2], 1, [6, 6, 6], out=out)
    # V1: S3 + curve equations (x1, x2, y1, y2), y large
    if not which or "v1" in which:
        f = box_support(4, {0: 2, 1: 2})
        g1 = curve_support(4, 0, 2)
        g2 = curve_support(4, 1, 3)
        run_system("V1 S3 + 1 curve eq", 4, [False, False, True, False], [f, curve_support(4, 0, 2)], 1, [7, 6, 3, 0], out=out)
        run_system("V1 S3 + 2 curve eqs", 4, [False, False, True, True], [f, g1, g2], 1, [7, 7, 3, 3], out=out)
    if "v1l2" in which:
        f = box_support(4, {0: 2, 1: 2})
        run_system("V1 S3 + 1 curve eq", 4, [False, False, True, False], [f, curve_support(4, 0, 2)], 2, [10, 9, 4, 0], out=out)
    # V2: S4 chain S3(x1,x2,t), S3(x3,t; z) with t large  (x1, x2, x3, t)
    if not which or "v2" in which:
        f1 = box_support(4, {0: 2, 1: 2, 3: 2})
        f2 = box_support(4, {2: 2, 3: 2})
        run_system("V2 S4 chain (t large)", 4, [False, False, False, True], [f1, f2], 1, [6, 6, 6, 4], out=out)
    json.dump(out, open(sys.argv[1], "w"), indent=1)


if __name__ == "__main__":
    main()

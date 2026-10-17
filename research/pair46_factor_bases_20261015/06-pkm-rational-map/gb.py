"""Minimal Buchberger over F_p (grevlex, sugar normal strategy, Gebauer-Moeller update).

Instrument for item 06: records D_max (max sugar degree of a processed pair) and D_corr
(max sugar degree of a pair whose reduction was nonzero), the analogue of PKM's D_av and
D_av^corr. Polynomials are dicts {exponent tuple: coefficient mod p}.
"""
from __future__ import annotations

import time


def gkey(e):
    return (sum(e), tuple(-x for x in reversed(e)))


def lm(f):
    return max(f, key=gkey)


def divides(a, b):
    return all(x <= y for x, y in zip(a, b))


def lcm(a, b):
    return tuple(max(x, y) for x, y in zip(a, b))


def coprime(a, b):
    return all(x == 0 or y == 0 for x, y in zip(a, b))


def monic(f, p):
    m = lm(f)
    inv = pow(f[m], -1, p)
    return {k: v * inv % p for k, v in f.items()}


def sub_mul(h, g, c, t, p):
    """h -= c * x^t * g (in place)."""
    for e, v in g.items():
        k = tuple(x + y for x, y in zip(e, t))
        nv = (h.get(k, 0) - c * v) % p
        if nv:
            h[k] = nv
        else:
            h.pop(k, None)


def reduce_full(f, G, p):
    """Full reduction of f by the list of monic polys G (with cached leading monomials)."""
    h = dict(f)
    r = {}
    while h:
        m = max(h, key=gkey)
        c = h[m]
        for g, gm in G:
            if divides(gm, m):
                t = tuple(x - y for x, y in zip(m, gm))
                sub_mul(h, g, c, t, p)
                break
        else:
            r[m] = c
            del h[m]
    return r


def spoly(f, fm, g, gmn, p):
    L = lcm(fm, gmn)
    s = {}
    sub_mul(s, f, p - 1, tuple(x - y for x, y in zip(L, fm)), p)  # +f*t
    sub_mul(s, g, 1, tuple(x - y for x, y in zip(L, gmn)), p)      # -g*t
    return s


def groebner(F, p, timeout=60.0):
    t0 = time.process_time()
    polys = []   # (poly, lm, sugar)
    active = []  # indices into polys
    B = []       # (sugar, key(lcm), i, j, lcm)
    stats = {"D_max": 0, "D_corr": 0, "Dl_max": 0, "Dl_corr": 0, "nonzero": 0, "pairs": 0, "timeout": False}

    def update(hidx):
        nonlocal B, active
        h, hm, hs = polys[hidx]
        C = [(g, lcm(polys[g][1], hm)) for g in active]
        D = []
        while C:
            g, L = C.pop(0)
            if coprime(polys[g][1], hm) or (not any(divides(L2, L) for _, L2 in C)
                                             and not any(divides(L2, L) for _, L2 in D)):
                D.append((g, L))
        E = [(g, L) for g, L in D if not coprime(polys[g][1], hm)]
        Bn = []
        for item in B:
            s, k, i, j, L = item
            if divides(hm, L) and lcm(polys[i][1], hm) != L and lcm(polys[j][1], hm) != L:
                continue
            Bn.append(item)
        for g, L in E:
            gs = polys[g][2] - sum(polys[g][1])
            sug = sum(L) + max(gs, hs - sum(hm))
            Bn.append((sug, gkey(L), g, hidx, L))
        B = Bn
        active = [g for g in active if not divides(hm, polys[g][1])] + [hidx]

    def basis():
        return [(polys[g][0], polys[g][1]) for g in active]

    for f in F:
        r = reduce_full(f, basis(), p)
        if r:
            r = monic(r, p)
            polys.append((r, lm(r), max(sum(e) for e in f)))
            update(len(polys) - 1)
    while B:
        if time.process_time() - t0 > timeout:
            stats["timeout"] = True
            break
        B.sort(key=lambda it: (it[0], it[1]))
        s, k, i, j, L = B.pop(0)
        stats["pairs"] += 1
        stats["D_max"] = max(stats["D_max"], s)
        stats["Dl_max"] = max(stats["Dl_max"], sum(L))
        S = spoly(polys[i][0], polys[i][1], polys[j][0], polys[j][1], p)
        r = reduce_full(S, basis(), p)
        if r:
            stats["D_corr"] = max(stats["D_corr"], s)
            stats["Dl_corr"] = max(stats["Dl_corr"], sum(L))
            stats["nonzero"] += 1
            r = monic(r, p)
            polys.append((r, lm(r), s))
            update(len(polys) - 1)
            if len(r) == 1 and sum(lm(r)) == 0:
                B = []
                active = [len(polys) - 1]
                break
    stats["time"] = time.process_time() - t0
    G = [polys[g][0] for g in active]
    return G, stats


def reduced(G, p):
    """Minimal, interreduced, monic basis."""
    G = sorted((monic(g, p) for g in G), key=lambda g: gkey(lm(g)))
    keep = []
    for g in G:
        if not any(divides(lm(h), lm(g)) for h in keep):
            keep.append(g)
    out = []
    for g in keep:
        others = [(h, lm(h)) for h in keep if h is not g]
        m = lm(g)
        tail = {k: v for k, v in g.items() if k != m}
        r = reduce_full(tail, others, p)
        r[m] = 1
        out.append(r)
    return out


def evaluate(f, pt, p):
    s = 0
    for e, c in f.items():
        t = c
        for x, k in zip(pt, e):
            if k:
                t = t * pow(x, k, p) % p
        s += t
    return s % p


def std_monomial_count(G, nvars, cap=200000):
    """Dimension of F_p[x]/I for zero-dim I (None if not zero-dim or too large)."""
    lms = [lm(g) for g in G]
    if any(sum(m) == 0 for m in lms):
        return 0
    bounds = []
    for i in range(nvars):
        pure = [m[i] for m in lms if sum(m) == m[i] and m[i] > 0]
        if not pure:
            return None
        bounds.append(min(pure))
    total = 1
    for b in bounds:
        total *= b
    if total > cap:
        return None
    count = 0

    def rec(i, e):
        nonlocal count
        if i == nvars:
            if not any(divides(m, e) for m in lms):
                count += 1
            return
        for k in range(bounds[i]):
            rec(i + 1, e + (k,))
    rec(0, ())
    return count

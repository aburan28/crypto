"""
Minimal multivariate Groebner-basis engine over F_p (Buchberger, grevlex),
used as a pure-Python stand-in for an F4/Macaulay "degree of regularity" probe
on small Semaev relation systems.

Polynomials are dicts {exponent_tuple: coeff mod p}.  Intended for tiny systems
(2-3 variables, low degree); it instruments the computation to report the
*solving degree* (max total degree of any reduced/ S-poly remainder reached),
the number of reductions, and the final basis size -- the practical proxies for
F4 behaviour requested in Phase 7.3.
"""

from __future__ import annotations
from itertools import combinations


def _deg(mono):
    return sum(mono)


def grevlex_key(mono):
    # total degree, then reverse of exponent vector negated (grevlex)
    return (_deg(mono),) + tuple(-x for x in reversed(mono))


def lead(poly):
    """Leading monomial (grevlex) and its coefficient."""
    lm = max(poly, key=grevlex_key)
    return lm, poly[lm]


def p_add(A, B, p):
    R = dict(A)
    for m, c in B.items():
        R[m] = (R.get(m, 0) + c) % p
        if R[m] == 0:
            del R[m]
    return R


def p_scale_mono(A, mono, coeff, p):
    return {tuple(e + mono[i] for i, e in enumerate(m)): (c * coeff) % p
            for m, c in A.items()}


def _mono_div(a, b):
    """a / b if divisible (exponentwise), else None."""
    q = []
    for ai, bi in zip(a, b):
        if ai < bi:
            return None
        q.append(ai - bi)
    return tuple(q)


def reduce_poly(f, G, p, stats):
    """Multivariate division of f by the list G; returns the remainder."""
    r = {}
    f = dict(f)
    while f:
        lm, lc = lead(f)
        stats["max_deg"] = max(stats["max_deg"], _deg(lm))
        divided = False
        for g in G:
            glm, glc = g["lead"]
            q = _mono_div(lm, glm)
            if q is not None:
                factor = (lc * pow(glc, -1, p)) % p
                f = p_add(f, p_scale_mono(g["poly"], q, (-factor) % p, p), p)
                stats["reductions"] += 1
                divided = True
                break
        if not divided:
            r[lm] = lc
            del f[lm]
    return r


def spoly(f, g, p):
    flm, flc = f["lead"]
    glm, glc = g["lead"]
    lcm = tuple(max(a, b) for a, b in zip(flm, glm))
    cf = _mono_div(lcm, flm)
    cg = _mono_div(lcm, glm)
    A = p_scale_mono(f["poly"], cf, pow(flc, -1, p), p)
    B = p_scale_mono(g["poly"], cg, pow(glc, -1, p), p)
    return p_add(A, B if False else {m: (-c) % p for m, c in B.items()}, p)


def _wrap(poly):
    return {"poly": poly, "lead": lead(poly)}


def buchberger(F, p, max_basis=400):
    """
    Compute a Groebner basis (grevlex). Returns (basis_polys, stats) where stats
    has solving_degree (max remainder/S-poly degree), reductions, basis_size.
    """
    G = [_wrap(f) for f in F if f]
    stats = {"max_deg": 0, "reductions": 0}
    pairs = list(combinations(range(len(G)), 2))
    while pairs:
        i, j = pairs.pop()
        s = spoly(G[i], G[j], p)
        if not s:
            continue
        r = reduce_poly(s, G, p, stats)
        if r:
            G.append(_wrap(r))
            ni = len(G) - 1
            pairs.extend((k, ni) for k in range(ni))
            if len(G) > max_basis:
                stats["aborted"] = True
                break
    stats["basis_size"] = len(G)
    stats["solving_degree"] = stats["max_deg"]
    return [g["poly"] for g in G], stats


def substitute_last(poly, value, p, nvars):
    """Specialise the last variable of an nvars-poly to a constant; return an
    (nvars-1)-variable poly dict."""
    R = {}
    for m, c in poly.items():
        val = (c * pow(value, m[-1], p)) % p
        key = m[:-1]
        R[key] = (R.get(key, 0) + val) % p
        if R[key] == 0:
            del R[key]
    return R

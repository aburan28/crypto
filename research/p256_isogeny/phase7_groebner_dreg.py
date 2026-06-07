#!/usr/bin/env python3
"""
Phase 7.3 - Groebner "degree of regularity" probe on Semaev relation systems.

For a target point with x-coordinate x_R, the 2-point factor-base decomposition
system is

    S_3(x1, x2, x_R) = 0          (collinearity / summation polynomial)
    g(x1) = prod_i (x1 - f_i) = 0  (x1 in the factor base {f_i})
    g(x2) = 0                      (x2 in the factor base)

a 0-dimensional system whose solutions are the FB pairs decomposing the target.
We run a pure-Python Buchberger (grevlex) and record the *solving degree* (max
remainder/S-poly degree reached), number of reductions and basis size -- the
practical F4 proxies -- then compare:

    root (P-256 mod q)  vs  same-order isogenous neighbour
    vs same-order random control  vs random curve
    vs GENERIC polynomial of identical support (H3 control)

H3 is confirmed-absent if the real Semaev system does not solve at a lower degree
(or with materially less work) than the generic same-shape control.
"""

from __future__ import annotations
import json
import random
import sys

from ecfp import Curve, INF, sqrt_mod, isogenous_neighbors
import semaev
import groebner as gb


def brute_order(E):
    p, a, b = E.p, E.a, E.b
    c = 1
    for x in range(p):
        r = (x * x * x + a * x + b) % p
        if r == 0:
            c += 1
        elif pow(r, (p - 1) // 2, p) == 1:
            c += 2
    return c


def fb_points(E, K):
    pts, x = [], 0
    while len(pts) < K and x < E.p:
        r = (x * x * x + E.a * x + E.b) % E.p
        y = sqrt_mod(r, E.p)
        if y is not None:
            pts.append((x, y))
        x += 1
    return pts


def build_system(E, fb, x_R, p):
    """[S3(x1,x2,x_R), g(x1), g(x2)] as 2-variable polys."""
    S3 = semaev.semaev3(E.a, E.b, p)                 # vars (x1,x2,x3)
    S3sub = gb.substitute_last(S3, x_R, p, 3)        # -> (x1,x2)
    gx1 = {(0, 0): 1}
    gx2 = {(0, 0): 1}
    for (fx, _) in fb:
        gx1 = _polymul2(gx1, {(1, 0): 1, (0, 0): (-fx) % p}, p)
        gx2 = _polymul2(gx2, {(0, 1): 1, (0, 0): (-fx) % p}, p)
    return [S3sub, gx1, gx2]


def _polymul2(A, B, p):
    R = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            m = (ma[0] + mb[0], ma[1] + mb[1])
            R[m] = (R.get(m, 0) + ca * cb) % p
            if R[m] == 0:
                del R[m]
    return R


def generic_like(poly, p, rng):
    """Random polynomial with the same monomial support as `poly`."""
    return {m: rng.randrange(1, p) for m in poly}


def solve_metrics(system, p):
    _, stats = gb.buchberger(system, p)
    return stats["solving_degree"], stats["reductions"], stats["basis_size"]


def decomposable_targets(E, fb, want, rng):
    """x-coords of sums of two FB points (guaranteed-consistent targets)."""
    pts = [(x, y) for (x, y) in fb]
    pts += [(x, (-y) % E.p) for (x, y) in fb]
    out = []
    while len(out) < want:
        U = rng.choice(pts)
        V = rng.choice(pts)
        S = E.add(U, V)
        if S is not INF:
            out.append(S[0])
    return out


def audit_curve(E, K, n_targets, p, rng, label):
    fb = fb_points(E, K)
    if len(fb) < K:
        return None
    targets = decomposable_targets(E, fb, n_targets, rng)
    sds, reds, gen_sds = [], [], []
    for x_R in targets:
        system = build_system(E, fb, x_R, p)
        sd, rd, _ = solve_metrics(system, p)
        sds.append(sd)
        reds.append(rd)
        # generic same-support control replaces the Semaev poly
        gsystem = [generic_like(system[0], p, rng)] + system[1:]
        gsd, _, _ = solve_metrics(gsystem, p)
        gen_sds.append(gsd)
    return {
        "label": label, "a": E.a, "b": E.b, "fb_size": K, "targets": n_targets,
        "semaev_solving_degree_mean": round(sum(sds) / len(sds), 3),
        "semaev_solving_degree_max": max(sds),
        "reductions_mean": round(sum(reds) / len(reds), 1),
        "generic_solving_degree_mean": round(sum(gen_sds) / len(gen_sds), 3),
    }


def main():
    import params as PP
    q = int(sys.argv[1]) if len(sys.argv) > 1 else 12007
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    n_targets = 10
    n_ctrl = 3
    rng = random.Random(7)

    a0, b0 = (-3) % q, PP.B % q
    E0 = Curve(a0, b0, q)
    n_q = brute_order(E0)
    t_q = q + 1 - n_q

    # same-order neighbour (smallest working odd ell)
    neighbor = None
    for ell in [3, 5, 7, 11, 13, 17, 19, 23]:
        if ell * ell > 4 * q:
            continue
        nb = isogenous_neighbors(E0, ell, t_q)
        if nb:
            neighbor = nb[0][2]
            break

    # controls
    same_order, random_curves = [], []
    tries = 0
    while len(same_order) < n_ctrl and tries < 6000:
        tries += 1
        a, b = rng.randrange(q), rng.randrange(q)
        if (4 * a ** 3 + 27 * b ** 2) % q == 0:
            continue
        E = Curve(a, b, q)
        if brute_order(E) == n_q:
            same_order.append(E)
    while len(random_curves) < n_ctrl:
        a, b = rng.randrange(q), rng.randrange(q)
        if (4 * a ** 3 + 27 * b ** 2) % q == 0:
            continue
        random_curves.append(Curve(a, b, q))

    rows = [audit_curve(E0, K, n_targets, q, rng, "root")]
    if neighbor:
        rows.append(audit_curve(neighbor, K, n_targets, q, rng, "neighbor"))
    for i, E in enumerate(same_order):
        rows.append(audit_curve(E, K, n_targets, q, rng, f"same_order_ctrl_{i}"))
    for i, E in enumerate(random_curves):
        rows.append(audit_curve(E, K, n_targets, q, rng, f"random_ctrl_{i}"))
    rows = [r for r in rows if r]

    sd_root = rows[0]["semaev_solving_degree_mean"]
    out = {
        "experiment": "Groebner degree-of-regularity on Semaev 2-decomposition systems",
        "q": q, "n_q": n_q, "fb_size_K": K, "targets_per_curve": n_targets,
        "rows": rows,
        "verdict": {
            "root_semaev_solving_degree": sd_root,
            "neighbor_semaev_solving_degree":
                next((r["semaev_solving_degree_mean"] for r in rows
                      if r["label"] == "neighbor"), None),
            "all_semaev_solving_degrees":
                {r["label"]: r["semaev_solving_degree_mean"] for r in rows},
            "semaev_vs_generic_same_shape":
                {r["label"]: [r["semaev_solving_degree_mean"],
                              r["generic_solving_degree_mean"]] for r in rows},
            "any_neighbor_or_ctrl_below_root":
                any(r["semaev_solving_degree_mean"] < sd_root for r in rows[1:]),
        },
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()

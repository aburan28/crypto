#!/usr/bin/env python3
"""Item 06: PKM (Petit-Kosters-Messeng, PKC 2016) rational-map factor bases at toy sizes.

Arms A1-A8 and metrics M1-M3 as in PROTOCOL.md. Single core, pure Python.
Usage: python3 run_pkm.py [max_n1] > run.log ; writes results.json
"""
from __future__ import annotations

import hashlib
import json
import math
import random
import sys
import time

sys.path.insert(0, "/home/user/crypto-autoresearcher/src")
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))

from crypto_autoresearcher.index_calculus.curve import Curve, is_probable_prime, factorize  # noqa: E402
from crypto_autoresearcher.index_calculus.semaev import s3  # noqa: E402
from crypto_autoresearcher.index_calculus.msolve import s3_terms  # noqa: E402
import gb  # noqa: E402

T_START = time.time()
CAP = 35 * 60
TIMEOUT = 60.0
NT = 12
NRHO = 12


def rng_for(*labels):
    h = hashlib.sha256("|".join(map(str, ("pkm06",) + labels)).encode()).digest()
    return random.Random(int.from_bytes(h, "big"))


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def count_points(p, a, b):
    return p + 1 + sum(legendre(x * x * x + a * x + b, p) for x in range(p))


def prime_order_curve(p, label):
    r = rng_for("curve", label, p)
    while True:
        a, b = r.randrange(p), r.randrange(p)
        if (4 * a ** 3 + 27 * b * b) % p == 0:
            continue
        N = count_points(p, a, b)
        if is_probable_prime(N):
            return Curve(p, a, b, N)


def prime_mult(n1):
    d = 2 ** n1
    p = 2 ** (2 * n1)
    while True:
        if p % d == 1 and is_probable_prime(p):
            return p
        p += 1


def prime_dickson(n1):
    d = 2 ** n1
    p = 2 ** (2 * n1)
    while True:
        if (p + 1) % d == 0 and (p + 1) // d >= 3 and is_probable_prime(p):
            return p
        p += 1


def dickson2_iter(x, n, p):
    for _ in range(n):
        x = (x * x - 2) % p
    return x


# ---------------- factor-base sets and their systems -------------------------
# A system is (nvars, list of polys); variables: block i (i=0,1) of size k, var 0 of block = x_i.

def s3_poly2(E, nv, i1, i2, X):
    p = E.p
    out = {}
    for (e1, e2, e3), c in s3_terms(E.a, E.b):
        ex = [0] * nv
        ex[i1] += e1
        ex[i2] += e2
        k = tuple(ex)
        out[k] = (out.get(k, 0) + c * pow(X, e3, p)) % p
    return {k: v for k, v in out.items() if v}


def var(nv, i, k=1):
    e = [0] * nv
    e[i] = k
    return tuple(e)


def poly_from_univ(coeffs, nv, i, p):
    out = {}
    for k, c in coeffs.items():
        if c % p:
            out[var(nv, i, k)] = c % p
    return out


def poly_mul_lin(coeffs, r, p):
    """coeffs (dict deg->c) times (x - r)."""
    out = {}
    for k, c in coeffs.items():
        out[k + 1] = (out.get(k + 1, 0) + c) % p
        out[k] = (out.get(k, 0) - r * c) % p
    return out


def flat_system(E, V, X, univ):
    p = E.p
    nv = 2
    return nv, [s3_poly2(E, nv, 0, 1, X), poly_from_univ(univ, nv, 0, p), poly_from_univ(univ, nv, 1, p)]


def prod_univ(V, p):
    c = {0: 1}
    for v in V:
        c = poly_mul_lin(c, v, p)
    return c


def chain_system(E, X, k, step_eqs, final_eq):
    """step_eqs(nv, a, b) -> poly for relation x_b = L_j(x_a) (j = 0..k-2);
    final_eq(nv, a) -> poly L_k(x_a) = 0."""
    nv = 2 * k
    F = [s3_poly2(E, nv, 0, k, X)]
    for blk in (0, k):
        for j in range(k - 1):
            F.append(step_eqs(nv, blk + j, blk + j + 1, j))
        F.append(final_eq(nv, blk + k - 1))
    return nv, F


def mult_chain(E, X, n1):
    p = E.p

    def step(nv, a, b, j):  # x_b - x_a^2
        return {var(nv, b): 1, var(nv, a, 2): p - 1}

    def final(nv, a):  # x_a^2 - 1
        return {var(nv, a, 2): 1, var(nv, a, 0): p - 1}
    return chain_system(E, X, n1, step, final)


def dickson_chain(E, X, n1, lam):
    p = E.p

    def step(nv, a, b, j):  # x_b - x_a^2 + 2
        return {var(nv, b): 1, var(nv, a, 2): p - 1, var(nv, a, 0): 2}

    def final(nv, a):  # x_a^2 - 2 - lam
        return {var(nv, a, 2): 1, var(nv, a, 0): (-2 - lam) % p}
    return chain_system(E, X, n1, step, final)


def isogeny_chain(E, X, n1, steps, chi):
    """steps[j] = (x0_j, v_j): L_j(x) = x + v/(x - x0)."""
    p = E.p

    def step(nv, a, b, j):  # x_b (x_a - x0) - (x_a^2 - x0 x_a + v)
        x0, v = steps[j]
        e_ab = [0] * nv
        e_ab[a] += 1
        e_ab[b] += 1
        f = {tuple(e_ab): 1, var(nv, b): (-x0) % p, var(nv, a, 2): p - 1,
             var(nv, a): x0 % p, var(nv, a, 0): (-v) % p}
        return {k: c for k, c in f.items() if c}

    def final(nv, a):  # x^2 - x0 x + v - chi (x - x0)
        x0, v = steps[n1 - 1]
        f = {var(nv, a, 2): 1, var(nv, a): (-x0 - chi) % p, var(nv, a, 0): (v + chi * x0) % p}
        return {k: c for k, c in f.items() if c}
    return chain_system(E, X, n1, step, final)


# ---------------- auxiliary curve for the isogeny instance -------------------

def find_aux_isogeny(p, n1):
    d = 2 ** n1
    r = rng_for("aux", p, n1)
    for _ in range(20000):
        a, b = r.randrange(p), r.randrange(p)
        if (4 * a ** 3 + 27 * b * b) % p == 0:
            continue
        N = count_points(p, a, b)
        if N % d:
            continue
        Ep = Curve(p, a, b, N)
        v2 = 0
        while N % (2 ** (v2 + 1)) == 0:
            v2 += 1
        odd = N >> v2
        for _t in range(30):
            Rr = Ep.random_point(r)
            h = Ep.mul(odd, Rr)
            e = 0
            Q = h
            while Q is not None:
                Q = Ep.add(Q, Q)
                e += 1
            if e >= n1:
                h = Ep.mul(2 ** (e - n1), h)
                res = isogeny_steps(Ep, h, n1, r)
                if res:
                    return Ep, h, res
    raise RuntimeError("no auxiliary curve")


def isogeny_steps(Ep, h, n1, r):
    p = Ep.p
    a, b = Ep.a, Ep.b
    R0 = Ep.random_point(r)
    pts_h, pt_r = h, R0
    steps = []
    for k in range(n1):
        C = Curve(p, a, b)
        T = C.mul(2 ** (n1 - k - 1), pts_h)
        if T is None or T[1] != 0:
            return None
        x0 = T[0]
        v = (3 * x0 * x0 + a) % p
        w = x0 * v % p

        def phi(P):
            if P is None or P[0] == x0:
                return None
            inv = pow(P[0] - x0, -1, p)
            return ((P[0] + v * inv) % p, P[1] * (1 - v * inv * inv) % p)
        pts_h, pt_r = phi(pts_h), phi(pt_r)
        a, b = (a - 5 * v) % p, (b - 7 * w) % p
        C2 = Curve(p, a, b)
        for P in (pts_h, pt_r):
            if P is not None and not C2.is_on_curve(P):
                raise AssertionError("Velu image off curve")
        steps.append((x0, v))
    if pts_h is not None:
        raise AssertionError("kernel not killed")
    if pt_r is None or pt_r[1] == 0:
        return None
    chi = pt_r[0]
    # V = x(R0 + <h>)
    V = set()
    P = R0
    for _ in range(2 ** n1):
        V.add(P[0])
        P = Ep.add(P, h)
    if len(V) != 2 ** n1:
        return None
    return steps, chi, sorted(V)


# ---------------- relation enumeration (M1) ---------------------------------

def relations(E, V, R):
    p = E.p
    X = R[0]
    Vl = sorted(V)
    lift = {x: E.lift_x(x) for x in Vl}
    nbsols = 0
    valid = 0
    for i, x1 in enumerate(Vl):
        for x2 in Vl[i:]:
            if s3(E, x1, x2, X) == 0:
                nbsols += 1
                P1, P2 = lift[x1], lift[x2]
                if P1 is None or P2 is None:
                    continue
                for s1 in (1, -1):
                    for s2 in (1, -1):
                        Q1 = P1 if s1 == 1 else E.neg(P1)
                        Q2 = P2 if s2 == 1 else E.neg(P2)
                        if E.add(Q1, Q2) == R:
                            valid += 1
    return nbsols, valid


def solutions(E, V, R):
    X = R[0]
    Vl = sorted(V)
    return [(x1, x2) for x1 in Vl for x2 in Vl if s3(E, x1, x2, X) == 0]


# ---------------- rho (M3) ---------------------------------------------------

def rho_solve(E, P, Q, r):
    N = E.order
    E.ops.group_ops = 0
    M = [(r.randrange(N), r.randrange(N)) for _ in range(16)]
    Ms = [E.add(E.mul(c, P), E.mul(d, Q)) for c, d in M]
    E.ops.group_ops = 0
    t0 = time.process_time()

    def f(S, c, d):
        j = (S[0] if S else 0) % 16
        return E.add(S, Ms[j]), (c + M[j][0]) % N, (d + M[j][1]) % N
    c0, d0 = r.randrange(N), r.randrange(N)
    S0 = E.add(E.mul(c0, P), E.mul(d0, Q))
    E.ops.group_ops = 0
    t0 = time.process_time()
    T, ct, dt = S0, c0, d0
    H, ch, dh = S0, c0, d0
    while True:
        T, ct, dt = f(T, ct, dt)
        H, ch, dh = f(*f(H, ch, dh))
        if T == H:
            break
    el = time.process_time() - t0
    ops = E.ops.group_ops
    if (dh - dt) % N == 0:
        return None, ops, el
    k = (ct - ch) * pow(dh - dt, -1, N) % N
    return k, ops, el


# ---------------- driver -----------------------------------------------------

def run_gb(nv, F, E, V, R, check_sympy):
    G, st = gb.groebner(F, E.p, timeout=TIMEOUT)
    sols = solutions(E, V, R)
    ok = True
    if not st["timeout"]:
        # every enumerated (x1,x2) solution extends to a full solution; check GB vanishes on it
        dim = gb.std_monomial_count(G, nv)
        st["dim"] = dim
        st["nsols_ordered"] = len(sols)
        if dim is not None and dim < len(sols):
            ok = False
        if check_sympy:
            st["sympy_match"] = sympy_check(F, G, nv, E.p)
            ok = ok and st["sympy_match"]
    st["ok"] = ok
    return st


def sympy_check(F, G, nv, p):
    import sympy
    xs = sympy.symbols(f"z0:{nv}")

    def to_expr(f):
        return sum(c * sympy.Mul(*[x ** k for x, k in zip(xs, e)]) for e, c in f.items())
    try:
        GS = sympy.groebner([to_expr(f) for f in F], *xs, modulus=p, order="grevlex")
    except Exception:
        return None
    mine = gb.reduced(G, p)
    ms = set()
    for g in mine:
        ms.add(tuple(sorted((e, c % p) for e, c in g.items())))
    ss = set()
    for g in GS.exprs:
        P = sympy.Poly(g, *xs, modulus=p)
        lc = int(P.LC(order="grevlex")) % p
        inv = pow(lc, -1, p)
        ss.add(tuple(sorted((tuple(m), int(c) * inv % p) for m, c in zip(P.monoms(), P.coeffs()))))
    return ms == ss


def main():
    max_n1 = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    out = {"protocol_sha256": open("PROTOCOL.sha256").read().split()[0], "levels": [],
           "stopped": {}, "cut": []}
    alive = {a: True for a in ["A1_pkm_mult_chain", "A2_subgroup_flat", "A3_pkm_isogeny_chain",
                               "A4_dickson_chain", "A5_dickson_flat", "A6_small_x", "A7_random_x"]}
    # exact pair-46 facts
    p46 = 5947698030651343
    out["pair46"] = {"p": p46, "p_minus_1": {str(k): v for k, v in factorize(p46 - 1).items()},
                     "p_plus_1": {str(k): v for k, v in factorize(p46 + 1).items()}}
    print("pair46", out["pair46"], flush=True)
    for n1 in range(2, max_n1 + 1):
        d = 2 ** n1
        p = prime_mult(n1)
        E = prime_order_curve(p, "mult")
        pD = prime_dickson(n1)
        ED = prime_order_curve(pD, "dickson")
        lv = {"n1": n1, "d": d, "p": p, "N": E.order, "a": E.a, "b": E.b,
              "pD": pD, "ND": ED.order, "aD": ED.a, "bD": ED.b, "arms": {}}
        # sets
        Vmult = [x for x in range(1, p) if pow(x, d, p) == 1]
        assert len(Vmult) == d
        Ep, h, (steps, chi, Viso) = find_aux_isogeny(p, n1)
        lv["aux"] = {"a": Ep.a, "b": Ep.b, "N": Ep.order, "h": h, "chi": chi, "steps": steps}
        r = rng_for("dick", pD)
        while True:
            x0 = r.randrange(3, pD)
            if legendre(x0 * x0 - 4, pD) == -1:
                lam = dickson2_iter(x0, n1, pD)
                if lam not in (2, pD - 2):
                    break
        Vd = [x for x in range(pD) if dickson2_iter(x, n1, pD) == lam]
        assert len(Vd) == d, (len(Vd), d)
        lv["dickson_lambda"] = lam
        Vsmall = list(range(d))
        Vrand = sorted(rng_for("rand", p).sample(range(p), d))
        sets = {"A1_pkm_mult_chain": (E, Vmult), "A2_subgroup_flat": (E, Vmult),
                "A3_pkm_isogeny_chain": (E, Viso), "A4_dickson_chain": (ED, Vd),
                "A5_dickson_flat": (ED, Vd), "A6_small_x": (E, Vsmall), "A7_random_x": (E, Vrand)}
        # X1/X2 exact check: chained membership L(x)=0 roots equal V
        targets = {}
        for key, EE in (("E", E), ("ED", ED)):
            rr = rng_for("targets", key, EE.p, n1)
            targets[key] = [EE.random_point(rr) for _ in range(NT)]
        # planted control
        rr = rng_for("planted", p, n1)
        P0 = E.random_point(rr)
        mults = [E.mul(j, P0) for j in range(1, d + 1)]
        Vpl = sorted({Q[0] for Q in mults if Q is not None})
        pl_targets = [E.add(E.mul(rr.randint(1, d), P0), E.mul(rr.randint(1, d), P0)) for _ in range(NT)]
        Fpl = len({x for x in Vpl if E.lift_x(x)})
        pl = [relations(E, Vpl, R) for R in pl_targets]
        lv["arms"]["A8_planted"] = {"F": Fpl, "valid": [v for _, v in pl], "nbsols": [s for s, _ in pl],
                                    "pred_valid": 2 * Fpl ** 2 / E.order}
        for arm, (EE, V) in sets.items():
            tg = targets["ED" if EE is ED else "E"]
            Fsz = sum(1 for x in V if EE.lift_x(x) is not None)
            rel = [relations(EE, V, R) for R in tg]
            rec = {"F": Fsz, "V": len(V), "valid": [v for _, v in rel], "nbsols": [s for s, _ in rel],
                   "pred_valid": 2 * Fsz ** 2 / EE.order, "gb": []}
            if alive[arm] and time.time() - T_START < CAP:
                for R in tg:
                    X = R[0]
                    if arm == "A1_pkm_mult_chain":
                        nv, F = mult_chain(EE, X, n1)
                    elif arm == "A4_dickson_chain":
                        nv, F = dickson_chain(EE, X, n1, lam)
                    elif arm == "A3_pkm_isogeny_chain":
                        nv, F = isogeny_chain(EE, X, n1, steps, chi)
                    elif arm == "A2_subgroup_flat":
                        nv, F = flat_system(EE, V, X, {d: 1, 0: EE.p - 1})
                    elif arm == "A5_dickson_flat":
                        u = {0: 1}
                        # D_{2^n1}(x) - lam as a univariate polynomial by iterated squaring
                        cur = {1: 1}
                        for _ in range(n1):
                            sq = {}
                            for i1, c1 in cur.items():
                                for i2, c2 in cur.items():
                                    sq[i1 + i2] = (sq.get(i1 + i2, 0) + c1 * c2) % EE.p
                            sq[0] = (sq.get(0, 0) - 2) % EE.p
                            cur = sq
                        cur[0] = (cur.get(0, 0) - lam) % EE.p
                        nv, F = flat_system(EE, V, X, cur)
                    else:
                        nv, F = flat_system(EE, V, X, prod_univ(V, EE.p))
                    st = run_gb(nv, F, EE, V, R, check_sympy=(n1 <= 3))
                    rec["gb"].append(st)
                    if st["timeout"] or time.time() - T_START > CAP:
                        break
                Ts = [s["time"] for s in rec["gb"]]
                rec["T_mean"] = sum(Ts) / len(Ts)
                if any(s["timeout"] for s in rec["gb"]) or rec["T_mean"] > 15 or len(Ts) < NT:
                    alive[arm] = False
                    out["stopped"][arm] = n1
            else:
                out["cut"].append(f"{arm} n1={n1} (stopped or cap)")
            lv["arms"][arm] = rec
            g = rec["gb"]
            print(f"n1={n1} p={EE.p} {arm}: |V|={len(V)} |F|={Fsz} valid={sum(rec['valid'])} "
                  f"nbsols={sum(rec['nbsols'])} pred={rec['pred_valid']*NT:.2f} "
                  + (f"T={rec['T_mean']:.4f}s Dmax={max(s['D_max'] for s in g)} "
                     f"Dcorr={[s['D_corr'] for s in g]} Dlcorr={[s['Dl_corr'] for s in g]} ok={all(s['ok'] for s in g)} "
                     f"sympy={[s.get('sympy_match') for s in g][:2]}" if g else "no GB"),
                  flush=True)
        # rho on E
        rr = rng_for("rho", p, n1)
        rho = []
        for _ in range(NRHO):
            P = E.random_point(rr)
            k = rr.randrange(1, E.order)
            Q = E.mul(k, P)
            kk, ops, el = None, 0, 0.0
            tot_ops, tot_el = 0, 0.0
            while kk is None:
                kk, ops, el = rho_solve(E, P, Q, rr)
                tot_ops += ops
                tot_el += el
            assert kk == k
            rho.append({"ops": tot_ops, "time": tot_el})
        lv["rho"] = rho
        print(f"n1={n1} rho mean ops={sum(r_['ops'] for r_ in rho)/NRHO:.1f} "
              f"time={sum(r_['time'] for r_ in rho)/NRHO:.5f}s  elapsed={time.time()-T_START:.0f}s",
              flush=True)
        out["levels"].append(lv)
        json.dump(out, open("results.json", "w"), indent=1, default=str)
        if time.time() - T_START > CAP or not any(alive.values()):
            out["cut"].append(f"levels > {n1} not run (cap or all arms stopped)")
            break
    out["wall_seconds"] = time.time() - T_START
    json.dump(out, open("results.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()

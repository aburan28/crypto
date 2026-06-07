#!/usr/bin/env python3
"""
Phase 7 - Semaev relation-yield harness at toy scale, with strict controls.

Question (the plan's highest-value bucket): does an isogenous same-order
neighbour of P-256 generate index-calculus relations any more easily than the
root or than a random same-order curve?

Experiment (Red-team Checks 1-4 baked in):
  * reduce P-256 coefficients mod a toy prime q  ->  root curve E_q (order n_q)
  * build its rational ell-isogenous SAME-ORDER neighbour(s) (our machinery)
  * controls:
        - SAME-ORDER random curves (brute scan for order == n_q)
        - RANDOM curves (any order, generic baseline)
  * factor base FB = the K on-curve points of smallest x-coordinate (so |FB| is
    IDENTICAL across every curve -> same-work comparison)
  * metric: relation yield = fraction of random target points decomposable as a
    sum of m factor-base points (with signs), via Semaev-style decomposition.
  * repeat across two field sizes to check the signal scales (not FF noise).

A Semaev S_3 correctness cross-check confirms the 2-point decomposition oracle.
"""

from __future__ import annotations
import json
import random
import sys

from ecfp import Curve, INF, isogenous_neighbors
import semaev


# --------------------------------------------------------------------------
# toy point-counting and factor bases
# --------------------------------------------------------------------------

def brute_order(E: Curve) -> int:
    p, a, b = E.p, E.a, E.b
    cnt = 1
    for x in range(p):
        r = (x * x * x + a * x + b) % p
        if r == 0:
            cnt += 1
        elif pow(r, (p - 1) // 2, p) == 1:
            cnt += 2
    return cnt


def smallest_x_factor_base(E: Curve, K: int):
    """The K on-curve points of smallest x-coordinate (one per x). |FB| = K."""
    from ecfp import sqrt_mod
    fb = []
    x = 0
    while len(fb) < K and x < E.p:
        r = (x * x * x + E.a * x + E.b) % E.p
        y = sqrt_mod(r, E.p)
        if y is not None:
            fb.append((x, y))
        x += 1
    return fb


def decompose_yield(E: Curve, fb, m: int, trials: int, rng: random.Random):
    """Fraction of random target points that split into m FB points (with +-)."""
    pts = set()
    for (x, y) in fb:
        pts.add((x, y))
        pts.add((x, (-y) % E.p))
    pts = list(pts)
    pts_set = set(pts)

    def two_set():           # all pairwise sums U+V (U,V in pts)
        s = {}
        for U in pts:
            for V in pts:
                s.setdefault(_key(E.add(U, V)), True)
        return s

    twosums = two_set() if m >= 3 else None

    def decomposable(T):
        if m == 2:
            for U in pts:
                if _key(E.add(T, E.neg(U))) in pts_set_key:
                    return True
            return False
        if m == 3:
            for U in pts:
                if _key(E.add(T, E.neg(U))) in twosums:
                    return True
            return False
        raise ValueError("m in {2,3}")

    pts_set_key = {_key(P) for P in pts}
    hit = 0
    for _ in range(trials):
        T = random_point(E, rng)
        if decomposable(T):
            hit += 1
    return hit / trials


def _key(P):
    return "INF" if P is INF else P


def random_point(E: Curve, rng: random.Random):
    from ecfp import sqrt_mod
    while True:
        x = rng.randrange(E.p)
        r = (x * x * x + E.a * x + E.b) % E.p
        y = sqrt_mod(r, E.p)
        if y is not None:
            return (x, y if rng.getrandbits(1) else (-y) % E.p)


# --------------------------------------------------------------------------
# control sampling
# --------------------------------------------------------------------------

def find_same_order_controls(p, n_q, want, exclude_ab, rng, max_tries=4000):
    out = []
    tries = 0
    while len(out) < want and tries < max_tries:
        tries += 1
        a = rng.randrange(p)
        b = rng.randrange(p)
        if (4 * a ** 3 + 27 * b ** 2) % p == 0 or (a, b) in exclude_ab:
            continue
        E = Curve(a, b, p)
        if brute_order(E) == n_q:
            out.append(E)
    return out


def random_controls(p, want, rng):
    out = []
    while len(out) < want:
        a = rng.randrange(p)
        b = rng.randrange(p)
        if (4 * a ** 3 + 27 * b ** 2) % p == 0:
            continue
        out.append(Curve(a, b, p))
    return out


# --------------------------------------------------------------------------
# semaev S_3 cross-check of the 2-decomposition oracle
# --------------------------------------------------------------------------

def semaev3_check(E: Curve, rng: random.Random, trials=200):
    """For random P1,P2 with P3=-(P1+P2), S_3(x1,x2,x3) must vanish."""
    a, b, p = E.a, E.b, E.p
    S3 = semaev.semaev3(a, b, p)
    ok = 0
    for _ in range(trials):
        P1 = random_point(E, rng)
        P2 = random_point(E, rng)
        S = E.add(P1, P2)
        if S is INF:
            continue
        x1, x2, x3 = P1[0], P2[0], S[0]
        val = 0
        for (e1, e2, e3), c in S3.items():
            val = (val + c * pow(x1, e1, p) * pow(x2, e2, p) * pow(x3, e3, p)) % p
        if val == 0:
            ok += 1
        else:
            return False, 0
    return True, ok


# --------------------------------------------------------------------------

def run_one_field(q, ell, K, m, trials, n_ctrl, seed):
    rng = random.Random(seed)
    import params as PP
    a0, b0 = (-3) % q, PP.B % q
    if (4 * a0 ** 3 + 27 * b0 ** 2) % q == 0:
        return None
    E0 = Curve(a0, b0, q)
    n_q = brute_order(E0)
    t_q = q + 1 - n_q

    # isogenous same-order neighbour: use the smallest odd ell that actually has
    # a rational ell-isogeny for this reduced curve.
    neighbor, used_ell = None, None
    for cand in [ell, 3, 5, 7, 11, 13, 17, 19, 23]:
        if cand * cand > 4 * q:
            continue
        neigh = isogenous_neighbors(E0, cand, t_q)
        if neigh:
            neighbor, used_ell = neigh[0][2], cand
            break
    ell = used_ell if used_ell else ell

    same_ctrls = find_same_order_controls(q, n_q, n_ctrl, {(a0, b0)}, rng)
    rand_ctrls = random_controls(q, n_ctrl, rng)

    ok3, _ = semaev3_check(E0, rng)

    def yld(E):
        fb = smallest_x_factor_base(E, K)
        if len(fb) < K:
            return None
        return round(decompose_yield(E, fb, m, trials, random.Random(seed + 17)), 4)

    res = {
        "q": q, "n_q": n_q, "trace_q": t_q, "ell": ell,
        "fb_size_K": K, "m": m, "targets": trials,
        "semaev3_oracle_consistent": ok3,
        "root_yield": yld(E0),
        "neighbor_yield": yld(neighbor) if neighbor else None,
        "neighbor_order_eq_root": (brute_order(neighbor) == n_q) if neighbor else None,
        "same_order_control_yields": [yld(E) for E in same_ctrls],
        "random_control_yields": [yld(E) for E in rand_ctrls],
    }
    so = [y for y in res["same_order_control_yields"] if y is not None]
    res["same_order_control_mean"] = round(sum(so) / len(so), 4) if so else None
    rc = [y for y in res["random_control_yields"] if y is not None]
    res["random_control_mean"] = round(sum(rc) / len(rc), 4) if rc else None
    return res


def main():
    # two field sizes (scaling check), m=3 decomposition, identical |FB|=K
    configs = [
        dict(q=12007, ell=3, K=14, m=3, trials=400, n_ctrl=4, seed=11),
        dict(q=40009, ell=3, K=14, m=3, trials=400, n_ctrl=4, seed=22),
    ]
    out = []
    for cfg in configs:
        print(f"running q={cfg['q']} ...", file=sys.stderr)
        r = run_one_field(**cfg)
        out.append(r)
    summary = {
        "experiment": "Semaev relation-yield, root vs same-order neighbour vs controls",
        "results": out,
        "verdict": _verdict(out),
    }
    print(json.dumps(summary, indent=2))


def _verdict(out):
    lines = []
    for r in out:
        if not r:
            continue
        ry, ny = r["root_yield"], r["neighbor_yield"]
        som, rcm = r["same_order_control_mean"], r["random_control_mean"]
        lines.append({
            "q": r["q"],
            "root": ry, "neighbor": ny,
            "same_order_ctrl_mean": som, "random_ctrl_mean": rcm,
            "neighbor_vs_root_ratio": round(ny / ry, 3) if (ry and ny) else None,
            "neighbor_vs_same_order_ratio": round(ny / som, 3) if (som and ny) else None,
        })
    return lines


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Item 03 floor-pullback: push crater factor bases through a descending
v-isogeny (Velu) on toy analogues of pair 46 and measure 2-term decomposition
yield on the floor against random-x controls.  See PROTOCOL.md (frozen)."""
import json, math, random, sys, time
import numpy as np

J0 = -884736
OUT = sys.argv[1] if len(sys.argv) > 1 else "results.json"
T0 = time.time()


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def factor(n):
    f, q = [], 2
    while q * q <= n:
        if n % q == 0:
            f.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        f.append(n)
    return f


def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def sqrt_mod(a, p):
    a %= p
    if a == 0:
        return 0
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    # Tonelli-Shanks
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2; s += 1
    z = 2
    while legendre(z, p) != -1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


class Curve:
    def __init__(self, a, b, p):
        self.a, self.b, self.p = a % p, b % p, p

    def f(self, x):
        p = self.p
        return (x * x * x + self.a * x + self.b) % p

    def add(self, P, Q):
        p = self.p
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P; x2, y2 = Q
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            l = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, p) % p
        else:
            l = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (l * l - x1 - x2) % p
        return (x3, (l * (x1 - x3) - y1) % p)

    def neg(self, P):
        return None if P is None else (P[0], (-P[1]) % self.p)

    def mul(self, k, P):
        R = None
        if k < 0:
            k, P = -k, self.neg(P)
        while k:
            if k & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            k >>= 1
        return R

    def rand_point(self, rng):
        p = self.p
        while True:
            x = rng.randrange(p)
            y2 = self.f(x)
            if legendre(y2, p) == 1:
                y = sqrt_mod(y2, p)
                return (x, y if rng.random() < .5 else (-y) % p)

    def j(self):
        p = self.p
        num = 1728 * 4 * pow(self.a, 3, p)
        den = (4 * pow(self.a, 3, p) + 27 * self.b * self.b) % p
        return num * pow(den, -1, p) % p

    def twist(self, d):
        p = self.p
        return Curve(self.a * d * d, self.b * d * d * d, p)


def velu(E, kx):
    """kx: x-coords of a half-kernel (odd prime degree), in F_p."""
    p = E.p
    vs, us = [], []
    for xq in kx:
        gx = (3 * xq * xq + E.a) % p
        vs.append(2 * gx % p)
        us.append(4 * E.f(xq) % p)
    V = sum(vs) % p
    W = sum(u + xq * vq for u, xq, vq in zip(us, kx, vs)) % p
    E2 = Curve(E.a - 5 * V, E.b - 7 * W, p)

    def phi(P):
        if P is None:
            return None
        x, y = P
        X, D = x, 1
        for xq, vq, uq in zip(kx, vs, us):
            if x == xq:
                return None
        # sums with batched inversion avoided: direct (toy sizes)
        X, Dx = x, 1
        for xq, vq, uq in zip(kx, vs, us):
            iv = pow(x - xq, -1, p)
            iv2 = iv * iv % p
            X += vq * iv + uq * iv2
            Dx -= vq * iv2 + 2 * uq * iv2 * iv
        return (X % p, y * Dx % p)
    return E2, phi


def log_table(E, G, N):
    p = E.p
    L = np.zeros(p, dtype=np.int32)
    YY = np.full(p, -1, dtype=np.int32)
    P = G
    for i in range(1, (N - 1) // 2 + 1):
        L[P[0]] = i; YY[P[0]] = P[1]
        P = E.add(P, G)
    return L, YY


def dlog(P, L, YY, N):
    i = int(L[P[0]])
    assert i > 0
    return i if int(YY[P[0]]) == P[1] else N - i


def yield2(logs, N):
    a = np.asarray(logs, dtype=np.int64)
    s = (a[:, None] + a[None, :]) % N
    hit = np.zeros(N, dtype=bool)
    hit[s.ravel()] = True
    return float(hit.sum()) / N


def smallx_points(E, m):
    pts, x = [], 0
    while len(pts) < m:
        y2 = E.f(x)
        if legendre(y2, E.p) == 1:
            y = sqrt_mod(y2, E.p)
            pts += [(x, y), (x, (-y) % E.p)]
        x += 1
    return pts[:m]


def primitive_root(p):
    fs = factor(p - 1)
    g = 2
    while any(pow(g, (p - 1) // q, p) == 1 for q in fs):
        g += 1
    return g


def coset_points(E, m):
    p = E.p
    g = primitive_root(p)
    divs = sorted(d for d in range(1, p) if (p - 1) % d == 0)
    for d in divs:
        if d < m // 2:
            continue
        zeta = pow(g, (p - 1) // d, p)
        pts, x = [], g
        for k in range(d):
            y2 = E.f(x)
            if legendre(y2, p) == 1:
                y = sqrt_mod(y2, p)
                pts += [(x, y), (x, (-y) % p)]
            x = x * zeta % p
        if len(pts) >= m:
            return pts[:m], {"d": d, "g": g}
    raise RuntimeError("no coset")


def randomx_points(E, m, rng):
    pts, seen = [], set()
    while len(pts) < m:
        x = rng.randrange(E.p)
        if x in seen:
            continue
        seen.add(x)
        y2 = E.f(x)
        if legendre(y2, E.p) == 1:
            y = sqrt_mod(y2, E.p)
            pts += [(x, y), (x, (-y) % E.p)]
    return pts[:m]


def nonresidue(p):
    d = 2
    while legendre(d, p) != -1:
        d += 1
    return d


def run_instance(v):
    t_start = time.time()
    p = 5 * v * v - v + 1
    N = 5 * v * v - 2 * v + 4
    rng = random.Random(1000 + v)
    j = J0 % p
    E = Curve(3 * j * (1728 - j), 2 * j * (1728 - j) ** 2, p)
    dn = nonresidue(p)
    if not all(E.mul(N, E.rand_point(rng)) is None for _ in range(8)):
        E = E.twist(dn)
    assert all(E.mul(N, E.rand_point(rng)) is None for _ in range(8))
    assert E.j() == j
    Et = E.twist(dn)
    assert all(Et.mul(5 * v * v, Et.rand_point(rng)) is None for _ in range(8))
    # two independent order-v points on the twist
    Q1 = None
    while Q1 is None:
        Q1 = Et.mul(5, Et.rand_point(rng))  # PROTOCOL said [5v]; twist exponent is 5v so [5] is needed (deviation D1)
    sub1 = set()
    P = Q1
    for _ in range(v):
        if P is not None:
            sub1.add(P)
        P = Et.add(P, Q1)
    Q2 = None
    while Q2 is None or Q2 in sub1:
        Q2 = Et.mul(5, Et.rand_point(rng))
    idn = pow(dn, -1, p)
    chosen, tried = None, []
    gens = [Q2] + [Et.add(Q1, Et.mul(k, Q2)) for k in range(v)]
    for gi, K in enumerate(gens):
        kx, P = [], K
        for _ in range((v - 1) // 2):
            kx.append(P[0] * idn % p)  # twist X -> E-x via x = X/d
            P = Et.add(P, K)
        E2, phi = velu(E, kx)
        j2 = E2.j()
        tried.append({"subgroup": gi, "j_codomain": j2, "horizontal": j2 == j})
        if j2 != j:
            chosen = (E2, phi, gi)
            break
    E2, phi, gi = chosen
    # certify floor
    assert all(E2.mul(N, E2.rand_point(rng)) is None for _ in range(8))
    E2t = E2.twist(nonresidue(p))
    cyc = any(E2t.mul(5 * v, E2t.rand_point(rng)) is not None for _ in range(8))
    assert cyc, "floor twist not cyclic"
    # homomorphism check
    for _ in range(20):
        A, B = E.rand_point(rng), E.rand_point(rng)
        lhs = phi(E.add(A, B)); rhs = E2.add(phi(A), phi(B))
        assert lhs == rhs
        assert E2.f(lhs[0]) == lhs[1] * lhs[1] % p
    G = E.rand_point(rng)
    G2 = phi(G)
    t1 = time.time()
    Lc, Yc = log_table(E, G, N)
    Lf, Yf = log_table(E2, G2, N)
    t_tables = time.time() - t1
    m = 2 * round(math.sqrt(N) / 2)
    res = {"v": v, "p": p, "N": N, "m": m, "kernel_subgroup_index": gi,
           "subgroups_tried": tried, "j_floor": E2.j(), "floor_twist_cyclic": cyc,
           "floor_a": E2.a, "floor_b": E2.b, "crater_a": E.a, "crater_b": E.b}

    def floor_yield(pts):
        return yield2([dlog(P, Lf, Yf, N) for P in pts], N)

    def xfrac(pts):
        return sum(1 for P in pts if P[0] < p / 16) / len(pts)

    S1 = smallx_points(E, m)
    S2, cinfo = coset_points(E, m)
    arms = {}
    for name, S in (("A1_pull_smallx", S1), ("A2_pull_coset", S2)):
        img = [phi(P) for P in S]
        y_floor = floor_yield(img)
        y_crater = yield2([dlog(P, Lc, Yc, N) for P in S], N)
        arms[name] = {"Y": y_floor, "Y_crater_X0": y_crater,
                      "X0_exact_equal": y_floor == y_crater,
                      "xfrac_lt_p16": xfrac(img), "n_points": len(set(img))}
    C1 = smallx_points(E2, m)
    arms["C1_floor_smallx"] = {"Y": floor_yield(C1), "xfrac_lt_p16": xfrac(C1)}
    c0 = []
    for s in range(10):
        pts = randomx_points(E2, m, random.Random(10 * v + s))
        c0.append({"seed": s, "Y": floor_yield(pts), "xfrac_lt_p16": xfrac(pts)})
    arms["C0_random_x"] = c0
    pc = []
    P = G2
    for i in range(m // 2):
        pc += [P, E2.neg(P)]
        P = E2.add(P, G2)
    arms["Pplus_interval"] = {"Y": floor_yield(pc)}
    ys = np.array([c["Y"] for c in c0])
    mu, sd = float(ys.mean()), float(ys.std(ddof=1))
    half = 2.262 * sd / math.sqrt(len(ys))  # t_{0.975,9}
    stats = {"mu0": mu, "sd0": sd, "mu0_ci95": [mu - half, mu + half],
             "random_set_heuristic": 1 - math.exp(-(m * m / 2) / N)}
    for name in ("A1_pull_smallx", "A2_pull_coset", "C1_floor_smallx", "Pplus_interval"):
        Y = arms[name]["Y"]
        stats[name] = {"z": (Y - mu) / sd, "ratio": Y / mu,
                       "ratio_ci95": [Y / (mu + half), Y / (mu - half)]}
    res.update({"arms": arms, "stats": stats, "coset": cinfo,
                "seconds": time.time() - t_start, "seconds_tables": t_tables})
    return res


def main():
    lo, hi = (920, 2000) if "--ext" in sys.argv else (29, 919)
    qual = [v for v in range(lo, hi + 1) if is_prime(v) and is_prime(5 * v * v - v + 1)
            and is_prime(5 * v * v - 2 * v + 4)]
    L = len(qual)
    idx = sorted(set(round(k * (L - 1) / 5) for k in range(6)))
    chosen = [qual[i] for i in idx] if "--ext" not in sys.argv else qual
    print("qualifying v:", qual, flush=True)
    print("chosen:", chosen, flush=True)
    out = {"qualifying_v": qual, "chosen_v": chosen, "instances": [], "cut": []}
    for v in chosen:
        if time.time() - T0 > 30 * 60:
            out["cut"].append(v); continue
        r = run_instance(v)
        st = r["stats"]
        print(f"v={v} p={r['p']} m={r['m']} mu0={st['mu0']:.4f} sd0={st['sd0']:.5f} "
              + " ".join(f"{k}:z={st[k]['z']:+.2f}" for k in
                         ("A1_pull_smallx", "A2_pull_coset", "C1_floor_smallx", "Pplus_interval"))
              + f" X0={r['arms']['A1_pull_smallx']['X0_exact_equal']},{r['arms']['A2_pull_coset']['X0_exact_equal']}"
              + f" t={r['seconds']:.1f}s", flush=True)
        out["instances"].append(r)
        json.dump(out, open(OUT, "w"), indent=1)
    out["wall_seconds"] = time.time() - T0
    json.dump(out, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()

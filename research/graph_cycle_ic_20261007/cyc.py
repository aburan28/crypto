"""Graph-cycle index calculus with the linear (A, P) pair oracle: no linear algebra.

Curve y^2 + xy = x^3 + 1 over F_{2^n} with #E = 4r, r prime (n in {13, 19, 23}).
P generates the order-r subgroup, Q = kP.

Walk: W_{j+1} = W_j + R_{h(W_j)} with R_i = a_i P + b_i Q, tracking W_j = a_j P + b_j Q.
Test: one linear (A, P) solve (../geometric_v_linear_20261006/glin.py) asks whether
W_j = s1 A + s2 B with A, B lifts of abscissae in a geometric V.  Each success is an
edge between nodes A and B with label
    s1 lam_A + s2 lam_B = 4 (a_j + b_j k)   (mod r),   lam_X := log_P(4 X).
Union-find with affine potentials over Z/r (node value = sigma * root + alpha + beta k)
closes the first usable cycle:
  - an even cycle gives alpha + beta k = 0, so k = -alpha / beta;
  - an odd cycle fixes the root's value; a second constraint on that component
    then gives k.
Every recovered k is verified by k P == Q.

Two modes.  'walk' drives the tests with an r-adding walk and also watches that
walk for a rho collision, recording which solves first: a walk longer than ~sqrt(r)
collides, so a method needing ~2^(n-l) > sqrt(r) tests is dominated by the rho in
its own walk.  'fresh' draws independent uniform W = aP + bQ per test (scalar
multiplications counted separately) to measure the graph law in isolation.

Reference: a separate plain Pollard rho, same r-adding walk (no negation or
Frobenius speedup), collisions found by storing every visited point.
"""
import json
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "geometric_v_linear_20261006"))
sys.path.insert(0, os.path.join(HERE, "..", "linearization_reach_20260930"))
import glin  # noqa: E402
import lr    # noqa: E402

R_PRIME = {13: 2003, 19: 130873, 23: 2095853}   # r = #E / 4, prime (trace recurrence)
NSTEPS = 20


def smul(F, k, P):
    R, A = None, P
    while k:
        if k & 1:
            R = lr.add(F, R, A)
        A = lr.add(F, A, A)
        k >>= 1
    return R


def generator(F, r, rng):
    while True:
        T = lr.lift(F, rng.getrandbits(F.n))
        if T is None:
            continue
        P = smul(F, 4, T)
        if P is not None and smul(F, r, P) is None:
            return P


class Graph:
    """Union-find; value(v) = sigma * value(root) + alpha + beta * k (mod r)."""

    def __init__(self, r):
        self.r = r
        self.parent, self.pot, self.known = {}, {}, {}

    def find(self, v):
        if v not in self.parent:
            self.parent[v], self.pot[v] = v, (1, 0, 0)
            return v, (1, 0, 0)
        p = self.parent[v]
        if p == v:
            return v, (1, 0, 0)
        root, (s2, a2, b2) = self.find(p)
        s1, a1, b1 = self.pot[v]            # v = s1 p + a1 + b1 k ; p = s2 root + a2 + b2 k
        comp = (s1 * s2, (s1 * a2 + a1) % self.r, (s1 * b2 + b1) % self.r)
        self.parent[v], self.pot[v] = root, comp
        return root, comp

    def _solve_k(self, a, b):
        """a + b k = 0 (mod r)."""
        if b % self.r == 0:
            return None
        return (-a * pow(b, -1, self.r)) % self.r

    def add_edge(self, u, su, v, sv, ca, cb):
        """su val(u) + sv val(v) = ca + cb k.  Returns k if determined, else None."""
        r = self.r
        ru, (pu, au, bu) = self.find(u)
        rv, (pv, av, bv) = self.find(v)
        # su(pu X_ru + fu) + sv(pv X_rv + fv) = c
        cu, cv = su * pu, sv * pv
        ra = (ca - su * au - sv * av) % r
        rb = (cb - su * bu - sv * bv) % r
        # cu X_ru + cv X_rv = ra + rb k
        if ru != rv:
            # X_rv = cv * (ra + rb k - cu X_ru)  (cv = +-1 is its own inverse)
            self.parent[rv] = ru
            self.pot[rv] = (-cv * cu, (cv * ra) % r, (cv * rb) % r)
            if rv in self.known:
                ka, kb = self.known.pop(rv)    # X_rv = ka + kb k
                # ka + kb k = s X_ru + cv ra + cv rb k, s = -cv cu = +-1
                #   ->  X_ru = s (ka - cv ra) + s (kb - cv rb) k
                s = -cv * cu
                xa, xb = (s * (ka - cv * ra)) % r, (s * (kb - cv * rb)) % r
                return self._fix_root(ru, xa, xb)
            return None
        coef = (cu + cv)
        if coef == 0:
            return self._solve_k(-ra, -rb)     # 0 = ra + rb k
        inv = pow(coef % r, -1, r)             # coef = +-2
        return self._fix_root(ru, (ra * inv) % r, (rb * inv) % r)

    def _fix_root(self, root, xa, xb):
        if root in self.known:
            ka, kb = self.known[root]
            return self._solve_k(xa - ka, xb - kb)
        self.known[root] = (xa, xb)
        return None


def cycle_ic(n, l, k, P, Q, rng, cap, mode):
    """mode 'walk': r-adding walk, with rho collision detection on the same walk;
    returns which of the two solves first.  mode 'fresh': independent uniform
    W = aP + bQ per test (two scalar multiplications each, counted separately)."""
    F, r = lr.Field(n), R_PRIME[n]
    base = glin.Base(F, l, "geometric", rng)
    nodes = sum(1 for x in glin_span(base) if x not in (0, 1) and lr.lift(F, x))
    steps = []
    for _ in range(NSTEPS):
        sa, sb = rng.randrange(r), rng.randrange(r)
        steps.append((sa, sb, lr.add(F, smul(F, sa, P), smul(F, sb, Q))))
    a, b = rng.randrange(1, r), rng.randrange(1, r)
    W = lr.add(F, smul(F, a, P), smul(F, b, Q))
    G = Graph(r)
    seen = {}
    tests = edges = scalar_mults = 0
    while tests < cap:
        if mode == "walk" and W in seen:
            a2, b2 = seen[W]
            if (b2 - b) % r:
                kk = ((a - a2) * pow((b2 - b) % r, -1, r)) % r
                return dict(solved_by="rho", k_ok=(smul(F, kk, P) == Q), tests=tests,
                            edges=edges, nodes=nodes, scalar_mults=scalar_mults)
        if mode == "walk":
            seen[W] = (a, b)
        if W is not None:
            tests += 1
            _, cand = glin.solve(base, W[0])
            for A_, P_ in (cand or []):
                hit = edge_from(F, base, W, A_, P_)
                if hit is None:
                    continue
                (x2, s2), (x3, s3) = hit
                edges += 1
                kk = G.add_edge(x2, s2, x3, s3, (4 * a) % r, (4 * b) % r)
                if kk is not None:
                    return dict(solved_by="cycle", k_ok=(smul(F, kk, P) == Q), tests=tests,
                                edges=edges, nodes=nodes, scalar_mults=scalar_mults)
                break
        if mode == "walk":
            i = (W[0] % NSTEPS) if W is not None else 0
            sa, sb, SR = steps[i]
            W = lr.add(F, W, SR)
            a, b = (a + sa) % r, (b + sb) % r
        else:
            a, b = rng.randrange(1, r), rng.randrange(1, r)
            W = lr.add(F, smul(F, a, P), smul(F, b, Q))
            scalar_mults += 2
    return dict(solved_by=None, k_ok=False, tests=tests, edges=edges, nodes=nodes,
                scalar_mults=scalar_mults)


def glin_span(base):
    span = [0]
    for v in base.V:
        span = span + [s ^ v for s in span]
    return span


def edge_from(F, base, W, A, Pp):
    """If W = s2 A2 + s3 A3 for canonical lifts A2, A3 of the roots, return ((x2,s2),(x3,s3))."""
    if A == 0:
        return None
    roots = lr.quad_roots(F, 1, A, Pp)
    if len(roots) != 2 or not all(base.in_V(z) and z not in (0, 1) for z in roots):
        return None
    pts = [lr.lift(F, z) for z in roots]
    if None in pts:
        return None
    for s2 in (1, -1):
        for s3 in (1, -1):
            q2 = pts[0] if s2 == 1 else lr.neg(pts[0])
            q3 = pts[1] if s3 == 1 else lr.neg(pts[1])
            if lr.add(F, q2, q3) == W:
                return (roots[0], s2), (roots[1], s3)
    return None


def rho(n, k, P, Q, rng):
    """Plain r-adding-walk rho; stores visited points; returns iterations to collision-solve."""
    F, r = lr.Field(n), R_PRIME[n]
    steps = []
    for _ in range(NSTEPS):
        a, b = rng.randrange(r), rng.randrange(r)
        steps.append((a, b, lr.add(F, smul(F, a, P), smul(F, b, Q))))
    a, b = rng.randrange(1, r), rng.randrange(1, r)
    W = lr.add(F, smul(F, a, P), smul(F, b, Q))
    seen, it = {}, 0
    while True:
        it += 1
        if W in seen:
            a2, b2 = seen[W]
            if (b2 - b) % r:
                kk = ((a - a2) * pow((b2 - b) % r, -1, r)) % r
                return dict(iterations=it, k_ok=(smul(F, kk, P) == Q))
        seen[W] = (a, b)
        i = (W[0] % NSTEPS) if W is not None else 0
        sa, sb, SR = steps[i]
        W = lr.add(F, W, SR)
        a, b = (a + sa) % r, (b + sb) % r


def run(n, ls, instances, seed):
    rng = random.Random(seed)
    F, r = lr.Field(n), R_PRIME[n]
    out = []
    for it in range(instances):
        P = generator(F, r, rng)
        k = rng.randrange(1, r)
        Q = smul(F, k, P)
        rr = rho(n, k, P, Q, rng)
        row = dict(n=n, instance=it, r=r, rho_iterations=rr["iterations"], rho_ok=rr["k_ok"],
                   S_rho=round(rr["iterations"] / math.sqrt(r), 3))
        for l in ls:
            cap = 2 ** (n - l + 5)
            for mode in ("walk", "fresh"):
                c = cycle_ic(n, l, k, P, Q, rng, cap, mode)
                row[f"l{l}_{mode}"] = dict(c, S_tests=round(c["tests"] / math.sqrt(r), 3),
                                           edges_per_node=round(c["edges"] / c["nodes"], 3) if c["nodes"] else None)
        out.append(row)
        print(json.dumps(row), flush=True)
    return out


if __name__ == "__main__":
    n, instances, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    ls = [int(x) for x in sys.argv[4].split(",")]
    res = run(n, ls, instances, seed)
    if len(sys.argv) > 5:
        json.dump(res, open(sys.argv[5], "w"), indent=1)

#!/usr/bin/env python3
"""Independent check of the checkable claims in the 2026-09-26 Frobenius
quotient follow-up note (see README.md).  Pure Python, no solver, no inputs.

Checks, on E: y^2 + xy = x^3 + 1 over F_2^n:
  1. curve / odd-subgroup orders from the trace recurrence (n = 13, 19, 29, 131);
  2. the line-formulation identities, including tangent cases;
  3. the quadratic subgroup certificate against true 4E membership;
  4. the uniform-target compatibility mean (2n)^4/|H| on an independently
     chosen factor space, by exhaustive pair matching.
Writes a JSON summary to stdout.  Seeded; deterministic.
"""
import itertools, json, math, random, sys

MODULI = {  # low-weight irreducibles; 29 and 131 match the note / ECC2K-130
    13: (1 << 13) | 0b11011,             # X^13+X^4+X^3+X+1
    19: (1 << 19) | 0b100111,            # X^19+X^5+X^2+X+1
    29: (1 << 29) | 0b101,               # X^29+X^2+1
    131: (1 << 131) | (1 << 13) | 0b111, # X^131+X^13+X^2+X+1
}


class F:
    def __init__(self, n):
        self.n, self.m = n, MODULI[n]

    def mul(self, a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a >> self.n:
                a ^= self.m
        return r

    def sq(self, a):
        return self.mul(a, a)

    def pow(self, a, e):
        r = 1
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.sq(a)
            e >>= 1
        return r

    def inv(self, a):
        assert a
        return self.pow(a, (1 << self.n) - 2)

    def sqrt(self, a):
        return self.pow(a, 1 << (self.n - 1))

    def tr(self, a):
        t, x = 0, a
        for _ in range(self.n):
            t ^= x
            x = self.sq(x)
        assert t in (0, 1)
        return t

    def half_trace(self, c):  # solves z^2+z=c when Tr(c)=0, n odd
        z, x = 0, c
        for _ in range((self.n + 1) // 2):
            z ^= x
            x = self.sq(self.sq(x))
        assert self.sq(z) ^ z == c
        return z


def orders(n):
    t0, t1 = 2, -1  # #E(F_2) = 4 for y^2+xy=x^3+1
    for _ in range(n - 1):
        t0, t1 = t1, -t1 - 2 * t0
    N = 2 ** n + 1 - t1
    return N, N // 4


def is_probable_prime(m):
    if m < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if m % p == 0:
            return m == p
    d, s = m - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, m)
        if x in (1, m - 1):
            continue
        for _ in range(s - 1):
            x = x * x % m
            if x == m - 1:
                break
        else:
            return False
    return True


class Curve:
    O = None

    def __init__(self, n):
        self.f = F(n)
        self.N, self.h = orders(n)

    def on(self, P):
        if P is None:
            return True
        f, (x, y) = self.f, P
        return f.sq(y) ^ f.mul(x, y) == f.mul(f.sq(x), x) ^ 1

    def neg(self, P):
        return None if P is None else (P[0], P[0] ^ P[1])

    def add(self, P, Q):
        f = self.f
        if P is None:
            return Q
        if Q is None:
            return P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if y1 ^ y2 == x2 or x1 == 0:  # Q = -P
                return None
            lam = x1 ^ f.mul(y1, f.inv(x1))
            x3 = f.sq(lam) ^ lam
            y3 = f.sq(x1) ^ f.mul(lam ^ 1, x3)
            return (x3, y3)
        lam = f.mul(y1 ^ y2, f.inv(x1 ^ x2))
        x3 = f.sq(lam) ^ lam ^ x1 ^ x2
        y3 = f.mul(lam, x1 ^ x3) ^ x3 ^ y1
        return (x3, y3)

    def mul(self, k, P):
        R = None
        while k:
            if k & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            k >>= 1
        return R

    def frob(self, P):
        return None if P is None else (self.f.sq(P[0]), self.f.sq(P[1]))

    def lift(self, x):
        """A point with this x-coordinate, or None if x is not rational."""
        f = self.f
        if x == 0:
            return (0, 1)
        c = x ^ f.inv(f.sq(x))  # y = x z, z^2+z = x + 1/x^2
        if f.tr(c):
            return None
        return (x, f.mul(x, f.half_trace(c)))

    def random_point(self, rng):
        while True:
            P = self.lift(rng.getrandbits(self.f.n))
            if P is not None:
                return P if rng.getrandbits(1) else self.neg(P)


# ---------------------------------------------------------------- 1. orders
def check_orders():
    out = {}
    for n in (13, 19, 29, 131):
        N, h = orders(n)
        out[n] = {"curve_order": N, "N_mod_8": N % 8, "odd_subgroup_order": h,
                  "odd_subgroup_is_odd": h % 2 == 1,
                  "odd_subgroup_probable_prime": is_probable_prime(h),
                  "log2_odd_subgroup": math.log2(h)}
    return out


# ----------------------------------------------------- 2. line identities
def check_line(n, trials, rng):
    E, f = Curve(n), F(n)
    ok = tangents = 0
    for i in range(trials):
        A = E.random_point(rng)
        B = A if i % 5 == 0 else E.random_point(rng)  # every 5th: tangent
        S = E.add(A, B)
        if S is None or A[0] == B[0] and A != B:
            continue
        R = S  # A + B + (-R) = 0: A, B, -R collinear
        a, b = R
        t, u = A[0], B[0]
        if A == B:
            tangents += 1
            ell = t ^ f.mul(A[1], f.inv(t))  # tangent slope
        else:
            ell = f.mul(A[1] ^ B[1], f.inv(t ^ u))
        nu = b ^ a ^ f.mul(a, ell)
        e1 = u == t ^ f.sq(ell) ^ ell ^ a
        e2 = f.mul(t, u) == f.mul(a, f.sq(ell)) ^ f.sq(a) ^ b ^ a
        # reconstruct the intermediate points from (t, ell, nu)
        e3 = E.on((t, f.mul(ell, t) ^ nu)) and E.on((u, f.mul(ell, u) ^ nu))
        e4 = (t, f.mul(ell, t) ^ nu) == A and (u, f.mul(ell, u) ^ nu) == B
        ok += e1 and e2 and e3 and e4
    return {"trials": trials, "passed": ok, "tangent_cases": tangents}


# ------------------------------------------------- 3. subgroup certificate
def certificate(f, x):
    """Exists q: q^2 + sqrt(x) q + 1 = 0, Tr(q) = 0 (with Tr(x) = 0, x != 0)."""
    if x == 0 or f.tr(x):
        return False
    c = f.sqrt(x)
    rhs = f.inv(f.sq(c))  # q = c z, z^2 + z = 1/c^2
    if f.tr(rhs):
        return False
    z = f.half_trace(rhs)
    for q in (f.mul(c, z), f.mul(c, z ^ 1)):
        assert f.sq(q) ^ f.mul(c, q) ^ 1 == 0
        if f.tr(q) == 0:
            return True
    return False


def check_certificate(n, xs):
    E = Curve(n)
    accepted = agree = 0
    for x in xs:
        P = E.lift(x)
        member = P is not None and x != 0 and E.mul(E.h, P) is None
        cert = certificate(E.f, x)
        accepted += cert
        agree += cert == member
    return {"checked": len(xs), "accepted": accepted, "agree": agree}



def _orbits(E, span):
    seen, orbits = set(), []
    for v in sorted(span):
        if v == 0 or E.f.tr(v):
            continue
        P = E.lift(v)
        if P is None or E.mul(E.h, P) is not None or P in seen:
            continue
        orb, Q = [], P
        for _ in range(E.f.n):
            orb += [Q, E.neg(Q)]
            Q = E.frob(Q)
        assert Q == P and len(set(orb)) == 2 * E.f.n
        seen.update(orb)
        orbits.append(orb)
    return orbits


# ------------------------------------------------ 4. compatibility density
def check_density(n, s, n_tuples, n_targets, rng):
    E, f = Curve(n), F(n)
    # independent factor space: first seeded random s-dim subspace (intersected
    # with Tr = 0) that yields at least four valid orbits in H = 4E
    while True:
        basis = [rng.getrandbits(n) for _ in range(s)]
        span = {0}
        for b in basis:
            span |= {v ^ b for v in span}
        if len(span) == 1 << s and len(_orbits(E, span)) >= 4:
            break
    orbits = _orbits(E, span)
    all_tuples = list(itertools.combinations(range(len(orbits)), 4))
    tuples = all_tuples if len(all_tuples) <= n_tuples else rng.sample(all_tuples, n_tuples)
    gen = None
    while gen is None:
        gen = E.mul(4, E.random_point(rng))
    targets = []
    while len(targets) < n_targets:
        T = E.mul(rng.randrange(1, E.h), gen)
        if T is not None:
            targets.append(T)
    total = positive = cases = 0
    for (i, j, k, l) in tuples:
        left = {}
        for P in orbits[i]:
            for Q in orbits[j]:
                S = E.add(P, Q)
                left[S] = left.get(S, 0) + 1
        right = [E.add(P, Q) for P in orbits[k] for Q in orbits[l]]
        for T in targets:
            cnt = sum(left.get(E.add(T, E.neg(S)), 0) for S in right)
            cases += 1
            total += cnt
            positive += cnt > 0
    return {"n": n, "s": s, "basis": [hex(b) for b in basis],
            "valid_orbits": len(orbits),
            "tuples": len(tuples), "targets": n_targets,
            "positive_cases": positive, "cases": cases,
            "observed_mean": total / cases,
            "uniform_target_mean": (2 * n) ** 4 / E.h}


def main():
    rng = random.Random(20260926)
    res = {"orders": check_orders()}
    res["log2_full_size_mean"] = math.log2(262 ** 4) - math.log2(orders(131)[1])
    res["line"] = {n: check_line(n, 100, rng) for n in (13, 19)}
    res["certificate"] = {
        13: check_certificate(13, range(1, 1 << 13)),
        19: check_certificate(19, [rng.getrandbits(19) or 1 for _ in range(4000)]),
    }
    res["density"] = [check_density(13, 4, 5, 8, rng),
                      check_density(19, 6, 32, 8, rng),
                      check_density(29, 7, 16, 8, rng)]
    json.dump(res, sys.stdout, indent=2, default=str)
    print()


if __name__ == "__main__":
    main()

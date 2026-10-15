"""Frobenius characteristic polynomials over F_2 from point counts, and targets."""
from fractions import Fraction

Q = 2


def charpoly_from_counts(N):
    """N[k-1] = #C(F_{2^k}), k = 1..g.  Returns P_C(T) descending, length 2g+1,
    or None if Newton's identities give a non-integer (invalid input)."""
    g = len(N)
    S = [Q ** k + 1 - N[k - 1] for k in range(1, g + 1)]  # power sums of roots
    c = [1]
    for k in range(1, g + 1):
        acc = S[k - 1] + sum(c[i] * S[k - i - 1] for i in range(1, k))
        if acc % k:
            return None
        c.append(-acc // k)
    full = c + [Q ** (g - j) * c[j] for j in range(g - 1, -1, -1)]
    return tuple(full)


def counts_from_charpoly(P, kmax):
    """Inverse: P descending (deg 2g); returns #C(F_{2^k}) for k=1..kmax."""
    n = len(P) - 1
    c = list(P)
    S = []
    for k in range(1, kmax + 1):
        acc = -k * c[k] if k <= n else 0
        acc -= sum(c[i] * S[k - i - 1] for i in range(1, min(k, n + 1)) if k - i - 1 >= 0)
        S.append(acc)
    return [Q ** k + 1 - S[k - 1] for k in range(1, kmax + 1)]


def polymul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return tuple(r)


def divides(d, p):
    """Exact division test of monic integer polys (descending)."""
    p = list(p)
    n, m = len(p) - 1, len(d) - 1
    if m > n:
        return False
    for i in range(n - m + 1):
        coef = p[i]
        if coef:
            for j in range(m + 1):
                p[i + j] -= coef * d[j]
    return all(x == 0 for x in p[n - m + 1:])


def weil_restriction(t, n):
    s0, s1 = 2, t
    for _ in range(1, n):
        s0, s1 = s1, s1 * t - 2 * s0
    P = [0] * (2 * n + 1)
    P[0], P[n], P[2 * n] = 1, -s1, 2 ** n
    return tuple(P)


def polydiv(p, d):
    p = list(p)
    n, m = len(p) - 1, len(d) - 1
    qd = []
    for i in range(n - m + 1):
        coef = p[i]
        qd.append(coef)
        for j in range(m + 1):
            p[i + j] -= coef * d[j]
    assert all(x == 0 for x in p[n - m + 1:])
    return tuple(qd)


def targets():
    T = {}
    for t, name in ((-1, "E0"), (1, "E1")):
        E = (1, -t, 2)
        for n in (3, 5, 7):
            W = weil_restriction(t, n)
            A = polydiv(W, E)
            T[f"W{n}({name})"] = W
            if n != 7:
                T[f"A{n}({name})"] = A
            else:
                T[f"A7({name})"] = A
    # factors of A7 (computed in PARI/GP; verified below)
    T["B7a(E0)"] = (1, -4, 9, -15, 18, -16, 8)
    T["B7b(E0)"] = (1, 3, 2, -1, 4, 12, 8)
    T["B7a(E1)"] = (1, -3, 2, 1, 4, -12, 8)
    T["B7b(E1)"] = (1, 4, 9, 15, 18, 16, 8)
    for s in ("E0", "E1"):
        assert polymul(T[f"B7a({s})"], T[f"B7b({s})"]) == T[f"A7({s})"]
    return T


def dim(P):
    return (len(P) - 1) // 2


def is_weil_valid(P, g):
    """Necessary condition: place counts from N_k nonnegative integers for k <= 2g."""
    N = counts_from_charpoly(P, 2 * g)
    from sympy import mobius, divisors
    for k in range(1, 2 * g + 1):
        s = sum(mobius(k // d) * N[d - 1] for d in divisors(k))
        if s < 0 or s % k:
            return False
    return True


def weil_polys(g):
    """All monic integer Weil 2-polynomials of degree 2g (roots |z| = sqrt 2)."""
    import itertools
    import numpy as np
    if g == 0:
        return [(1,)]
    from math import comb
    bounds = [int(comb(2 * g, i) * 2 ** (i / 2)) for i in range(1, g + 1)]
    out = []
    for a in itertools.product(*[range(-b, b + 1) for b in bounds]):
        c = (1,) + a
        P = list(c) + [2 ** (g - j) * c[j] for j in range(g - 1, -1, -1)]
        r = np.roots(P)
        if np.all(np.abs(np.abs(r) - 2 ** 0.5) < 1e-6):
            out.append(tuple(P))
    return out

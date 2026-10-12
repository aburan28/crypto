#!/usr/bin/env python3
"""End-to-end control: recover a planted logarithm from four-point relations.

A negative sweep is evidence only if the same pipeline finds a logarithm
when one is findable.  This module plants `d`, sets `Q = [d]P` on a small
analogue of the challenge curve, and recovers `d` from four-point relations
alone -- three support points and one known target `[a]P + [b]Q` per
relation.

**No ground truth enters the solve.**  The only structural constant used is
the Frobenius scalar `s`, and it is obtained from the characteristic
equation `sigma^2 + sigma + 2 = 0` of the Koblitz curve, with the correct
root of the two selected by testing `sigma(P) = [s]P` as a point identity.
No discrete logarithm is computed anywhere.

**Why the unknowns are orbits, not points.**  `log(sigma^k R) = s^k log(R)`,
so the `m` points of a `sigma`-orbit share one unknown.  A support of `T`
orbits gives a system in `T + 1` unknowns -- the orbit logarithms and `d` --
whatever its point count.  That is the same collapse that leaves the
ECC2K-130 support carrying 28 unknowns rather than 3668, exercised here
where the system can actually be closed.
"""

from __future__ import annotations

import random
from itertools import product


def _sqrt_mod(a: int, p: int):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, rr = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, rr = i, b * b % p, t * b % b % p if False else t * b * b % p, rr * b % p
    return rr


def frobenius_scalar(E, P, r: int) -> int:
    """`s` with `sigma(R) = [s]R` on the order-`r` subgroup.

    `#E(F_2) = 4` gives base-field trace `t = -1`, so `sigma^2 + sigma + 2 = 0`
    and `s = (-1 +- sqrt(-7)) / 2 mod r`.  The root is chosen by a point
    identity, not by a logarithm.
    """
    disc = _sqrt_mod(-7 % r, r)
    if disc is None:
        raise RuntimeError("-7 is not a square mod r; no rational Frobenius scalar")
    inv2 = pow(2, r - 2, r)
    for s in ((-1 + disc) * inv2 % r, (-1 - disc) * inv2 % r):
        if E.mul(P, s) == E.frobenius(P):
            return s
    raise RuntimeError("neither root of the characteristic equation acts as sigma")


def solve_mod_p(rows, rhs, p: int):
    """Gaussian elimination of `rows . x = rhs` over `F_p`, `p` prime.

    Returns a solution, or `None` when the system is underdetermined.
    """
    n = len(rows[0])
    M = [list(row) + [b % p] for row, b in zip(rows, rhs)]
    where = [-1] * n
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, len(M)) if M[i][c] % p), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [v * inv % p for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        where[c] = r
        r += 1
        if r == len(M):
            break
    if any(w < 0 for w in where):
        return None
    for i in range(r, len(M)):
        if all(v % p == 0 for v in M[i][:n]) and M[i][n] % p:
            return None                      # inconsistent
    return [M[where[c]][n] % p for c in range(n)]


def target_relations(E, support, reps, target, four_torsion, limit):
    """Relations `target + eps1 R + eps2 R' + eps3 R'' in E[4]`.

    Returns `(indices, signs)` per relation, with the three support indices
    distinct.  Found by tabulating pair sums and probing each `E[4]` target,
    the same shape as `fourpoint.exhaustive_four_point`.
    """
    B = len(reps)
    table: dict = {}
    for k in range(B):
        for l in range(k + 1, B):
            for sk, sl in product((1, -1), repeat=2):
                S = E.add(reps[k] if sk > 0 else E.neg(reps[k]),
                          reps[l] if sl > 0 else E.neg(reps[l]))
                table.setdefault(S, []).append((k, l, sk, sl))

    out = []
    for i in range(B):
        for si in (1, -1):
            Ri = reps[i] if si > 0 else E.neg(reps[i])
            head = E.add(target, Ri)
            for T in four_torsion:
                want = E.add(T, E.neg(head)) if head is not None else T
                for (k, l, sk, sl) in table.get(want, ()):
                    if k == i or l == i:
                        continue
                    out.append(((i, k, l), (si, sk, sl)))
                    if len(out) >= limit:
                        return out
    return out


def recover_planted(m: int, *, seed: int = 0, extra_targets: int = 40,
                    log=print):
    """Plant `d`, recover it from four-point relations, and verify `[d]P = Q`."""
    from smallcurve import small_curve, four_torsion
    from normalbasis import find_normal_elements, NormalSupport

    F, E, order = small_curve(m)
    e4 = four_torsion(E, order)
    r = order
    while r % 2 == 0:
        r //= 2

    rng = random.Random(seed)

    # A point of order r.
    P = None
    while P is None:
        x = rng.getrandbits(m) & F.mask
        pts = E.points_over(x)
        if pts:
            cand = E.mul(pts[0], order // r)
            if cand is not None:
                P = cand
    assert E.mul(P, r) is None

    d = rng.randrange(2, r)
    Q = E.mul(P, d)
    s = frobenius_scalar(E, P, r)

    # The largest support available from a handful of normal elements.
    best = None
    for a in find_normal_elements(F, 30, rng):
        sup = NormalSupport(F, E, a)
        if best is None or sup.size > best.size:
            best = sup
    support = best
    reps = support.orbit_points(E)
    T = support.orbit_count
    log(f"m={m} r={r} B={support.size} orbits={T} unknowns={T + 1} "
        f"planted d={d}")

    # Relations against random known targets [a]P + [b]Q.
    rows, rhs = [], []
    for _ in range(extra_targets):
        a, b = rng.randrange(r), rng.randrange(1, r)
        tgt = E.add(E.mul(P, a), E.mul(Q, b))
        if tgt is None:
            continue
        for (idx, signs) in target_relations(E, support, reps, tgt, e4, limit=4):
            row = [0] * (T + 1)
            for i, sg in zip(idx, signs):
                t, k = divmod(i, support.m)
                row[t] = (row[t] + sg * pow(s, k, r)) % r
            row[T] = 4 * b % r                     # the coefficient of d
            rows.append(row)
            rhs.append((-4 * a) % r)
        if len(rows) >= 4 * (T + 1):
            break

    log(f"  collected {len(rows)} relations for {T + 1} unknowns")
    if len(rows) < T + 1:
        return {"m": m, "recovered": False, "reason": "too few relations",
                "relations": len(rows), "unknowns": T + 1}

    sol = solve_mod_p(rows, rhs, r)
    if sol is None:
        return {"m": m, "recovered": False, "reason": "system underdetermined",
                "relations": len(rows), "unknowns": T + 1}
    d_hat = sol[T]
    ok = E.mul(P, d_hat) == Q
    log(f"  recovered d={d_hat}  [d]P == Q: {ok}")
    return {"m": m, "support_size": support.size, "orbits": T,
            "unknowns": T + 1, "relations": len(rows),
            "planted_d": d, "recovered_d": d_hat,
            "recovered": bool(ok), "verified_by_point_identity": bool(ok)}

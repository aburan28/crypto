#!/usr/bin/env python3
"""B4's F0 candidates, checked properly: Koblitz K_a / GF(2^n), 127 <= n <= 190,
with a prime r | #K_a(GF(2^n)), 2^25 <= r < 2^52, such that the signed Frobenius
group of order 2n acts freely on the subgroup of order r:

  - pi acts on that subgroup as lambda mod r, a root of x^2 - mu x + 2, with
    lambda^n = 1;
  - the action is free when lambda has order exactly n and, for even n,
    lambda^(n/2) != -1 (else pi^(n/2) = -1 and every orbit halves).

Every degree's order is searched for any such factor, not only the largest,
by trial division and then ECM under a time budget per order.

    python3 instances.py instances.json 2.5     # sympy 1.14.0, Python 3.11.15

ECM under a time budget is not deterministic, so `instances.json` is the
record of what this search found, on 2026-10-01. B4 uses six of its rows,
the prime degrees named in the design (§4); `../../conformance/v2-b4/
make_cases.py` checks each of them again exactly, with the standard
library only.
"""
import json
import math
import multiprocessing as mp
import sys

import sympy

LO, HI = 25, 52


def trace(n, a):
    mu = 1 if a == 1 else -1
    t0, t1 = 2, mu
    for _ in range(n - 1):
        t0, t1 = t1, mu * t1 - 2 * t0
    return t1


def order(n, a):
    return (1 << n) + 1 - trace(n, a)


def small_factors(N, q):
    """Prime factors of N below 2^HI that sympy finds; put on queue q."""
    found = set()
    f = sympy.factorint(N, limit=1 << 24)
    for p in f:
        if sympy.isprime(p) and p < (1 << HI):
            found.add(p)
    rest = [p for p in f if not sympy.isprime(p)]
    for c in rest:
        while True:
            ds = sympy.ntheory.ecm(c, B1=200000, B2=20000000, max_curve=200)
            for d in ds:
                if sympy.isprime(d) and d < (1 << HI):
                    found.add(d)
            q.put(sorted(found))
            if all(sympy.isprime(d) for d in ds):
                break
            c = max(d for d in ds if not sympy.isprime(d))
    q.put(sorted(found))


def factors_within(N, seconds):
    q = mp.Queue()
    proc = mp.Process(target=small_factors, args=(N, q))
    proc.start()
    proc.join(seconds)
    if proc.is_alive():
        proc.terminate()
    best = []
    while not q.empty():
        best = q.get()
    return best


def eigenvalue(n, a, r):
    """lambda mod r with lambda^2 - mu lambda + 2 = 0 and lambda^n = 1, or None."""
    mu = 1 if a == 1 else -1
    disc = (mu * mu - 8) % r
    for s in sympy.sqrt_mod(disc, r, all_roots=True) or []:
        lam = (mu + s) * pow(2, -1, r) % r
        if pow(lam, n, r) == 1:
            return lam
    return None


def free(n, a, r):
    lam = eigenvalue(n, a, r)
    if lam is None:
        return False, None
    order_lam = sympy.n_order(lam, r)
    if order_lam != n:
        return False, lam
    if n % 2 == 0 and pow(lam, n // 2, r) == r - 1:
        return False, lam
    return True, lam


def main():
    out_path, seconds = sys.argv[1], float(sys.argv[2])
    rows = []
    for n in range(127, 191):
        for a in (0, 1):
            N = order(n, a)
            candidates = [p for p in factors_within(N, seconds) if p >= (1 << LO)]
            for r in sorted(candidates):
                ok, lam = free(n, a, r)
                if not ok:
                    continue
                row = {"n": n, "a": a, "trace": trace(n, a), "order": str(N), "r": str(r),
                       "log2_r": round(math.log2(r), 2), "cofactor": str(N // r),
                       "r_squared_divides": N % (r * r) == 0,
                       "eigenvalue": str(lam),
                       "proper_intermediate_subfields": [d for d in sympy.divisors(n) if 1 < d < n],
                       "prime_degree": sympy.isprime(n)}
                rows.append(row)
                print(json.dumps({k: row[k] for k in ("n", "a", "log2_r", "r", "prime_degree")}), flush=True)
    json.dump(rows, open(out_path, "w"), indent=1)


if __name__ == "__main__":
    main()

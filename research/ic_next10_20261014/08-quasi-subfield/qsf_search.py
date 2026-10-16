#!/usr/bin/env python3
"""Arms A1 (Frobenius-stable linearised QSPs) and A3 (multiplicative QSPs).

A1: an F_2-subspace V of F_{2^n} is Frobenius-stable iff V = ker g(sigma) for a
divisor g of t^n - 1 over F_2; its subspace polynomial is sum g_i X^{2^i}, so
n0 = deg g and l = log2 deg lambda = second-highest exponent of g.  Every divisor
is enumerated (products of the irreducible factors, which PARI/GP supplies).

A3: L = X^{2^{n'}} - X^a = X^a (X^{2^{n'}-a} - 1) has D+1 distinct roots in
F_{2^n} when D = 2^{n'} - a divides 2^n - 1.  Every divisor D of 2^n - 1 is
enumerated (factorisation from PARI/GP).

beta = l * n / n'^2 (Euler-Petit); HKPYY generic-beating threshold beta < 1/(2*4.876).

Usage: python3 -I qsf_search.py <outdir>
"""
import csv
import json
import math
import subprocess
import sys
from itertools import product

KAPPA = 4.876
BETA_STAR = 1 / (2 * KAPPA)


def gp(cmd):
    out = subprocess.run(["gp", "-q", "-f"], input=cmd, capture_output=True, text=True, check=True)
    return out.stdout.strip()


def primes_upto(n):
    return [p for p in range(2, n + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]


def f2_factors(n):
    """Irreducible factors of t^n - 1 over F_2 (n odd), as int bitmasks."""
    s = gp(f"f=factormod(x^{n}-1,2)[,1]; for(i=1,#f, print(Vec(lift(f[i]))))")
    facs = []
    for line in s.splitlines():
        coeffs = json.loads(line)  # leading first
        mask = 0
        for c in coeffs:
            mask = (mask << 1) | (int(c) & 1)
        facs.append(mask)
    return facs


def pmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def stable_frontier(n):
    """min l per n0 over all divisors of t^n - 1 (n odd)."""
    facs = f2_factors(n)
    prods = [1]
    for f in facs:
        prods = prods + [pmul(p, f) for p in prods]
    best = {}
    for g in prods:
        d = g.bit_length() - 1
        if d <= 0 or d >= n:
            continue
        low = g ^ (1 << d)
        ell = low.bit_length() - 1  # second-highest exponent
        if d not in best or ell < best[d][0]:
            best[d] = (ell, g)
    return facs, best


def divisors_of_mersenne(n):
    s = gp(f"print(divisors(2^{n}-1))")
    return [int(x) for x in s.strip("[]").split(",")]


def mult_rows(n):
    rows = []
    for D in divisors_of_mersenne(n):
        if D < 2:
            continue
        base = max(1, math.ceil(math.log2(D)))
        for nprime in (base, base + 1):
            a = (1 << nprime) - D
            if a < 1 or nprime >= n:
                continue  # n' >= n is not a proper factor base (degenerate)
            frac = (D + 1) / (1 << nprime)
            if frac < 0.5:
                continue
            ell = math.log2(a)
            beta = ell * n / nprime ** 2
            rows.append(dict(n=n, nprime=nprime, D=D, a=a, root_fraction=frac,
                             log2_deg_lambda=ell, ell_over_nprime=ell / nprime, beta=beta,
                             subfield=(a == 1)))
    return rows


def main():
    out = sys.argv[1]
    primes = [p for p in primes_upto(127) if p >= 3]
    # ---------- A1 ----------
    a1 = []
    for n in primes:
        facs, best = stable_frontier(n)
        for n0 in sorted(best):
            ell, g = best[n0]
            beta = ell * n / n0 ** 2
            a1.append(dict(n=n, n0=n0, min_ell=ell, beta=beta, ell_over_n0=ell / n0,
                           g_hex=hex(g), qsp=beta < 1, in_scope=n <= 61,
                           trivial=(n0 in (1, n - 1))))
    # positive control at composite odd n is not needed; n=7 and n=31 are prime.
    with open(f"{out}/a1_stable_frontier.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(a1[0].keys()))
        w.writeheader()
        w.writerows(a1)
    # ---------- A3 ----------
    a3 = []
    for n in [p for p in primes if p <= 61]:
        a3 += mult_rows(n)
    controls = []
    for n in (4, 6, 8, 10, 12):
        controls += mult_rows(n)
    with open(f"{out}/a3_multiplicative.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(a3[0].keys()))
        w.writeheader()
        w.writerows(a3 + controls)

    # ---------- summary ----------
    summ = {"beta_star": BETA_STAR, "kappa": KAPPA}
    # A1 positive controls
    def a1get(n, n0):
        return next((r for r in a1 if r["n"] == n and r["n0"] == n0), None)
    summ["a1_controls"] = {f"n{n}_n0{n0}": a1get(n, n0) for n, n0 in [(7, 3), (7, 4), (31, 15), (31, 16), (73, 9)]}
    nontriv = [r for r in a1 if r["in_scope"] and not r["trivial"] and r["min_ell"] >= 1]
    summ["a1_best_nontrivial_in_scope"] = min(nontriv, key=lambda r: r["beta"])
    summ["a1_qsp_nontrivial_in_scope"] = [r for r in nontriv if r["qsp"]]
    summ["a1_best_nontrivial_all"] = min([r for r in a1 if not r["trivial"] and r["min_ell"] >= 1], key=lambda r: r["beta"])
    summ["a1_below_threshold"] = [r for r in a1 if r["min_ell"] >= 1 and r["beta"] < BETA_STAR]
    nonsub = [r for r in a3 if not r["subfield"]]
    summ["a3_best_prime"] = min(nonsub, key=lambda r: r["beta"]) if nonsub else None
    summ["a3_qsp_prime"] = sorted([r for r in nonsub if r["beta"] < 1], key=lambda r: r["beta"])
    summ["a3_below_threshold_prime"] = [r for r in nonsub if r["beta"] < BETA_STAR]
    summ["a3_rows_prime"] = len(a3)
    def ctl(n, npr, a):
        return next((r for r in controls if r["n"] == n and r["nprime"] == npr and r["a"] == a), None)
    summ["a3_controls"] = {"n4_X8-X3": ctl(4, 3, 3), "n6_np5_a11": ctl(6, 5, 11),
                           "n8_np7_a43": ctl(8, 7, 43), "n8_np6_a13": ctl(8, 6, 13)}
    with open(f"{out}/search_summary.json", "w") as fh:
        json.dump(summ, fh, indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k != "a3_qsp_prime"}, indent=1))
    print("A3 QSP rows at prime n:", len(summ["a3_qsp_prime"]))
    for r in summ["a3_qsp_prime"][:20]:
        print(r)


if __name__ == "__main__":
    main()

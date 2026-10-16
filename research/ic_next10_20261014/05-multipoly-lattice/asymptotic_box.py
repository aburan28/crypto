"""Asymptotic (ℓ → ∞) reach of the per-group box family for shared-target Semaev systems.

System: k targets, each S_{m+1}(x_S, x_{A_i}, z_i) with |S| = s shared summands and |A_i| = q = m − s
private ones, all |x| < p^δ, degree d = 2^{m−1} per variable, top-corner leading monomials.
Scaled coordinates u = β/(dℓ). Gain fraction g(u) = min(1, min_S u, Σ_i min_{A_i} u).
Box R = [0,Λ_S]^s × [0,Λ_A]^{kq}. Condition ∫_R g > δ d ∫_R Σu, so
δ_box(Λ) = E[g] / (d · E[Σu]),  E[Σu] = (sΛ_S + kqΛ_A)/2.
E[g] = ∫_0^1 P(min_S > t) · P(Σ_i min_{A_i} > t) dt  (independence).

Also prints the relation economics: n = s + kq small unknowns, r = k equations; IC beats rho iff
(r − 1/2)/(n − 1) < δ < 1/4 (cost exponent max(r − (n−1)δ, 2δ) < 1/2).

Usage: python3 -I asymptotic_box.py out.json
"""
import json
import sys
from fractions import Fraction

import numpy as np
from scipy.optimize import minimize
from scipy.signal import fftconvolve

N = 20_000  # grid points (search); the reported value is re-evaluated at N_FINE
N_FINE = 200_000


def set_grid(n):
    global N
    N = n


def surv_min(q, lam, t):
    """P(min of q iid U[0,lam] > t)."""
    if q == 0:
        return np.ones_like(t)
    x = np.clip(1.0 - t / lam, 0.0, 1.0)
    return x ** q


def expected_g(s, q, k, lam_s, lam_a):
    """E[min(1, min_S, Σ_i min_{A_i})] = ∫_0^tmax P(min_S > t) P(Σ min_{A_i} > t) dt."""
    tmax = min(1.0, lam_s, k * lam_a if q else 1.0)
    T = np.linspace(0.0, tmax, N + 1)
    H = T[1] - T[0]
    a = surv_min(s, lam_s, T)
    if q == 0 or k == 0:
        b = np.ones_like(T)
    else:
        cdf = 1.0 - surv_min(q, lam_a, T)
        mass = np.diff(cdf)  # cell j ~ value (j + 1/2) H
        tot = mass.copy()
        for _ in range(k - 1):
            tot = fftconvolve(tot, mass)[:N]
        c = np.cumsum(tot)
        idx = np.floor(T / H - k / 2.0 + 1e-9).astype(np.int64)
        b = 1.0 - np.where(idx >= 0, c[np.clip(idx, 0, N - 1)], 0.0)
    f = a * b
    return float(np.sum((f[1:] + f[:-1]) * 0.5) * H)


def delta_box(s, q, k, d, lam_s, lam_a):
    eg = expected_g(s, q, k, lam_s, lam_a)
    es = (s * lam_s + k * q * lam_a) / 2.0
    return eg / (d * es)


def best(s, q, k, d):
    if q == 0 or k == 0:
        f = lambda x: -delta_box(s, 0, 0, d, abs(x[0]) + 1e-6, 1.0)
        r = minimize(f, [0.7], method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-12})
        return -r.fun, (abs(float(r.x[0])) + 1e-6, None)
    bestv, bestx = -1, None
    for x0 in ([1.0, 0.5], [0.8, 0.8]):
        f = lambda x: -delta_box(s, q, k, d, abs(x[0]) + 1e-6, abs(x[1]) + 1e-6)
        r = minimize(f, x0, method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-12, "maxiter": 4000})
        if -r.fun > bestv:
            bestv, bestx = -r.fun, (abs(float(r.x[0])) + 1e-6, abs(float(r.x[1])) + 1e-6)
    return bestv, bestx


def economics(n, r):
    lo = Fraction(2 * r - 1, 2 * (n - 1)) if n > 1 else None
    beats = lo is not None and lo < Fraction(1, 4)
    best_delta = Fraction(r, n + 1)
    best_exp = Fraction(2 * r, n + 1)
    return lo, beats, best_delta, best_exp


def main():
    out = []
    rows = []
    # (label, m, s, k): k = number of targets, s = shared summands
    cases = []
    for m in (2, 3, 4):
        cases.append((f"single S{m+1}", m, m, 1))
        for s in range(1, m):
            for k in (2, 3):
                cases.append((f"S{m+1} x{k} targets, s={s}", m, s, k))
    print("| system | m | s | k | n | r | d | δ*_box (ℓ→∞) | Λ_S, Λ_A | vs single | ceiling 1/(md) | needed δ > | beats rho possible? | best exponent |")
    print("|:--|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|:--|--:|")
    single = {}
    for label, m, s, k in cases:
        d = 2 ** (m - 1)
        q = m - s
        set_grid(20_000)
        if s == m:
            v, x = best(m, 0, 0, d)
            set_grid(N_FINE)
            v = delta_box(m, 0, 0, d, x[0], 1.0)
            single[m] = v
            n, r = m, 1
        else:
            v, x = best(s, q, k, d)
            set_grid(N_FINE)
            v = delta_box(s, q, k, d, x[0], x[1])
            n, r = s + k * q, k
        lo, beats, bd, be = economics(n, r)
        ratio = v / single[m]
        ceil = 1.0 / (m * d)
        print(
            f"| {label} | {m} | {s} | {k} | {n} | {r} | {d} | {v:.6f} | {x[0]:.3f}, {x[1] if x[1] is None else round(x[1], 3)} | {ratio:.4f} | {ceil:.6f} | {float(lo):.4f} | {'yes' if beats else 'no'} (needs δ*≥{float(lo):.4f}) | {be} at δ={bd} |"
        )
        out.append({"system": label, "m": m, "s": s, "k": k, "n": n, "r": r, "d": d, "delta_box": v,
                    "lam_s": x[0], "lam_a": x[1], "ratio_vs_single": ratio, "ceiling": ceil,
                    "needed_lo": float(lo), "needed_lo_exact": str(lo), "window_nonempty": beats,
                    "best_exponent": str(be)})
    # closed-form check for the single polynomial
    for m in (2, 3, 4):
        exact = 2 / (m * (m + 1) * 2 ** (m - 1))
        print(f"check single m={m}: numeric {single[m]:.8f} vs 2/(m(m+1)d) = {exact:.8f} (diff {abs(single[m]-exact):.1e})")
        out.append({"check": f"single m={m}", "numeric": single[m], "closed_form": exact})
    json.dump(out, open(sys.argv[1], "w"), indent=1)


if __name__ == "__main__":
    main()

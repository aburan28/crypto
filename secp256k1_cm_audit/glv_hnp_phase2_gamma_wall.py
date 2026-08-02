"""
Thread 23b — calibrate the GLV-HNP Phase-2 K1 wall against the Gaussian
heuristic, and test nu_hat as a predictor of the *wall location*.

Context.  EXP P5 in glv_hnp_phase2_projected.py measured, on 8 fresh 13-bit
curves, the largest K1 at which the Phase-2 attack still recovers d 5/5, and
converted it to

    gamma_wall = eff * n^(1/m) / (3/(2*pi*e)),   eff = K1*K2/n

the Gaussian-heuristic margin at the wall (gamma == 1 <=> ||v_planted|| == GH
of the BDD lattice).  It found spearman(nu_hat, gamma_wall) = -0.762 (n=8)
versus +0.190 for lam*.  n=8 is too small to conclude, so this script
  (a) enlarges the sample (up to 24 curves per bit-size),
  (b) replicates at a second bit-size,
  (c) adds a permutation p-value and a log-linear fit,
  (d) reports the same correlation for the 8 previously-falsified invariants'
      closest available stand-in (lam*) as the negative control.

nu_hat = lambda_1(L2)/sqrt(det L2) for L2 = <(n*S_K1,0), (-lam*S_K1, S_K2)>,
the 2026-07-29 separator (AUC 0.935 against the June C1/C2 classes).

Run: python3 glv_hnp_phase2_gamma_wall.py
"""

import math
import random

import sympy

from glv_hnp_phase2_projected import (
    ec_mul, tonelli_shanks, find_generator, j0_traces, glv_roots,
    lam_star, nu_hat, scales, gh_margin, sweep, build_curve,
)

GH_CONST = 3.0 / (2.0 * math.pi * math.e)      # 0.175664...


def eisenstein_decompose_fast(p):
    """p = a^2 - a*b + b^2 for p = 1 mod 3, in O(sqrt(p)).

    Completing the square: (2b-a)^2 = 4p - 3a^2, so a runs only to
    sqrt(4p/3) and b is read off when 4p-3a^2 is a perfect square."""
    if p % 3 != 1:
        return None
    a_max = math.isqrt(4 * p // 3) + 1
    for a in range(0, a_max + 1):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        r = math.isqrt(disc)
        if r * r != disc:
            continue
        for num in (a + r, a - r):
            if num % 2 == 0:
                b = num // 2
                if a * a - a * b + b * b == p:
                    return (a, b)
    return None


def search_curves_fast(lo, hi, per_bin=3, nbins=8, max_primes=6000):
    """j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3, bucketed by
    lam* into nbins bins over [0, 0.5]."""
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    seen = 0
    while p < hi and seen < max_primes:
        seen += 1
        if p % 3 == 1:
            eis = eisenstein_decompose_fast(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1:
                        continue
                    if not sympy.isprime(n_cand):
                        continue
                    roots = glv_roots(n_cand)
                    if roots is None:
                        continue
                    lam = roots[0]
                    idx = min(nbins - 1,
                              int(lam_star(lam, n_cand) / (0.5 / nbins)))
                    if len(bins[idx]) >= per_bin:
                        continue
                    cur = build_curve(p, n_cand)
                    if cur is None:
                        continue
                    bins[idx].append((p, cur[1], n_cand, lam, cur[4]))
        if all(len(v) >= per_bin for v in bins.values()):
            break
        p = int(sympy.nextprime(p))
    out = []
    for i in range(nbins):
        out.extend(bins[i])
    return out


def find_wall(curve, m, seeds, k1_max=48):
    """Largest K1 with full recovery over `seeds` (method KAN).  Scans upward
    and stops one step after the first failure, so a non-monotone blip does
    not extend the wall arbitrarily."""
    k1_wall = 0
    for k1 in range(2, k1_max + 1):
        agg, used, _ = sweep(curve, m, k1, seeds, methods=("KAN",))
        if used and agg["KAN"] == used:
            k1_wall = k1
        elif k1_wall:
            break
    return k1_wall


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def _rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r

def spearman(x, y):
    rx, ry = _rank(x), _rank(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx)
                    * sum((b - my) ** 2 for b in ry))
    return num / den if den else float('nan')

def perm_p(x, y, trials=20000, seed=7):
    """Two-sided permutation p-value for |spearman|."""
    obs = abs(spearman(x, y))
    rng = random.Random(seed)
    ys = list(y)
    hits = 0
    for _ in range(trials):
        rng.shuffle(ys)
        if abs(spearman(x, ys)) >= obs - 1e-12:
            hits += 1
    return (hits + 1) / (trials + 1)

def linfit(x, y):
    """Least squares y = a + b*x; returns (a, b, r2)."""
    nn = len(x)
    mx, my = sum(x) / nn, sum(y) / nn
    sxx = sum((a - mx) ** 2 for a in x)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    b = sxy / sxx if sxx else float('nan')
    a = my - b * mx
    ss_tot = sum((c - my) ** 2 for c in y)
    ss_res = sum((c - (a + b * d)) ** 2 for c, d in zip(y, x))
    r2 = 1 - ss_res / ss_tot if ss_tot else float('nan')
    return a, b, r2


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 23b — gamma_wall calibration and nu_hat as a wall predictor")
    print("=" * 78)
    print(f"GH constant 3/(2*pi*e) = {GH_CONST:.6f}")
    print("gamma = eff * n^(1/m) / GH_CONST ;  gamma<1 <=> ||v_planted|| < GH(L0)")

    SEEDS = [42, 1234, 9999, 555, 31337]
    M_SIGS = 10

    all_rows = []
    for lo, hi, per_bin, tag in ((2 ** 13, 2 ** 14, 3, "13-bit"),
                                 (2 ** 15, 2 ** 16, 3, "15-bit")):
        print("\n" + "-" * 78)
        print(f"{tag}: curve search")
        print("-" * 78)
        curves = search_curves_fast(lo, hi, per_bin=per_bin, nbins=8)
        print(f"  {len(curves)} curves found")
        print(f"  {'p':>7}{'n':>7}{'K2':>5}{'lam*':>7}{'nu@2':>7}"
              f"{'K1wall':>8}{'eff_w':>8}{'g_wall':>8}")
        rows = []
        for cur in curves:
            p, b, n, lam, G = cur
            k2 = math.isqrt(n) + 1
            k1_wall = find_wall(cur, M_SIGS, SEEDS)
            if not k1_wall:
                print(f"  {p:>7}{n:>7}{k2:>5}{lam_star(lam,n):>7.3f}"
                      f"{'-':>7}{'none':>8}")
                continue
            g_w = gh_margin(n, M_SIGS, k1_wall, k2)
            s = scales(n, 2, k2)
            nu2 = nu_hat(n, lam, s[0], s[2])
            rows.append({'tag': tag, 'p': p, 'n': n, 'lam_star': lam_star(lam, n),
                         'nu2': nu2, 'k1w': k1_wall, 'eff': k1_wall * k2 / n,
                         'gw': g_w})
            print(f"  {p:>7}{n:>7}{k2:>5}{rows[-1]['lam_star']:>7.3f}"
                  f"{nu2:>7.3f}{k1_wall:>8}{rows[-1]['eff']:>8.4f}{g_w:>8.2f}")
        all_rows.extend(rows)

        if len(rows) >= 5:
            gws = [r['gw'] for r in rows]
            nus = [r['nu2'] for r in rows]
            lss = [r['lam_star'] for r in rows]
            print(f"\n  {tag}: gamma_wall min={min(gws):.2f} "
                  f"median={sorted(gws)[len(gws)//2]:.2f} max={max(gws):.2f}")
            r_nu = spearman(nus, gws)
            r_ls = spearman(lss, gws)
            print(f"  spearman(nu_hat, gamma_wall) = {r_nu:+.3f}  "
                  f"perm p = {perm_p(nus, gws):.4f}   [predictor]")
            print(f"  spearman(lam*,   gamma_wall) = {r_ls:+.3f}  "
                  f"perm p = {perm_p(lss, gws):.4f}   [negative control]")

    # -----------------------------------------------------------------
    print("\n" + "-" * 78)
    print("POOLED (both bit-sizes)")
    print("-" * 78)
    if len(all_rows) >= 8:
        gws = [r['gw'] for r in all_rows]
        nus = [r['nu2'] for r in all_rows]
        lss = [r['lam_star'] for r in all_rows]
        effs = [r['eff'] for r in all_rows]
        print(f"  N = {len(all_rows)} curves")
        print(f"  gamma_wall : min={min(gws):.2f} "
              f"median={sorted(gws)[len(gws)//2]:.2f} max={max(gws):.2f}")
        print(f"  eff_wall   : min={min(effs):.4f} "
              f"median={sorted(effs)[len(effs)//2]:.4f} max={max(effs):.4f}")
        print(f"  spearman(nu_hat, gamma_wall) = {spearman(nus, gws):+.3f}  "
              f"perm p = {perm_p(nus, gws):.4f}")
        print(f"  spearman(lam*,   gamma_wall) = {spearman(lss, gws):+.3f}  "
              f"perm p = {perm_p(lss, gws):.4f}")
        a, b, r2 = linfit(nus, [math.log(g) for g in gws])
        print(f"  fit: log(gamma_wall) = {a:+.3f} {b:+.3f}*nu_hat   R^2={r2:.3f}")
        print(f"       => gamma_wall ~ {math.exp(a):.2f} * exp({b:.3f}*nu_hat)")
        # residual scatter of the corrected quantity
        corr = [math.log(g) - (a + b * nu) for g, nu in zip(gws, nus)]
        sd0 = math.sqrt(sum((math.log(g) - sum(math.log(z) for z in gws)/len(gws))**2
                            for g in gws) / len(gws))
        sd1 = math.sqrt(sum(c * c for c in corr) / len(corr))
        print(f"  sd log(gamma_wall) = {sd0:.3f}  ->  after nu_hat correction "
              f"{sd1:.3f}  ({100*(1-sd1/sd0):.0f}% reduction)")

    print("\nDone.")

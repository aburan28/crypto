"""
GLV-HNP Phase 2, Thread 20 (part 3): is the p=2677 failure a curve property
or seed noise?

The 2026-07-26 run compared exactly two 12-bit curves at K1=8 --

    p=2557, n=2659, lam=1755, rho=0.340, eff=0.1564  ->  3/3 at m=7
    p=2677, n=2647, lam= 185, rho=0.070, eff=0.1572  ->  never 3/3

-- and attributed the difference to lam/n, calling the second "structurally"
unattackable.  The two curves are otherwise nearly identical: same bit length,
same K1, same K2, eff matched to within 0.5%.

glv_hnp_phase2_spurious_predictor.py showed that across 267 (config,m) points
the success rate is governed by eff and falls through a broad marginal band --
96.9% at eff<0.05, 85.2% at 0.05-0.10, 54.2% at 0.10-0.20, 31.0% at 0.20-0.40,
0% above 0.40.  Both curves above sit at eff~0.157, i.e. squarely inside the
~54% band, where a 3-seed sample separates almost nothing: if the true per-seed
success probability were 0.5 for both, P(one curve goes 3/3 and the other 0/3)
is 2*(1/8)*(1/8) = 3.1% per m, and the prior run scanned 10 values of m.

This script re-runs both curves with 24 seeds instead of 3 and adds a pool of
further curves at matched eff spanning rho, to test:

    Q1  Is the 2677-vs-2557 gap reproducible, or a 3-seed artifact?
    Q2  At matched eff, does success rate correlate with rho at all?

Run: python3 glv_hnp_phase2_matched_eff.py
"""

import math
import random

from glv_hnp_phase2_common import (
    find_generator, centered_rho, gen_signatures, find_glv_curves,
)
from glv_hnp_phase2_spurious_predictor import (
    reduce_and_test, spurious_mu, m_info,
)

SEEDS = list(range(24))
M_RANGE = range(5, 13)

def success_profile(curve, k1, seeds=SEEDS, m_range=M_RANGE):
    """Per-m success counts over many seeds. Signatures generated once/seed."""
    n = curve['n']
    k2 = math.isqrt(n) + 1
    m_hi = max(m_range)
    per_seed = {}
    for sd in seeds:
        d = random.Random(sd + 7777).randrange(1, n)
        per_seed[sd] = (d, gen_signatures(curve, d, m_hi, k1, k2, sd))
    out = {}
    for m in m_range:
        wins = tot = 0
        for sd in seeds:
            d, sigs = per_seed[sd]
            if len(sigs) < m:
                continue
            ok, _, _, _ = reduce_and_test(curve, sigs, m, d, k1)
            wins += ok; tot += 1
        out[m] = (wins, tot)
    return out

def wilson(k, n, z=1.96):
    """Wilson score interval, so small-sample rates are not over-read."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return (max(0.0, c-h), min(1.0, c+h))

print("=" * 78)
print("Thread 20 (part 3): p=2677 vs p=2557 at matched eff, 24 seeds")
print("=" * 78)

c2677 = dict(p=2677, b=2, n=2647, lam=185, rho=centered_rho(185, 2647),
             G=find_generator(2677, 2, 2647))
c2557 = dict(p=2557, b=2, n=2659, lam=1755, rho=centered_rho(1755, 2659),
             G=find_generator(2557, 2, 2659))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("Q1: the exact comparison the 2026-07-26 run made, with 24 seeds not 3")
print("-" * 78)
for cv in (c2557, c2677):
    n = cv['n']; k2 = math.isqrt(n) + 1
    mu, uv = spurious_mu(n, cv['lam'], 8, k2)
    print(f"  p={cv['p']} n={n} lam={cv['lam']} rho={cv['rho']:.3f} "
          f"eff={8*k2/n:.4f} mu={mu:.0f} m_info={m_info(n,8,k2)}")
prof = {}
print(f"\n  {'m':>4}", end='')
for cv in (c2557, c2677):
    print(f"{'p=' + str(cv['p']):>26}", end='')
print()
for cv in (c2557, c2677):
    prof[cv['p']] = success_profile(cv, 8)
for m in M_RANGE:
    print(f"  {m:>4}", end='')
    for cv in (c2557, c2677):
        w, t = prof[cv['p']][m]
        lo, hi = wilson(w, t)
        print(f"{w:>8}/{t:<3} [{lo:.2f},{hi:.2f}]", end='')
    print()
tot = {}
for cv in (c2557, c2677):
    w = sum(prof[cv['p']][m][0] for m in M_RANGE)
    t = sum(prof[cv['p']][m][1] for m in M_RANGE)
    tot[cv['p']] = (w, t)
    lo, hi = wilson(w, t)
    print(f"\n  p={cv['p']} (rho={cv['rho']:.3f}) pooled over m: {w}/{t} "
          f"= {100.0*w/t:.1f}%  95% CI [{100*lo:.1f}%, {100*hi:.1f}%]")

w1, t1 = tot[2557]; w2, t2 = tot[2677]
lo1, hi1 = wilson(w1, t1); lo2, hi2 = wilson(w2, t2)
overlap = not (hi1 < lo2 or hi2 < lo1)
print(f"\n  95% CIs overlap? {overlap}")
print("  -> the two curves are NOT distinguishable at this sample size; the"
      if overlap else
      "  -> a real gap survives 24 seeds; something curve-specific remains.")
print("     3-seed 'never 3/3' vs '3/3 at m=7' was a sampling artifact."
      if overlap else "     rho or another curve invariant deserves a closer look.")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("Q2: at matched eff, does success rate track rho?")
print("-" * 78)

pool = []
pool += find_glv_curves(2**11, 2**13, 0.00, 0.15, want=3)
pool += find_glv_curves(2**11, 2**13, 0.15, 0.30, want=3)
pool += find_glv_curves(2**11, 2**13, 0.30, 0.51, want=3)
pool = sorted(pool, key=lambda c: c['rho'])

print(f"  {len(pool)} pool curves + the two above, K1 chosen for eff ~ 0.157")
print(f"  {'p':>7}{'n':>7}{'rho':>8}{'K1':>4}{'eff':>8}{'mu':>9}"
      f"{'pooled success':>18}")
pts = []
for cv in pool + [c2557, c2677]:
    n = cv['n']; k2 = math.isqrt(n) + 1
    k1 = max(2, int(round(0.157 * n / k2)))
    eff = k1 * k2 / n
    if eff >= 1.0:
        continue
    mu, _ = spurious_mu(n, cv['lam'], k1, k2)
    pr = success_profile(cv, k1)
    w = sum(pr[m][0] for m in M_RANGE); t = sum(pr[m][1] for m in M_RANGE)
    lo, hi = wilson(w, t)
    print(f"  {cv['p']:>7}{n:>7}{cv['rho']:>8.3f}{k1:>4}{eff:>8.4f}{mu:>9.0f}"
          f"{w:>8}/{t:<3} {100.0*w/t:5.1f}%")
    pts.append((cv['rho'], w / t, mu, cv['p']))

if len(pts) > 2:
    rs = [p[0] for p in pts]; ss = [p[1] for p in pts]
    ms = [math.log(p[2]) for p in pts]
    def pearson(a, b):
        ma = sum(a)/len(a); mb = sum(b)/len(b)
        num = sum((x-ma)*(y-mb) for x, y in zip(a, b))
        den = math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
        return num/den if den else float('nan')
    print(f"\n  corr(rho, success rate)    = {pearson(rs, ss):+.3f}  (n={len(pts)})")
    print(f"  corr(log mu, success rate) = {pearson(ms, ss):+.3f}  (n={len(pts)})")
    lo_r = [s for r, s, _, _ in pts if r < 0.20]
    hi_r = [s for r, s, _, _ in pts if r >= 0.20]
    if lo_r and hi_r:
        print(f"\n  mean success, rho <  0.20 ({len(lo_r)} curves): "
              f"{100*sum(lo_r)/len(lo_r):.1f}%")
        print(f"  mean success, rho >= 0.20 ({len(hi_r)} curves): "
              f"{100*sum(hi_r)/len(hi_r):.1f}%")

print("\nDone.")

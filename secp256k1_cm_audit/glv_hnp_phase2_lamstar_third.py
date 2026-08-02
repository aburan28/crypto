"""
GLV-HNP Phase 2, Thread 23 addendum (U6): the lam* ~ 1/3 anomaly.

Thread 23's U4 swept 6 fresh 17-bit curves at m=12 over eff = K1*K2/n and found
that at eff >= 0.10 the curves split BIMODALLY, not monotonically, in
lam* = min(lam, n-lam)/n:

    lam*     eff=0.10  0.15  0.20  0.25
    0.0273     1/5     1/5   0/5   0/5
    0.1372     1/5     0/5   0/5   0/5
    0.2115     0/5     0/5   0/5   0/5
    0.3184     5/5     5/5   5/5   3/5     <-- lam* ~ 1/3
    0.3346     5/5     5/5   4/5   3/5     <-- lam* ~ 1/3
    0.4456     0/5     0/5   0/5   0/5

The two survivors both sit at lam* ~ 1/3, and lam*=0.4456 (the LARGEST lam*)
fails while lam*=0.3184 succeeds - so this is not "bigger lam* is better".
The two historical 12-bit curves agree: 2557 has lam*=0.3400 (far K1 wall,
~12-16) and 2677 has lam*=0.0699 (near K1 wall, ~4-6).  8/8 so far.

This is in tension with neither 2026-07-26 (which claimed a *monotone* lam*
threshold; still falsified - lam*=0.0068 recovers 5/5 at eff=0.05) nor
2026-07-29 T3 (which showed lam* is not a *viability* threshold; still true).
It is a new, weaker claim: at moderate eff, lam* modulates the success rate
with a PEAK at 1/3.

H23a:  |lam* - 1/3| separates success from failure at moderate eff,
       and beats lam* itself, the lambda-block quantity rho (falsified as a
       standalone predictor in 2026-07-29 T2), and distance to other small
       rationals (1/4, 1/2, 0).

Falsifier: run ~40 fresh curves with lam* spread over (0, 0.5) at three eff
levels.  If |lam* - 1/3| does not beat the majority baseline, H23a is dead and
Phase-2 predictor hunting should stop.  If it does, the next question is
mechanistic: lam ~ n/3 means |3*lam - n| is small, i.e. lam/n has continued
fraction [0; 3, large, ...], so the test also records the CF partial quotients
and delta3 = |3*lam - n|/n.

RESULTS (2026-08-02 run; 40 curves, 17-18 bit, m=12, 5 seeds, 3 eff levels):

  H23a FALSIFIED.  |lam* - 1/3| pooled AUC 0.567, best accuracy 0.758 against
  a 0.733 majority baseline - no signal.  delta3 is identical to it (0.567),
  |lam*-1/4| is 0.436, and lam* itself is 0.328.  The bimodal pattern seen on
  U4's six curves was a small-sample artifact; curves at lam* = 0.2976, 0.3424,
  0.3581, 0.3661 (all within 0.04 of 1/3) fail at every eff level, while
  lam* = 0.4449, 0.4585, 0.4816 succeed at all three.

  UNEXPECTED: rho - the lambda-block quantity from 2026-07-29 T2 - is a strong
  separator, but with the OPPOSITE SIGN to the one H20 predicted.
  SMALL rho => success:

      eff    AUC(rho)   best acc   tau      baseline
      0.10   0.945      0.925      0.6057   0.650
      0.15   0.977      0.975      0.4181   0.750
      0.20   0.926      0.925      0.3893   0.800
      pooled 0.891      0.892      0.4325   0.733

  At eff=0.15 the rule "rho <= 0.42 => success" is perfect on the positive
  side: all 9 curves with rho <= 0.418 recover, no curve with rho >= 0.52
  recovers, and the single error is one false negative (rho=0.607, 3/5).

  This is consistent with, not contrary to, 2026-07-29 T2, whose own three data
  points already had the failure curve at the LARGEST rho (0.686 vs 0.610 and
  0.399).  T2 recorded that as a falsification because H20 predicted
  "rho >= ~1 needed for success"; the sign was simply backwards.  T2's 17-bit
  sweep could not see it because it was run at eff=0.05, where 19/20 curves
  succeed and the 0.95 majority baseline is unbeatable.

  Note rho is NOT a curve-level invariant - it depends on (n, lam, K1, K2) -
  so it does not contradict 2026-07-29 T5's argument that no curve-level
  invariant can predict a BDD outcome.  It is an instance-level invariant,
  computable before any signature is seen.

  maxq (largest continued-fraction partial quotient of lam/n) carries the same
  signal weakly and in the same direction (AUC 0.220 as printed, i.e. 0.780
  under "large maxq => success"); large partial quotient => small rho.  rho
  dominates it because rho also accounts for the S_K1/S_K2 scaling.

Run: python3 glv_hnp_phase2_lamstar_third.py
"""

import math
import random

import sympy
from fpylll import IntegerMatrix, LLL

from glv_hnp_phase2_projected import (
    modinv, ec_mul, tonelli_shanks, find_generator,
    eisenstein_decompose, j0_traces, glv_roots, lam_star,
    scales, gen_signatures, build_proj, planted_proj, recover_proj,
    norm, gamma_gap, build_curve,
)

# ---------------------------------------------------------------------------
# Candidate separators
# ---------------------------------------------------------------------------

def cf_partial_quotients(a, b, k=6):
    """First k partial quotients of a/b."""
    out = []
    while b and len(out) < k:
        q, r = divmod(a, b)
        out.append(q)
        a, b = b, r
    return out

def gauss_reduce_2d(u, v):
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    def dot(w, z):
        return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v):
        u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0:
            break
        q = (2 * num + den) // (2 * den) if num >= 0 else \
            -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u

def rho_lambda_block(n, lam, k1_bound, k2_bound, m):
    """rho = mu / E||v'||, mu = shortest vector of the 2D lambda block."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    mu = math.sqrt(w[0] ** 2 + w[1] ** 2)
    pv = math.sqrt(m * (k1_bound * S_K1) ** 2 / 3.0
                   + m * (k2_bound * S_K2) ** 2 / 3.0 + S_KAN ** 2)
    return mu / pv

def separators(p, n, lam, m, k1_bound, k2_bound):
    ls = lam_star(lam, n)
    lam_lo = min(lam % n, n - (lam % n))
    return {
        'lam*': ls,
        '|lam*-1/3|': abs(ls - 1.0 / 3.0),
        '|lam*-1/4|': abs(ls - 0.25),
        '|lam*-1/2|': abs(ls - 0.5),
        'rho': rho_lambda_block(n, lam, k1_bound, k2_bound, m),
        'delta3': abs(3 * lam_lo - n) / n,
        'maxq': max(cf_partial_quotients(lam_lo, n, 6)[1:] or [0]),
    }

# ---------------------------------------------------------------------------
# Success measurement (PROJ arm; PROJ == ORIG on every cell tested in U2/U3)
# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 555, 31337]

def proj_wins(curve, m, k1_bound, seeds=SEEDS):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    wins = 0
    for seed in seeds:
        d_secret = random.Random(seed + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
        if len(sigs) < m:
            continue
        M = build_proj(sigs, n, lam, k1_bound, k2_bound)
        A = IntegerMatrix.from_matrix(M)
        LLL.reduction(A)
        red = [[A[i][j] for j in range(2 * m + 1)] for i in range(2 * m + 2)]
        if recover_proj(red, sigs, n, lam, k1_bound, k2_bound, d_secret):
            wins += 1
    return wins

def search_curves_binned(lo, hi, nbins=20, per_bin=2):
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1 or not sympy.isprime(n_cand):
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

def auc(rows, key, label_key):
    """AUC of `key` as a score for P(label=1); reported for -key too, so the
    better-oriented value is the one to read."""
    pos = [r[key] for r in rows if r[label_key]]
    neg = [r[key] for r in rows if not r[label_key]]
    if not pos or not neg:
        return float('nan')
    tot = sum(1.0 if a < b else 0.5 if a == b else 0.0 for a in pos for b in neg)
    return tot / (len(pos) * len(neg))

def best_threshold(rows, key, label_key):
    """Best single threshold `key <= tau` -> predict success.  Returns
    (accuracy, tau, baseline)."""
    vals = sorted({r[key] for r in rows})
    n_tot = len(rows)
    base = max(sum(1 for r in rows if r[label_key]),
               sum(1 for r in rows if not r[label_key])) / n_tot
    best = (0.0, None)
    for tau in vals:
        acc = sum(1 for r in rows
                  if (r[key] <= tau) == bool(r[label_key])) / n_tot
        if acc > best[0]:
            best = (acc, tau)
    return best[0], best[1], base


# ===========================================================================
print("=" * 78)
print("Thread 23 / U6 - is |lam* - 1/3| the Phase-2 separator at moderate eff?")
print("=" * 78)

M = 12
EFFS = [0.10, 0.15, 0.20]

print("\nsearching 17-18 bit j=0 GLV curves, lam* binned over (0, 0.5) ...")
curves = search_curves_binned(2 ** 16, 2 ** 18, nbins=20, per_bin=2)
print(f"{len(curves)} curves found\n")

KEYS = ['lam*', '|lam*-1/3|', '|lam*-1/4|', '|lam*-1/2|', 'rho', 'delta3', 'maxq']

all_rows = {e: [] for e in EFFS}

print(f"{'p':>7} {'n':>7} {'lam*':>7} {'|l*-1/3|':>9} {'rho':>6} "
      f"{'delta3':>7} {'maxq':>6} | " +
      " ".join(f"{'eff=' + str(e):>9}" for e in EFFS))
for (p, b, n, lam, G) in curves:
    k2b = math.isqrt(n) + 1
    cells = []
    for e in EFFS:
        k1b = max(2, int(round(e * n / k2b)))
        w = proj_wins((p, b, n, lam, G), M, k1b)
        s = separators(p, n, lam, M, k1b, k2b)
        s['wins'] = w
        s['ok'] = w >= 3
        s['gamma'] = gamma_gap(M, n, k1b, k2b)
        all_rows[e].append(s)
        cells.append(f"{w}/5")
    s0 = separators(p, n, lam, M, max(2, int(round(0.15 * n / k2b))), k2b)
    print(f"{p:>7} {n:>7} {s0['lam*']:>7.4f} {s0['|lam*-1/3|']:>9.4f} "
          f"{s0['rho']:>6.3f} {s0['delta3']:>7.4f} {s0['maxq']:>6} | " +
          " ".join(f"{c:>9}" for c in cells))

print("\n" + "-" * 78)
print("Separator quality (label: >=3/5 seeds recover).  AUC is for the")
print("orientation 'smaller key => more likely success'; 0.5 = no signal.")
print("-" * 78)
for e in EFFS:
    rows = all_rows[e]
    npos = sum(1 for r in rows if r['ok'])
    print(f"\neff={e}   {npos}/{len(rows)} curves succeed   "
          f"mean gamma={sum(r['gamma'] for r in rows)/len(rows):.3f}")
    print(f"  {'separator':<12} {'AUC':>6} {'best acc':>9} {'tau':>9} "
          f"{'baseline':>9}")
    for k in KEYS:
        a = auc(rows, k, 'ok')
        acc, tau, base = best_threshold(rows, k, 'ok')
        tau_s = f"{tau:.4f}" if tau is not None and tau == tau else "-"
        print(f"  {k:<12} {a:>6.3f} {acc:>9.3f} {tau_s:>9} {base:>9.3f}")

print("\n" + "-" * 78)
print("Pooled over all three eff levels")
print("-" * 78)
pooled = [r for e in EFFS for r in all_rows[e]]
npos = sum(1 for r in pooled if r['ok'])
print(f"{npos}/{len(pooled)} cells succeed")
print(f"  {'separator':<12} {'AUC':>6} {'best acc':>9} {'tau':>9} {'baseline':>9}")
for k in KEYS:
    a = auc(pooled, k, 'ok')
    acc, tau, base = best_threshold(pooled, k, 'ok')
    tau_s = f"{tau:.4f}" if tau is not None and tau == tau else "-"
    print(f"  {k:<12} {a:>6.3f} {acc:>9.3f} {tau_s:>9} {base:>9.3f}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)

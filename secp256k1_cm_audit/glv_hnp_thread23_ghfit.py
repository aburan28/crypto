"""
GLV-HNP Phase 2, Thread 23 part 2: is the Gaussian-heuristic ratio THE predictor?

Thread 23 part 1 (glv_hnp_thread23_projection.py) found that quotienting out the
trivial vector n*S_D*e_m changes nothing: variants A (Kannan, dim 2m+2) and B
(projected Kannan, dim 2m+1) agree on 100% of cells, and even an exact-CVP
oracle on the projected BDD instance moves the K1 wall by at most one seed.

Its E0 table suggested why: all three embeddings have essentially the same
Gaussian-heuristic ratio

    r  =  ||v_planted||  /  GH(L),      GH(L) = sqrt(dim/2*pi*e) * det^(1/dim)

and r crosses 1 right where the empirical wall sits.  If r is the predictor,
then (i) it is a *lattice-level* invariant, which is exactly why the six
*curve-level* invariants tried 2026-06-21..06-29 (delta/n, kappa(M), q_cf,
max_q_cf, max_a, a_corn/n) all failed, and (ii) it subsumes nu_hat, the
lambda-block invariant that scored AUC 0.935 on 2026-07-29.

This script pools a large grid (curve x m x K1 x seed), all variant A, and:
  F1  reports the success rate binned by r
  F2  computes AUC and the best single-threshold accuracy for r, for nu_hat,
      and for eff = K1*K2/n, on the same pooled data
  F3  locates the empirical r at which the success rate crosses 50%
  F4  reports how far the wall sits below the information-theoretic bound

Run: python3 glv_hnp_thread23_ghfit.py
"""

import math
import os
import random
import sys

from fpylll import IntegerMatrix, LLL

HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_PATH = os.path.join(HERE, "glv_hnp_phase2_lambda_threshold.py")
with open(_SRC_PATH) as _f:
    _src = _f.read()
_helpers = _src.split("# EXPERIMENT T1")[0]
exec(compile(_helpers, _SRC_PATH + " (helpers)", "exec"), globals())


def gh_ratio(n, m, K1):
    """||v_planted|| / GH(L) for the dim-(2m+2) Kannan lattice (variant A)."""
    K2 = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
    dim = 2 * m + 2
    logdet = (m * math.log(n * S_K1) + math.log(S_D)
              + m * math.log(S_K2) + math.log(S_KAN))
    gh = math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(logdet / dim)
    tgt = planted_norm_expected(m, n, K1, K2, S_K1, S_D, S_K2, S_KAN)
    return tgt / gh


def nu_hat(n, lam, K1):
    """lambda_1(L2)/sqrt(det L2) for L2 = <(n*S_K1,0), (-lam*S_K1, S_K2)>,
    the 2026-07-29 separator."""
    K2 = math.isqrt(n) + 1
    S_K1, _S_D, S_K2, _S_KAN = scales(n, K1, K2)
    mu, _w = lambda_block_mu(n, lam, S_K1, S_K2)
    det = n * S_K1 * S_K2
    return mu / math.sqrt(det)


def it_slack(n, m, K1):
    """eff / n^(-1/m).  < 1 means the instance is information-theoretically
    solvable (enough signature bits to pin down k1, k2 and d)."""
    K2 = math.isqrt(n) + 1
    return (K1 * K2 / n) / (n ** (-1.0 / m))


def trial(curve, m, K1, seed):
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    d_secret = random.Random(seed ^ 0xC0FFEE).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, K1, K2, seed)
    if len(sigs) < m:
        return None
    _S_K1, _S_D, _S_K2, S_KAN = scales(n, K1, K2)
    M = build_glv_lattice(sigs, n, lam, K1, K2)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    dim = 2 * m + 2
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(rows, m, n, S_KAN, d_secret) is not None


def auc(rows, key, positive_is_low=True):
    """Mann-Whitney AUC of `key` as a predictor of success."""
    pos = [r[key] for r in rows if r['ok']]
    neg = [r[key] for r in rows if not r['ok']]
    if not pos or not neg:
        return float('nan')
    wins = ties = 0
    for a in pos:
        for b in neg:
            if (a < b) if positive_is_low else (a > b):
                wins += 1
            elif a == b:
                ties += 1
    return (wins + 0.5 * ties) / (len(pos) * len(neg))


def best_threshold(rows, key, positive_is_low=True):
    vals = sorted({r[key] for r in rows})
    best = (0.0, None)
    for i in range(len(vals)):
        t = vals[i]
        correct = 0
        for r in rows:
            pred = (r[key] <= t) if positive_is_low else (r[key] >= t)
            if pred == r['ok']:
                correct += 1
        acc = correct / len(rows)
        if acc > best[0]:
            best = (acc, t)
    return best


print("=" * 78)
print("Thread 23 part 2 — the Gaussian-heuristic ratio as the Phase-2 predictor")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

curves = []
# the two historical 12-bit curves
for p, b, n, lam in ((2557, 2, 2659, 1755), (2677, 2, 2647, 185)):
    G = find_generator(p, b, n)
    assert G is not None
    curves.append((p, b, n, lam, G))
# fresh 14-bit and 17-bit curves, lam* spread across [0, 0.5]
curves += search_curves(2 ** 13, 2 ** 14, per_bin=1, nbins=6)
curves += search_curves(2 ** 16, 2 ** 17, per_bin=1, nbins=8)
print(f"\ncurves: {len(curves)}  "
      f"(bit sizes {sorted({c[2].bit_length() for c in curves})})")

M_GRID = [6, 8, 10, 12, 16]
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]

rows = []
for ci, curve in enumerate(curves):
    p, b, n, lam, G = curve
    for m in M_GRID:
        for K1 in K1_GRID:
            r = gh_ratio(n, m, K1)
            if r > 3.0:                    # hopeless, do not waste LLL calls
                continue
            for s in SEEDS:
                ok = trial(curve, m, K1, s)
                if ok is None:
                    continue
                rows.append({
                    'n': n, 'nb': n.bit_length(), 'lam': lam, 'm': m, 'K1': K1,
                    'r': r, 'nu': nu_hat(n, lam, K1),
                    'eff': K1 * (math.isqrt(n) + 1) / n,
                    'its': it_slack(n, m, K1), 'ok': ok,
                })
    sys.stdout.write(f"\r  curve {ci+1}/{len(curves)}  trials so far: {len(rows)}")
    sys.stdout.flush()
print(f"\n\ntotal trials: {len(rows)}   successes: {sum(r['ok'] for r in rows)}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("F1: success rate binned by the Gaussian-heuristic ratio r")
print("-" * 78)
EDGES = [0, .4, .6, .7, .8, .9, .95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.7, 3.01]
print(f"{'r bin':>14} {'trials':>7} {'succ':>6} {'rate':>7}")
for i in range(len(EDGES) - 1):
    lo, hi = EDGES[i], EDGES[i + 1]
    sel = [r for r in rows if lo <= r['r'] < hi]
    if not sel:
        continue
    s = sum(r['ok'] for r in sel)
    print(f"[{lo:>5.2f},{hi:>5.2f}) {len(sel):>7} {s:>6} {s/len(sel):>7.3f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("F2: AUC and best-threshold accuracy, same pooled data")
print("-" * 78)
base = max(sum(r['ok'] for r in rows), len(rows) - sum(r['ok'] for r in rows)) / len(rows)
print(f"majority baseline accuracy: {base:.3f}\n")
print("AUC is reported in the orientation each predictor was proposed in")
print("(low r / low eff / low nu_hat = success; nu_hat's sign is the one found")
print("on 2026-07-29: 'a short rival vector makes the attack easier').\n")
print(f"{'predictor':<28} {'AUC':>7} {'best acc':>9} {'threshold':>11}")
for name, key, low in (("r  = |t|/GH(L)  [this run]", 'r', True),
                       ("nu_hat  [2026-07-29]", 'nu', True),
                       ("eff = K1*K2/n", 'eff', True),
                       ("eff / IT bound", 'its', True)):
    a = auc(rows, key, low)
    acc, thr = best_threshold(rows, key, low)
    print(f"{name:<28} {a:>7.3f} {acc:>9.3f} {thr:>11.4f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("F3: where does the success rate cross 50%?")
print("-" * 78)
fine = []
step = 0.05
x = 0.3
while x < 2.0:
    sel = [r for r in rows if x <= r['r'] < x + step]
    if len(sel) >= 15:
        fine.append((x + step / 2, sum(r['ok'] for r in sel) / len(sel), len(sel)))
    x += step
cross = None
for i in range(1, len(fine)):
    if fine[i - 1][1] >= 0.5 > fine[i][1]:
        cross = fine[i - 1][0], fine[i][0]
for c, rate_, k in fine:
    bar = "#" * int(round(rate_ * 40))
    print(f"  r={c:>5.2f}  n={k:>4}  {rate_:>5.3f}  {bar}")
print(f"\n50% crossing between r = {cross[0]:.2f} and {cross[1]:.2f}"
      if cross else "\nno clean 50% crossing")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("F4: how far below the information-theoretic bound is the wall?")
print("-" * 78)
succ = [r for r in rows if r['ok']]
if succ:
    print(f"max eff/IT among SUCCESSES : {max(r['its'] for r in succ):.4f}")
    print(f"90th pct eff/IT (successes): "
          f"{sorted(r['its'] for r in succ)[int(0.9*len(succ))]:.4f}")
print(f"eff/IT = 1 is the information-theoretic limit; the gap is the headroom")
print(f"a non-lattice method (e.g. Bleichenbacher FFT) could in principle use.")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("F5: does nu_hat retain any power AFTER controlling for r?")
print("-" * 78)
print("r is a function of (n, m, K1) only, so within a fixed (n-bits, m, K1)")
print("stratum r is (nearly) constant and any residual separation across curves")
print("is genuine curve-level signal.  This is the fair test of the 2026-07-29")
print("claim, which was measured at a single (m, K1).\n")
strata = {}
for r in rows:
    strata.setdefault((r['nb'], r['m'], r['K1']), []).append(r)
mixed = {k: v for k, v in strata.items()
         if 0 < sum(x['ok'] for x in v) < len(v) and len(v) >= 20}
print(f"mixed strata (both outcomes present, >=20 trials): {len(mixed)}")
if mixed:
    aucs = []
    print(f"\n{'nb':>3} {'m':>3} {'K1':>3} {'trials':>7} {'rate':>6} {'nu AUC':>7}")
    for k in sorted(mixed):
        v = mixed[k]
        a = auc(v, 'nu', True)
        aucs.append(a)
        print(f"{k[0]:>3} {k[1]:>3} {k[2]:>3} {len(v):>7} "
              f"{sum(x['ok'] for x in v)/len(v):>6.3f} {a:>7.3f}")
    mean_auc = sum(aucs) / len(aucs)
    above = sum(1 for a in aucs if a > 0.5)
    print(f"\nmean within-stratum nu_hat AUC: {mean_auc:.3f} "
          f"({above}/{len(aucs)} strata above 0.5)")
    print("0.5 = no signal.  Compare the 0.935 reported 2026-07-29 on the C1/C2 slice.")
    print("\nbreakdown by n bit-length (does the signal survive scaling?):")
    by_nb = {}
    for k in sorted(mixed):
        by_nb.setdefault(k[0], []).append(auc(mixed[k], 'nu', True))
    for nb in sorted(by_nb):
        v = by_nb[nb]
        print(f"  {nb} bits: {len(v):>2} strata, mean nu_hat AUC "
              f"{sum(v)/len(v):>5.3f}, {sum(1 for a in v if a > 0.5)}/{len(v)} above 0.5")

print("\n" + "=" * 78)
print("done")
print("=" * 78)

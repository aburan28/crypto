"""
GLV-HNP Phase 2, Thread 23 part 4: calibrate the Psi threshold and extrapolate.

Part 3 established that the controlling quantity is

    Psi(m) = gamma(m) * delta^(2m+2),
    Psi_min = sqrt(2*pi*e/3) * sqrt(eff) * exp( 2*sqrt( log(n/eff)*log(delta) ) )
    attained at   dim* = 2m*+2 = sqrt( log(n/eff) / log(delta) ).

This script asks the only question that matters for the paper: is the empirical
threshold Psi* (below which recovery succeeds) STABLE across n?  If it is, the
Phase-2 attack's reach on a real curve follows from the closed form.

Inputs are the recorded outcomes of three independent measurement campaigns:
  - 2026-07-29 T4    (12-bit, m=12, 5 seeds)          log line ~6060
  - 2026-08-03 T23-E (16-bit, m=12, 5 seeds x 10 curves)
  - 2026-08-03 T23-C (12-bit, m=12, 5 seeds, re-measured)

No new lattice reduction is run here; this is a pure post-hoc calibration of the
closed form against data already in the repository.

Run: python3 glv_hnp_phase2_thread23_calibration.py
"""

import math

DELTA = 1.021
C = math.sqrt(2 * math.pi * math.e / 3.0)


def psi_min(n, eff, delta=DELTA):
    return C * math.sqrt(eff) * math.exp(2 * math.sqrt(math.log(n / eff)
                                                       * math.log(delta)))


def dim_star(n, eff, delta=DELTA):
    return math.sqrt(math.log(n / eff) / math.log(delta))


# (label, n, K1, K2, wins, trials)   -- all at m = 12
DATA = [
    # --- 12-bit, curve 2557 (lam* = 0.340), 2026-07-29 T4 / 2026-08-03 T23-C
    ("2557", 2659, 2,  52, 5, 5), ("2557", 2659, 3,  52, 5, 5),
    ("2557", 2659, 4,  52, 5, 5), ("2557", 2659, 6,  52, 5, 5),
    ("2557", 2659, 8,  52, 5, 5), ("2557", 2659, 12, 52, 4, 5),
    ("2557", 2659, 16, 52, 1, 5), ("2557", 2659, 24, 52, 0, 5),
    # --- 12-bit, curve 2677 (lam* = 0.070)
    ("2677", 2647, 2,  52, 5, 5), ("2677", 2647, 3,  52, 5, 5),
    ("2677", 2647, 4,  52, 5, 5), ("2677", 2647, 6,  52, 2, 5),
    ("2677", 2647, 8,  52, 0, 5), ("2677", 2647, 12, 52, 0, 5),
    ("2677", 2647, 16, 52, 0, 5), ("2677", 2647, 24, 52, 0, 5),
    # --- 16-bit sweep, 2026-08-03 T23-E (K2 = isqrt(n)+1)
    ("17b/65719", 65719, 12, 257, 5, 5), ("17b/66499", 66499, 12, 258, 5, 5),
    ("17b/65203", 65203, 12, 256, 5, 5), ("17b/66751", 66751, 12, 259, 4, 5),
    ("17b/66037", 66037, 12, 257, 5, 5), ("17b/66301", 66301, 12, 258, 5, 5),
    ("17b/65287", 65287, 12, 256, 5, 5), ("17b/65053", 65053, 12, 256, 5, 5),
    ("17b/66109", 66109, 12, 258, 5, 5), ("17b/67369", 67369, 12, 260, 5, 5),
    ("17b/65719", 65719, 38, 257, 1, 5), ("17b/66499", 66499, 38, 258, 0, 5),
    ("17b/65203", 65203, 38, 256, 0, 5), ("17b/66751", 66751, 38, 259, 0, 5),
    ("17b/66037", 66037, 38, 257, 0, 5), ("17b/66301", 66301, 38, 258, 0, 5),
    ("17b/65287", 65287, 38, 256, 5, 5), ("17b/65053", 65053, 38, 256, 0, 5),
    ("17b/66109", 66109, 38, 258, 0, 5), ("17b/67369", 67369, 38, 260, 5, 5),
    ("17b/65719", 65719, 63, 257, 0, 5), ("17b/66499", 66499, 64, 258, 0, 5),
    ("17b/65203", 65203, 63, 256, 0, 5), ("17b/66751", 66751, 64, 259, 0, 5),
    ("17b/66037", 66037, 64, 257, 0, 5), ("17b/66301", 66301, 64, 258, 0, 5),
    ("17b/65287", 65287, 63, 256, 1, 5), ("17b/65053", 65053, 63, 256, 0, 5),
    ("17b/66109", 66109, 64, 258, 0, 5), ("17b/67369", 67369, 64, 260, 3, 5),
]

print("=" * 84)
print("Thread 23 part 4 — is the Psi threshold stable across n?")
print("=" * 84)
print(f"delta = {DELTA}, C = sqrt(2*pi*e/3) = {C:.4f}\n")

rows = []
for label, n, k1, k2, w, t in DATA:
    eff = k1 * k2 / n
    rows.append({'label': label, 'n': n, 'k1': k1, 'eff': eff,
                 'psi': psi_min(n, eff), 'dim*': dim_star(n, eff),
                 'w': w, 't': t, 'ok': w == t})

rows.sort(key=lambda r: r['psi'])
print(f"{'curve':>11} {'n':>7} {'K1':>4} {'eff':>7} {'dim*':>6} {'Psi_min':>8} "
      f"{'obs':>6}")
for r in rows:
    print(f"{r['label']:>11} {r['n']:>7} {r['k1']:>4} {r['eff']:>7.4f} "
          f"{r['dim*']:>6.1f} {r['psi']:>8.3f} "
          f"{str(r['w'])+'/'+str(r['t']):>6}")

# --- best single threshold ------------------------------------------------
print("\n" + "-" * 84)
print("Best single threshold classifier (predict full recovery iff Psi < thr)")
print("-" * 84)
best = (0, None)
cands = sorted(r['psi'] for r in rows)
for i in range(len(cands) + 1):
    thr = (cands[0] - 0.01 if i == 0 else
           cands[-1] + 0.01 if i == len(cands) else
           (cands[i - 1] + cands[i]) / 2)
    acc = sum(1 for r in rows if (r['psi'] < thr) == r['ok'])
    if acc > best[0]:
        best = (acc, thr)
print(f"  best accuracy {best[0]}/{len(rows)} = {best[0]/len(rows):.3f} "
      f"at Psi* = {best[1]:.3f}")
maj = max(sum(1 for r in rows if r['ok']), sum(1 for r in rows if not r['ok']))
print(f"  majority baseline {maj}/{len(rows)} = {maj/len(rows):.3f}")

print("\n  per-bit-length breakdown at the global Psi*:")
for grp, sel in (("12-bit", lambda r: r['n'] < 10000),
                 ("16-bit", lambda r: r['n'] >= 10000)):
    sub = [r for r in rows if sel(r)]
    acc = sum(1 for r in sub if (r['psi'] < best[1]) == r['ok'])
    lo = max([r['psi'] for r in sub if r['ok']], default=float('nan'))
    hi = min([r['psi'] for r in sub if not r['ok']], default=float('nan'))
    print(f"    {grp}: {acc}/{len(sub)}   highest success Psi={lo:.3f}, "
          f"lowest failure Psi={hi:.3f}")

# --- extrapolation --------------------------------------------------------
print("\n" + "-" * 84)
print("Extrapolation: what bias does Phase 2 need on a real curve?")
print(f"Solve Psi_min(n, eff) = Psi* = {best[1]:.3f} for eff; report K1 with")
print("K2 = sqrt(n) (balanced GLV decomposition, the realistic nonce model).")
print("-" * 84)
print(f"{'bits':>6} {'log2 eff':>10} {'dim*':>7} {'m*':>6} "
      f"{'k1 bits':>9} {'bias':>7}")
for bits in (32, 48, 64, 96, 128, 160, 192, 256, 384, 521):
    n = 2.0 ** bits
    lo, hi = 1e-30, 1.0
    for _ in range(300):
        mid = math.sqrt(lo * hi)
        if psi_min(n, mid) < best[1]:
            lo = mid
        else:
            hi = mid
    eff = lo
    ds = dim_star(n, eff)
    k1_bits = math.log2(eff) + bits / 2.0     # K1 = eff*n/K2 = eff*sqrt(n)
    print(f"{bits:>6} {math.log2(eff):>10.2f} {ds:>7.1f} {(ds-2)/2:>6.1f} "
          f"{k1_bits:>9.1f} {bits/2 - k1_bits:>7.1f}")

print("\n  'bias' = bits of k1 that must be zero, relative to the balanced")
print("  GLV nonce k1 ~ sqrt(n).  Balanced GLV (bias = 0) gives eff = 1 and")
print("  Psi_min = inf: NOT attackable at any m.")

# ===========================================================================
# Combine Psi with nu_hat, the lam-geometry predictor from 2026-07-29
# (commit e845207 / glv_hnp_nuhat_vs_c1c2.py; AUC 0.935 on the C1/C2 sample).
#   L2      = < (n*S_K1, 0), (-lam*S_K1, S_K2) >   -- one of m identical copies
#   nu_hat  = lambda_1(L2) / sqrt(det L2),   det L2 = n*S_K1*S_K2
# Psi captures the scale/dimension budget; nu_hat captures the lam-dependent
# density of NON-planted short vectors.  They should be complementary: all four
# Psi misclassifications above are high-lam* curves.
# ===========================================================================
print("\n" + "-" * 84)
print("Psi x nu_hat: does the 2026-07-29 lam-geometry predictor close the gap?")
print("-" * 84)


def lagrange_gauss(u, v):
    """Exact 2-D Lagrange-Gauss reduction; returns the shortest vector norm."""
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    if nrm2(u) > nrm2(v):
        u, v = v, u
    while True:
        q = round((u[0] * v[0] + u[1] * v[1]) / nrm2(u))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            return math.sqrt(nrm2(u))
        u, v = v, u


def nu_hat(n, lam, k1, k2):
    S_K1 = n // k1
    S_K2 = max(1, n // k2)
    l1 = lagrange_gauss((n * S_K1, 0), (-lam * S_K1, S_K2))
    return l1 / math.sqrt(n * S_K1 * S_K2)


LAMS = {2659: 1755, 2647: 185, 65719: None}   # 12-bit lams are known exactly
# lam for the 16-bit curves: recomputed from n via the GLV root (both roots
# give the same nu_hat, since L2(lam) and L2(n-1-lam) are unimodularly equal).
import sys, os                                                    # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from glv_hnp_phase2_thread23_bdd import glv_roots                 # noqa: E402

for r in rows:
    n = r['n']
    if n in (2659, 2647):
        lam = LAMS[n]
    else:
        gr = glv_roots(n)
        lam = gr[0] if gr else None
    r['lam'] = lam
    r['lam_star'] = min(lam, n - lam) / n if lam else float('nan')
    r['nu'] = nu_hat(n, lam, r['k1'], round(math.isqrt(n)) + 1) if lam else None

usable = [r for r in rows if r['nu'] is not None]
print(f"{'curve':>11} {'K1':>4} {'lam*':>6} {'Psi':>7} {'nu_hat':>8} "
      f"{'Psi/nu':>8} {'obs':>6} {'Psi-only':>9}")
for r in sorted(usable, key=lambda r: r['psi'] / r['nu']):
    pred_psi = r['psi'] < best[1]
    print(f"{r['label']:>11} {r['k1']:>4} {r['lam_star']:>6.3f} "
          f"{r['psi']:>7.3f} {r['nu']:>8.4f} {r['psi']/r['nu']:>8.3f} "
          f"{str(r['w'])+'/'+str(r['t']):>6} "
          f"{('ok' if pred_psi == r['ok'] else 'MISS'):>9}")


def best_thr(items, key, direction=-1):
    """direction=-1: predict success iff value < thr."""
    vals = sorted(key(r) for r in items)
    bb = (0, None)
    for i in range(len(vals) + 1):
        thr = (vals[0] - 1e-9 if i == 0 else vals[-1] + 1e-9 if i == len(vals)
               else (vals[i - 1] + vals[i]) / 2)
        acc = sum(1 for r in items
                  if ((key(r) < thr) if direction < 0 else (key(r) >= thr))
                  == r['ok'])
        if acc > bb[0]:
            bb = (acc, thr)
    return bb


N = len(usable)
a_psi = best_thr(usable, lambda r: r['psi'])
a_nu = best_thr(usable, lambda r: r['nu'], direction=+1)
a_rat = best_thr(usable, lambda r: r['psi'] / r['nu'])
print(f"\n  Psi alone      : {a_psi[0]}/{N} = {a_psi[0]/N:.3f} at thr {a_psi[1]:.3f}")
print(f"  nu_hat alone   : {a_nu[0]}/{N} = {a_nu[0]/N:.3f} at thr {a_nu[1]:.4f}")
print(f"  Psi/nu_hat     : {a_rat[0]}/{N} = {a_rat[0]/N:.3f} at thr {a_rat[1]:.3f}")

bb = (0, None, None)
for j in range(-40, 41):
    c = j / 10.0
    acc, thr = best_thr(usable, lambda r, c=c: math.log(r['psi'])
                        - c * math.log(r['nu']))
    if acc > bb[0]:
        bb = (acc, c, thr)
print(f"  log Psi - c*log nu_hat : {bb[0]}/{N} = {bb[0]/N:.3f} "
      f"at c = {bb[1]:.1f}, thr = {bb[2]:.3f}")
print(f"  majority baseline      : {maj}/{N} = {maj/N:.3f}")

# --- leave-one-out cross-validation ---------------------------------------
# The combined score fits 2 parameters (c and the threshold) on 46 points, so
# its in-sample 45/46 is not evidence on its own.  LOOCV is.
print("\n  leave-one-out cross-validation (fit on N-1, predict the held-out):")


def loocv_psi(items):
    hits = 0
    for i in range(len(items)):
        tr = items[:i] + items[i+1:]
        _, thr = best_thr(tr, lambda r: r['psi'])
        hits += ((items[i]['psi'] < thr) == items[i]['ok'])
    return hits


def loocv_combo(items):
    hits = 0
    for i in range(len(items)):
        tr = items[:i] + items[i+1:]
        bb2 = (0, None, None)
        for j in range(-40, 41):
            c = j / 10.0
            acc, thr = best_thr(tr, lambda r, c=c: math.log(r['psi'])
                                - c * math.log(r['nu']))
            if acc > bb2[0]:
                bb2 = (acc, c, thr)
        s = math.log(items[i]['psi']) - bb2[1] * math.log(items[i]['nu'])
        hits += ((s < bb2[2]) == items[i]['ok'])
    return hits


lp = loocv_psi(usable)
lc = loocv_combo(usable)
print(f"    Psi alone (1 param)          : {lp}/{N} = {lp/N:.3f}")
print(f"    log Psi - c*log nu (2 params): {lc}/{N} = {lc/N:.3f}")
print(f"    majority baseline            : {maj}/{N} = {maj/N:.3f}")
print("=" * 84)

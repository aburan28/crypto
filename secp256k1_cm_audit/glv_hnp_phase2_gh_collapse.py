"""
GLV-HNP Phase 2, Thread 23c: does ||pv||/GH collapse the (m, eff, n) grid?

Where this comes from
---------------------
23a (glv_hnp_phase2_projected.py): quotienting out the trivial n*S_D*e_m
vector changed nothing -- 96/100 vs 96/100, 18/100 vs 18/100, 9/100 vs 9/100.
Masking is not the obstruction.

23b (glv_hnp_phase2_density_wall.py): the obstruction is lattice density, but
the closed form proposed there -- a fixed threshold eff* = 3/(2*pi*e) = 0.176
independent of m -- was FALSIFIED at the m used in practice:

  * ||pv||/GH is strongly m-dependent (m=6/12/24/48 at eff=0.05 gives
    1.637/0.976/0.730/0.626); 3/(2*pi*e) is only the m -> infinity limit and
    convergence is slow (n^(1/D) decays only like exp(log n / 2m)).
  * success DOES rise with m at fixed eff (eff=0.14: 0.00 0.07 0.03 0.27 0.33
    0.37 for m = 6 8 12 16 24 32), so the wall is not a fixed eff.
    This also corrects T4b of 2026-07-29, which read "more data does not
    rescue it" off a single curve x 5 seeds (0,0,1,0,1) -- 30 trials per cell
    here show the rise.

So the quantity to test is the scalar itself, not a threshold in eff.  For the
baseline (2m+2)-dim lattice, with D = 2m+2 and eff = K1*K2/n:

    ||pv||/GH  =  n^(1/D) * eff^(m/D) * sqrt( 2*pi*e * (2m+4) / (3*D) )

This script asks whether that one number predicts recovery across the whole
(m, eff, n) grid, and where its empirical LLL threshold tau sits.

Grid: bits in {14,17,20} x 8 curves, m in {6,8,12,16,24,32},
      eff in {0.04,0.06,0.08,0.12,0.16,0.22,0.30}, 3 seeds.

Competing predictors scored by AUC on the same trials: eff alone, m alone,
and nu_hat = lambda_1(L2)/sqrt(det L2) for the 2D lambda-block
L2 = <(n*S_K1,0), (-lam*S_K1,S_K2)> (the 2026-07-29 separator).

Run: python3 glv_hnp_phase2_gh_collapse.py
"""

import math
import sys

from glv_hnp_phase2_projected import (  # noqa: E402
    scales, gen_signatures, lam_star, search_curves, trial,
)

TWO_PI_E = 2 * math.pi * math.e


def gh_ratio(n, k1_bound, k2_bound, m):
    """||pv||/GH for the baseline lattice, using the actual integer scales."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    dim = 2 * m + 2
    log_det = (m * math.log(n * S_K1) + math.log(S_D)
               + m * math.log(S_K2) + math.log(S_KAN))
    pv = math.sqrt(m * (k1_bound * S_K1) ** 2 / 3.0 + (n * S_D) ** 2 / 3.0
                   + m * (k2_bound * S_K2) ** 2 / 3.0 + S_KAN ** 2)
    g = math.sqrt(dim / TWO_PI_E) * math.exp(log_det / dim)
    return pv / g


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
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u


def nu_hat(n, lam, k1_bound, k2_bound):
    """lambda_1(L2)/sqrt(det L2) for the 2D lambda-block (2026-07-29)."""
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    l1 = math.sqrt(w[0] * w[0] + w[1] * w[1])
    return l1 / math.sqrt(n * S_K1 * S_K2)


def auc(scores, labels):
    """AUC for 'lower score => success'.  Ties get 0.5."""
    pos = [s for s, y in zip(scores, labels) if y]
    neg = [s for s, y in zip(scores, labels) if not y]
    if not pos or not neg:
        return float('nan')
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks, i = [0.0] * len(scores), 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        r = (i + j) / 2.0 + 1.0
        for t in range(i, j + 1):
            ranks[order[t]] = r
        i = j + 1
    rsum = sum(r for r, y in zip(ranks, labels) if y)
    np_, nn = len(pos), len(neg)
    # low score = success, so invert
    return 1.0 - (rsum - np_ * (np_ + 1) / 2.0) / (np_ * nn)


print("=" * 78)
print("Thread 23c -- does ||pv||/GH collapse the Phase-2 (m, eff, n) grid?")
print("=" * 78)
sys.stdout.flush()

BITS = [14, 17, 20]
M_GRID = [6, 8, 12, 16, 24, 32]
EFF_GRID = [0.04, 0.06, 0.08, 0.12, 0.16, 0.22, 0.30]
SEEDS = [42, 1234, 9999]

curves = {}
for bits in BITS:
    cs = search_curves(2 ** (bits - 1), 2 ** bits, per_bin=1, nbins=8)[:8]
    curves[bits] = cs
    print(f"{bits:>2} bits: {len(cs)} curves, lam* in "
          f"[{min(lam_star(c[3], c[2]) for c in cs):.4f}, "
          f"{max(lam_star(c[3], c[2]) for c in cs):.4f}]")
sys.stdout.flush()

rows = []   # (ratio, eff, m, bits, nu, ok)
total = sum(len(curves[b]) for b in BITS) * len(M_GRID) * len(EFF_GRID) * len(SEEDS)
print(f"\nrunning {total} trials ...")
sys.stdout.flush()

done = 0
for bits in BITS:
    for cur in curves[bits]:
        p, b, n, lam, G = cur
        k2 = math.isqrt(n) + 1
        d_secret = n // 3 + 7
        for eff in EFF_GRID:
            k1 = max(2, round(eff * n / k2))
            eff_act = k1 * k2 / n
            nu = nu_hat(n, lam, k1, k2)
            for m in M_GRID:
                r = gh_ratio(n, k1, k2, m)
                for s in SEEDS:
                    t = trial(cur, m, d_secret, k1, s, 'kannan')
                    done += 1
                    if t is None:
                        continue
                    rows.append((r, eff_act, m, bits, nu, bool(t['ok'])))
        print(f"  {bits}b p={p} done ({done}/{total})")
        sys.stdout.flush()

print(f"\ncollected {len(rows)} trials, "
      f"{sum(1 for r in rows if r[5])} successes")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("C1: success rate binned by ||pv||/GH, split by m and by n")
print("-" * 78)
print("if the scalar collapses the grid, the columns agree within noise.\n")

EDGES = [0.0, 0.7, 0.85, 1.0, 1.15, 1.3, 1.5, 1.75, 2.0, 2.5, 99.0]

def binidx(x):
    for i in range(len(EDGES) - 1):
        if EDGES[i] <= x < EDGES[i + 1]:
            return i
    return len(EDGES) - 2

hdr = f"{'pv/GH bin':>14}{'all':>10}"
for m in M_GRID:
    hdr += f"{'m=' + str(m):>8}"
for bits in BITS:
    hdr += f"{str(bits) + 'b':>7}"
print(hdr)
for i in range(len(EDGES) - 1):
    sel = [r for r in rows if binidx(r[0]) == i]
    if not sel:
        continue
    line = f"{EDGES[i]:>6.2f}-{EDGES[i+1]:<7.2f}"
    line += f"{sum(1 for r in sel if r[5])/len(sel):>6.2f}({len(sel):>3})"
    for m in M_GRID:
        s2 = [r for r in sel if r[2] == m]
        line += f"{(sum(1 for r in s2 if r[5])/len(s2)) if s2 else float('nan'):>8.2f}"
    for bits in BITS:
        s2 = [r for r in sel if r[3] == bits]
        line += f"{(sum(1 for r in s2 if r[5])/len(s2)) if s2 else float('nan'):>7.2f}"
    print(line)

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("C2: predictor comparison (AUC, lower score = predicted success)")
print("-" * 78)
labels = [r[5] for r in rows]
print(f"  ||pv||/GH      {auc([r[0] for r in rows], labels):.4f}")
print(f"  eff alone      {auc([r[1] for r in rows], labels):.4f}")
print(f"  -m (more sigs) {auc([-r[2] for r in rows], labels):.4f}")
print(f"  nu_hat         {auc([r[4] for r in rows], labels):.4f}")
print(f"  baseline       0.5000")

print("\n  within-cell AUC of nu_hat (fixed m and eff, so pv/GH is constant):")
cells, num, den = {}, 0.0, 0
for r in rows:
    cells.setdefault((r[2], round(r[1], 4), r[3]), []).append(r)
for k, v in cells.items():
    ys = [x[5] for x in v]
    if 0 < sum(ys) < len(ys):
        a = auc([x[4] for x in v], ys)
        if not math.isnan(a):
            num += a * len(v); den += len(v)
print(f"    {num/den:.4f}  (weighted over {den} trials in mixed cells)"
      if den else "    n/a")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("C3: empirical threshold tau  (50% success in pv/GH)")
print("-" * 78)
pts = sorted(rows, key=lambda r: r[0])
W = max(50, len(pts) // 25)
tau = None
prev = None
for i in range(0, len(pts) - W, max(1, W // 4)):
    win = pts[i:i + W]
    x = sum(r[0] for r in win) / W
    y = sum(1 for r in win if r[5]) / W
    if prev and prev[1] >= 0.5 > y:
        (x0, y0) = prev
        tau = x0 + (0.5 - y0) * (x - x0) / (y - y0)
        break
    prev = (x, y)
print(f"  sliding-window ({W} trials) 50% crossing:  tau = "
      f"{tau:.3f}" if tau else "  no crossing")
print(f"  (pv/GH = 1 is the Gaussian-heuristic point; LLL keeps working past it)")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("C4: what tau implies -- max attackable eff vs m, extrapolated in n")
print("-" * 78)
if tau:
    print(f"solving  n^(1/D) * eff^(m/D) * sqrt(2*pi*e*(2m+4)/(3D)) = tau = {tau:.3f}")
    print(f"for eff, with D = 2m+2.  (Heuristic: assumes the collapse in C1")
    print(f"holds at cryptographic size; it is verified here only to 20 bits.)\n")
    print(f"{'m':>6}" + "".join(f"{str(bb) + '-bit':>11}" for bb in (20, 64, 128, 256)))
    for m in (12, 32, 64, 128, 256, 512, 1024):
        D = 2 * m + 2
        s = math.sqrt(TWO_PI_E * (2 * m + 4) / (3 * D))
        cells = []
        for bb in (20, 64, 128, 256):
            log2_eff = (math.log2(tau) - bb / D - math.log2(s)) * D / m
            e = 2 ** log2_eff
            cells.append(f"{e:>11.2e}" if e < 0.01 else f"{e:>11.4f}")
        print(f"{m:>6}" + "".join(cells))
    print(f"\nasymptote (m -> inf, any n): eff -> (tau/sqrt(2*pi*e/3))^2 = "
          f"{(tau / math.sqrt(TWO_PI_E / 3.0)) ** 2:.3f}")
    print("\nRead for secp256k1 (n ~ 2^256, K2 = sqrt(n) = 2^128, K1 = eff*n/K2):")
    for m in (12, 128, 1024):
        D = 2 * m + 2
        s = math.sqrt(TWO_PI_E * (2 * m + 4) / (3 * D))
        log2_eff = (math.log2(tau) - 256 / D - math.log2(s)) * D / m
        print(f"  m={m:>5}: eff = 2^{log2_eff:>7.2f}  ->  k1 < 2^{128 + log2_eff:.1f}"
              f"  ({128 - (128 + log2_eff):.0f} bits of k1 must be known)")

print("\n" + "=" * 78)
print("done")
print("=" * 78)

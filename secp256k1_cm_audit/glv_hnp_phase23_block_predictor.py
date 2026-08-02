"""
GLV-HNP Phase 2, Thread 23 (part 3): a curve-level predictor that actually
works -- the geometry of the 2-D lambda-block.

Threads 20/23 established:
  * eff = K1*K2/n is the dominant variable, but at FIXED eff two curves of
    the same size can differ by 5/5 vs 0/5 (12-bit/2557 vs 12-bit/2677);
  * lam*/n is not a predictor (Thread 20 T3);
  * mu = |lambda_1(B)| alone is not a predictor (Thread 20 T2, hypothesis
    rho = mu/||v_planted|| >= 1, falsified);
  * gamma = ||v||/lambda_1^GH is not a predictor either (Thread 23 U1/U3):
    the lattice is not GH-random, it is m identical copies of one 2-D block.

The block is the whole story.  Write

    B = < (n*S_K1, 0), (-lam*S_K1, S_K2) >,     det B = n*S_K1*S_K2,

Gauss-reduced to (b1, b2) with |b1| = mu1 <= mu2 = |b2|.  det B is the SAME
for all curves at a given (n, eff) -- only the SHAPE of B varies with lam.
Per signature, the true coset contributes an energy

    Eplant = (K1*S_K1)^2/3 + (K2*S_K2)^2/3   ~  2n^2/3      (curve-free)

while a wrong coset contributes  Erand = E[dist^2(u, B)]  for u uniform mod
B -- and THAT depends on the shape: a skew block (mu2/mu1 large) has a long
thin Voronoi cell and a large Erand, so wrong cosets are pushed further away.
Define

    nu = sqrt( Eplant / Erand )     (per-signature signal-to-noise)

H23b: nu -- not lam*, not mu, not gamma -- is the curve-level predictor.
      nu < 1 is the per-signature condition; the union bound over the n
      candidate cosets then needs  m * ln(1/nu^2) > 2*ln n, i.e.

          m_star(nu, n) = ceil( 2*ln n / ln(1/nu^2) ).

Falsifier: hold eff fixed (so det B, Eplant, and gamma are all fixed) and
vary the curve.  If nu has no more predictive power than lam* under those
conditions, H23b is dead.

RESULT (2026-08-02 run; see glv_hnp_phase23_block_predictor_output.txt):
  H23b SURVIVES the falsifier, cleanly.  14 j=0 GLV curves near 2^16, m=12,
  5 seeds, at eff = 0.098 and eff = 0.149 (28 cells, eff/gamma/det B/Eplant
  all held constant within each group):

      classifier                      eff=0.098   eff=0.149   pooled
      majority baseline                 0.714       0.786      0.750
      best threshold on nu              1.000       1.000      1.000
      best threshold on skew=mu2/mu1    1.000       1.000      1.000
      best threshold on lam*            0.714       0.786      0.750   (= baseline)
      parameter-free  m*(nu,n) <= m     0.929       0.786      0.857

      Spearman rho vs #successes:  nu -0.698,  skew +0.497,  lam* -0.119

  So the curve-level predictor is the SHAPE of the 2-D lambda-block, and
  lam*/n carries no information at all once eff is fixed.  Note the curve
  p=65713, lam*=0.0068 -- the one Thread 20 T3 used to falsify the lam*
  threshold -- succeeds 5/5 at eff=0.098 (skew 3.57) and fails 0/5 at
  eff=0.149 (skew 2.33), exactly as skew predicts and with no reference to
  lam* at all.  This closes the run of failed curve-level separators logged
  2026-06-21..2026-06-29: every one of them tested lam, mu1, or a global
  lattice invariant.  The right invariant is the RATIO mu2/mu1.

Run: python3 glv_hnp_phase23_block_predictor.py
"""

import math

from glv_hnp_phase23_gh_law import (
    find_generator, lam_star, scales, rate, gamma_orig, search_curves,
)
from glv_hnp_phase23_coset_exact import Block2D

SEEDS = [42, 1234, 9999, 555, 31337]
M_SIGS = 12


def block_stats(n, lam, k1, k2):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
    B = Block2D(n, lam, S_K1, S_K2)
    erand = B.mean_random_dist2(trials=3000)
    eplant = (k1 * S_K1) ** 2 / 3.0 + (k2 * S_K2) ** 2 / 3.0
    nu = math.sqrt(eplant / erand) if erand > 0 else float('inf')
    if nu < 1:
        m_star = math.ceil(2 * math.log(n) / math.log(1.0 / (nu * nu)))
    else:
        m_star = -1
    return B, erand, eplant, nu, m_star


def best_acc(rows, key, greater_is_success):
    vals = sorted({r[key] for r in rows})
    best = (0.0, None)
    for i in range(len(vals) + 1):
        thr = vals[0] - 1 if i == 0 else vals[i - 1]
        c = sum(1 for r in rows
                if (r[key] > thr if greater_is_success else r[key] <= thr) == r['ok'])
        if c / len(rows) > best[0]:
            best = (c / len(rows), thr)
    return best


if __name__ == '__main__':

    print("=" * 78)
    print("Thread 23 part 3 -- the lambda-block shape as curve-level predictor")
    print("=" * 78)

    # --- W1: the two historical 12-bit curves, same n-size, same eff -------
    print("\n" + "-" * 78)
    print("W1: the 2557 / 2677 discrepancy explained")
    print("-" * 78)
    HIST = [("8-bit/199", 211, 2, 199, 106),
            ("12-bit/2557", 2557, 2, 2659, 1755),
            ("12-bit/2677", 2677, 2, 2647, 185)]
    hist = []
    for label, p, b, n, lam in HIST:
        G = find_generator(p, b, n)
        assert G is not None and (lam * lam + lam + 1) % n == 0, label
        hist.append((label, (p, b, n, lam, G)))
    print(f"{'curve':>13} {'lam*':>6} {'K1':>4} {'eff':>6} {'mu1':>10} {'mu2':>10} "
          f"{'skew':>6} {'nu':>7} {'m*':>5} {'LLL(m=12)':>10}")
    for label, curve in hist:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        for k1 in (2, 4, 6, 8, 12, 16):
            B, er, ep, nu, ms = block_stats(n, lam, k1, k2)
            w, t = rate(curve, M_SIGS, k1, SEEDS)
            print(f"{label:>13} {lam_star(lam,n):>6.3f} {k1:>4} {k1*k2/n:>6.3f} "
                  f"{B.mu1:>10.4g} {B.mu2:>10.4g} {B.mu2/B.mu1:>6.2f} {nu:>7.3f} "
                  f"{ms:>5} {str(w)+'/'+str(t):>10}")

    # --- W2: 17-bit sweep at FIXED eff -- the decisive test ----------------
    print("\n" + "-" * 78)
    print("W2: 17-bit curves at FIXED eff (det B, Eplant, gamma all constant)")
    print("-" * 78)
    curves = search_curves(1 << 16, 1 << 17, 14)
    print(f"{len(curves)} curves")
    rows = []
    for eff_t in (0.098, 0.149):
        print(f"\n=== eff ~ {eff_t}   (m={M_SIGS}, {len(SEEDS)} seeds) ===")
        print(f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>4} {'eff':>6} {'gamma':>7} "
              f"{'skew':>7} {'nu':>7} {'m*':>5} {'LLL':>6}")
        for (p, b, n, lam, G) in curves:
            k2 = math.isqrt(n) + 1
            k1 = max(2, int(eff_t * n / k2))
            eff = k1 * k2 / n
            B, er, ep, nu, ms = block_stats(n, lam, k1, k2)
            g = gamma_orig(M_SIGS, n, k1, k2)
            w, t = rate((p, b, n, lam, G), M_SIGS, k1, SEEDS)
            rows.append({'p': p, 'n': n, 'lam_star': lam_star(lam, n),
                         'eff': eff, 'gamma': g, 'skew': B.mu2 / B.mu1,
                         'nu': nu, 'm_star': ms, 'ok': w == t, 'w': w})
            print(f"{p:>7} {n:>7} {lam_star(lam,n):>6.3f} {k1:>4} {eff:>6.3f} "
                  f"{g:>7.3f} {B.mu2/B.mu1:>7.2f} {nu:>7.3f} {ms:>5} "
                  f"{str(w)+'/'+str(t):>6}")

    print("\n" + "-" * 78)
    print("classifier accuracy at FIXED eff (eff and gamma are constants here,")
    print("so they carry no information -- only curve-level statistics can)")
    print("-" * 78)
    for eff_t, sub in (("0.098", rows[:len(rows) // 2]),
                       ("0.149", rows[len(rows) // 2:]),
                       ("both", rows)):
        if not sub:
            continue
        base = max(sum(1 for r in sub if r['ok']),
                   sum(1 for r in sub if not r['ok'])) / len(sub)
        print(f"\n  eff={eff_t}  n={len(sub)} cells   majority baseline {base:.3f}")
        for key, gis in (('nu', False), ('skew', True), ('lam_star', True)):
            a, thr = best_acc(sub, key, gis)
            print(f"    best threshold on {key:<9}: {a:.3f}  (thr={thr:.4f})")
        pf = sum(1 for r in sub
                 if (r['m_star'] > 0 and r['m_star'] <= M_SIGS) == r['ok']) / len(sub)
        print(f"    PARAMETER-FREE m* <= m  : {pf:.3f}")

    # rank correlation nu vs success count
    def spearman(xs, ys):
        def rk(v):
            s = sorted(range(len(v)), key=lambda i: v[i])
            r = [0.0] * len(v)
            for pos, i in enumerate(s): r[i] = pos
            return r
        a, bb = rk(xs), rk(ys)
        n_ = len(xs)
        ma, mb = sum(a) / n_, sum(bb) / n_
        num = sum((a[i] - ma) * (bb[i] - mb) for i in range(n_))
        da = math.sqrt(sum((x - ma) ** 2 for x in a))
        db = math.sqrt(sum((x - mb) ** 2 for x in bb))
        return num / (da * db) if da and db else 0.0

    print("\n  Spearman rho vs LLL success count (5 seeds), pooled:")
    for key in ('nu', 'skew', 'lam_star'):
        print(f"    {key:<9}: {spearman([r[key] for r in rows], [r['w'] for r in rows]):+.3f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)

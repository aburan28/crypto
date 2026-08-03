"""
GLV-HNP Phase 2, Thread 23 addendum: what IS lambda_1 of the projected lattice?

E1-E4 (glv_hnp_thread23_projected.py) established:
  * the projected lattice P is the exact quotient of the original lattice O by
    the trivial direction n*e_m (HNF-verified), and
  * removing that direction changes NOTHING: O and P agree on 16/16 K1-wall
    cells and 36/36 eff-sweep cells, and
  * sv/pv still sits BELOW 1 in P (0.475 - 0.790).

So a second, non-trivial rival vector is shorter than the planted vector.
This script identifies it.

Hypothesis H23a: the rival is a "lambda-block" vector, i.e. it lives in one of
the m two-dimensional sublattices

    L2_i = < (n*S_K1, 0), (-lam*S_K1, S_K2) >   in coords (i, m+i)

whose exact shortest vector mu_i is computable in O(log n) by Lagrange-Gauss
reduction.  These are the same blocks whose normalised invariant
nu_hat = lambda_1(L2)/sqrt(det L2) was found on 2026-07-29 (commit e845207) to
separate the June C1/C2 classes with AUC 0.935.  If H23a holds, then

    lambda_1(P) = min_i mu_i

and the ratio  min_i mu_i / ||v_planted||  is the exact, closed-form,
O(m log n) predictor of whether the planted vector can ever be lambda_1 --
computable without running LLL at all.

Falsifier for H23a: if the LLL shortest vector has energy outside a single
(i, m+i) coordinate pair, or if ||sv|| < min_i mu_i (a shorter rival exists
that is not a lambda-block vector), H23a is wrong.

Run:  python3 secp256k1_cm_audit/glv_hnp_thread23_sv_anatomy.py
"""

import math
import random
import sys

from fpylll import IntegerMatrix, LLL
import sympy

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from glv_hnp_thread23_projected import (          # noqa: E402
    find_generator, build_curve, eisenstein_decompose, j0_traces, glv_roots,
    lam_star, scales, gen_signatures, norm,
    build_projected_lattice, planted_vector_projected,
)

# ---------------------------------------------------------------------------
# Exact 2D Lagrange-Gauss reduction (verbatim from glv_hnp_phase2_lambda_threshold.py)
# ---------------------------------------------------------------------------

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

def lambda_block_mu(n, lam, S_K1, S_K2):
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] * w[0] + w[1] * w[1]), w


def anatomy(curve, m, k1_bound, seed):
    """LLL the projected lattice; dissect its shortest vector."""
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2b)
    d = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d, m, n, lam, p, k1_bound, k2b, seed)
    if len(sigs) < m:
        return None
    M, _ = build_projected_lattice(sigs, n, lam, k1_bound, k2b)
    A = IntegerMatrix.from_matrix([list(r) for r in M])
    LLL.reduction(A)
    red = [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]
    rows = [r for r in red if any(r)]
    sv = min(rows, key=norm)
    pv = planted_vector_projected(sigs, n, k1_bound, k2b)

    e_k1 = sum(x * x for x in sv[:m])
    e_k2 = sum(x * x for x in sv[m:2 * m])
    e_kan = sv[2 * m] ** 2
    tot = e_k1 + e_k2 + e_kan

    # is sv supported on a single (i, m+i) coordinate pair?
    live = [i for i in range(m) if sv[i] or sv[m + i]]
    single_block = (len(live) == 1 and sv[2 * m] == 0)

    mu, _w = lambda_block_mu(n, lam, S_K1, S_K2)
    return {
        'n': n, 'lam': lam, 'lam*': lam_star(lam, n), 'm': m, 'K1': k1_bound,
        'sv': norm(sv), 'pv': norm(pv), 'mu': mu,
        'frac_k1': e_k1 / tot, 'frac_k2': e_k2 / tot, 'frac_kan': e_kan / tot,
        'nblocks': len(live), 'single_block': single_block,
        'kan0': sv[2 * m] == 0,
        'sv_over_mu': norm(sv) / mu if mu else float('nan'),
        'mu_over_pv': mu / norm(pv),
    }


print("=" * 92)
print("Thread 23 addendum -- anatomy of lambda_1 in the projected lattice P")
print("=" * 92)

HIST = [
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]
cases = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None
    cases.append((label, (p, b, n, lam, G), k1, m))

print("\n" + "-" * 92)
print("A1: is the LLL shortest vector a single lambda-block vector?")
print("-" * 92)
print(f"{'curve':<20} {'K1':>3} {'seed':>6} {'sv/pv':>7} {'mu/pv':>7} {'sv/mu':>7} "
      f"{'#blocks':>8} {'kan=0':>6} {'E_k1':>6} {'E_k2':>6}")
agree = tot = 0
for label, curve, k1, m in cases:
    for seed in (42, 1234, 9999):
        r = anatomy(curve, m, k1, seed)
        tot += 1
        agree += bool(r['single_block'] and abs(r['sv_over_mu'] - 1.0) < 1e-9)
        print(f"{label:<20} {k1:>3} {seed:>6} {r['sv']/r['pv']:>7.3f} "
              f"{r['mu_over_pv']:>7.3f} {r['sv_over_mu']:>7.4f} "
              f"{r['nblocks']:>8} {str(r['kan0']):>6} "
              f"{r['frac_k1']:>6.3f} {r['frac_k2']:>6.3f}")
print(f"\nH23a holds (sv is exactly one lambda-block vector) in {agree}/{tot} cases.")

print("\n" + "-" * 92)
print("A2: does mu/pv predict recovery?  Sweep the K1 grid on both 12-bit curves.")
print("    (recovery data from E3: identical for O and P)")
print("-" * 92)

E3_TRUTH = {
    ("12-bit/2557", 2): 5, ("12-bit/2557", 3): 5, ("12-bit/2557", 4): 5,
    ("12-bit/2557", 6): 5, ("12-bit/2557", 8): 5, ("12-bit/2557", 12): 1,
    ("12-bit/2557", 16): 0, ("12-bit/2557", 24): 0,
    ("12-bit/2677 FAIL", 2): 5, ("12-bit/2677 FAIL", 3): 5,
    ("12-bit/2677 FAIL", 4): 5, ("12-bit/2677 FAIL", 6): 0,
    ("12-bit/2677 FAIL", 8): 0, ("12-bit/2677 FAIL", 12): 0,
    ("12-bit/2677 FAIL", 16): 0, ("12-bit/2677 FAIL", 24): 0,
}
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
print(f"{'curve':<20} {'K1':>3} {'mu/pv':>7} {'sv/pv':>7} {'wins(E3)':>9}")
data = []
for label, curve, _k1, m in cases:
    if not label.startswith("12-bit"):
        continue
    for k1 in K1_GRID:
        r = anatomy(curve, m, k1, 42)
        w = E3_TRUTH[(label, k1)]
        data.append((label, k1, r['mu_over_pv'], r['sv'] / r['pv'], w))
        print(f"{label:<20} {k1:>3} {r['mu_over_pv']:>7.3f} "
              f"{r['sv']/r['pv']:>7.3f} {w:>9}")

succ = [x[2] for x in data if x[4] == 5]
fail = [x[2] for x in data if x[4] == 0]
print(f"\n  mu/pv over full successes (5/5): min={min(succ):.3f} max={max(succ):.3f}")
print(f"  mu/pv over full failures  (0/5): min={min(fail):.3f} max={max(fail):.3f}")
sep = max(succ) < min(fail) or min(succ) > max(fail)
print(f"  separable by a single threshold: {sep}")
if sep:
    lo, hi = (max(succ), min(fail)) if max(succ) < min(fail) else (max(fail), min(succ))
    print(f"  threshold band: ({lo:.3f}, {hi:.3f})")

print("\n" + "-" * 92)
print("A3: out-of-sample check on fresh 17-bit curves at eff=0.15 (E4 truth).")
print("-" * 92)

def search_curves(lo, hi, want=12, nbins=6):
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_c = p + 1 - t
                    if n_c < 2 or n_c % 3 != 1 or not sympy.isprime(n_c):
                        continue
                    rr = glv_roots(n_c)
                    if rr is None:
                        continue
                    idx = min(nbins - 1, int(lam_star(rr[0], n_c) / (0.5 / nbins)))
                    if len(bins[idx]) >= want // nbins:
                        continue
                    cur = build_curve(p, n_c)
                    if cur is None:
                        continue
                    bins[idx].append((p, cur[1], n_c, rr[0], cur[4]))
        if all(len(v) >= want // nbins for v in bins.values()):
            break
        p = int(sympy.nextprime(p))
    out = []
    for i in range(nbins):
        out.extend(bins[i])
    return out

# E4 truth at eff=0.15, m=12, seeds [42,1234,9999,555,31337]  (O column)
E4_TRUTH_015 = {
    65719: 1, 65269: 0, 65203: 0, 66751: 0, 66037: 0, 65647: 0,
    65899: 5, 65167: 5, 65287: 5, 65053: 0, 66943: 0, 67369: 5,
}
curves17 = search_curves(2 ** 16, 2 ** 17, want=12, nbins=6)
print(f"{'p':>7} {'n':>7} {'lam*':>7} {'K1':>4} {'mu/pv':>7} {'sv/pv':>7} {'wins(E4)':>9}")
rows = []
for cur in curves17:
    p, b, n, lam, G = cur
    if n not in E4_TRUTH_015:
        continue
    k2b = math.isqrt(n) + 1
    k1b = max(2, int(0.15 * n / k2b))
    r = anatomy(cur, 12, k1b, 42)
    w = E4_TRUTH_015[n]
    rows.append((n, r['mu_over_pv'], r['sv'] / r['pv'], w))
    print(f"{p:>7} {n:>7} {lam_star(lam, n):>7.4f} {k1b:>4} "
          f"{r['mu_over_pv']:>7.3f} {r['sv']/r['pv']:>7.3f} {w:>9}")

s = [x[1] for x in rows if x[3] == 5]
f = [x[1] for x in rows if x[3] == 0]
if s and f:
    print(f"\n  mu/pv successes: min={min(s):.3f} max={max(s):.3f}  (n={len(s)})")
    print(f"  mu/pv failures : min={min(f):.3f} max={max(f):.3f}  (n={len(f)})")
    # best single threshold, BOTH directions.  The 2026-07-29 nu_hat result
    # (commit e845207) found the sign to be opposite to the naive guess -- a
    # SHORT rival makes the attack easier -- so testing only ">=" would have
    # reported a spurious null.
    cands = sorted(set(x[1] for x in rows))
    best = (0.0, None, None)
    for c in cands:
        for direction in ('>=', '<='):
            pred = (lambda v, c=c: v >= c) if direction == '>=' else (lambda v, c=c: v <= c)
            acc = sum((pred(x[1]) == (x[3] == 5)) for x in rows) / len(rows)
            if acc > best[0]:
                best = (acc, c, direction)
    base = max(len(s), len(f)) / len(rows)
    print(f"  best threshold  mu/pv {best[2]} {best[1]:.3f}: accuracy {best[0]:.3f} "
          f"(majority baseline {base:.3f})")
    print(f"  direction: {'LOW mu/pv predicts success' if best[2] == '<=' else 'HIGH mu/pv predicts success'}"
          f"  -- 2026-07-29 nu_hat sign was 'short rival => easier'")

    # Is lambda_1(P) actually attained by the lambda-block vector at m=12?
    tight = [x for x in rows if abs(x[2] - x[1]) < 1e-9]
    print(f"\n  lambda_1(P) == mu (single lambda-block) in {len(tight)}/{len(rows)} "
          f"of the m=12 curves")
    loose = [x for x in rows if abs(x[2] - x[1]) >= 1e-9]
    if loose:
        print("  cases where a MULTI-block combination beats mu:")
        for n_, mp, sp, w in loose:
            print(f"    n={n_}  mu/pv={mp:.3f}  sv/pv={sp:.3f}  wins={w}")
        print(f"  all of them have mu/pv >= {min(x[1] for x in loose):.3f} "
              f"-- mu is the binding constraint only while it is small")

print("\n" + "=" * 92)
print("done")
print("=" * 92)


print("\n" + "-" * 92)
print("A4: pooled model across all three slices (12-bit m=8, 12-bit m=10, 17-bit m=12).")
print("    Two regimes, only ONE fitted parameter:")
print("      R_A  mu/pv > 1   <=>  the planted vector IS lambda_1(P)  ->  predict success")
print("      R_B  mu/pv <= 1  ->   success iff nu_hat = mu/sqrt(det L2) <= tau")
print("    R_A is derived, not fitted.  tau is the single fitted parameter.")
print("    Direction of tau was pre-registered by the 2026-07-29 nu_hat result (e845207).")
print("-" * 92)

POOLED = [
    # (label, n, lam, m, K1, wins/5)     -- A2 rows
    ("12b/2659", 2659, 1755, 8, 2, 5), ("12b/2659", 2659, 1755, 8, 3, 5),
    ("12b/2659", 2659, 1755, 8, 4, 5), ("12b/2659", 2659, 1755, 8, 6, 5),
    ("12b/2659", 2659, 1755, 8, 8, 5), ("12b/2659", 2659, 1755, 8, 12, 1),
    ("12b/2659", 2659, 1755, 8, 16, 0), ("12b/2659", 2659, 1755, 8, 24, 0),
    ("12b/2647", 2647, 185, 10, 2, 5), ("12b/2647", 2647, 185, 10, 3, 5),
    ("12b/2647", 2647, 185, 10, 4, 5), ("12b/2647", 2647, 185, 10, 6, 0),
    ("12b/2647", 2647, 185, 10, 8, 0), ("12b/2647", 2647, 185, 10, 12, 0),
    ("12b/2647", 2647, 185, 10, 16, 0), ("12b/2647", 2647, 185, 10, 24, 0),
]
for n_, w in [(65719, 1), (65269, 0), (65203, 0), (66751, 0), (66037, 0), (65647, 0),
              (65899, 5), (65167, 5), (65287, 5), (65053, 0), (66943, 0), (67369, 5)]:
    lam_ = glv_roots(n_)[0]
    k2b_ = math.isqrt(n_) + 1
    POOLED.append((f"17b/{n_}", n_, lam_, 12, max(2, int(0.15 * n_ / k2b_)), w))

feat = []
print(f"{'curve':<12} {'m':>3} {'K1':>3} {'mu/pv':>7} {'nu_hat':>8} {'regime':>7} {'wins':>5}")
for lbl, n_, lam_, m_, K1_, w in POOLED:
    k2b_ = math.isqrt(n_) + 1
    S_K1_, _, S_K2_, _ = scales(n_, K1_, k2b_)
    mu_, _ = lambda_block_mu(n_, lam_, S_K1_, S_K2_)
    nuhat_ = mu_ / math.sqrt(n_ * S_K1_ * S_K2_)       # det L2 = n*S_K1*S_K2
    pv_ = n_ * math.sqrt(2 * m_ / 3 + 1)               # E||v_planted|| in P
    mp_ = mu_ / pv_
    feat.append((lbl, mp_, nuhat_, w))
    print(f"{lbl:<12} {m_:>3} {K1_:>3} {mp_:>7.3f} {nuhat_:>8.4f} "
          f"{'A' if mp_ > 1 else 'B':>7} {w:>5}")

def score(tau):
    hit = 0
    for _lbl, mp_, nu_, w in feat:
        pred = (mp_ > 1.0) or (nu_ <= tau)
        hit += (pred == (w == 5))
    return hit

taus = sorted(set(x[2] for x in feat))
best_tau, best_hit = max(((t, score(t)) for t in taus), key=lambda z: z[1])
nsucc = sum(1 for x in feat if x[3] == 5)
base = max(nsucc, len(feat) - nsucc) / len(feat)
print(f"\n  best tau = {best_tau:.4f}:  {best_hit}/{len(feat)} = {best_hit/len(feat):.3f} "
      f"(majority baseline {base:.3f})")

# nu_hat alone, for comparison
def score_nu_only(tau):
    return sum(((nu_ <= tau) == (w == 5)) for _l, _m, nu_, w in feat)
bt2, bh2 = max(((t, score_nu_only(t)) for t in taus), key=lambda z: z[1])
print(f"  nu_hat alone (no regime A): best tau = {bt2:.4f}: {bh2}/{len(feat)} = "
      f"{bh2/len(feat):.3f}")
print(f"  regime A alone (mu/pv > 1): "
      f"{sum(((mp_ > 1) == (w == 5)) for _l, mp_, _n, w in feat)}/{len(feat)}")

print("\n  misclassified by the two-regime rule at best tau:")
for lbl, mp_, nu_, w in feat:
    pred = (mp_ > 1.0) or (nu_ <= best_tau)
    if pred != (w == 5):
        print(f"    {lbl:<12} mu/pv={mp_:.3f} nu_hat={nu_:.4f} "
              f"predicted={'succ' if pred else 'fail'} actual={w}/5")
print("\n  CAVEAT: tau is fitted in-sample on these 28 configs.  The DIRECTION and")
print("  regime A are not fitted.  Out-of-sample validation is the next step.")

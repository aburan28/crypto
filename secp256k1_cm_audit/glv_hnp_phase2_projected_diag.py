"""
Thread 23, follow-up diagnostic.

`glv_hnp_phase2_projected.py` showed that projecting out the d-column (arm P)
removes the trivial shortest vector n*S_D*e_m and moves sv/pv from ~0.35 to
0.43-0.86, yet the recovery COUNTS are identical to the original lattice at
every eff and every K1.  Two questions follow:

  D1  Is the agreement instance-level, or only aggregate?  (Identical counts
      could be coincidence; per-(curve, seed) agreement could not.)

  D2  If the d-vector is gone, what is the remaining rival that keeps the
      planted vector off lambda_1?  Hypothesis: the m independent copies of the
      2D lambda-block  L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >  sitting in
      column pairs (i, m+i).  These are the same L2 whose shortest vector mu
      defines the nu_hat separator found on 2026-07-29 (AUC 0.935).
      Prediction: sv(arm P) == mu exactly whenever mu < ||v_planted||, and the
      shortest vector's support is a single column pair.

Run: python3 glv_hnp_phase2_projected_diag.py
"""

import math
import random

import sympy

import glv_hnp_phase2_projected as T23
from glv_hnp_phase2_projected import (
    scales, gen_signatures, norm, lam_star, find_generator, build_curve,
    eisenstein_decompose, j0_traces, glv_roots,
    build_lattice_O, planted_O, recover_O,
    build_lattice_P, planted_P, recover_P, lll_rows,
)

# ---------------------------------------------------------------------------
# Exact 2D Gauss/Lagrange reduction (verbatim from
# glv_hnp_phase2_lambda_threshold.py:188 so mu matches the 2026-07-29 numbers)
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

def nu_hat(n, lam, S_K1, S_K2):
    """mu / sqrt(det L2), the 2026-07-29 separator.  det L2 = n*S_K1*S_K2."""
    mu, _ = lambda_block_mu(n, lam, S_K1, S_K2)
    return mu / math.sqrt(n * S_K1 * S_K2)

# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 555, 31337]
M_SIGS = 12

print("=" * 78)
print("Thread 23 diagnostic -- what still blocks lambda_1 after projection?")
print("=" * 78)

HIST = [
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]
hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    hist.append((label, (p, b, n, lam, G), k1, m))

def search_curves(lo, hi, want=10):
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    nc = p + 1 - t
                    if nc < 2 or nc % 3 != 1 or not sympy.isprime(nc):
                        continue
                    roots = glv_roots(nc)
                    if roots is None:
                        continue
                    cur = build_curve(p, nc)
                    if cur is None:
                        continue
                    out.append((p, cur[1], nc, roots[0], cur[4]))
                    break
        p = int(sympy.nextprime(p))
    return out

CURVES17 = search_curves(2**16, 2**17, want=10)
print(f"\n{len(CURVES17)} fresh 17-bit j=0 GLV curves "
      f"(same search as EXP C, so the same curves)")

# ===========================================================================
print("\n" + "-" * 78)
print("D1 -- instance-level agreement between arm O and arm P")
print("-" * 78)

agree = disagree = 0
disagreements = []
both_ok = neither = 0
for eff in (0.05, 0.10, 0.15, 0.20, 0.25, 0.35):
    for (p, b, n, lam, G) in CURVES17:
        k2b = math.isqrt(n) + 1
        k1b = max(2, int(eff * n / k2b))
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            res = T23.trial((p, b, n, lam, G), M_SIGS, d_trial, k1b, seed,
                            arms=('O', 'P'))
            if res is None:
                continue
            o, pp = res['O'][0], res['P'][0]
            if o == pp:
                agree += 1
                both_ok += (o and pp)
                neither += (not o and not pp)
            else:
                disagree += 1
                disagreements.append((eff, n, seed, o, pp))

tot = agree + disagree
print(f"trials: {tot}   agree: {agree} ({100.0*agree/tot:.2f}%)   "
      f"disagree: {disagree}")
print(f"  both recovered: {both_ok}   neither recovered: {neither}")
if disagreements:
    print("  disagreements (eff, n, seed, O, P):")
    for row in disagreements[:20]:
        print("   ", row)
else:
    print("  -> arm P succeeds on EXACTLY the instance set arm O succeeds on.")

# ===========================================================================
print("\n" + "-" * 78)
print("D2 -- identity of the arm-P shortest vector")
print("-" * 78)
print("For each instance: is sv(P) equal to mu = lambda_1(L2), and is its")
print("support a single column pair (i, m+i)?")
print()
print(f"{'n':>7} {'eff':>5} {'sd':>6} | {'sv/pv':>7} {'mu/pv':>7} "
      f"{'sv==mu':>7} {'1-pair':>7} {'kan=0':>6} {'nu_hat':>7} {'rec':>4}")
print("-" * 78)

stats = {'sv_eq_mu': 0, 'one_pair': 0, 'kan_zero': 0, 'n': 0,
         'planted_is_sv': 0}
rows_for_pred = []

for eff in (0.05, 0.15, 0.25):
    for (p, b, n, lam, G) in CURVES17[:5]:
        k2b = math.isqrt(n) + 1
        k1b = max(2, int(eff * n / k2b))
        S_K1, S_D, S_K2, S_KAN = scales(n, k1b, k2b)
        mu, wvec = lambda_block_mu(n, lam, S_K1, S_K2)
        nh = nu_hat(n, lam, S_K1, S_K2)
        for seed in SEEDS[:3]:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, d_trial, M_SIGS, n, lam, p, k1b, k2b, seed)
            if len(sigs) < M_SIGS:
                continue
            m = M_SIGS
            rows = build_lattice_P(sigs, n, lam, k1b, k2b)
            red = lll_rows(rows, 2 * m + 1)
            nzrows = [r for r in red if any(r)]
            sv_row = min(nzrows, key=norm)
            sv = norm(sv_row)
            pv = norm(planted_P(sigs, n, k1b, k2b))
            rec = recover_P(red, sigs, n, lam, k1b, k2b, d_trial) is not None

            sv_eq_mu = abs(sv - mu) < 1e-6 * max(sv, mu)
            # support test: nonzero only in columns i and m+i for one i
            nzc = [j for j, x in enumerate(sv_row) if x != 0]
            pairs = {j % m for j in nzc if j < 2 * m}
            one_pair = (len(pairs) == 1 and all(j < 2 * m for j in nzc))
            kan_zero = (sv_row[2 * m] == 0)

            stats['n'] += 1
            stats['sv_eq_mu'] += sv_eq_mu
            stats['one_pair'] += one_pair
            stats['kan_zero'] += kan_zero
            stats['planted_is_sv'] += (abs(sv - pv) < 1e-9 * pv)
            rows_for_pred.append((eff, n, mu, pv, nh, rec))

            print(f"{n:>7} {eff:>5.2f} {seed:>6} | {sv/pv:>7.4f} {mu/pv:>7.4f} "
                  f"{str(sv_eq_mu):>7} {str(one_pair):>7} {str(kan_zero):>6} "
                  f"{nh:>7.4f} {str(rec):>4}")

N = stats['n']
print()
print(f"summary over {N} instances:")
print(f"  sv(P) == mu exactly          : {stats['sv_eq_mu']}/{N}")
print(f"  sv(P) supported on one pair  : {stats['one_pair']}/{N}")
print(f"  sv(P) has Kannan coord 0     : {stats['kan_zero']}/{N}")
print(f"  planted vector IS lambda_1   : {stats['planted_is_sv']}/{N}")

# ===========================================================================
print("\n" + "-" * 78)
print("D3 -- does mu/pv predict recovery?  STRATIFIED BY eff")
print("-" * 78)
print("The unstratified split is confounded: at eff=0.05 essentially every")
print("instance recovers regardless of mu, and mu/pv > 1 happens to be common")
print("there because K1 is small.  Within a fixed eff the direction reverses.")
print()

def pct(rs):
    return f"{sum(1 for r in rs if r[-1]):>2}/{len(rs):<2}" if rs else " 0/0 "

print("  unstratified (CONFOUNDED -- reported only to be retracted):")
gt = [r for r in rows_for_pred if r[2] > r[3]]
lt = [r for r in rows_for_pred if r[2] <= r[3]]
print(f"    mu >  pv : {len(gt):>3} instances, recovered {pct(gt)}")
print(f"    mu <= pv : {len(lt):>3} instances, recovered {pct(lt)}")
print()
print("  stratified by eff:")
print(f"    {'eff':>5} | {'lo mu/pv half':>16} {'hi mu/pv half':>16} | "
      f"{'mean mu/pv rec':>15} {'mean mu/pv fail':>16}")
for eff in (0.05, 0.15, 0.25):
    rs = [r for r in rows_for_pred if r[0] == eff]
    if not rs:
        continue
    rs_sorted = sorted(rs, key=lambda r: r[2] / r[3])
    h = len(rs_sorted) // 2
    lo, hi = rs_sorted[:h], rs_sorted[h:]
    rec = [r[2] / r[3] for r in rs if r[-1]]
    fail = [r[2] / r[3] for r in rs if not r[-1]]
    mrec = sum(rec) / len(rec) if rec else float('nan')
    mfail = sum(fail) / len(fail) if fail else float('nan')
    print(f"    {eff:>5.2f} | {pct(lo):>16} {pct(hi):>16} | "
          f"{mrec:>15.4f} {mfail:>16.4f}")
print()
print("  same, using nu_hat = mu/sqrt(det L2) (the 2026-07-29 separator):")
print(f"    {'eff':>5} | {'lo nu_hat half':>16} {'hi nu_hat half':>16} | "
      f"{'mean nu rec':>15} {'mean nu fail':>16}")
for eff in (0.05, 0.15, 0.25):
    rs = [r for r in rows_for_pred if r[0] == eff]
    if not rs:
        continue
    rs_sorted = sorted(rs, key=lambda r: r[4])
    h = len(rs_sorted) // 2
    lo, hi = rs_sorted[:h], rs_sorted[h:]
    rec = [r[4] for r in rs if r[-1]]
    fail = [r[4] for r in rs if not r[-1]]
    mrec = sum(rec) / len(rec) if rec else float('nan')
    mfail = sum(fail) / len(fail) if fail else float('nan')
    print(f"    {eff:>5.2f} | {pct(lo):>16} {pct(hi):>16} | "
          f"{mrec:>15.4f} {mfail:>16.4f}")
print()
print("Note: mu/pv is a function of (n, lam, K1, K2, m) only -- it does not")
print("depend on the signatures or on d.  ||v_planted|| ~ n*sqrt(2m/3+1) grows")
print("with m while mu does not, so mu/pv DECREASES with m: adding signatures")
print("makes the planted vector relatively longer, never shorter.")

print("\n" + "=" * 78)
print("done")
print("=" * 78)

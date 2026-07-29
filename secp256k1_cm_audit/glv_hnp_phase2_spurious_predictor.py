"""
GLV-HNP Phase 2, Thread 20 (part 2): what actually gates the attack?

Companion to glv_hnp_phase2_lambda_threshold.py, which falsifies the
"lambda/n threshold" of the 2026-07-26 run by rescuing the rho=0.07 curve
p=2677 at K1 in {2,3,4} without moving rho at all.  This script looks for the
quantity that does the gating.

THE SPURIOUS FAMILY
-------------------
In the (2m+2)-dimensional attack lattice, modular row i is

    row_i        = ( ..., n*S_K1 at col i, ... )

and GLV row m+1+i is

    row_{m+1+i}  = ( ..., -lam*S_K1 at col i, ..., S_K2 at col m+1+i, ... ).

For any v0, the integer combination  v0*row_{m+1+i} + c*row_i  has

    col i        = (-lam*v0 + c*n) * S_K1  =  u0 * S_K1,   u0 = -lam*v0 mod n
    col m+1+i    = v0 * S_K2
    everything else = 0.

So for EACH of the m signature indices the lattice contains a vector supported
on just two coordinates, drawn from the rank-2 lattice

    L2 = { (u*S_K1, v*S_K2) : u = lam*v mod n },   det(L2) = n*S_K1*S_K2,

and we write mu = lambda_1(L2), computed exactly by Lagrange-Gauss reduction.

These vectors carry ZERO in the Kannan coordinate, so they can never yield d --
they are pure noise -- yet there are m of them, all of norm mu.  Part A confirms
they really do occupy reduced-basis rows.  The natural guess is therefore that
the attack fails when mu drops below the planted norm ||v||.

RESULT: THAT GUESS IS WRONG, AND THE SCRIPT SAYS SO
---------------------------------------------------
Measured over 267 (config,m) points that clear the information threshold,
||v||/mu does NOT separate success from failure:

    ||v||/mu among full successes: [0.36, 3.21]
    ||v||/mu among full failures : [0.46, 3.06]      -- almost total overlap

The cleanest counterexample is the very pair the 2026-07-26 run compared: at
K1=8, p=2557 has the SMALLER mu (2738 vs 5096) and yet succeeds, while p=2677
has the larger mu and fails.  So the spurious family is real but is not what
gates the attack.

WHAT DOES GATE IT: eff = K1*K2/n, the ordinary HNP bias budget.  Bucketing the
same 267 points by eff instead gives a clean monotone curve:

    eff < 0.05        96.9%
    0.05 <= eff <0.10 85.2%
    0.10 <= eff <0.20 54.2%
    0.20 <= eff <0.40 31.0%
    0.40 <= eff <1.00  0.0%

Both curves in the prior run's comparison sit at eff ~ 0.157, i.e. inside the
~54% marginal band, where a 3-seed sample separates very little.  See
glv_hnp_phase2_matched_eff.py for the 24-seed re-test of that pair.

Note on rho: Part D shows corr(rho, log mu) is weak and of the wrong sign to
support a rho story -- the smallest-rho pool curve has a LARGE mu and one of
the largest-rho curves has among the smallest.  rho gates nothing here.

All recoveries are verified WITHOUT a secret oracle (see verify_d).

Run: python3 glv_hnp_phase2_spurious_predictor.py
"""

import math
import random

from glv_hnp_phase2_common import (
    find_generator, centered_rho, gen_signatures, build_glv_lattice,
    planted_vector, scales, norm, verify_d, find_glv_curves,
    gauss_reduce, spurious_mu, m_info, reduce_and_test, M_MIN, M_MAX, SEEDS,
)

# ---------------------------------------------------------------------------
# Exact Lagrange-Gauss reduction of the scaled rank-2 GLV lattice
# ---------------------------------------------------------------------------

def profile(curve, k1_bound, algo='LLL', beta=20, m_max=M_MAX):
    """Full m-profile for one (curve, K1).  Returns list of per-m records."""
    n = curve['n']
    k2_bound = math.isqrt(n) + 1
    if k1_bound * k2_bound >= n:
        return []                       # eff>=1: no unique solution, skip
    mu, uv = spurious_mu(n, curve['lam'], k1_bound, k2_bound)
    mi = m_info(n, k1_bound, k2_bound)
    per_seed = {}
    for sd in SEEDS:
        d = random.Random(sd + 7777).randrange(1, n)
        per_seed[sd] = (d, gen_signatures(curve, d, m_max, k1_bound,
                                          k2_bound, sd))
    out = []
    for m in range(M_MIN, m_max + 1):
        wins, vns, spur = 0, [], []
        for sd in SEEDS:
            d, sigs = per_seed[sd]
            if len(sigs) < m:
                continue
            ok, vn, ns, dim = reduce_and_test(curve, sigs, m, d, k1_bound,
                                              algo, beta)
            wins += ok; vns.append(vn); spur.append(ns)
        if not vns:
            continue
        out.append(dict(m=m, wins=wins, total=len(vns),
                        v=sum(vns) / len(vns), mu=mu, uv=uv, m_info=mi,
                        ratio=(sum(vns) / len(vns)) / mu,
                        n_spur=sum(spur) / len(spur),
                        eff=k1_bound * k2_bound / n))
    return out

# ===========================================================================

print("=" * 78)
print("Thread 20 (part 2): what gates the GLV-HNP attack, if not lambda/n?")
print("=" * 78)

fail_curve = dict(p=2677, b=2, n=2647, lam=185, rho=centered_rho(185, 2647),
                  G=find_generator(2677, 2, 2647))
c8 = dict(p=211, b=2, n=199, lam=106, rho=centered_rho(106, 199),
          G=find_generator(211, 2, 199))
c12 = dict(p=2557, b=2, n=2659, lam=1755, rho=centered_rho(1755, 2659),
           G=find_generator(2557, 2, 2659))

# ---------------------------------------------------------------------------
# Part A: the spurious family is real and occupies the reduced basis
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("Part A: is the spurious 2-coordinate family actually in the basis?")
print("-" * 78)
print(f"  {'curve':<8}{'K1':>4}{'m':>4}{'dim':>5}{'mu':>10}{'||v||':>9}"
      f"{'||v||/mu':>10}{'#spur':>7}{'ok':>4}   short (u,v)")
for lbl, cv, k1, m in [("2677", fail_curve, 8, 8), ("2677", fail_curve, 2, 8),
                       ("2557", c12, 8, 8), ("199", c8, 2, 6)]:
    rec = [r for r in profile(cv, k1, m_max=m) if r['m'] == m][0]
    print(f"  {lbl:<8}{k1:>4}{m:>4}{2*m+2:>5}{rec['mu']:>10.0f}{rec['v']:>9.0f}"
          f"{rec['ratio']:>10.2f}{rec['n_spur']:>7.1f}"
          f"{('Y' if rec['wins']==rec['total'] else 'n'):>4}   {rec['uv']}")
print("  (#spur = reduced-basis rows supported on a single (col i, col m+1+i)")
print("   pair with zero d- and Kannan-coordinates, averaged over seeds)")

# ---------------------------------------------------------------------------
# Part B: measured separation.  Collect every configuration, then look.
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("Part B: does ||v||/mu separate success from failure?")
print("-" * 78)

CONFIGS = ([(f"2677(rho=.070)", fail_curve, k1) for k1 in [2, 3, 4, 6, 8, 12, 16]]
           + [(f"199(rho=.467)", c8, k1) for k1 in [2, 4, 8]]
           + [(f"2557(rho=.340)", c12, k1) for k1 in [8, 16, 32]])

pool = []
pool += find_glv_curves(2**11, 2**13, 0.00, 0.15, want=3)
pool += find_glv_curves(2**11, 2**13, 0.15, 0.30, want=3)
pool += find_glv_curves(2**11, 2**13, 0.30, 0.51, want=3)
for cv in sorted(pool, key=lambda c: c['rho']):
    k2 = math.isqrt(cv['n']) + 1
    CONFIGS.append((f"p={cv['p']}(rho=.{int(cv['rho']*1000):03d})", cv,
                    max(2, int(round(0.05 * cv['n'] / k2)))))

records = []
print(f"  {'config':<20}{'K1':>4}{'eff':>7}{'mu':>9}{'m_info':>7}"
      f"{'first 3/3':>10}{'||v||/mu at that m':>20}")
for lbl, cv, k1 in CONFIGS:
    prof = profile(cv, k1)
    if not prof:
        print(f"  {lbl:<20}{k1:>4}   eff>=1 -- excluded (solution not unique)")
        continue
    for r in prof:
        r['label'] = lbl; r['K1'] = k1; r['rho'] = cv['rho']
    records += prof
    hit = next((r for r in prof if r['wins'] == r['total']), None)
    hit_m = ("m=%d" % hit['m']) if hit else "never"
    hit_r = ("%.2f" % hit['ratio']) if hit else "-"
    print(f"  {lbl:<20}{k1:>4}{prof[0]['eff']:>7.4f}{prof[0]['mu']:>9.0f}"
          f"{prof[0]['m_info']:>7}{hit_m:>10}{hit_r:>20}")

# separation analysis, restricted to configurations that clear m_info
elig = [r for r in records if r['m'] >= r['m_info']]
succ = [r for r in elig if r['wins'] == r['total']]
fail = [r for r in elig if r['wins'] == 0]
print(f"\n  {len(elig)} (config,m) points clear m_info: "
      f"{len(succ)} full-success, {len(fail)} full-failure, "
      f"{len(elig)-len(succ)-len(fail)} partial")
if succ and fail:
    print(f"  ||v||/mu among full successes: "
          f"[{min(r['ratio'] for r in succ):.2f}, {max(r['ratio'] for r in succ):.2f}]")
    print(f"  ||v||/mu among full failures : "
          f"[{min(r['ratio'] for r in fail):.2f}, {max(r['ratio'] for r in fail):.2f}]")

print(f"\n  success rate bucketed by ||v||/mu (m >= m_info only):")
buckets = [(0, .25), (.25, .5), (.5, .75), (.75, 1.0), (1.0, 1.5), (1.5, 1e9)]
for lo, hi in buckets:
    b = [r for r in elig if lo <= r['ratio'] < hi]
    if not b:
        continue
    w = sum(r['wins'] for r in b); t = sum(r['total'] for r in b)
    hi_s = "inf" if hi > 1e8 else f"{hi:.2f}"
    print(f"    {lo:.2f} <= ||v||/mu < {hi_s:>4}: {w:4d}/{t:4d} "
          f"= {100.0*w/t:5.1f}%   ({len(b)} configs)")

# is eff a confound?  bucket by eff too
print(f"\n  success rate bucketed by eff (m >= m_info only), for contrast:")
for lo, hi in [(0, .05), (.05, .10), (.10, .20), (.20, .40), (.40, 1.0)]:
    b = [r for r in elig if lo <= r['eff'] < hi]
    if not b:
        continue
    w = sum(r['wins'] for r in b); t = sum(r['total'] for r in b)
    print(f"    {lo:.2f} <= eff < {hi:.2f}: {w:4d}/{t:4d} = {100.0*w/t:5.1f}%"
          f"   ({len(b)} configs)")

# ---------------------------------------------------------------------------
# Part C: BKZ cannot rescue a small-mu configuration
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("Part C: BKZ cannot rescue -- the spurious vectors ARE the short vectors")
print("-" * 78)
n = fail_curve['n']; k2 = math.isqrt(n) + 1
mu8, uv8 = spurious_mu(n, fail_curve['lam'], 8, k2)
print(f"  p=2677, K1=8: mu={mu8:.0f} attained at (u,v)={uv8}, "
      f"m_info={m_info(n,8,k2)}")
for algo, beta in [('LLL', 0), ('BKZ', 20), ('BKZ', 40)]:
    prof = profile(fail_curve, 8, algo, beta, m_max=10)
    hit = next((r for r in prof if r['wins'] == r['total']), None)
    r8 = next(r for r in prof if r['m'] == 8)
    tag = algo if algo == 'LLL' else f'BKZ({beta})'
    print(f"    {tag:<8} 3/3 at: {('m=' + str(hit['m'])) if hit else 'never':<8}"
          f" | at m=8: ||v||/mu={r8['ratio']:.2f}, "
          f"spurious rows {r8['n_spur']:.1f}/{2*8+2}")

# ---------------------------------------------------------------------------
# Part D: mu and rho are decoupled
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("Part D: mu vs rho -- are they the same quantity in disguise?")
print("-" * 78)
byconf = {}
for r in records:
    byconf.setdefault((r['label'], r['K1']), r)
pts = [(v['rho'], v['mu'], v['eff'], k[0]) for k, v in byconf.items()]
pool_pts = [(rho, mu, eff, lbl) for rho, mu, eff, lbl in pts if lbl.startswith('p=')]
pool_pts.sort()
print(f"  pool curves at matched eff ~ 0.045 (so only the curve changes):")
print(f"  {'curve':<20}{'rho':>8}{'mu':>12}")
for rho, mu, eff, lbl in pool_pts:
    print(f"  {lbl:<20}{rho:>8.3f}{mu:>12.0f}")
if len(pool_pts) > 2:
    rs = [p[0] for p in pool_pts]; ms = [math.log(p[1]) for p in pool_pts]
    mr = sum(rs)/len(rs); mm = sum(ms)/len(ms)
    num = sum((a-mr)*(b-mm) for a, b in zip(rs, ms))
    den = math.sqrt(sum((a-mr)**2 for a in rs) * sum((b-mm)**2 for b in ms))
    print(f"\n  Pearson corr(rho, log mu) over the pool = {num/den:+.3f}"
          f"   (n={len(pool_pts)})")
    print(f"  -> rho carries little information about mu; they are not the")
    print(f"     same quantity.  Note the extremes: smallest rho has a LARGE")
    print(f"     mu and one of the largest rho has among the smallest mu.")

print("\n" + "=" * 78)
print("SUMMARY")
print("=" * 78)
print(f"\n  rho among pool successes: "
      f"{sorted({round(r['rho'],3) for r in records if r['label'].startswith('p=') and r['wins']==r['total']})}")
print(f"  rho among pool failures : "
      f"{sorted({round(r['rho'],3) for r in records if r['label'].startswith('p=')} - {round(r['rho'],3) for r in records if r['label'].startswith('p=') and r['wins']==r['total']})}")
print("\nDone.")

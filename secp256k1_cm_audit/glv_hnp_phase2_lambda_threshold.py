"""
GLV-HNP Phase 2, Thread 20: is the attack-viability threshold really lambda/n?

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-26 run #1):
  The 20-bit scaling run found that a 12-bit curve with lam/n = 0.07 (p=2677,
  n=2647, lam=185) resisted both LLL and BKZ(40) at every m tested, whereas
  curves with lam/n in {0.34, 0.53, 0.66} succeeded.  The run concluded a
  "lambda/n threshold for attack viability" and attributed it to the -lam rows
  failing to separate from the modular rows in Gram-Schmidt.

This script tests that conclusion against a competing explanation.

TWO OBSERVATIONS THE PRIOR RUN DID NOT ACCOUNT FOR
--------------------------------------------------
(O1) lam is only defined mod n.  Lattice row (m+1+i) carries -lam*S_K1 in
     column i, and modular row i carries n*S_K1 in the same column and is zero
     elsewhere.  So replacing lam by lam-n is an integer row operation: it
     generates the SAME lattice.  The attack therefore cannot depend on lam/n,
     only on the CENTERED ratio

         rho = min(lam, n-lam) / n   in (0, 1/2].

     Recomputing the prior table in rho: 0.53 -> 0.47, 0.66 -> 0.34,
     0.34 -> 0.34, 0.07 -> 0.07.  (Part 0 below verifies the lattice identity.)

(O2) rho is not a free parameter.  The two GLV eigenvalues satisfy
     lam_1 + lam_2 = n-1, so centered(lam_1) = centered(lam_2) +/- 1: rho is an
     invariant of n, not a choice.  To vary rho you must change curve, which
     changes n, K2 = isqrt(n)+1, and hence the whole lattice geometry.  In the
     prior run's data rho is confounded with n: the single failure had both the
     smallest rho AND (jointly with 2557) the smallest n.

COMPETING HYPOTHESIS (H_geo)
----------------------------
Success is governed by lattice geometry alone, not by rho.  The planted vector

    v = (k1_0*S_K1, ..., d*S_D, ..., k2_0*S_K2, ..., S_KANNAN)

is recoverable iff it is actually the shortest vector LLL can reach, i.e. iff
||v|| is below what LLL attains on a lattice of that determinant and dimension.
With S_K1 = n//K1, S_K2 = n//K2, S_D = 1, S_KANNAN = n the determinant is

    det = (n*S_K1)^m * S_D * S_K2^m * S_KANNAN  ~  n^(3m+1) / (K1*K2)^m

so the normalised target length is

    ||v|| / det^(1/dim)  ~  sqrt((2m+4)/3) * n^((1-m)/(2m+2)) * (K1*K2)^(m/(2m+2))

which tends to sqrt(eff) = sqrt(K1*K2/n) as m grows.  Crucially this depends on
K1 and m but NOT on lam.

THE DECISIVE FALSIFIER (Part 2)
-------------------------------
Hold the failure curve p=2677 fixed -- rho stays pinned at 0.07 -- and vary only
K1, which moves the geometry but cannot move rho.

    H_lambda (prior run) predicts: fails for every K1, at every m.
    H_geo    predicts: succeeds once K1 is small enough / m large enough.

Part 3 runs the reverse control: degrade a high-rho curve's geometry by raising
K1 and see whether it fails despite rho ~ 0.47.

Part 4 compares rho against the measured geometry as a predictor over a pool of
curves spanning both.

All recoveries are verified WITHOUT a secret oracle: a candidate d is accepted
only if every signature's implied nonce admits a GLV decomposition inside the
declared (K1, K2) box.  (The prior scripts compared d_cand against d_secret,
which cannot distinguish "recovered" from "coincidence"; it also cannot be run
by a real attacker.)

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

from glv_hnp_phase2_common import (
    find_generator, glv_eigenvalues, centered_rho, find_glv_curves,
    gen_signatures, scales, build_glv_lattice, planted_vector, norm,
    log2_det, verify_d,
)

# ---------------------------------------------------------------------------
# One instance: reduce, measure geometry, attempt recovery
# ---------------------------------------------------------------------------

def run_instance(curve, m, d_secret, k1_bound, seed=42, algo='LLL', beta=20):
    n, lam = curve['n'], curve['lam']
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(curve, d_secret, m, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    dim = 2 * m + 2
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    v = planted_vector(sigs, n, d_secret, k1_bound, k2_bound)

    A = IntegerMatrix.from_matrix(M)
    if algo == 'LLL':
        LLL.reduction(A)
    else:
        BKZ.reduction(A, BKZ.Param(beta))
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]

    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    ok = False
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN:
            continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if verify_d(d_cand, sigs, n, lam, k1_bound, k2_bound):
            ok = True
            break

    v_norm = norm(v)
    b1_norm = min(norm(r) for r in rows)
    ld = log2_det(m, n, k1_bound, k2_bound)
    gh_log2 = 0.5 * math.log2(dim / (2 * math.pi * math.e)) + ld / dim
    return dict(success=ok, v_norm=v_norm, b1_norm=b1_norm,
                gap=b1_norm / v_norm,
                v_over_gh=v_norm / (2 ** gh_log2),
                delta_req=2 ** ((math.log2(v_norm) - ld / dim) / dim),
                delta_lll=2 ** ((math.log2(b1_norm) - ld / dim) / dim),
                dim=dim)

def sweep(curve, k1_bound, m_range, seeds, algo='LLL', beta=20, verbose=True,
          label=''):
    """Returns (first_m_all_seeds, per_m_results)."""
    n = curve['n']
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    m_thresh = (math.ceil(math.log(n) / math.log(1.0 / eff))
                if eff < 1.0 else float('inf'))
    if verbose:
        tag = algo if algo == 'LLL' else f'{algo}({beta})'
        print(f"  {label}  p={n and curve['p']} n={n}({n.bit_length()}b) "
              f"lam={curve['lam']} rho={curve['rho']:.3f} "
              f"K1={k1_bound} K2={k2_bound} eff={eff:.4f} "
              f"m_thresh~{m_thresh:.1f} [{tag}]")
    per_m = {}
    first = None
    for m in m_range:
        wins, stats = 0, []
        for sd in seeds:
            d_trial = random.Random(sd + 7777).randrange(1, n)
            res = run_instance(curve, m, d_trial, k1_bound, sd, algo, beta)
            if res is None:
                continue
            wins += res['success']
            stats.append(res)
        if not stats:
            continue
        avg = lambda k: sum(s[k] for s in stats) / len(stats)
        per_m[m] = dict(wins=wins, total=len(stats), gap=avg('gap'),
                        v_over_gh=avg('v_over_gh'),
                        delta_req=avg('delta_req'), delta_lll=avg('delta_lll'),
                        dim=stats[0]['dim'])
        if verbose:
            print(f"    m={m:<3d} {wins}/{len(stats)}  "
                  f"||v||/GH={avg('v_over_gh'):6.3f}  "
                  f"||b1||/||v||={avg('gap'):6.3f}  "
                  f"d_req={avg('delta_req'):.4f} d_LLL={avg('delta_lll'):.4f}")
        if first is None and wins == len(stats):
            first = m
    return first, per_m


SEEDS = [42, 1234, 9999]

# ===========================================================================
# Part 0 -- lam and lam-n generate the same lattice (observation O1)
# ===========================================================================

print("=" * 78)
print("Thread 20: is the GLV-HNP viability threshold really lambda/n?")
print("=" * 78)
print("\n" + "-" * 78)
print("Part 0: lam vs lam-n  --  do they give the same attack?")
print("-" * 78)

fail_curve = dict(p=2677, b=2, n=2647, lam=185, lam2=2647 - 1 - 185,
                  rho=centered_rho(185, 2647), G=None)
fail_curve['G'] = find_generator(2677, 2, 2647)
print(f"  failure curve from 2026-07-26: p=2677 n=2647 lam=185")
print(f"    eigenvalues of x^2+x+1 mod n: {glv_eigenvalues(2647)}")
print(f"    lam/n     = {185/2647:.4f}   (quantity used by the prior run)")
print(f"    rho       = {centered_rho(185,2647):.4f}   (centered, basis-invariant)")
l1, l2 = glv_eigenvalues(2647)
print(f"    rho(lam_1)= {centered_rho(l1,2647):.4f}, "
      f"rho(lam_2)= {centered_rho(l2,2647):.4f}  -> equal, as lam_1+lam_2=n-1")

_sigs = gen_signatures(fail_curve, 12345, 4, 8, math.isqrt(2647) + 1, seed=42)
_M_a = build_glv_lattice(_sigs, 2647, 185, 8, math.isqrt(2647) + 1)
_M_b = build_glv_lattice(_sigs, 2647, 185 - 2647, 8, math.isqrt(2647) + 1)
_A = IntegerMatrix.from_matrix(_M_a); LLL.reduction(_A)
_B = IntegerMatrix.from_matrix(_M_b); LLL.reduction(_B)
_ra = sorted(tuple(_A[i][j] for j in range(_A.ncols)) for i in range(_A.nrows))
_rb = sorted(tuple(_B[i][j] for j in range(_B.ncols)) for i in range(_B.nrows))
print(f"    LLL(basis with lam) == LLL(basis with lam-n)? {_ra == _rb}")
print("    => the attack cannot depend on lam/n, only on rho. O1 confirmed."
      if _ra == _rb else "    => UNEXPECTED: bases differ.")

# ===========================================================================
# Part 1 -- reproduce the prior run's four data points, in rho
# ===========================================================================

print("\n" + "-" * 78)
print("Part 1: reproduce 2026-07-26 baseline, recomputed in rho")
print("-" * 78)

c8 = dict(p=211, b=2, n=199, lam=106, lam2=92, rho=centered_rho(106, 199),
          G=find_generator(211, 2, 199))
c12 = dict(p=2557, b=2, n=2659, lam=1755, lam2=903, rho=centered_rho(1755, 2659),
           G=find_generator(2557, 2, 2659))

base = []
for lbl, cv, k1 in [("8-bit/199", c8, 2), ("12-bit/2557", c12, 8),
                    ("12-bit/2677", fail_curve, 8)]:
    f, pm = sweep(cv, k1, range(3, 13), SEEDS, label=lbl)
    base.append((lbl, cv, k1, f, pm))
    print(f"    -> first 3/3 at m={f}" if f else "    -> never 3/3")

print("\n  prior run's table, corrected:")
print(f"  {'curve':<14}{'lam/n':>8}{'rho':>8}{'K1':>5}{'first 3/3':>11}")
for lbl, cv, k1, f, pm in base:
    print(f"  {lbl:<14}{cv['lam']/cv['n']:>8.3f}{cv['rho']:>8.3f}{k1:>5}"
          f"{(str(f) if f else 'never'):>11}")

# ===========================================================================
# Part 2 -- THE FALSIFIER: hold rho=0.07 fixed, vary only K1
# ===========================================================================

print("\n" + "-" * 78)
print("Part 2: FALSIFIER -- failure curve p=2677 (rho=0.070 pinned), sweep K1")
print("-" * 78)
print("  H_lambda predicts: no K1 works.   H_geo predicts: small K1 works.\n")

part2 = []
for k1 in [2, 3, 4, 6, 8, 12, 16]:
    f, pm = sweep(fail_curve, k1, range(3, 15), SEEDS,
                  label=f"K1={k1}", verbose=False)
    n = fail_curve['n']; k2 = math.isqrt(n) + 1
    eff = k1 * k2 / n
    best = pm.get(f) if f else None
    print(f"  K1={k1:<3d} eff={eff:.4f}  first 3/3: "
          f"{(('m=' + str(f)) if f else 'never'):<8}"
          + (f"  at m={f}: ||v||/GH={best['v_over_gh']:.3f} "
             f"||b1||/||v||={best['gap']:.3f}" if best else ""))
    part2.append((k1, eff, f, pm))

rescued = [k1 for k1, eff, f, pm in part2 if f is not None]
print(f"\n  K1 values that recover d at rho=0.07: {rescued if rescued else 'NONE'}")

# ===========================================================================
# Part 3 -- reverse control: degrade a high-rho curve's geometry
# ===========================================================================

print("\n" + "-" * 78)
print("Part 3: REVERSE CONTROL -- high-rho curve, degrade geometry via large K1")
print("-" * 78)
print("  If rho drove viability, 12-bit/2557 (rho=0.34) should keep working.\n")

part3 = []
for k1 in [8, 16, 32, 64, 128]:
    n = c12['n']; k2 = math.isqrt(n) + 1
    eff = k1 * k2 / n
    if eff >= 1.0:
        # The (k1,k2) box then covers more than n residues, so a WRONG d still
        # admits a decomposition for every signature and verify_d becomes
        # vacuous: the HNP instance has no unique solution.  Such rows carry no
        # information and are excluded rather than reported as "successes".
        print(f"  K1={k1:<4d} eff={eff:.4f}  EXCLUDED (eff>=1: solution not "
              f"unique, verification vacuous)")
        continue
    f, pm = sweep(c12, k1, range(3, 15), SEEDS, label=f"K1={k1}", verbose=False)
    print(f"  K1={k1:<4d} eff={eff:.4f}  first 3/3: "
          f"{(('m=' + str(f)) if f else 'never')}")
    part3.append((k1, eff, f, pm))

# ===========================================================================
# Part 4 -- rho vs geometry as predictors over a curve pool
# ===========================================================================

print("\n" + "-" * 78)
print("Part 4: curve pool -- does rho predict anything once eff is fixed?")
print("-" * 78)

pool = []
pool += find_glv_curves(2**11, 2**13, 0.00, 0.15, want=3)
pool += find_glv_curves(2**11, 2**13, 0.15, 0.30, want=3)
pool += find_glv_curves(2**11, 2**13, 0.30, 0.51, want=3)
print(f"  pool: {len(pool)} curves, 11-13 bit, rho spanning "
      f"[{min(c['rho'] for c in pool):.3f}, {max(c['rho'] for c in pool):.3f}]\n")

print(f"  {'p':>7}{'n':>8}{'rho':>8}{'K1':>5}{'eff':>8}{'first3/3':>10}"
      f"{'||v||/GH':>10}")
rows4 = []
for cv in sorted(pool, key=lambda c: c['rho']):
    n = cv['n']; k2 = math.isqrt(n) + 1
    # hold eff ~ 0.05 across the pool so geometry is matched and only rho varies
    k1 = max(2, int(round(0.05 * n / k2)))
    f, pm = sweep(cv, k1, range(3, 15), SEEDS, verbose=False)
    eff = k1 * k2 / n
    ref = pm.get(f) if f else (pm.get(max(pm)) if pm else None)
    vgh = "%.3f" % ref['v_over_gh'] if ref else "-"
    print(f"  {cv['p']:>7}{n:>8}{cv['rho']:>8.3f}{k1:>5}{eff:>8.4f}"
          f"{(str(f) if f else 'never'):>10}{vgh:>10}")
    rows4.append((cv, k1, eff, f, pm))

# ===========================================================================
# Summary
# ===========================================================================

print("\n" + "=" * 78)
print("SUMMARY")
print("=" * 78)

print(f"\n  Part 0: LLL(lam) == LLL(lam-n): {_ra == _rb}"
      f"  -> lam/n is not a basis-invariant quantity; rho is.")

print(f"\n  Part 2 (rho pinned at 0.070, only K1 varies):")
for k1, eff, f, pm in part2:
    print(f"    K1={k1:<3d} eff={eff:.4f} -> {('m=' + str(f)) if f else 'never'}")
if rescued:
    print(f"    VERDICT: H_lambda FALSIFIED. The rho=0.07 curve is attackable")
    print(f"             at K1 in {rescued}; rho never changed.")
else:
    print(f"    VERDICT: H_lambda survives this test -- no K1 rescued the curve.")

print(f"\n  Part 3 (rho pinned at {c12['rho']:.3f}, only K1 varies):")
for k1, eff, f, pm in part3:
    print(f"    K1={k1:<4d} eff={eff:.4f} -> {('m=' + str(f)) if f else 'never'}")
broke = [k1 for k1, eff, f, pm in part3 if f is None]
if broke:
    print(f"    VERDICT: high-rho curve FAILS at K1 in {broke} -- large rho is")
    print(f"             not sufficient for viability.")
else:
    print(f"    VERDICT: high-rho curve never broke in the tested K1 range.")

print(f"\n  Part 4 (eff held ~0.05, rho varies over the pool):")
succ = [(cv['rho'], f) for cv, k1, eff, f, pm in rows4 if f is not None]
fail = [cv['rho'] for cv, k1, eff, f, pm in rows4 if f is None]
print(f"    recovered: {len(succ)}/{len(rows4)} curves")
if succ:
    print(f"      rho range among successes: "
          f"[{min(r for r, _ in succ):.3f}, {max(r for r, _ in succ):.3f}]")
if fail:
    print(f"      rho among failures: {[f'{r:.3f}' for r in fail]}")
else:
    print(f"      no failures -> rho does not gate viability at fixed eff.")

print("\nDone.")

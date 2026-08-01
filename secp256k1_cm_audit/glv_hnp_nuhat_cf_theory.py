"""
GLV-HNP Thread 23: closed form for lambda_1(L2), and why six June invariants failed.

Context
-------
Thread 20c/20d (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run #2) found that

    nu_hat = lambda_1(L2) / sqrt(det L2),
    L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >          (`rival_sublattice_nu`)

separates the June C1/C2 classes with AUC 0.935 where six scale-free curve
invariants (delta/n, kappa(M), q_cf, max_q_cf, max_a, a_corn/n) and two
lambda-ratio invariants (lam/n, mu) all sat at the trivial baseline.  But
nu_hat was *computed*, not *understood*: a Lagrange-Gauss reduction is a black
box.  The 2026-07-29 next-step proposal was to derive it.

This script does that.  Every vector of L2 is

    a*(n*S_K1, 0) + b*(-lam*S_K1, S_K2) = ( (a*n - b*lam)*S_K1 , b*S_K2 ),

so for fixed b the best a makes the first coordinate S_K1*||b*lam||_n, where
||x||_n = min_{c in Z} |x - c*n|.  Hence

    lambda_1(L2)^2 = min_b [ (S_K1*||b*lam||_n)^2 + (b*S_K2)^2 ].         (*)

Only b with b*S_K2 < n*S_K1 can win, and among those the b that lower
||b*lam||_n below every smaller b -- the "records" -- are exactly the
best approximations of the second kind to lam/n, i.e. the continued-fraction
CONVERGENT DENOMINATORS q_j (semiconvergents are best approximations of the
first kind and are NOT records for ||b*lam||_n).  So (*) reduces to a search
over the O(log n) convergents, and the minimiser is the convergent nearest the
balance point

    b* = sqrt(n * S_K1 / S_K2) = sqrt(n * K2 / K1).                       (**)

Experiments
-----------
  E1  Structure theorem.  Exact Lagrange-Gauss lambda_1 vs the convergent
      minimum of (*).  Predicted: identical on every instance.
  E2  Single-convergent predictor.  Use ONLY the convergent nearest b* in
      (**).  Falsifier (2026-07-29): if corr(nu_hat_cf1, nu_hat) < 0.9 the
      balance-point story is wrong.
  E3  Why the June invariants failed.  The argmin index j in (*) is a function
      of the SCALE K2/K1, not of lam alone: sweeping K1 over a fixed curve
      walks the argmin down the convergent sequence and moves nu_hat over most
      of its range.  Any scale-free function of the CF of lam/n (max_a, q_cf,
      ...) is therefore constant along a direction in which nu_hat varies, so
      it cannot predict recovery.  Quantified as a rank correlation that is
      ~0 by construction.

Run: python3 glv_hnp_nuhat_cf_theory.py
"""

import math
import random
import sys
import time

import sympy

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_t20a", __file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_nuhat_base.py")
_t20a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20a)

glv_eigenvalues = _t20a.glv_eigenvalues
mu_of = _t20a.mu_of
scales = _t20a.scales
gauss_reduce_2d = _t20a.gauss_reduce_2d
rival_sublattice_nu = _t20a.rival_sublattice_nu

BITS_LO, BITS_HI = 2 ** 19, 2 ** 20
N_CURVES = 200
K1_FIXED = 72                      # Exp S protocol (2026-06-29)

# ---------------------------------------------------------------------------
# Continued fractions
# ---------------------------------------------------------------------------

def cf_convergents(a, b):
    """
    Convergent denominators q_j of a/b, together with the residual
    r_j = |q_j*a - p_j*b| = ||q_j*a||_b.  Returns [(q_j, r_j, a_j)] where a_j
    is the j-th partial quotient.  q_0 = 1 always (the j=0 convergent).
    """
    out = []
    x, y = a, b
    p_prev, p_cur = 1, 0        # p_{-1}, p_{-2} bookkeeping (numerators)
    q_prev, q_cur = 0, 1
    while y:
        quo, rem = divmod(x, y)
        p_prev, p_cur = quo * p_prev + p_cur, p_prev
        q_prev, q_cur = quo * q_prev + q_cur, q_prev
        out.append((q_prev, abs(q_prev * a - p_prev * b), quo))
        x, y = y, rem
    return out

def max_partial_quotient(a, b):
    """max_a, falsified by Exp S 2026-06-29.  Scale-free by construction."""
    return max(q for _, _, q in cf_convergents(a, b)) if b else 0

def lambda1_from_cf(n, lam, k1_bound, k2_bound):
    """
    lambda_1(L2) via the convergent search (*), plus diagnostics:
      j_star   index of the winning convergent (-1 = the b=0 vector)
      q_star   its denominator
      b_star   the balance point sqrt(n*S_K1/S_K2) of (**)
      n_cands  number of convergents that could win (q_j*S_K2 < n*S_K1)
    """
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    best = float(n * s_k1)                      # b = 0
    j_star, q_star = -1, 0
    n_cands = 0
    for j, (q, r, _) in enumerate(cf_convergents(lam, n)):
        if q * s_k2 >= n * s_k1:
            break
        n_cands += 1
        cand = math.hypot(r * s_k1, q * s_k2)
        if cand < best:
            best, j_star, q_star = cand, j, q
    return {
        'lam1': best, 'j_star': j_star, 'q_star': q_star,
        'b_star': math.sqrt(n * s_k1 / s_k2), 'n_cands': n_cands,
    }

def _anchor(cands, b_star, corrected):
    """
    Index of the convergent to anchor the window on.

    corrected=False -- the 2026-07-29 rule: q_j closest to b* in log scale.
    corrected=True  -- the sharpened rule.  The naive rule implicitly assumes
        r_j = ||q_j*lam||_n ~ n/q_j, but the true CF bound is
            1/(q_{j+1} + q_j) < r_j/n < 1/q_{j+1},
        i.e. r_j ~ n/q_{j+1}.  Balancing r_j*S_K1 against q_j*S_K2 therefore
        gives q_j*q_{j+1} ~ n*S_K1/S_K2 = b*^2, so the anchor is the j whose
        GEOMETRIC MEAN sqrt(q_j*q_{j+1}) -- not q_j itself -- sits at b*.
        Because q_{j+1} > q_j this always anchors at or below the naive index,
        which is exactly the one-sided offset measured in E2b.
    """
    if not corrected:
        return min(range(len(cands)),
                   key=lambda i: abs(math.log(cands[i][0]) - math.log(b_star)))
    def key(i):
        q_next = cands[i + 1][0] if i + 1 < len(cands) else cands[i][0] ** 2
        return abs(0.5 * (math.log(cands[i][0]) + math.log(q_next))
                   - math.log(b_star))
    return min(range(len(cands)), key=key)

def nu_hat_window(n, lam, k1_bound, k2_bound, width=0, corrected=False):
    """
    Window predictor: locate an anchor convergent near the balance point b* and
    minimise (*) over the +-`width` convergents around it.  width=0 with
    corrected=False is the single-convergent predictor proposed on 2026-07-29;
    width -> inf recovers the exact lambda_1(L2) of E1.
    """
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    b_star = math.sqrt(n * s_k1 / s_k2)
    cands = [(q, r) for q, r, _ in cf_convergents(lam, n) if q * s_k2 < n * s_k1]
    if not cands:
        return (n * s_k1) / math.sqrt(n * s_k1 * s_k2)
    j0 = _anchor(cands, b_star, corrected)
    lo, hi = max(0, j0 - width), min(len(cands), j0 + width + 1)
    best = float(n * s_k1)
    for q, r in cands[lo:hi]:
        best = min(best, math.hypot(r * s_k1, q * s_k2))
    return best / math.sqrt(n * s_k1 * s_k2)

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def pearson(xs, ys):
    k = len(xs)
    mx, my = sum(xs) / k, sum(ys) / k
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else float('nan')

def rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r

def spearman(xs, ys):
    return pearson(rank(xs), rank(ys))

# ---------------------------------------------------------------------------
# Sample
# ---------------------------------------------------------------------------

def harvest_lambdas(n_curves):
    """
    (n, lam) pairs with n prime, n = 1 mod 3, lam^2+lam+1 = 0 mod n.  The
    C1/C2 experiments need the full curve (p, b, G); Thread 23 is pure number
    theory on (n, lam, K1, K2), so only n and lam are needed here.
    """
    out, seen = [], set()
    n = int(sympy.nextprime(BITS_LO))
    while n < BITS_HI and len(out) < n_curves:
        if n % 3 == 1 and n not in seen:
            roots = glv_eigenvalues(n)
            if roots is not None:
                seen.add(n)
                out.append({'n': n, 'lam': roots[0]})
        n = int(sympy.nextprime(n))
    return out

# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("GLV-HNP Thread 23: closed form for lambda_1(L2)")
    print("=" * 78)

    t0 = time.time()
    sample = harvest_lambdas(N_CURVES)
    for c in sample:
        c['k2'] = math.isqrt(c['n']) + 1
    print(f"sample: {len(sample)} 20-bit (n, lam) pairs in {time.time()-t0:.1f}s")

    # ---------------- E1: structure theorem ----------------
    print("\n" + "-" * 78)
    print("E1  exact Lagrange-Gauss lambda_1  ==  convergent minimum (*)")
    print("-" * 78)
    K1_GRID = [8, 24, 72, 216, 648]
    print(f"{'K1':>6} {'instances':>10} {'exact match':>12} {'max rel err':>12}"
          f" {'mean n_cands':>13}")
    e1_total = e1_match = 0
    for k1 in K1_GRID:
        worst, match, cands = 0.0, 0, 0
        for c in sample:
            exact, _, _ = rival_sublattice_nu(c['n'], c['lam'], k1, c['k2'])
            cf = lambda1_from_cf(c['n'], c['lam'], k1, c['k2'])
            cands += cf['n_cands']
            rel = abs(exact - cf['lam1']) / exact
            worst = max(worst, rel)
            if rel < 1e-9:
                match += 1
        e1_total += len(sample)
        e1_match += match
        print(f"{k1:>6} {len(sample):>10} {match:>7}/{len(sample):<4}"
              f" {worst:>12.2e} {cands/len(sample):>13.1f}")
    print(f"\n  overall {e1_match}/{e1_total} exact "
          f"({100.0*e1_match/e1_total:.2f}%)")
    print("  -> lambda_1(L2) is the min over CF convergent denominators of")
    print("     sqrt( (S_K1*||q_j*lam||_n)^2 + (q_j*S_K2)^2 ),  O(log n) work.")

    # random (non-CM) lambda, to show the theorem is about the lattice not GLV
    print("\n  control: 500 uniformly random lam mod n (no CM structure)")
    rng = random.Random(31337)
    ctrl_match = 0
    for _ in range(500):
        c = rng.choice(sample)
        lam = rng.randrange(1, c['n'])
        exact, _, _ = rival_sublattice_nu(c['n'], lam, K1_FIXED, c['k2'])
        cf = lambda1_from_cf(c['n'], lam, K1_FIXED, c['k2'])
        if abs(exact - cf['lam1']) / exact < 1e-9:
            ctrl_match += 1
    print(f"    {ctrl_match}/500 exact")

    # ---------------- E2: single-convergent predictor ----------------
    print("\n" + "-" * 78)
    print("E2  single convergent at the balance point b* = sqrt(n*K2/K1)")
    print("-" * 78)
    print("falsifier (2026-07-29): corr(nu_hat_cf1, nu_hat) < 0.9 kills the story")
    print(f"{'K1':>6} {'pearson':>9} {'spearman':>9} {'argmin==nearest':>16}"
          f" {'|dj|<=1':>9} {'mean q*/b*':>11}")
    for k1 in K1_GRID:
        exact_nu, cf1_nu = [], []
        hit, near1, ratios = 0, 0, []
        for c in sample:
            _, _, nh = rival_sublattice_nu(c['n'], c['lam'], k1, c['k2'])
            exact_nu.append(nh)
            cf1_nu.append(nu_hat_window(c['n'], c['lam'], k1, c['k2'], 0))
            info = lambda1_from_cf(c['n'], c['lam'], k1, c['k2'])
            s_k1, _, s_k2, _ = scales(c['n'], k1, c['k2'])
            cands = [(q, r) for q, r, _ in cf_convergents(c['lam'], c['n'])
                     if q * s_k2 < c['n'] * s_k1]
            if cands and info['j_star'] >= 0:
                jn = min(range(len(cands)),
                         key=lambda i: abs(math.log(cands[i][0])
                                           - math.log(info['b_star'])))
                hit += (jn == info['j_star'])
                near1 += (abs(jn - info['j_star']) <= 1)
                ratios.append(info['q_star'] / info['b_star'])
        pr = pearson(cf1_nu, exact_nu)
        sp = spearman(cf1_nu, exact_nu)
        mr = sum(ratios) / len(ratios) if ratios else float('nan')
        print(f"{k1:>6} {pr:>9.4f} {sp:>9.4f} {hit:>10}/{len(sample):<5}"
              f" {near1:>5}/{len(sample):<3} {mr:>11.3f}")

    # ---------------- E2b: how wide a window is needed ----------------
    print("\n" + "-" * 78)
    print("E2b  minimise (*) over a +-w window of convergents around b*")
    print("-" * 78)
    print("w=0 is E2; w=inf is the exact lambda_1 of E1.  How big must w be?")
    print(f"{'K1':>6} " + " ".join(f"{'w='+str(w):>9}" for w in (0, 1, 2, 3))
          + f" {'exact @w=1':>11} {'exact @w=2':>11}")
    for k1 in K1_GRID:
        exact_nu = [rival_sublattice_nu(c['n'], c['lam'], k1, c['k2'])[2]
                    for c in sample]
        row, ex1, ex2 = [], 0, 0
        for w in (0, 1, 2, 3):
            pred = [nu_hat_window(c['n'], c['lam'], k1, c['k2'], w)
                    for c in sample]
            row.append(pearson(pred, exact_nu))
            if w == 1:
                ex1 = sum(1 for a, b in zip(pred, exact_nu) if abs(a - b) < 1e-12)
            if w == 2:
                ex2 = sum(1 for a, b in zip(pred, exact_nu) if abs(a - b) < 1e-12)
        print(f"{k1:>6} " + " ".join(f"{v:>9.4f}" for v in row)
              + f" {ex1:>7}/{len(sample):<3} {ex2:>7}/{len(sample):<3}")

    # ---------------- E2c: corrected balance point ----------------
    print("\n" + "-" * 78)
    print("E2c  corrected anchor: sqrt(q_j*q_{j+1}) at b*, not q_j at b*")
    print("-" * 78)
    print("r_j ~ n/q_{j+1} (not n/q_j), so the balance is q_j*q_{j+1} ~ b*^2.")
    print(f"{'K1':>6} {'naive w=0':>10} {'corr. w=0':>10} {'corr. w=1':>10}"
          f" {'exact @w=0':>11} {'exact @w=1':>11}")
    for k1 in K1_GRID:
        exact_nu = [rival_sublattice_nu(c['n'], c['lam'], k1, c['k2'])[2]
                    for c in sample]
        naive = [nu_hat_window(c['n'], c['lam'], k1, c['k2'], 0, False)
                 for c in sample]
        cor0 = [nu_hat_window(c['n'], c['lam'], k1, c['k2'], 0, True)
                for c in sample]
        cor1 = [nu_hat_window(c['n'], c['lam'], k1, c['k2'], 1, True)
                for c in sample]
        ex0 = sum(1 for a, b in zip(cor0, exact_nu) if abs(a - b) < 1e-12)
        ex1 = sum(1 for a, b in zip(cor1, exact_nu) if abs(a - b) < 1e-12)
        print(f"{k1:>6} {pearson(naive, exact_nu):>10.4f}"
              f" {pearson(cor0, exact_nu):>10.4f}"
              f" {pearson(cor1, exact_nu):>10.4f}"
              f" {ex0:>7}/{len(sample):<3} {ex1:>7}/{len(sample):<3}")

    print("\n  distribution of the offset j* - j_nearest (K1=72):")
    offs = {}
    for c in sample:
        s_k1, _, s_k2, _ = scales(c['n'], K1_FIXED, c['k2'])
        cands = [(q, r) for q, r, _ in cf_convergents(c['lam'], c['n'])
                 if q * s_k2 < c['n'] * s_k1]
        info = lambda1_from_cf(c['n'], c['lam'], K1_FIXED, c['k2'])
        if not cands or info['j_star'] < 0:
            continue
        b_star = info['b_star']
        jn = min(range(len(cands)),
                 key=lambda i: abs(math.log(cands[i][0]) - math.log(b_star)))
        offs[info['j_star'] - jn] = offs.get(info['j_star'] - jn, 0) + 1
    for k in sorted(offs):
        print(f"    {k:>+3} : {offs[k]:>4}")

    # ---------------- E3: why the June invariants failed ----------------
    print("\n" + "-" * 78)
    print("E3  scale-dependence: the argmin index walks with K1")
    print("-" * 78)
    print("For a FIXED curve, lam and its whole continued fraction are fixed,")
    print("so every scale-free CF invariant (max_a, q_cf, ...) is constant.")
    print("nu_hat is not:")
    demo = sample[0]
    print(f"\n  n = {demo['n']}  lam = {demo['lam']}  "
          f"max_a = {max_partial_quotient(demo['lam'], demo['n'])} (constant)")
    print(f"{'K1':>7} {'b*':>10} {'j*':>4} {'q*':>10} {'nu_hat':>8}")
    for k1 in (2, 8, 24, 72, 216, 648, 2048):
        info = lambda1_from_cf(demo['n'], demo['lam'], k1, demo['k2'])
        _, _, nh = rival_sublattice_nu(demo['n'], demo['lam'], k1, demo['k2'])
        print(f"{k1:>7} {info['b_star']:>10.1f} {info['j_star']:>4}"
              f" {info['q_star']:>10} {nh:>8.4f}")

    print("\n  across the sample at K1=72: correlation of nu_hat with the")
    print("  scale-free CF invariants (these are the June falsifications):")
    nus = [rival_sublattice_nu(c['n'], c['lam'], K1_FIXED, c['k2'])[2]
           for c in sample]
    maxa = [max_partial_quotient(c['lam'], c['n']) for c in sample]
    # cf_len stands in for the June "q_cf" family; q_0 is identically 1 so it
    # has zero variance and no correlation is even defined for it.
    cflen = [len(cf_convergents(c['lam'], c['n'])) for c in sample]
    mus = [mu_of(c['lam'], c['n']) for c in sample]
    for name, v in (('max_a', maxa), ('cf_len', cflen), ('mu', mus)):
        print(f"    spearman(nu_hat, {name:>6}) = {spearman(nus, v):+.4f}")

    print("\n  and the scale-DEPENDENT quantity that does carry the signal:")
    print("    q*/b*  (how well a convergent lands on the balance point)")
    r_over = [lambda1_from_cf(c['n'], c['lam'], K1_FIXED, c['k2'])['q_star']
              / lambda1_from_cf(c['n'], c['lam'], K1_FIXED, c['k2'])['b_star']
              for c in sample]
    lr = [abs(math.log(x)) for x in r_over]
    print(f"    spearman(nu_hat, |log(q*/b*)|) = {spearman(nus, lr):+.4f}")

    # ---------------- E4: does the closed form survive to 256 bits? ----------
    print("\n" + "-" * 78)
    print("E4  closed form at 256 bits (secp256k1 and random 256-bit lam)")
    print("-" * 78)
    n_s = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
    lam_s = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
    assert (lam_s * lam_s + lam_s + 1) % n_s == 0
    k2_s = math.isqrt(n_s) + 1
    print(f"{'eff':>7} {'K1':>26} {'nu_hat':>9} {'cf w=1':>9} {'exact':>7}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1_s = max(2, int(eff * math.isqrt(n_s)))
        _, _, nh = rival_sublattice_nu(n_s, lam_s, k1_s, k2_s)
        cf = nu_hat_window(n_s, lam_s, k1_s, k2_s, 1, True)
        print(f"{eff:>7.2f} {k1_s:>26} {nh:>9.4f} {cf:>9.4f}"
              f" {'yes' if abs(nh-cf) < 1e-12 else 'NO':>7}")

    rng2 = random.Random(20260801)
    ok = 0
    for _ in range(300):
        lam = rng2.randrange(1, n_s)
        k1 = max(2, int(rng2.choice([0.02, 0.05, 0.1, 0.3, 0.7])
                        * math.isqrt(n_s)))
        _, _, nh = rival_sublattice_nu(n_s, lam, k1, k2_s)
        if abs(nh - nu_hat_window(n_s, lam, k1, k2_s, 1, True)) < 1e-9:
            ok += 1
    print(f"\n  300 random 256-bit lam, eff in [0.02, 0.7]: {ok}/300 exact")

    print("\nDone.")

if __name__ == '__main__':
    sys.exit(main())

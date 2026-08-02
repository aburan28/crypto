"""
GLV-HNP — Thread 23: closed form for lambda_1(L2), i.e. for the nu_hat separator.

Context
-------
2026-07-29 (`e845207`) established that

    nu_hat = lambda_1(L2) / sqrt(det L2),
    L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >,   det L2 = n*S_K1*S_K2

separates the June C1/C2 classes with AUC 0.935, after eight scale-free curve
invariants had failed (delta/n, kappa(M), q_cf, max_q_cf, max_a, a_corn/n,
lam/n, mu).  Its next-step proposal:

    "nu_hat is currently computed, not understood. [...] minimising
     |a*lam mod n|*S_K1 against a*S_K2 is best rational approximation to lam/n
     *at the scale* S_K2/S_K1 = K1/K2. That is why the raw CF invariants all
     failed in June -- they are scale-free."

    Falsifier: "if predicted-from-CF nu_hat correlates <0.9 with computed
    nu_hat, the convergent-scale story is wrong."

This script answers it.  A general vector of L2 is

    a*(n*S_K1, 0) + b*(-lam*S_K1, S_K2) = ( (a*n - b*lam)*S_K1 , b*S_K2 )

so for fixed b >= 0 the optimal a minimises |b*lam - a*n| =: theta(b), and

    ||v||^2 = theta(b)^2 * S_K1^2 + b^2 * S_K2^2.                        (*)

Both terms are monotone in their own argument, so if b is NOT a best rational
approximation denominator of the second kind for lam/n, there is a smaller b'
with theta(b') <= theta(b) and the corresponding vector is shorter in BOTH
coordinates.  Best approximations of the second kind are exactly the continued
fraction convergent denominators of lam/n.  Hence

    lambda_1(L2)^2 = min over CF convergent denominators q_j of lam/n
                     (including q = 0, theta = n) of
                     theta_j^2 * S_K1^2 + q_j^2 * S_K2^2                 (T)

which is a closed form, not a heuristic.  E1 verifies (T) exactly against
Lagrange-Gauss reduction; E2-E4 turn it into the scale-dependent predictor;
E5 places secp256k1.

Normalisation used throughout.  With  A = sqrt(n*S_K2/S_K1),
B = sqrt(n*S_K1/S_K2) = sqrt(n*K2/K1)  (so A*B = n),

    nu_hat^2 = min_j [ (theta_j/A)^2 + (q_j/B)^2 ]

and since theta_j ~ n/q_{j+1} = A*B/q_{j+1}, writing x_j = q_j/B,

    nu_hat^2 ~ min_j [ 1/x_{j+1}^2 + x_j^2 ].                            (**)

So nu_hat depends only on where the convergent ladder of lam/n sits relative to
the single scale B -- which is exactly what a scale-free CF invariant cannot see.

Run: python3 glv_hnp_nuhat_cf_closed_form.py
"""

import math
import random
import sys
import time

import sympy

import importlib.util

_spec = importlib.util.spec_from_file_location(
    "_nb", __file__.rsplit("/", 1)[0] + "/glv_hnp_nuhat_base.py")
_nb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_nb)

glv_eigenvalues = _nb.glv_eigenvalues
mu_of = _nb.mu_of
scales = _nb.scales
rival_sublattice_nu = _nb.rival_sublattice_nu          # float Gauss reduction
gauss_reduce_2d_float = _nb.gauss_reduce_2d

SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

# ---------------------------------------------------------------------------
# Exact integer primitives
# ---------------------------------------------------------------------------

def n2(w):
    return w[0] * w[0] + w[1] * w[1]


def gauss_reduce_2d_exact(u, v):
    """Lagrange-Gauss reduction in exact integer arithmetic.

    Returns the shortest nonzero vector (as a tuple).  Unlike the float version
    in glv_hnp_nuhat_base.py this is provably correct at any bit size.
    """
    if n2(u) < n2(v):
        u, v = v, u
    while True:
        d = u[0] * v[0] + u[1] * v[1]
        nv = n2(v)
        q = (2 * d + nv) // (2 * nv)          # nearest integer to d/nv
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if n2(u) <= n2(v):
            return u


def cf_convergents(num, den):
    """Continued fraction convergents h_j/k_j of num/den.

    Returns (partial_quotients, [(h_j, k_j), ...]) with the j = -1 seed omitted.
    """
    a_list = []
    conv = []
    hm1, hm2 = 1, 0
    km1, km2 = 0, 1
    x, y = num, den
    while y:
        a = x // y
        a_list.append(a)
        h = a * hm1 + hm2
        k = a * km1 + km2
        conv.append((h, k))
        hm2, hm1 = hm1, h
        km2, km1 = km1, k
        x, y = y, x - a * y
    return a_list, conv


def lambda1_via_cf(n, lam, s_k1, s_k2):
    """lambda_1(L2)^2 by formula (T), plus the winning convergent index.

    Candidate j = -1 is the vector (n*S_K1, 0) (b = 0, theta = n).
    """
    a_list, conv = cf_convergents(lam, n)
    best = (n * s_k1) ** 2          # b = 0
    best_j = -1
    best_q = 0
    best_theta = n
    for j, (h, k) in enumerate(conv):
        theta = abs(k * lam - h * n)
        val = (theta * s_k1) ** 2 + (k * s_k2) ** 2
        if val < best:
            best, best_j, best_q, best_theta = val, j, k, theta
    return best, best_j, best_q, best_theta, a_list, conv


def k1k2(n, eff):
    """Same convention as glv_hnp_phase2_nuhat_control.py:k1k2 (eff = K1*K2/n)."""
    return max(2, int(eff * math.sqrt(n))), math.isqrt(n) + 1


def nuhat_exact(n, lam, eff):
    """nu_hat computed exactly (integer Gauss reduction), plus CF diagnostics."""
    k1, k2 = k1k2(n, eff)
    s_k1, _, s_k2, _ = scales(n, k1, k2)
    sv = gauss_reduce_2d_exact((n * s_k1, 0), (-lam * s_k1, s_k2))
    l1sq_gauss = n2(sv)
    l1sq_cf, j_star, q_star, theta_star, a_list, conv = lambda1_via_cf(
        n, lam, s_k1, s_k2)
    det = n * s_k1 * s_k2
    B = math.sqrt(n * s_k1 / s_k2)
    return {
        'n': n, 'lam': lam, 'eff': eff, 'k1': k1, 'k2': k2,
        's_k1': s_k1, 's_k2': s_k2, 'det': det, 'B': B,
        'l1sq_gauss': l1sq_gauss, 'l1sq_cf': l1sq_cf,
        'match': l1sq_gauss == l1sq_cf,
        'nu_hat': math.sqrt(l1sq_gauss / det),
        'j_star': j_star, 'q_star': q_star, 'theta_star': theta_star,
        'a_list': a_list, 'conv': conv,
        'x_star': q_star / B if B else 0.0,
        'a_next': a_list[j_star + 1] if 0 <= j_star + 1 < len(a_list) else None,
        'max_a': max(a_list) if a_list else 0,
    }


def nuhat_cf_approx(rec):
    """nu_hat by the two-term approximation (**): theta_j ~ n/q_{j+1}."""
    n, B = rec['n'], rec['B']
    conv, s_k1, s_k2, det = rec['conv'], rec['s_k1'], rec['s_k2'], rec['det']
    best = (n * s_k1) ** 2 / det          # j = -1 term, exact
    qs = [0] + [k for (_, k) in conv]
    for j in range(len(qs) - 1):
        q_j, q_next = qs[j], qs[j + 1]
        if q_next == 0:
            continue
        theta_approx = n / q_next
        val = ((theta_approx * s_k1) ** 2 + (q_j * s_k2) ** 2) / det
        best = min(best, val)
    return math.sqrt(best)


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def spearman(xs, ys):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(xs), ranks(ys)
    mx = sum(rx) / len(rx)
    my = sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else 0.0


def pearson(xs, ys):
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    dy = math.sqrt(sum((b - my) ** 2 for b in ys))
    return num / (dx * dy) if dx and dy else 0.0


# ---------------------------------------------------------------------------
# Curve sampling.  L2 depends only on (n, lam, K1, K2), never on p or on the
# curve equation, so no point counting is needed here.
# ---------------------------------------------------------------------------

def sample_glv_moduli(bits, count, rng):
    """Primes n = 1 mod 3 of the given size, with a GLV eigenvalue lam."""
    out = []
    while len(out) < count:
        cand = rng.randrange(1 << (bits - 1), 1 << bits)
        n = int(sympy.nextprime(cand))
        if n % 3 != 1:
            continue
        roots = glv_eigenvalues(n)
        if roots is None:
            continue
        out.append((n, roots[0]))
    return out


# ---------------------------------------------------------------------------
# E1 -- the closed form (T) is exact
# ---------------------------------------------------------------------------

def e1_exactness(rng):
    print("=" * 78)
    print("E1 -- closed form (T): lambda_1(L2) is attained at a CF convergent")
    print("=" * 78)
    print("For each sample, exact integer Lagrange-Gauss reduction is compared")
    print("against min over CF convergent denominators of lam/n.  A single")
    print("mismatch falsifies (T).\n")
    print(f"{'bits':>5} {'eff':>6} {'N':>5} {'CF == Gauss':>12} "
          f"{'float-Gauss == exact':>21}")
    total = bad = float_bad = 0
    detail = []
    for bits in (20, 24, 32, 64, 128, 256):
        moduli = sample_glv_moduli(bits, 40, rng)
        for eff in (0.05, 0.10, 0.25, 0.50):
            ok = fok = 0
            for n, lam in moduli:
                rec = nuhat_exact(n, lam, eff)
                ok += rec['match']
                # cross-check the float Gauss reduction actually used in July
                f = gauss_reduce_2d_float((n * rec['s_k1'], 0),
                                          (-lam * rec['s_k1'], rec['s_k2']))
                fok += abs(f - math.sqrt(rec['l1sq_gauss'])) <= \
                    1e-9 * math.sqrt(rec['l1sq_gauss'])
                detail.append(rec)
            total += len(moduli)
            bad += len(moduli) - ok
            float_bad += len(moduli) - fok
            print(f"{bits:>5} {eff:>6.2f} {len(moduli):>5} "
                  f"{ok}/{len(moduli):<10} {fok}/{len(moduli)}")
    print(f"\n  CF closed form (T): {total - bad}/{total} exact matches"
          f"   -> {'VERIFIED' if bad == 0 else 'FALSIFIED'}")
    print(f"  float Lagrange-Gauss vs exact: {total - float_bad}/{total} agree"
          f"   -> {'float version is safe at 256 bits' if float_bad == 0 else 'FLOAT VERSION UNRELIABLE'}")
    return detail


# ---------------------------------------------------------------------------
# E2 -- which convergent wins, and where it sits relative to B
# ---------------------------------------------------------------------------

def e2_which_convergent(detail):
    print()
    print("=" * 78)
    print("E2 -- the winning convergent sits at the scale B = sqrt(n*K2/K1)")
    print("=" * 78)
    xs = [r['x_star'] for r in detail if r['q_star'] > 0]
    print(f"  samples with q* > 0: {len(xs)}/{len(detail)}")
    xs_sorted = sorted(xs)
    def q(f):
        return xs_sorted[min(len(xs_sorted) - 1, int(f * len(xs_sorted)))]
    print(f"  x* = q*/B quantiles: q05={q(.05):.3f} q25={q(.25):.3f} "
          f"q50={q(.50):.3f} q75={q(.75):.3f} q95={q(.95):.3f}")
    within = sum(1 for x in xs if 0.1 <= x <= 10.0) / len(xs)
    print(f"  fraction with x* in [0.1, 10]: {within*100:.1f}%")
    print("  (x* clustered around 1 == the winning convergent is the one whose")
    print("   denominator is nearest B; this is the 'scale' in the proposal.)")
    # index j* as a fraction of the CF length
    fr = [r['j_star'] / max(1, len(r['a_list']) - 1) for r in detail
          if r['q_star'] > 0]
    print(f"  j*/len(CF): mean={sum(fr)/len(fr):.3f} "
          f"min={min(fr):.3f} max={max(fr):.3f}")


# ---------------------------------------------------------------------------
# E3 -- the stated falsifier: CF-predicted nu_hat vs computed nu_hat
# ---------------------------------------------------------------------------

def e3_falsifier(detail):
    print()
    print("=" * 78)
    print("E3 -- FALSIFIER: corr(predicted-from-CF nu_hat, computed nu_hat)")
    print("=" * 78)
    print("Threshold set 2026-07-29: correlation < 0.9 kills the story.\n")
    print(f"{'bits':>5} {'eff':>6} {'N':>5} {'pearson':>9} {'spearman':>9} "
          f"{'max rel err':>12}")
    by = {}
    for r in detail:
        by.setdefault((r['n'].bit_length() // 4 * 4, r['eff']), []).append(r)
    allp = []
    for key in sorted(by):
        grp = by[key]
        exact = [r['nu_hat'] for r in grp]
        approx = [nuhat_cf_approx(r) for r in grp]
        p = pearson(exact, approx)
        s = spearman(exact, approx)
        err = max(abs(a - e) / e for a, e in zip(approx, exact))
        allp.append(p)
        print(f"{key[0]:>5} {key[1]:>6.2f} {len(grp):>5} {p:>9.4f} "
              f"{s:>9.4f} {err:>12.4f}")
    print(f"\n  worst-cell pearson = {min(allp):.4f}  -> "
          f"{'SURVIVES (all >= 0.9)' if min(allp) >= 0.9 else 'FALSIFIED'}")


# ---------------------------------------------------------------------------
# E4 -- scale-dependent vs scale-free partial quotients
# ---------------------------------------------------------------------------

def e4_partial_quotient(detail):
    print()
    print("=" * 78)
    print("E4 -- why the June scale-free CF invariants had to fail")
    print("=" * 78)
    print("Claim: nu_hat is governed by the partial quotient of lam/n *at the")
    print("scale B* (a_{j*+1}), not by any global CF statistic such as max_a.\n")
    print(f"{'bits':>5} {'eff':>6} {'N':>5} {'sp(nu,1/sqrt(a_next))':>22} "
          f"{'sp(nu,max_a) [June]':>21}")
    by = {}
    for r in detail:
        if r['a_next']:
            by.setdefault((r['n'].bit_length() // 4 * 4, r['eff']), []).append(r)
    rows = []
    for key in sorted(by):
        grp = by[key]
        nus = [r['nu_hat'] for r in grp]
        an = [1.0 / math.sqrt(r['a_next']) for r in grp]
        ma = [float(r['max_a']) for r in grp]
        s1 = spearman(nus, an)
        s2 = spearman(nus, ma)
        rows.append((s1, s2))
        print(f"{key[0]:>5} {key[1]:>6.2f} {len(grp):>5} {s1:>22.3f} "
              f"{s2:>21.3f}")
    m1 = sum(abs(a) for a, _ in rows) / len(rows)
    m2 = sum(abs(b) for _, b in rows) / len(rows)
    print(f"\n  mean |spearman|:  at-scale a_next = {m1:.3f}   "
          f"scale-free max_a = {m2:.3f}")
    # binned response of nu_hat to a_next
    pooled = [r for grp in by.values() for r in grp]
    print("\n  nu_hat response to the at-scale partial quotient a_{j*+1}:")
    print(f"    {'a_next':>10} {'N':>6} {'mean nu_hat':>12} {'2/sqrt(a)':>11}")
    buckets = [(1, 1), (2, 2), (3, 4), (5, 8), (9, 16), (17, 10 ** 9)]
    for lo, hi in buckets:
        grp = [r for r in pooled if lo <= r['a_next'] <= hi]
        if not grp:
            continue
        mid = math.sqrt(lo * min(hi, 32))
        label = f"{lo}" if lo == hi else f"{lo}-{hi if hi < 10**9 else '+'}"
        print(f"    {label:>10} {len(grp):>6} "
              f"{sum(r['nu_hat'] for r in grp)/len(grp):>12.3f} "
              f"{2/math.sqrt(mid):>11.3f}")


# ---------------------------------------------------------------------------
# E5 -- secp256k1 under the closed form
# ---------------------------------------------------------------------------

def e5_secp256k1(rng):
    print()
    print("=" * 78)
    print("E5 -- secp256k1 under the closed form")
    print("=" * 78)
    n = SECP_N
    roots = glv_eigenvalues(n)
    assert roots is not None
    lam = roots[0]
    assert (lam * lam + lam + 1) % n == 0
    print(f"  n   = {n}")
    print(f"  lam = {lam}")
    print(f"  mu  = {mu_of(lam, n):.6f}")
    a_list, _ = cf_convergents(lam, n)
    print(f"  CF of lam/n: length {len(a_list)}, max partial quotient "
          f"{max(a_list)} (scale-free -- the June invariant)")
    print()
    print(f"  {'eff':>6} {'nu_hat':>9} {'j*':>4} {'q*/B':>7} {'a_next':>7} "
          f"{'CF-approx':>10}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        rec = nuhat_exact(n, lam, eff)
        assert rec['match'], "closed form failed on secp256k1"
        print(f"  {eff:>6.2f} {rec['nu_hat']:>9.4f} {rec['j_star']:>4} "
              f"{rec['x_star']:>7.3f} {str(rec['a_next']):>7} "
              f"{nuhat_cf_approx(rec):>10.4f}")
    print()
    print("  Reference: 2026-07-29 reported nu_hat = 0.8709 / 0.6624 / 0.5852 /")
    print("  0.6851 at eff = 0.05 / 0.10 / 0.25 / 0.50 (float Gauss reduction).")
    print()
    print("  Null distribution of the at-scale partial quotient a_{j*+1} over")
    print("  random 256-bit lam (eff = 0.05):")
    sample = []
    for _ in range(2000):
        lam_r = rng.randrange(1, n)
        rec = nuhat_exact(n, lam_r, 0.05)
        if rec['a_next']:
            sample.append((rec['a_next'], rec['nu_hat']))
    sample.sort()
    cut = [a for a, _ in sample]
    def qq(f):
        return cut[min(len(cut) - 1, int(f * len(cut)))]
    print(f"    a_next quantiles: q25={qq(.25)} q50={qq(.50)} q75={qq(.75)} "
          f"q90={qq(.90)} q99={qq(.99)}")
    lo = [nu for a, nu in sample if a >= 8]
    hi = [nu for a, nu in sample if a == 1]
    print(f"    a_next >= 8 (N={len(lo)}): mean nu_hat = "
          f"{sum(lo)/len(lo):.3f}")
    print(f"    a_next == 1 (N={len(hi)}): mean nu_hat = "
          f"{sum(hi)/len(hi):.3f}")
    frac = sum(1 for _, nu in sample if nu <= 0.645) / len(sample)
    print(f"    fraction of random lam below the C2 ceiling nu_hat = 0.645: "
          f"{frac*100:.1f}%")


def main():
    t0 = time.time()
    rng = random.Random(20260802)
    print("GLV-HNP Thread 23 -- closed form for lambda_1(L2) / nu_hat")
    print("Continuation of 2026-07-29 (e845207).  No LLL is used anywhere in")
    print("this script: every number below is exact integer arithmetic or a")
    print("Lagrange-Gauss reduction of a 2-dimensional lattice.\n")
    detail = e1_exactness(rng)
    e2_which_convergent(detail)
    e3_falsifier(detail)
    e4_partial_quotient(detail)
    e5_secp256k1(rng)
    print(f"\nDone in {time.time() - t0:.1f}s.")
    return 0


if __name__ == '__main__':
    sys.exit(main())

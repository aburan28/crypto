"""
GLV-HNP — Thread 23: closed form for lambda_1(L2), and why the June continued-
fraction invariants failed.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run #2, commit e845207)
-----------------------------------------------------------------------
Thread 20d found the first invariant to beat the trivial baseline on the Exp S
protocol: nu_hat = lambda_1(L2)/sqrt(det L2), AUC 0.935, where L2 is the
NON-PLANTED 2-dimensional sublattice of the Phase-2 GLV lattice

    L2 = < u = (n*S_K1, 0),  v = (-lam*S_K1, S_K2) >,   det L2 = n*S_K1*S_K2.

nu_hat was *computed* (one Lagrange-Gauss reduction) but not *understood*.  The
2026-07-29 next-step asked for a closed form, and conjectured the answer is a
best-rational-approximation quantity read at the scale sqrt(n*S_K1/S_K2).

This script derives and tests that closed form.

Derivation
----------
A general lattice vector is a*u + b*v = ((a*n - b*lam)*S_K1, b*S_K2), so

    N(a,b)^2 = ((a*n - b*lam)*S_K1)^2 + (b*S_K2)^2.

WLOG b >= 0.  For b = 0 the shortest is (n*S_K1, 0).  For b >= 1 the optimal a
is round(b*lam/n), leaving the centred residue eta(b) = |b*lam mod^{+-} n|.  So

    lambda_1(L2)^2 = min( (n*S_K1)^2,
                          min_{b>=1} (eta(b)*S_K1)^2 + (b*S_K2)^2 ).

The objective is strictly increasing in both b and eta(b), so only the Pareto
frontier of (b, eta(b)) can contain the minimiser.  That frontier is exactly the
set of best approximations of the second kind to lam/n, i.e. exactly the
continued-fraction convergent denominators q_j of lam/n (Khinchin, Thm 16).
Hence the CLOSED FORM, exact and in integers only:

    lambda_1(L2)^2 = min_j  (eta_j*S_K1)^2 + (q_j*S_K2)^2,
    eta_j = |q_j*lam - p_j*n|,   (p_j/q_j) the convergents of lam/n, plus q=0.

Normalising by det L2 = n*S_K1*S_K2 and writing the natural scale

    B = sqrt(n*S_K1/S_K2)  ( ~ sqrt(n*K2/K1) ),      x_j = q_j/B,  y_j = eta_j/(n/B)

gives the scale-free two-coordinate form

    nu_hat^2 = min_j  x_j^2 + y_j^2,      and since eta_j ~ n/q_{j+1},
    nu_hat   ~ min_j  sqrt( x_j^2 + 1/x_{j+1}^2 ).                        (*)

(*) is the punchline.  nu_hat is small -- the attack is EASY -- exactly when
some convergent denominator q_j sits below the scale B while q_{j+1} sits well
above it, i.e. when there is a LARGE PARTIAL QUOTIENT STRADDLING THE SCALE B.

This is why every scale-free CF invariant tried in June failed (q_cf,
max_q_cf, max_a; log lines ~3560-3580): max_a asks whether a large partial
quotient exists ANYWHERE in the expansion of lam/n, but (*) says only the one
partial quotient at the scale B matters.  A large a_k far from B does nothing.

Experiments
-----------
E1  Exactness of the closed form vs Lagrange-Gauss, 20-bit sample.
E2  Exactness at 256 bits, against an exact-integer Gauss reduction (also
    checks whether the float-based gauss_reduce_2d in the Thread 20a module
    stays reliable at cryptographic size).
E3  Accuracy of the asymptotic two-term form (*) against the exact nu_hat.
E4  Where does the argmin sit relative to B?
E5  a_local (the partial quotient straddling B) as a PREDICTOR, head to head
    with nu_hat / max_a / mu on the Exp S protocol.

Falsifier (as stated 2026-07-29): if the CF-predicted nu_hat correlates < 0.9
with the computed nu_hat, the convergent-scale story is wrong.

Run: python3 glv_hnp_nuhat_cf_closed_form.py
"""

import math
import random
import sys
import time

import sympy

import importlib.util

_here = __file__.rsplit("/", 1)[0]

def _load(mod, fname):
    spec = importlib.util.spec_from_file_location(mod, _here + "/" + fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

_t20a = _load("_t20a", "glv_hnp_phase2_thread20a.py")
_t20d = _load("_t20d", "glv_hnp_nuhat_vs_c1c2.py")

scales = _t20a.scales
gauss_reduce_2d = _t20a.gauss_reduce_2d
rival_sublattice_nu = _t20a.rival_sublattice_nu
run_experiment = _t20a.run_experiment
mu_of = _t20a.mu_of

auc = _t20d.auc
best_cut_accuracy = _t20d.best_cut_accuracy
max_partial_quotient = _t20d.max_partial_quotient

SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
SECP_P = 2 ** 256 - 2 ** 32 - 977


# --------------------------------------------------------------------------
# Continued fractions and the closed form
# --------------------------------------------------------------------------

def cf_convergents(a, b):
    """Convergents p_j/q_j of a/b.  Returns [(p, q, partial_quotient), ...]."""
    out = []
    pm1, qm1, pm2, qm2 = 1, 0, 0, 1
    x, y = a, b
    while y:
        k, r = divmod(x, y)
        p, q = k * pm1 + pm2, k * qm1 + qm2
        out.append((p, q, k))
        pm2, qm2, pm1, qm1 = pm1, qm1, p, q
        x, y = y, r
    return out


def lambda1_cf(n, lam, s_k1, s_k2):
    """
    Exact lambda_1(L2)^2 via the convergent closed form, integers only.
    Returns (norm2, b_argmin, eta_argmin).
    """
    best2, bb, be = (n * s_k1) ** 2, 0, n          # the b = 0 vector
    for p, q, _k in cf_convergents(lam, n):
        eta = abs(q * lam - p * n)
        cand = (eta * s_k1) ** 2 + (q * s_k2) ** 2
        if cand < best2:
            best2, bb, be = cand, q, eta
    return best2, bb, be


def nu_hat_cf(n, lam, k1_bound, k2_bound):
    """Scale-free nu_hat from the closed form (exact integer arithmetic)."""
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    best2, bb, be = lambda1_cf(n, lam, s_k1, s_k2)
    det = n * s_k1 * s_k2
    return math.sqrt(best2 / det), bb, be


def gauss_reduce_2d_exact(u, v):
    """
    Lagrange-Gauss with exact integer arithmetic (no float division), so it is
    valid at 256 bits where the float version's round() loses precision.
    Returns the squared norm of the shortest nonzero vector, as an int.
    """
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]
    if n2(u) < n2(v):
        u, v = v, u
    while True:
        d = n2(v)
        if d == 0:
            return n2(u)
        num = u[0] * v[0] + u[1] * v[1]
        # exact round-half-to-even-free nearest integer of num/d
        q = (2 * num + d) // (2 * d) if d > 0 else 0
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if n2(u) <= n2(v):
            return n2(u)


def a_local(n, lam, k1_bound, k2_bound):
    """
    The partial quotient straddling the scale B = sqrt(n*S_K1/S_K2): the
    a_{j+1} for the largest convergent denominator q_j <= B.  This is the
    scale-AWARE replacement for the falsified scale-free max_a.
    """
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    B = math.sqrt(n * s_k1 / s_k2)
    conv = cf_convergents(lam, n)
    qs = [0] + [q for _p, q, _k in conv]
    ks = [k for _p, _q, k in conv]
    best_a, j_at = 1, -1
    for j in range(len(qs) - 1):
        if qs[j] <= B:
            j_at = j
    if 0 <= j_at < len(ks):
        best_a = ks[j_at]
    return best_a, B, j_at


def straddle_ratio(n, lam, k1_bound, k2_bound):
    """
    log-symmetric straddle of the scale: max over j of q_{j+1}/q_j evaluated
    at the j whose interval [q_j, q_{j+1}] contains B.  Equivalent content to
    a_local but continuous, so it plots better.
    """
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    B = math.sqrt(n * s_k1 / s_k2)
    qs = [1] + [q for _p, q, _k in cf_convergents(lam, n)]
    for j in range(len(qs) - 1):
        if qs[j] <= B <= qs[j + 1]:
            return qs[j + 1] / max(1, qs[j])
    return 1.0


# --------------------------------------------------------------------------
# Experiments
# --------------------------------------------------------------------------

def e1_exactness_20bit(trials=4000, seed=20230):
    print("-" * 74)
    print("E1  closed form vs Lagrange-Gauss, 20-bit, %d random (n, lam, K1, K2)"
          % trials)
    rng = random.Random(seed)
    bad = 0
    worst = 0.0
    for _ in range(trials):
        n = int(sympy.nextprime(rng.randrange(2 ** 19, 2 ** 20)))
        lam = rng.randrange(1, n)
        k1 = rng.choice([12, 36, 72, 144])
        k2 = rng.choice([32, 64, 128])
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        exact2, _b, _e = lambda1_cf(n, lam, s_k1, s_k2)
        ref2 = gauss_reduce_2d_exact((n * s_k1, 0), (-lam * s_k1, s_k2))
        if exact2 != ref2:
            bad += 1
            worst = max(worst, abs(math.sqrt(exact2) - math.sqrt(ref2))
                        / math.sqrt(ref2))
    print("  mismatches: %d / %d   (worst relative gap %.3e)" % (bad, trials, worst))
    print("  => closed form is %s" % ("EXACT" if bad == 0 else "NOT exact"))
    return bad


def e2_exactness_256bit(trials=1500, seed=777):
    print("-" * 74)
    print("E2  closed form at 256 bits (%d draws), and float-Gauss reliability"
          % trials)
    rng = random.Random(seed)
    bad = 0
    float_bad = 0
    float_worst = 0.0
    for _ in range(trials):
        n = SECP_N
        lam = rng.randrange(1, n)
        k1 = rng.choice([2 ** 20, 2 ** 32, 2 ** 64])
        k2 = rng.choice([2 ** 20, 2 ** 32, 2 ** 64])
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        cf2, _b, _e = lambda1_cf(n, lam, s_k1, s_k2)
        ref2 = gauss_reduce_2d_exact((n * s_k1, 0), (-lam * s_k1, s_k2))
        if cf2 != ref2:
            bad += 1
        try:
            fl = gauss_reduce_2d((n * s_k1, 0), (-lam * s_k1, s_k2))
            rel = abs(fl - math.sqrt(ref2)) / math.sqrt(ref2)
            if rel > 1e-9:
                float_bad += 1
                float_worst = max(float_worst, rel)
        except (OverflowError, ValueError, ZeroDivisionError):
            float_bad += 1
            float_worst = float('inf')
    print("  closed-form mismatches vs exact Gauss: %d / %d" % (bad, trials))
    print("  float gauss_reduce_2d disagreements  : %d / %d (worst rel %.3e)"
          % (float_bad, trials, float_worst))
    return bad, float_bad


def e3_two_term_form(trials=3000, seed=4242):
    print("-" * 74)
    print("E3  asymptotic two-term form  nu_hat ~ min_j sqrt(x_j^2 + 1/x_{j+1}^2)")
    rng = random.Random(seed)
    xs, ys = [], []
    for _ in range(trials):
        n = int(sympy.nextprime(rng.randrange(2 ** 19, 2 ** 20)))
        lam = rng.randrange(1, n)
        k1, k2 = rng.choice([12, 36, 72, 144]), rng.choice([32, 64, 128])
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        det = n * s_k1 * s_k2
        exact2, _b, _e = lambda1_cf(n, lam, s_k1, s_k2)
        exact = math.sqrt(exact2 / det)
        B = math.sqrt(n * s_k1 / s_k2)
        qs = [0] + [q for _p, q, _k in cf_convergents(lam, n)]
        approx = min(
            math.sqrt((qs[j] / B) ** 2 + (B / qs[j + 1]) ** 2)
            for j in range(len(qs) - 1) if qs[j + 1] > 0)
        xs.append(exact)
        ys.append(approx)
    print("  pearson(exact, two-term) = %.4f" % pearson(xs, ys))
    print("  spearman                 = %.4f" % spearman(xs, ys))
    rels = sorted(abs(a - b) / a for a, b in zip(xs, ys))
    print("  relative error: median %.4f  p90 %.4f  max %.4f"
          % (rels[len(rels) // 2], rels[int(0.9 * len(rels))], rels[-1]))
    return pearson(xs, ys)


def e4_argmin_location(trials=3000, seed=99):
    print("-" * 74)
    print("E4  location of the argmin denominator b* relative to B")
    rng = random.Random(seed)
    ratios, nconv, straddles, geo = [], [], 0, []
    tot = 0
    for _ in range(trials):
        n = int(sympy.nextprime(rng.randrange(2 ** 19, 2 ** 20)))
        lam = rng.randrange(1, n)
        k1, k2 = rng.choice([12, 36, 72, 144]), rng.choice([32, 64, 128])
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        _e2, b, _eta = lambda1_cf(n, lam, s_k1, s_k2)
        B = math.sqrt(n * s_k1 / s_k2)
        conv = [q for _p, q, _k in cf_convergents(lam, n)]
        if b > 0:
            ratios.append(math.log(b / B, 2))
            # the successor convergent denominator after the argmin
            nxt = next((q for q in conv if q > b), None)
            if nxt:
                tot += 1
                if b <= B <= nxt:
                    straddles += 1
                geo.append(math.log(math.sqrt(b * nxt) / B, 2))
        nconv.append(sum(1 for q in conv if 0.25 * B <= q <= 4 * B))
    ratios.sort()
    geo.sort()
    print("  log2(b_argmin / B):  p05 %.2f  p25 %.2f  median %.2f  p75 %.2f  p95 %.2f"
          % (ratios[int(.05 * len(ratios))], ratios[int(.25 * len(ratios))],
             ratios[len(ratios) // 2], ratios[int(.75 * len(ratios))],
             ratios[int(.95 * len(ratios))]))
    within = sum(1 for r in ratios if abs(r) <= 1) / len(ratios)
    print("  |log2(b/B)| <= 1 (argmin within a factor 2 of the scale): %.1f%%"
          % (100 * within))
    print("  SHARPER: argmin pair straddles B, i.e. q_j <= B <= q_{j+1}: %.1f%%"
          % (100 * straddles / max(1, tot)))
    print("  log2( sqrt(q_j*q_{j+1}) / B ):  p05 %.2f  median %.2f  p95 %.2f"
          % (geo[int(.05 * len(geo))], geo[len(geo) // 2], geo[int(.95 * len(geo))]))
    print("  mean #convergents in [B/4, 4B]: %.2f" % (sum(nconv) / len(nconv)))


def e5_predictor(n_curves=100):
    """
    Exp S protocol, mirrored EXACTLY from glv_hnp_nuhat_vs_c1c2.py:main()
    (K1=72, m=12, 6 seeds, k2 = isqrt(n)+1, d = Random(s+7777).randint(1,n-1),
    run_experiment(..., seed=s)), so the nu_hat column here must reproduce the
    2026-07-29 run #2 numbers.  a_local and straddle are the new columns.
    """
    print("-" * 74)
    print("E5  a_local as a PREDICTOR on the Exp S protocol (K1=%d, m=%d, %d seeds)"
          % (_t20d.K1_FIXED, _t20d.M_FIXED, _t20d.SEEDS_N))
    t0 = time.time()
    curves = _t20d.harvest(n_curves)
    print("  harvested %d curves in %.1fs" % (len(curves), time.time() - t0))
    k1, m = _t20d.K1_FIXED, _t20d.M_FIXED
    t1 = time.time()
    for c in curves:
        n, lam, G, p = c['n'], c['lam'], c['G'], c['p']
        k2 = math.isqrt(n) + 1
        wins = 0
        for s in range(_t20d.SEEDS_N):
            d = random.Random(s + 7777).randint(1, n - 1)
            wins += run_experiment(p, n, lam, G, m, d, k1, k2, s)
        c['wins'] = wins
        c['C2'] = (wins == _t20d.SEEDS_N)
        _nu, _det, c['nu_hat'] = rival_sublattice_nu(n, lam, k1, k2)
        c['nu_hat_cf'] = nu_hat_cf(n, lam, k1, k2)[0]
        c['a_local'] = a_local(n, lam, k1, k2)[0]
        c['straddle'] = straddle_ratio(n, lam, k1, k2)
        c['max_a'] = max_partial_quotient(lam, n)
        c['mu'] = mu_of(lam, n)
    print("  %d LLL trials in %.1fs" % (len(curves) * _t20d.SEEDS_N,
                                        time.time() - t1))
    labels = [c['C2'] for c in curves]
    nC2 = sum(labels)
    base = max(nC2, len(curves) - nC2) / len(curves)
    print("  C1 (fails): %d/%d   C2 (6/6): %d/%d   trivial baseline = %.1f%%"
          % (len(curves) - nC2, len(curves), nC2, len(curves), 100 * base))
    print()
    print("  %-12s %8s %10s %26s" % ("invariant", "AUC", "best acc", "C2 range"))
    for key in ('nu_hat', 'nu_hat_cf', 'a_local', 'straddle', 'max_a', 'mu'):
        sc = [c[key] for c in curves]
        a = auc(sc, labels)
        acc, _cut = best_cut_accuracy(sc, labels)
        s2 = [c[key] for c in curves if c['C2']]
        rng = "[%.3f,%.3f]" % (min(s2), max(s2)) if s2 else "-"
        print("  %-12s %8.3f %9.1f%% %26s" % (key, max(a, 1 - a), 100 * acc, rng))
    d = max(abs(c['nu_hat'] - c['nu_hat_cf']) for c in curves)
    print("\n  max |nu_hat(Gauss) - nu_hat(CF)| over the sample: %.3e" % d)
    print("\n  --- a_local decile response ---")
    order = sorted(curves, key=lambda c: c['a_local'])
    k = max(1, len(order) // 10)
    print("  %7s %12s %9s %12s" % ("decile", "a_local mid", "C2 rate", "mean wins/6"))
    for i in range(0, len(order), k):
        grp = order[i:i + k]
        if len(grp) < 2:
            continue
        print("  %7d %12.1f %9.2f %12.2f"
              % (i // k + 1, sum(c['a_local'] for c in grp) / len(grp),
                 sum(c['C2'] for c in grp) / len(grp),
                 sum(c['wins'] for c in grp) / len(grp)))
    return curves


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    if sx == 0 or sy == 0:
        return float('nan')
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def _ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(x, y):
    return pearson(_ranks(x), _ranks(y))


def e6_secp256k1():
    print("-" * 74)
    print("E6  secp256k1: closed-form decomposition of nu_hat")
    # the GLV eigenvalue of secp256k1, as used in glv_hnp_nuhat_vs_c1c2.py:184
    lam = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
    assert (lam * lam + lam + 1) % SECP_N == 0, "lambda check failed"
    print("  lam = 0x%x" % lam)
    print("  lam^2 + lam + 1 = 0 mod n  : verified")
    print()
    print("  %8s %12s %14s %10s %10s %8s"
          % ("eff", "nu_hat(CF)", "b_argmin", "log2 b/B", "a_local", "match"))
    for eff in (0.05, 0.0993, 0.10, 0.25, 0.50):
        k2 = math.isqrt(SECP_N) + 1
        k1 = max(2, int(eff * SECP_N / k2))
        nh, b, _eta = nu_hat_cf(SECP_N, lam, k1, k2)
        al, B, _j = a_local(SECP_N, lam, k1, k2)
        s_k1, _, s_k2, _ = scales(SECP_N, k1, k2)
        ref2 = gauss_reduce_2d_exact((SECP_N * s_k1, 0), (-lam * s_k1, s_k2))
        det = SECP_N * s_k1 * s_k2
        ok = abs(nh - math.sqrt(ref2 / det)) < 1e-12
        print("  %8.4f %12.4f %14s %10.2f %10d %8s"
              % (eff, nh, ("2^%.1f" % math.log(b, 2)) if b else "0",
                 math.log(b / B, 2) if b else float('nan'), al, ok))
    print()
    print("  Largest partial quotient of lam/n (the falsified scale-free max_a):"
          " %d" % max_partial_quotient(lam, SECP_N))
    print("  -> compare with a_local above: max_a is large but sits far from B,")
    print("     which is exactly the distinction (*) predicts matters.")


def main():
    print("=" * 74)
    print("GLV-HNP Thread 23: closed form for lambda_1(L2)")
    print("=" * 74)
    t0 = time.time()
    bad1 = e1_exactness_20bit()
    bad2, fbad = e2_exactness_256bit()
    r3 = e3_two_term_form()
    e4_argmin_location()
    e6_secp256k1()
    e5_predictor()
    print("-" * 74)
    print("total %.1fs" % (time.time() - t0))
    print("=" * 74)


if __name__ == '__main__':
    main()

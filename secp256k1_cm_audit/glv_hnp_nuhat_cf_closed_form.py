"""
GLV-HNP Thread 23 — closed form for nu_hat via continued fractions.

Background (2026-07-29, commit e845207).  The Phase 2 lattice contains m copies
of a non-planted 2D sublattice

    L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >,        det L2 = n*S_K1*S_K2

and the scale-free  nu_hat = lambda_1(L2)/sqrt(det L2)  was measured to separate
the June C1/C2 recovery classes at AUC 0.935, where six algebraic curve
invariants (delta/n, kappa(M), q_cf, max_q_cf, max_a, a_corn/n) and two
lambda-ratio invariants (lam/n, mu) had all failed.  nu_hat was COMPUTED
(Lagrange-Gauss), never UNDERSTOOD.  This script derives it.

THEORY.  A general vector of L2 is
    a*(n*S1, 0) + b*(-lam*S1, S2) = ( (a*n - b*lam)*S1,  b*S2 ),
so with  ||x||_n := min_a |x - a*n|  (centered residue),

    lambda_1(L2)^2 = min(  (n*S1)^2 ,  min_{b>=1} [ S1^2*||b*lam||_n^2 + S2^2*b^2 ]  ).

Any b attaining the inner min is a best approximation of the second kind to
lam/n: if some q < b had ||q*lam||_n <= ||b*lam||_n, the vector at q would be
strictly shorter in BOTH coordinates.  By Khinchin Thm 16 every best
approximation of the second kind is a convergent denominator.  Hence

    (H1)  lambda_1(L2) is attained at a CONTINUED-FRACTION CONVERGENT of lam/n.

Normalising by sqrt(det) = sqrt(n*S1*S2), write for the j-th convergent q_j
    t_j = q_j * sqrt(S2/(n*S1)),        u_j = ||q_j*lam||_n * sqrt(S1/(n*S2)),
so that  nu_hat = min_j sqrt(t_j^2 + u_j^2).  Since ||q_j*lam||_n ~ n/q_{j+1},

    (H2)  u_j ~ 1/t_{j+1},   hence   nu_hat ~ min_j sqrt(t_j^2 + t_{j+1}^{-2}).

t_j = 1 at the CRITICAL SCALE q* = sqrt(n*S1/S2) = sqrt(n*K2/K1): the active
convergent is the one straddling q*, NOT the largest or the global-max one.

Finally q_{j+1} ~ a_{j+1}*q_j gives t_{j+1} ~ a_{j+1}*t_j, and minimising
sqrt(t^2 + 1/(a*t)^2) over the continuum t gives t = 1/sqrt(a) and

    (H3)  nu_hat ~ sqrt(2/a*),   a* = the partial quotient AT THE CRITICAL SCALE.

(H3) is the payoff: it converts the empirical nu_hat into an algebraic quantity,
and it retro-explains the June falsifications — q_cf/max_q_cf/max_a are all
SCALE-FREE statistics of the whole expansion, whereas the geometry only ever
reads one LOCAL partial quotient, the one at q* = sqrt(n*K2/K1).

EXPERIMENTS
  C1  exactness of (H1): min-over-convergents vs Lagrange-Gauss, exact integers.
  C2  localisation of the active convergent around q*.
  C3  accuracy of the two-term form (H2).
  C4  accuracy of (H3), and local a* vs global max_a as predictors of nu_hat.
  C5  does a* predict ACTUAL LLL recovery?  (the June question, answered
      algebraically) — real 20-bit j=0 CM curves, LLL trials.
  C6  secp256k1 placement + consistency check against the 2026-07-29 numbers.

Run: python3 glv_hnp_nuhat_cf_closed_form.py
"""

import importlib.util
import math
import os
import random
import time

import sympy

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "_t20a", os.path.join(_HERE, "glv_hnp_phase2_lambda_threshold_20a.py"))
_t20a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20a)

scales = _t20a.scales
rival_sublattice_nu = _t20a.rival_sublattice_nu
glv_eigenvalues = _t20a.glv_eigenvalues
mu_of = _t20a.mu_of
eisenstein_decompose = _t20a.eisenstein_decompose
j0_traces = _t20a.j0_traces
identify_twist = _t20a.identify_twist
find_generator = _t20a.find_generator
run_experiment = _t20a.run_experiment

SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
SECP_LAM = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72


# ---------------------------------------------------------------------------
# Exact 2D shortest vector (integer arithmetic, no floats)
# ---------------------------------------------------------------------------

def gauss_reduce_2d_exact(u, v):
    """Lagrange-Gauss; returns the exact INTEGER squared norm of lambda_1."""
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]
    if n2(u) < n2(v):
        u, v = v, u
    while True:
        nv = n2(v)
        if nv == 0:
            return n2(u)
        # round-half-to-even on the exact rational <u,v>/<v,v>
        num = u[0] * v[0] + u[1] * v[1]
        q = (2 * num + nv) // (2 * nv)          # floor(num/nv + 1/2), exact
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if n2(u) <= n2(v):
            return n2(u)


def nu_exact(n, lam, k1_bound, k2_bound):
    """Exact (lambda_1^2, det) of L2."""
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    l1sq = gauss_reduce_2d_exact((n * s1, 0), (-lam * s1, s2))
    return l1sq, n * s1 * s2


# ---------------------------------------------------------------------------
# Continued fractions of lam/n
# ---------------------------------------------------------------------------

def cf_convergents(lam, n):
    """
    Convergent denominators q_j of lam/n together with the partial quotients.
    Returns list of (q_j, a_{j+1}) where a_{j+1} is the NEXT partial quotient
    (the one that governs how well q_j approximates), plus the raw a-list.
    """
    a_list = []
    x, y = lam, n
    while y:
        a_list.append(x // y)
        x, y = y, x % y
    # q_{-1} = 0, q_0 = 1, q_j = a_j*q_{j-1} + q_{j-2}   (a_0 = 0 since lam < n)
    qs = [1]
    prev, cur = 0, 1
    for a in a_list[1:]:
        prev, cur = cur, a * cur + prev
        qs.append(cur)
    # q_j is governed by a_{j+1}, since |q_j*lam - p_j*n| ~ n/q_{j+1}
    pairs = [(q, a_list[idx + 1] if idx + 1 < len(a_list) else None)
             for idx, q in enumerate(qs)]
    return qs, a_list, pairs


def centered_residue(x, n):
    """||x||_n = min_a |x - a*n|, as a nonnegative integer."""
    r = x % n
    return min(r, n - r)


def nu_from_cf(n, lam, k1_bound, k2_bound, use_semiconvergents=False):
    """
    Minimise S1^2*||b*lam||_n^2 + S2^2*b^2 over convergent denominators b,
    against the b=0 candidate (n*S1).  Returns (best_sq, best_q, det).
    """
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    qs, a_list, _ = cf_convergents(lam, n)
    best_sq = (n * s1) ** 2
    best_q = 0
    cands = list(qs)
    if use_semiconvergents:
        prev, cur = 0, 1
        for idx, a in enumerate(a_list[1:]):
            for t in range(1, a):
                cands.append(t * cur + prev)
            prev, cur = cur, a * cur + prev
    for q in cands:
        if q == 0 or q >= n:
            continue
        r = centered_residue(q * lam, n)
        val = (s1 * r) ** 2 + (s2 * q) ** 2
        if val < best_sq:
            best_sq, best_q = val, q
    return best_sq, best_q, n * s1 * s2


def critical_scale(n, k1_bound, k2_bound):
    """q* = sqrt(n*S1/S2), the denominator scale at which t_j = 1."""
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    return math.sqrt(n * s1 / s2)


def local_partial_quotient(n, lam, k1_bound, k2_bound):
    """
    a* = the partial quotient at the critical scale: the a_{j+1} of the last
    convergent with q_j <= q*.  Also returns that convergent index and q_j.
    """
    qstar = critical_scale(n, k1_bound, k2_bound)
    qs, a_list, pairs = cf_convergents(lam, n)
    a_star, j_star, q_j = 1, 0, 1
    for idx, (q, a_next) in enumerate(pairs):
        if q <= qstar:
            a_star = a_next if a_next is not None else 1
            j_star, q_j = idx, q
        else:
            break
    return a_star, j_star, q_j, qstar


def max_partial_quotient(lam, n):
    _, a_list, _ = cf_convergents(lam, n)
    return max(a_list[1:]) if len(a_list) > 1 else 1


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else 0.0


def pearson(xs, ys):
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    dy = math.sqrt(sum((b - my) ** 2 for b in ys))
    return num / (dx * dy) if dx and dy else 0.0


def auc(scores, labels):
    """AUC for `scores` predicting label==1 (higher score -> label 1)."""
    pos = [s for s, l in zip(scores, labels) if l == 1]
    neg = [s for s, l in zip(scores, labels) if l == 0]
    if not pos or not neg:
        return float('nan')
    wins = sum((1.0 if a > b else 0.5 if a == b else 0.0) for a in pos for b in neg)
    return wins / (len(pos) * len(neg))


# ---------------------------------------------------------------------------
# GLV (n, lam) pair sampling — no curve needed for C1-C4
# ---------------------------------------------------------------------------

def sample_glv_pairs(bits, count, seed=1):
    """Primes n = 1 mod 3 of the given bit size with lam^2+lam+1 = 0 mod n."""
    rng = random.Random(seed)
    out = []
    lo, hi = 1 << (bits - 1), (1 << bits) - 1
    tries = 0
    while len(out) < count and tries < count * 400:
        tries += 1
        n = int(sympy.nextprime(rng.randint(lo, hi)))
        if n > hi or n % 3 != 1:
            continue
        roots = glv_eigenvalues(n)
        if roots is None:
            continue
        out.append((n, roots[0]))
    return out


def k1k2_for(n, eff):
    k2 = math.isqrt(n) + 1
    k1 = max(2, int(eff * n / k2))
    return k1, k2


# ---------------------------------------------------------------------------
# C1 — exactness of (H1)
# ---------------------------------------------------------------------------

def exp_c1():
    print("\n" + "=" * 78)
    print("C1 — is lambda_1(L2) attained at a CF convergent of lam/n?  (H1)")
    print("=" * 78)
    print("Exact integer comparison: min-over-convergents  vs  Lagrange-Gauss.")
    print(f"{'bits':>5} {'eff':>6} {'pairs':>6} {'exact match':>12} {'pearson':>9} "
          f"{'spearman':>9} {'max rel err':>12}")
    total, agree = 0, 0
    for bits in (20, 24, 32, 64, 128, 256):
        pairs = sample_glv_pairs(bits, 60 if bits <= 32 else 40, seed=bits)
        for eff in (0.05, 0.10, 0.25):
            nh_g, nh_c, ok, worst = [], [], 0, 0.0
            for n, lam in pairs:
                k1, k2 = k1k2_for(n, eff)
                g_sq, det = nu_exact(n, lam, k1, k2)
                c_sq, _, det2 = nu_from_cf(n, lam, k1, k2)
                assert det == det2
                if g_sq == c_sq:
                    ok += 1
                else:
                    worst = max(worst, abs(c_sq - g_sq) / g_sq)
                nh_g.append(math.sqrt(g_sq / det))
                nh_c.append(math.sqrt(c_sq / det))
            total += len(pairs)
            agree += ok
            print(f"{bits:>5} {eff:>6.2f} {len(pairs):>6} "
                  f"{ok:>5}/{len(pairs):<6} {pearson(nh_g, nh_c):>9.6f} "
                  f"{spearman(nh_g, nh_c):>9.6f} {worst:>12.2e}")
    print(f"\n  TOTAL exact agreement: {agree}/{total} "
          f"({100.0*agree/total:.2f}%)")
    print("  (H1) is CONFIRMED exactly." if agree == total else
          "  (H1) FAILS on some inputs — see max rel err column.")
    return agree == total


# ---------------------------------------------------------------------------
# C2 — where does the active convergent sit?
# ---------------------------------------------------------------------------

def exp_c2():
    print("\n" + "=" * 78)
    print("C2 — localisation of the active convergent around q* = sqrt(n*K2/K1)")
    print("=" * 78)
    print("log2(q_active / q*) should concentrate on [-1, +1] if the critical")
    print("scale story is right; the active index should NOT be the last one.")
    print(f"{'bits':>5} {'eff':>6} {'|log2 ratio| med':>17} {'<=1':>7} {'<=2':>7} "
          f"{'j_act/j_max':>12} {'b=0 wins':>9}")
    for bits in (20, 32, 64, 128, 256):
        pairs = sample_glv_pairs(bits, 60 if bits <= 32 else 40, seed=1000 + bits)
        for eff in (0.05, 0.25):
            ratios, fracs, zeros = [], [], 0
            for n, lam in pairs:
                k1, k2 = k1k2_for(n, eff)
                _, q_act, _ = nu_from_cf(n, lam, k1, k2)
                qstar = critical_scale(n, k1, k2)
                if q_act == 0:
                    zeros += 1
                    continue
                ratios.append(abs(math.log2(q_act / qstar)))
                qs, _, _ = cf_convergents(lam, n)
                fracs.append(qs.index(q_act) / max(1, len(qs) - 1))
            ratios.sort()
            med = ratios[len(ratios) // 2] if ratios else float('nan')
            le1 = sum(1 for r in ratios if r <= 1) / max(1, len(ratios))
            le2 = sum(1 for r in ratios if r <= 2) / max(1, len(ratios))
            mf = sum(fracs) / max(1, len(fracs))
            print(f"{bits:>5} {eff:>6.2f} {med:>17.3f} {le1:>7.2f} {le2:>7.2f} "
                  f"{mf:>12.3f} {zeros:>9}")
    print("\n  j_act/j_max ~ 0.5 means the active convergent is in the MIDDLE of")
    print("  the expansion — this is exactly why scale-free CF statistics fail.")

    # C2b — the offset is not noise: (H3) minimises at t = 1/sqrt(a*), i.e.
    # predicts q_act = q*/sqrt(a*), a SIGNED offset of -log2(a*)/2 octaves.
    print("\n  C2b — signed offset test: is log2(q_act/q*) = -log2(a*)/2 ?")
    print(f"  {'bits':>5} {'eff':>6} {'med signed':>11} {'med pred':>9} "
          f"{'pearson':>8} {'med |resid|':>12}")
    for bits in (20, 64, 256):
        pairs = sample_glv_pairs(bits, 60 if bits <= 32 else 50, seed=3000 + bits)
        for eff in (0.05, 0.25):
            obs, pred, resid = [], [], []
            for n, lam in pairs:
                k1, k2 = k1k2_for(n, eff)
                _, q_act, _ = nu_from_cf(n, lam, k1, k2)
                if q_act == 0:
                    continue
                a_star, _, _, qstar = local_partial_quotient(n, lam, k1, k2)
                o = math.log2(q_act / qstar)
                p = -math.log2(max(1, a_star)) / 2.0
                obs.append(o)
                pred.append(p)
                resid.append(abs(o - p))
            obs_s, pred_s, res_s = sorted(obs), sorted(pred), sorted(resid)
            print(f"  {bits:>5} {eff:>6.2f} {obs_s[len(obs_s)//2]:>11.3f} "
                  f"{pred_s[len(pred_s)//2]:>9.3f} {pearson(obs, pred):>8.4f} "
                  f"{res_s[len(res_s)//2]:>12.3f}")


# ---------------------------------------------------------------------------
# C3 / C4 — two-term form and the sqrt(2/a*) closed form
# ---------------------------------------------------------------------------

def exp_c3_c4():
    print("\n" + "=" * 78)
    print("C3/C4 — two-term form (H2) and the closed form nu_hat ~ sqrt(2/a*) (H3)")
    print("=" * 78)
    print(f"{'bits':>5} {'eff':>6} {'r(H2)':>8} {'r(H3)':>8} {'sp(H3)':>8} "
          f"{'r(max_a)':>9} {'sp(max_a)':>10} {'med |H3/nu-1|':>14}")
    pooled = {'nu': [], 'a_star': [], 'max_a': []}
    for bits in (20, 32, 64, 128, 256):
        pairs = sample_glv_pairs(bits, 80 if bits <= 32 else 50, seed=2000 + bits)
        for eff in (0.05, 0.25):
            nu, h2, h3, inv_maxa, rel = [], [], [], [], []
            for n, lam in pairs:
                k1, k2 = k1k2_for(n, eff)
                g_sq, det = nu_exact(n, lam, k1, k2)
                nh = math.sqrt(g_sq / det)
                s1, _, s2, _ = scales(n, k1, k2)
                qs, _, _ = cf_convergents(lam, n)
                # (H2): min_j sqrt(t_j^2 + 1/t_{j+1}^2)
                sc = math.sqrt(s2 / (n * s1))
                best = float('inf')
                for idx in range(len(qs) - 1):
                    t_j, t_j1 = qs[idx] * sc, qs[idx + 1] * sc
                    if t_j1 == 0:
                        continue
                    best = min(best, math.sqrt(t_j ** 2 + 1.0 / t_j1 ** 2))
                a_star, _, _, _ = local_partial_quotient(n, lam, k1, k2)
                ma = max_partial_quotient(lam, n)
                pred3 = math.sqrt(2.0 / max(1, a_star))
                nu.append(nh)
                h2.append(best)
                h3.append(pred3)
                inv_maxa.append(math.sqrt(2.0 / ma))
                rel.append(abs(pred3 / nh - 1.0))
                pooled['nu'].append(nh)
                pooled['a_star'].append(a_star)
                pooled['max_a'].append(ma)
            rel.sort()
            print(f"{bits:>5} {eff:>6.2f} {pearson(nu, h2):>8.4f} "
                  f"{pearson(nu, h3):>8.4f} {spearman(nu, h3):>8.4f} "
                  f"{pearson(nu, inv_maxa):>9.4f} {spearman(nu, inv_maxa):>10.4f} "
                  f"{rel[len(rel)//2]:>14.3f}")
    print("\n  Pooled (all bit sizes / eff):")
    print(f"    spearman(nu_hat, a*)     = {spearman(pooled['nu'], pooled['a_star']):+.4f}"
          "   <- LOCAL partial quotient")
    print(f"    spearman(nu_hat, max_a)  = {spearman(pooled['nu'], pooled['max_a']):+.4f}"
          "   <- GLOBAL max (the June invariant)")
    print("  Both should be negative; |local| >> |global| is the retro-explanation")
    print("  of the 2026-06 max_a falsification.")


# ---------------------------------------------------------------------------
# C5 — does a* predict actual LLL recovery?
# ---------------------------------------------------------------------------

C5_BITS_LO, C5_BITS_HI = 1 << 19, 1 << 20
C5_EFF = 0.05
C5_M = 9
C5_SEEDS = 24
C5_CURVES = 150


def harvest_c5(count, bits_lo=C5_BITS_LO, bits_hi=C5_BITS_HI, max_primes=200000):
    """Real j=0 CM curves with a generator, spread over a*."""
    out, seen = [], set()
    scanned = 0
    p = int(sympy.nextprime(bits_lo))
    while p < bits_hi and scanned < max_primes and len(out) < count:
        scanned += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n = p + 1 - t
                    if n < 2 or n in seen or n % 3 != 1 or not sympy.isprime(n):
                        continue
                    roots = glv_eigenvalues(n)
                    if roots is None:
                        continue
                    b = identify_twist(p, n)
                    if b is None:
                        continue
                    G = find_generator(p, b, n)
                    if G is None:
                        continue
                    seen.add(n)
                    lam = roots[0]
                    k1, k2 = k1k2_for(n, C5_EFF)
                    g_sq, det = nu_exact(n, lam, k1, k2)
                    a_star, j_star, q_j, qstar = local_partial_quotient(n, lam, k1, k2)
                    out.append({'p': p, 'b': b, 'n': n, 'lam': lam, 'G': G,
                                'k1': k1, 'k2': k2,
                                'nu_hat': math.sqrt(g_sq / det),
                                'a_star': a_star, 'max_a': max_partial_quotient(lam, n),
                                'mu': mu_of(lam, n)})
                    break
        p = int(sympy.nextprime(p))
    return out


def exp_c5(bits=20, m_sigs=C5_M):
    print("\n" + "=" * 78)
    print(f"C5 — does the ALGEBRAIC a* predict LLL recovery?  ({bits}-bit curves)")
    print("=" * 78)
    t0 = time.time()
    curves = harvest_c5(C5_CURVES, 1 << (bits - 1), 1 << bits)
    print(f"  harvested {len(curves)} {bits}-bit j=0 CM curves in {time.time()-t0:.1f}s "
          f"(eff={C5_EFF}, m={m_sigs}, {C5_SEEDS} seeds each)")
    C5_M_LOCAL = m_sigs
    t0 = time.time()
    for c in curves:
        wins = 0
        for seed in range(C5_SEEDS):
            d = random.Random(seed + 7777).randint(1, c['n'] - 1)
            wins += run_experiment(c['p'], c['n'], c['lam'], c['G'], C5_M_LOCAL, d,
                                   c['k1'], c['k2'], seed)
        c['p_hat'] = wins / C5_SEEDS
    print(f"  {len(curves)*C5_SEEDS} LLL trials in {time.time()-t0:.1f}s")

    nu = [c['nu_hat'] for c in curves]
    ph = [c['p_hat'] for c in curves]
    print(f"\n  spearman(nu_hat, p_hat) = {spearman(nu, ph):+.4f}   (the 07-29 predictor)")
    print(f"  spearman(a*,     p_hat) = {spearman([c['a_star'] for c in curves], ph):+.4f}"
          "   <- ALGEBRAIC, local")
    print(f"  spearman(max_a,  p_hat) = {spearman([c['max_a'] for c in curves], ph):+.4f}"
          "   <- June invariant")
    print(f"  spearman(mu,     p_hat) = {spearman([c['mu'] for c in curves], ph):+.4f}"
          "   <- falsified 07-29")

    med = sorted(ph)[len(ph) // 2]
    labels = [1 if c['p_hat'] > med else 0 for c in curves]
    if 0 < sum(labels) < len(labels):
        print("\n  AUC for predicting above-median recovery:")
        print(f"    a*      : {auc([c['a_star'] for c in curves], labels):.3f}")
        print(f"    -nu_hat : {auc([-c['nu_hat'] for c in curves], labels):.3f}")
        print(f"    max_a   : {auc([c['max_a'] for c in curves], labels):.3f}")
        print(f"    mu      : {auc([c['mu'] for c in curves], labels):.3f}")

    print(f"\n  {'a*':>6} {'curves':>7} {'mean nu_hat':>12} {'mean p_hat':>11}")
    buckets = {}
    for c in curves:
        key = 1 if c['a_star'] <= 1 else 2 if c['a_star'] <= 2 else 4 if c['a_star'] <= 4 else 9
        buckets.setdefault(key, []).append(c)
    for key in sorted(buckets):
        g = buckets[key]
        lbl = {1: '1', 2: '2', 4: '3-4', 9: '>=5'}[key]
        print(f"  {lbl:>6} {len(g):>7} "
              f"{sum(x['nu_hat'] for x in g)/len(g):>12.3f} "
              f"{sum(x['p_hat'] for x in g)/len(g):>11.3f}")
    return curves


# ---------------------------------------------------------------------------
# C6 — secp256k1
# ---------------------------------------------------------------------------

def exp_c6():
    print("\n" + "=" * 78)
    print("C6 — secp256k1: active convergent, a*, and consistency with 2026-07-29")
    print("=" * 78)
    assert (SECP_LAM ** 2 + SECP_LAM + 1) % SECP_N == 0, "lam is not a GLV eigenvalue"
    print("  lam^2 + lam + 1 = 0 mod n  ... verified")
    print(f"\n  {'eff':>6} {'nu_hat':>9} {'07-29':>8} {'a*':>6} {'j*':>4} "
          f"{'log2 q_act':>11} {'log2 q*':>9} {'sqrt(2/a*)':>11}")
    prior = {0.05: 0.8709, 0.10: 0.6624, 0.25: 0.5852, 0.50: 0.6851}
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1, k2 = k1k2_for(SECP_N, eff)
        g_sq, det = nu_exact(SECP_N, SECP_LAM, k1, k2)
        nh = math.sqrt(g_sq / det)
        c_sq, q_act, _ = nu_from_cf(SECP_N, SECP_LAM, k1, k2)
        assert c_sq == g_sq, "CF form disagrees at 256 bits"
        a_star, j_star, _, qstar = local_partial_quotient(SECP_N, SECP_LAM, k1, k2)
        print(f"  {eff:>6.2f} {nh:>9.4f} {prior.get(eff, float('nan')):>8.4f} "
              f"{a_star:>6} {j_star:>4} {math.log2(q_act):>11.2f} "
              f"{math.log2(qstar):>9.2f} {math.sqrt(2.0/max(1,a_star)):>11.4f}")
    qs, a_list, _ = cf_convergents(SECP_LAM, SECP_N)
    print(f"\n  CF of lam/n: {len(a_list)-1} partial quotients, "
          f"max_a = {max(a_list[1:])}, "
          f"mean = {sum(a_list[1:])/len(a_list[1:]):.2f}")
    print(f"  first 20 a_j: {a_list[1:21]}")
    print("\n  Interpretation: max_a is a global outlier statistic; the geometry")
    print("  only reads a* at j*, which is mid-expansion and modest.")


def main():
    print("=" * 78)
    print("GLV-HNP Thread 23 — closed form for nu_hat via continued fractions")
    print("=" * 78)
    ok = exp_c1()
    exp_c2()
    exp_c3_c4()
    exp_c5(bits=20, m_sigs=9)
    exp_c5(bits=26, m_sigs=11)
    exp_c6()
    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print(f"  (H1) lambda_1(L2) at a CF convergent : {'CONFIRMED (exact)' if ok else 'FAILED'}")
    print("  (H2)/(H3): see C3/C4 correlations above.")


if __name__ == "__main__":
    main()

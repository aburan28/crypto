"""
GLV-HNP — Thread 23: closed form for nu_hat via continued fractions.

Background (2026-07-29, Thread 20b/20c)
---------------------------------------
nu_hat = lambda_1(L2) / sqrt(det L2) was found to be the first invariant that
separates the Phase-2 easy/hard classes (AUC 0.935 on the Exp-S replication),
after eight curve-level invariants had been falsified:

    delta/n, kappa(M), q_cf, max_q_cf, max_a, a_corn/n, lam/n, mu

L2 is the NON-planted 2-dimensional sublattice on coordinate pair (i, m+1+i):

    L2 = < u = (n*S1, 0),  v = (-lam*S1, S2) >,      det L2 = n*S1*S2

with S1 = S_K1 = n//K1 and S2 = S_K2 = max(1, n//K2).

nu_hat was *computed* (Lagrange-Gauss) but not *understood*.  The 2026-07-29
next-step proposal (Thread 23) conjectured that lambda_1(L2) is a classical
continued-fraction quantity, and that this would retroactively explain why the
scale-free CF invariants (q_cf, max_q_cf, max_a) all failed in June.

The claim proved and tested here
--------------------------------
A general element of L2 is

    a*u + b*v = ( (a*n - b*lam)*S1 , b*S2 ).

For fixed b the inner factor is minimised over a at the centred residue
||b*lam||_n := min_a |a*n - b*lam|.  Hence

    lambda_1(L2)^2 = min_{b >= 0}  ( ||b*lam||_n * S1 )^2 + ( b*S2 )^2 .   (*)

Because b*S2 is strictly increasing in b, any minimiser of (*) must be a strict
record of ||b*lam||_n: no smaller b' may have ||b'*lam||_n <= ||b*lam||_n.  The
records of ||b*lam||_n are exactly the *best approximations of the second kind*
to lam/n, which by the classical theorem are exactly the continued-fraction
convergent denominators q_j of lam/n.  Therefore

    lambda_1(L2) = min over CF convergents (h_j, q_j) of lam/n of
                   sqrt( (|q_j*lam - h_j*n| * S1)^2 + (q_j*S2)^2 )           (**)

which is EXACT, not approximate, and costs one Euclid run: O(log n).

T1 tests (**) against the Lagrange-Gauss ground truth.  The stated falsifier
was "predicted-from-CF nu_hat correlates < 0.9 with computed nu_hat"; if the
argument above is right the correlation is not 0.9 but exactly 1, with zero
mismatches, so T1 is run as an exact-equality test.

T2 locates the minimising convergent relative to the balance scale
b* = sqrt(n*S1/S2) = sqrt(n*K2/K1).

T3/T4 are the payoff: nu_hat is small (attack easy) exactly when lam/n has a
large partial quotient AT THE BALANCE SCALE, because ||q_j*lam||_n ~ n/q_{j+1}
= n/(a_{j+1}*q_j + q_{j-1}).  A scale-free statistic such as max_a picks the
largest partial quotient ANYWHERE in the expansion, which is almost always at
a different scale -- hence its June falsification.  T4 measures both.

T5 records the Hermite bound nu_hat <= (2/sqrt(3))^(1/2) = 1.074570 and the
secp256k1 placement, both recomputed through the closed form.

Run: python3 glv_hnp_nuhat_cf_closed_form.py
"""

import math
import random
import sys
import time

import sympy

import contextlib
import io
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_t20a", __file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_lambda_threshold.py")
_t20a = importlib.util.module_from_spec(_spec)
# NOTE: glv_hnp_phase2_lambda_threshold.py has NO `if __name__ == "__main__"`
# guard -- its whole T1..T5 experiment suite runs at import time (~4 s of
# output).  Suppress it so this script's output is readable.  The same defect
# silently affects glv_hnp_phase2_nuhat_control.py, which imports it the same
# way.  Left unfixed here to keep that file byte-identical to the 2026-07-29
# run it documents; see the log entry for 2026-08-02.
with contextlib.redirect_stdout(io.StringIO()):
    _spec.loader.exec_module(_t20a)

gauss_reduce_2d = _t20a.gauss_reduce_2d
lambda_block_mu = _t20a.lambda_block_mu
glv_roots = _t20a.glv_roots
lam_star = _t20a.lam_star
scales = _t20a.scales

HERMITE2 = math.sqrt(2.0 / math.sqrt(3.0))  # 1.0745699...


# ---------------------------------------------------------------------------
# Continued fractions
# ---------------------------------------------------------------------------

def cf_convergents(lam, n):
    """Convergents h_j/q_j of lam/n via Euclid.

    Returns list of (a_j, h_j, q_j, err_j) with err_j = |q_j*lam - h_j*n|,
    a_j the j-th partial quotient.  err_j = ||q_j*lam||_n for convergents.
    """
    lam %= n
    out = []
    hm1, hm2 = 1, 0   # h_{j-1}, h_{j-2}
    km1, km2 = 0, 1   # q_{j-1}, q_{j-2}
    a, b = lam, n
    while b:
        q = a // b
        a, b = b, a - q * b
        h = q * hm1 + hm2
        k = q * km1 + km2
        out.append((q, h, k, abs(k * lam - h * n)))
        hm1, hm2 = h, hm1
        km1, km2 = k, km1
        if k > n:
            break
    return out


def nu_hat_cf(n, lam, S1, S2):
    """lambda_1(L2)/sqrt(det L2) via the closed form (**).  Returns
    (nu_hat, lambda_1, b_star, j_min, q_min, a_next)."""
    lam %= n
    best2 = (n * S1) ** 2          # the b = 0 vector (n*S1, 0)
    j_min, q_min, a_next = -1, 0, 0
    conv = cf_convergents(lam, n)
    for j, (a_j, h_j, q_j, err_j) in enumerate(conv):
        c2 = (err_j * S1) ** 2 + (q_j * S2) ** 2
        if c2 < best2:
            best2 = c2
            j_min, q_min = j, q_j
            a_next = conv[j + 1][0] if j + 1 < len(conv) else 0
    lam1 = math.sqrt(best2)
    det = n * S1 * S2
    b_star = math.sqrt(n * S1 / S2)
    return lam1 / math.sqrt(det), lam1, b_star, j_min, q_min, a_next


def nu_hat_gauss(n, lam, S1, S2):
    """Ground truth via Lagrange-Gauss reduction."""
    lam1, _w = lambda_block_mu(n, lam, S1, S2)
    return lam1 / math.sqrt(n * S1 * S2), lam1


def max_partial_quotient(lam, n):
    """Scale-free June invariant max_a: largest partial quotient of lam/n."""
    return max((a for a, _h, _k, _e in cf_convergents(lam, n)), default=0)


def local_partial_quotient(lam, n, S1, S2):
    """The partial quotient a_{j+1} immediately after the convergent nearest
    (from below) the balance scale b* = sqrt(n*S1/S2).  This is the
    scale-LOCAL analogue of max_a."""
    b_star = math.sqrt(n * S1 / S2)
    conv = cf_convergents(lam, n)
    j_lo = -1
    for j, (_a, _h, q, _e) in enumerate(conv):
        if q <= b_star:
            j_lo = j
        else:
            break
    if j_lo + 1 < len(conv):
        return conv[j_lo + 1][0]
    return conv[-1][0] if conv else 0


# ---------------------------------------------------------------------------
# statistics
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
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx > 0 and dy > 0 else 0.0


def pearson(xs, ys):
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    dy = math.sqrt(sum((b - my) ** 2 for b in ys))
    return num / (dx * dy) if dx > 0 and dy > 0 else 0.0


def glv_pairs(bits, count, seed=20260802):
    """(n, lam) with n prime, n = 1 mod 3, lam^2+lam+1 = 0 mod n."""
    rng = random.Random(seed)
    out = []
    lo, hi = 1 << (bits - 1), 1 << bits
    guard = 0
    while len(out) < count and guard < count * 4000:
        guard += 1
        n = int(sympy.nextprime(rng.randrange(lo, hi)))
        if n >= hi or n % 3 != 1:
            continue
        roots = glv_roots(n)
        if roots is None:
            continue
        out.append((n, roots[0]))
    return out


# ---------------------------------------------------------------------------
# T1 -- exactness of the closed form
# ---------------------------------------------------------------------------

def T1_exactness():
    print("=" * 74)
    print("T1  closed form (**) vs Lagrange-Gauss ground truth")
    print("=" * 74)
    rng = random.Random(11)
    total = 0
    mismatch = 0
    worst = 0.0
    cf_nu, gs_nu = [], []
    rows = []
    for bits in (16, 20, 24, 32, 48, 64, 128, 256):
        bad = 0
        cnt = 0
        for _ in range(300):
            n = int(sympy.nextprime(rng.randrange(1 << (bits - 1), 1 << bits)))
            lam = rng.randrange(1, n)
            K1 = rng.choice([2, 4, 8, 16, 72, 144, 1 << (bits // 3)])
            K2 = rng.choice([2, 16, 52, 1 << (bits // 3)])
            S1, _SD, S2, _SK = scales(n, max(1, K1), max(1, K2))
            nu_c, l_c, _bs, _j, _q, _a = nu_hat_cf(n, lam, S1, S2)
            nu_g, l_g = nu_hat_gauss(n, lam, S1, S2)
            total += 1
            cnt += 1
            cf_nu.append(nu_c)
            gs_nu.append(nu_g)
            if l_c != l_g:
                rel = abs(l_c - l_g) / max(l_g, 1e-300)
                worst = max(worst, rel)
                if rel > 1e-12:
                    bad += 1
                    mismatch += 1
        rows.append((bits, cnt, bad))
    print(f"{'bits':>6} {'cases':>7} {'mismatches':>11}")
    for bits, cnt, bad in rows:
        print(f"{bits:>6} {cnt:>7} {bad:>11}")
    print(f"\ntotal cases    : {total}")
    print(f"mismatches     : {mismatch}")
    print(f"worst rel.err  : {worst:.3e}")
    print(f"pearson(cf,gs) : {pearson(cf_nu, gs_nu):.6f}")
    print(f"spearman(cf,gs): {spearman(cf_nu, gs_nu):.6f}")
    print(f"\nfalsifier was 'correlation < 0.9'  ->  "
          f"{'PASSED (exact)' if mismatch == 0 else 'FAILED'}")
    return mismatch == 0


# ---------------------------------------------------------------------------
# T2 -- where the minimiser sits relative to the balance scale
# ---------------------------------------------------------------------------

def T2_balance_scale():
    print()
    print("=" * 74)
    print("T2  minimising convergent vs balance scale b* = sqrt(n*S1/S2)")
    print("=" * 74)
    rng = random.Random(23)
    offs = {}
    ratios = []
    for _ in range(2000):
        bits = rng.choice([20, 24, 32, 64, 128, 256])
        n = int(sympy.nextprime(rng.randrange(1 << (bits - 1), 1 << bits)))
        lam = rng.randrange(1, n)
        K1 = rng.choice([4, 16, 72, 144])
        K2 = rng.choice([16, 52, 256])
        S1, _SD, S2, _SK = scales(n, K1, K2)
        _nu, _l1, b_star, j_min, q_min, _a = nu_hat_cf(n, lam, S1, S2)
        conv = cf_convergents(lam, n)
        j_lo = -1
        for j, (_a2, _h, q, _e) in enumerate(conv):
            if q <= b_star:
                j_lo = j
            else:
                break
        off = j_min - j_lo
        offs[off] = offs.get(off, 0) + 1
        if q_min > 0 and b_star > 0:
            ratios.append(math.log(q_min / b_star, 2))
    print("offset of argmin convergent index from the last convergent <= b*:")
    print(f"{'offset':>8} {'count':>8} {'share':>8}")
    tot = sum(offs.values())
    for off in sorted(offs):
        print(f"{off:>8} {offs[off]:>8} {offs[off] / tot:>7.1%}")
    ratios.sort()
    def qtl(v, x):
        return v[min(len(v) - 1, int(x * len(v)))]
    print(f"\nlog2(q_min / b*)  q05={qtl(ratios,0.05):+.2f} "
          f"q25={qtl(ratios,0.25):+.2f} q50={qtl(ratios,0.50):+.2f} "
          f"q75={qtl(ratios,0.75):+.2f} q95={qtl(ratios,0.95):+.2f}")
    within1 = sum(1 for r in ratios if abs(r) <= 1.0) / len(ratios)
    print(f"share with |log2(q_min/b*)| <= 1 (i.e. within a factor 2): {within1:.1%}")


# ---------------------------------------------------------------------------
# T3/T4 -- scale-local vs scale-free CF statistics on real GLV (n, lam)
# ---------------------------------------------------------------------------

def T3_T4_scale_locality():
    print()
    print("=" * 74)
    print("T3/T4  scale-free max_a vs scale-local a_{j+1} as predictors of nu_hat")
    print("=" * 74)
    print("(real j=0 GLV pairs: n prime, n=1 mod 3, lam^2+lam+1=0 mod n)")
    for bits, count in ((20, 400), (24, 300), (32, 200)):
        pairs = glv_pairs(bits, count, seed=20260802 + bits)
        for K1, K2 in ((72, 52), (144, 52), (36, 52)):
            nus, maxa, loca, mus = [], [], [], []
            for n, lam in pairs:
                S1, _SD, S2, _SK = scales(n, K1, K2)
                nu, _l1, _bs, _j, _q, _an = nu_hat_cf(n, lam, S1, S2)
                nus.append(nu)
                maxa.append(max_partial_quotient(lam, n))
                loca.append(local_partial_quotient(lam, n, S1, S2))
                mus.append(lam_star(lam, n))
            print(f"\nbits={bits}  K1={K1} K2={K2}  N={len(pairs)}")
            print(f"  spearman(max_a,   nu_hat) = {spearman(maxa, nus):+.3f}   "
                  f"[scale-free  -- the June invariant]")
            print(f"  spearman(a_local, nu_hat) = {spearman(loca, nus):+.3f}   "
                  f"[scale-local -- this thread]")
            print(f"  spearman(mu,      nu_hat) = {spearman(mus, nus):+.3f}   "
                  f"[falsified 2026-07-29]")


# ---------------------------------------------------------------------------
# T6 -- the two-parameter asymptotic formula
# ---------------------------------------------------------------------------

def nu_hat_formula(n, lam, S1, S2):
    """Asymptotic nu_hat from the two CF parameters at the balance scale.

    Put t = q_j / b* with b* = sqrt(n*S1/S2), and let a = a_{j+1}.  Using
    ||q_j*lam||_n ~ n/q_{j+1} ~ n/(a*q_j):

        (q_j*S2)^2 / (n*S1*S2)          = t^2
        (||q_j lam||_n*S1)^2/(n*S1*S2)  ~ 1/(a^2 * t^2)

    so   nu_hat^2  ~  t^2 + 1/(a^2 t^2).                                 (***)

    This is the point of T6: nu_hat is governed by TWO parameters, not one.
    a alone predicts nu_hat only when t happens to be near-constant across the
    sample -- which is exactly the K1-dependence seen in T3/T4.
    """
    b_star = math.sqrt(n * S1 / S2)
    conv = cf_convergents(lam, n)
    j_lo = -1
    for j, (_a, _h, q, _e) in enumerate(conv):
        if q <= b_star:
            j_lo = j
        else:
            break
    if j_lo < 0 or j_lo + 1 >= len(conv):
        return None, None, None
    q_j = conv[j_lo][2]
    a_next = conv[j_lo + 1][0]
    t = q_j / b_star
    if t <= 0 or a_next <= 0:
        return None, None, None
    return math.sqrt(t * t + 1.0 / (a_next * a_next * t * t)), t, a_next


def T6_two_parameter():
    print()
    print("=" * 74)
    print("T6  two-parameter formula  nu_hat^2 ~ t^2 + 1/(a^2 t^2),  t = q_j/b*")
    print("=" * 74)
    print("Columns: t2dom = share where the t^2 term carries the norm;")
    print("IQR(nu) = interquartile spread of nu_hat itself over the sample;")
    print("argmin@j_lo = share where the true SVP argmin is the convergent")
    print("immediately below the balance scale b* = sqrt(n*r).")
    print("\nNOTE: the first-guess explanation 'a_local works when the spread of")
    print("t is small' is NOT supported -- IQR(t) is ~0.5 at both r=0.722 and")
    print("r=0.361 while spearman(a,nu) goes +0.01 -> -0.94.  What does track")
    print("it is IQR(nu): the model is informative exactly where nu_hat is")
    print("dispersed (r <~ 0.36, IQR ~ 0.3) and uninformative where nu_hat is")
    print("compressed (r ~ 0.7-1.4, IQR ~ 0.11).  r ~ 0.72 remains anomalous:")
    print("t2dom there (27%) is close to r=0.361's (18%), yet spearman(a,nu)")
    print("collapses.  Recorded as open; see the 2026-08-02 log next step.\n")
    print(f"{'bits':>5} {'K1':>5} {'K2':>4} {'r=K2/K1':>8} {'IQR(t)':>7} "
          f"{'med t':>6} {'t2dom':>6} {'IQR(nu)':>8} "
          f"{'sp(a,nu)':>9} {'sp(fml,nu)':>11} {'argmin@j_lo':>12}")
    for bits, count in ((20, 400), (24, 300), (32, 200)):
        pairs = glv_pairs(bits, count, seed=20260802 + bits)
        for K1, K2 in ((36, 52), (72, 52), (144, 52), (288, 52)):
            nus, fml, ts, aa = [], [], [], []
            t2dom = 0
            at_jlo = 0
            for n, lam in pairs:
                S1, _SD, S2, _SK = scales(n, K1, K2)
                nu, _l1, b_star, j_min, q_min, _an = nu_hat_cf(n, lam, S1, S2)
                f, t, a = nu_hat_formula(n, lam, S1, S2)
                if f is None:
                    continue
                nus.append(nu); fml.append(f); ts.append(t); aa.append(a)
                # which of the two terms of (***) carries the norm?
                if t * t > 1.0 / (a * a * t * t):
                    t2dom += 1
                conv = cf_convergents(lam, n)
                j_lo = -1
                for j, (_a2, _h, q, _e) in enumerate(conv):
                    if q <= b_star:
                        j_lo = j
                    else:
                        break
                if j_min == j_lo:
                    at_jlo += 1
            if len(nus) < 30:
                continue
            st, sn = sorted(ts), sorted(nus)
            iqr = st[int(.75 * len(st))] - st[int(.25 * len(st))]
            iqrn = sn[int(.75 * len(sn))] - sn[int(.25 * len(sn))]
            print(f"{bits:>5} {K1:>5} {K2:>4} {K2/K1:>8.3f} {iqr:>7.3f} "
                  f"{st[len(st) // 2]:>6.3f} {t2dom / len(nus):>5.0%} "
                  f"{iqrn:>8.3f} "
                  f"{spearman(aa, nus):>+9.3f} {spearman(fml, nus):>+11.3f} "
                  f"{at_jlo / len(nus):>11.0%}")


# ---------------------------------------------------------------------------
# T7 -- nu_hat is scale-invariant: it depends on K2/K1 only, not on eff
# ---------------------------------------------------------------------------

def T7_scale_invariance():
    print()
    print("=" * 74)
    print("T7  nu_hat depends on (n, lam, r=K2/K1) ONLY -- not on eff=K1*K2/n")
    print("=" * 74)
    print("Proof: scaling (S1,S2) -> (c*S1,c*S2) scales lambda_1 by c and")
    print("sqrt(det)=sqrt(n*S1*S2) by c, so nu_hat is unchanged.  Since")
    print("S1=n//K1 and S2=n//K2, nu_hat is a function of K2/K1 alone.\n")
    n = int(sympy.nextprime((1 << 61) + 12345))
    lam = 1234567890123 % n
    print(f"{'K1':>8} {'K2':>8} {'eff=K1K2/n':>12} {'r=K2/K1':>9} {'nu_hat':>10}")
    for K1, K2 in ((72, 52), (144, 104), (720, 520), (7200, 5200)):
        S1, _SD, S2, _SK = scales(n, K1, K2)
        nu, _l1, _b, _j, _q, _a = nu_hat_cf(n, lam, S1, S2)
        print(f"{K1:>8} {K2:>8} {K1 * K2 / n:>12.3e} {K2 / K1:>9.4f} {nu:>10.6f}")
    print("\n=> identical nu_hat across four orders of magnitude of eff.")
    print("CONSEQUENCE: the 2026-07-29 secp256k1 table indexed by 'eff' was")
    print("mis-parameterised; the correct abscissa is r = K2/K1.  Restated below.")


# ---------------------------------------------------------------------------
# T5 -- Hermite bound and secp256k1 placement through the closed form
# ---------------------------------------------------------------------------

def T5_bound_and_secp():
    print()
    print("=" * 74)
    print("T5  Hermite bound and secp256k1 placement via the closed form")
    print("=" * 74)
    n_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
    roots = glv_roots(n_secp)
    assert roots is not None, "secp256k1 n must admit lam^2+lam+1=0"
    lam = roots[0]
    assert (lam * lam + lam + 1) % n_secp == 0
    print(f"lam (verified lam^2+lam+1 = 0 mod n): {hex(lam)}")
    print(f"mu = min(lam,n-lam)/n              : {lam_star(lam, n_secp):.6f}")

    print(f"\nHermite bound in dim 2: nu_hat <= (2/sqrt(3))^(1/2) = {HERMITE2:.6f}")
    print("(2026-07-29 measured max over the high group was 1.072 -- consistent)")

    print(f"\nsecp256k1 nu_hat as a function of the bias-shape ratio r = K2/K1")
    print("(percentile = share of 2000 random lam with smaller nu_hat, same r)\n")
    print(f"{'r=K2/K1':>9} {'nu_hat(cf)':>11} {'nu_hat(gauss)':>14} "
          f"{'q_min bits':>11} {'a_next':>8} {'pctile':>8}")
    rng = random.Random(2026)
    K2 = 1 << 40
    for r in (1 / 64, 1 / 8, 1.0, 8.0, 64.0):
        K1 = max(2, int(K2 / r))
        S1, _SD, S2, _SK = scales(n_secp, K1, K2)
        nu_c, _l1, _bs, _j, q_min, a_next = nu_hat_cf(n_secp, lam, S1, S2)
        nu_g, _lg = nu_hat_gauss(n_secp, lam, S1, S2)
        draws = sorted(nu_hat_cf(n_secp, rng.randrange(1, n_secp), S1, S2)[0]
                       for _ in range(2000))
        pct = sum(1 for d in draws if d < nu_c) / len(draws)
        qb = q_min.bit_length() if q_min else 0
        print(f"{r:>9.4f} {nu_c:>11.4f} {nu_g:>14.4f} {qb:>11} "
              f"{a_next:>8} {pct:>7.1%}")

    K1 = K2 = 1 << 40
    S1, _SD, S2, _SK = scales(n_secp, K1, K2)
    draws = sorted(nu_hat_cf(n_secp, rng.randrange(1, n_secp), S1, S2)[0]
                   for _ in range(4000))
    def qtl(v, x):
        return v[min(len(v) - 1, int(x * len(v)))]
    print(f"\nnull distribution of nu_hat over random lam (256-bit, r=1, 4000 draws):")
    print("  q05=%.3f q10=%.3f q25=%.3f q50=%.3f q75=%.3f q90=%.3f q95=%.3f"
          % tuple(qtl(draws, x) for x in (.05, .10, .25, .50, .75, .90, .95)))
    print("  share below 0.45: %.1f%%   share above 0.90: %.1f%%"
          % (100 * sum(1 for d in draws if d < 0.45) / len(draws),
             100 * sum(1 for d in draws if d > 0.90) / len(draws)))


def main():
    t0 = time.time()
    ok = T1_exactness()
    T2_balance_scale()
    T3_T4_scale_locality()
    T6_two_parameter()
    T7_scale_invariance()
    T5_bound_and_secp()
    print()
    print("=" * 74)
    print(f"done in {time.time() - t0:.1f}s   T1 exactness: {'PASS' if ok else 'FAIL'}")
    print("=" * 74)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

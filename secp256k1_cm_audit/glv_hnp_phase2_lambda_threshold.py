"""
GLV-HNP Phase 2 — Thread 20: is there a lambda/n threshold for LLL viability?

BACKGROUND
  The 2026-07-26 log entry concluded from 4 data points that the Phase-2 GLV
  lattice attack has a "lambda/n threshold":

      8-bit  n=199,    K1=2,  lam/n=0.533 -> LLL 3/3 at m=4   SUCCESS
      12-bit n=2659,   K1=8,  lam/n=0.660 -> LLL 3/3 at m=7   SUCCESS
      20-bit n=523969, K1=36, lam/n=0.340 -> LLL 3/3 at m=9   SUCCESS
      12-bit n=2647,   K1=8,  lam/n=0.070 -> never 3/3        FAILURE (BKZ-40 too)

  and proposed "LLL succeeds iff lam/n > threshold" with the threshold somewhere
  in (0.07, 0.34).

  This is suspicious for two reasons:
   (a) The 2026-06-27/06-28 sessions already FALSIFIED lam/n as a predictor of
       the Phase-2 K1_threshold (see log, Exp M/N/O/P/R).
   (b) lam and n-1-lam are BOTH roots of x^2+x+1 mod n, and the lattice entry
       -lam*S_K1 in row m+1+i can be reduced mod the modular row n*S_K1.  So the
       lattice only ever sees lam_eff = min(lam, n-lam).  A "threshold at 0.34"
       stated on lam/n in (0,1) is therefore not even well defined: lam/n=0.66
       and lam/n=0.34 give the SAME lam_eff/n=0.34, and the two "independent"
       successes at 0.66 and 0.34 may be the same measurement.

DESIGN
  Controlled cohort study.  Within a cohort we hold the bit-length and the
  bias efficiency eff = K1*K2/n fixed, and vary only lam.  We then score six
  candidate predictors against measured LLL recovery:

    P1  lam/n                      (the 2026-07-26 claim, as stated)
    P2  lam_eff/n = min(lam,n-lam)/n
    P3  mu/||v||   where mu = lambda_1 of the scaled 2-D sublattice
                   L_lam = {((c*lam mod +-n)*S_K1, c*S_K2)}, Lagrange-reduced,
                   and ||v|| the planted-vector norm
    P4  r_min/K1   where r_min = min_{0<c<K2} |c*lam mod +-n|   (ambiguity margin)
    P5  amb        = min_{0<c<K2} max(|r|/K1, |c|/K2)  (<1 => genuine ambiguity)
    P6  eff        = K1*K2/n   (information-theoretic control; held fixed)

  A predictor "survives" only if it separates the measured SUCCESS and FAILURE
  groups with no overlap.

  Exp 1  12-bit cohort, eff held at ~0.157 (the 2026-07-26 value), lam swept.
  Exp 2  20-bit cohort, eff held at ~0.157, including curves with lam_eff/n<0.10.
         Under P1/P2 these must FAIL; under P3/P4 they may succeed.  Decisive.
  Exp 3  Reproduce the 4 historical data points and score every predictor on them.
  Exp 4  Fixed curve (lam/n constant), sweep K1 -> moves P3/P4/P6 while P1/P2
         are pinned.  If recovery tracks K1, lam/n cannot be the controlling
         variable.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sys

import sympy
from fpylll import IntegerMatrix, LLL

SEEDS = [42, 1234, 9999]

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0: y^2 = x^3 + b)
# ---------------------------------------------------------------------------

def modinv(a, m):
    return pow(a, -1, m)


def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        s = 3 * x1 * x1 * modinv(2 * y1, p) % p
    else:
        s = (y2 - y1) * modinv(x2 - x1, p) % p
    x3 = (s * s - x1 - x2) % p
    y3 = (s * (x1 - x3) - y1) % p
    return (x3, y3)


def ec_mul(P, k, p):
    if k == 0: return None
    R, Q = None, P
    while k > 0:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R


def tonelli_shanks(a, p):
    a %= p
    if a == 0: return 0
    if pow(a, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p


def find_generator(p, b, n, rng_seed=12345):
    """Return P on y^2=x^3+b/F_p with n*P = O, P != O.

    n prime and n in the Hasse interval => ord(P)=n => n | #E and #E < 2n
    => #E = n.  So this both finds a generator and certifies the order.
    """
    rng = random.Random(rng_seed)
    for _ in range(4000):
        x = rng.randint(0, p - 1)
        y = tonelli_shanks((pow(x, 3, p) + b) % p, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            return P
    return None


# ---------------------------------------------------------------------------
# j=0 CM curve generation (Eisenstein)
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """(a,b) >= 0 with a^2 - a*b + b^2 = p, or None."""
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in (a + s, a - s):
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None


def j0_traces(a, b):
    """The 6 Frobenius traces of the 6 sextic twists of j=0 over F_p."""
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]


def glv_eigenvalues(n):
    """Both roots of x^2+x+1 = 0 mod n, or (None,None)."""
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    for r in (r1, r2):
        if (r * r + r + 1) % n != 0:
            return None, None
    return min(r1, r2), max(r1, r2)


def gen_curves(p_lo, p_hi, want, lam_ratio_filter=None, max_primes=400000):
    """Yield (p, b, n, lam, G) for prime-order j=0 curves with p in [p_lo,p_hi).

    lam is always the SMALLER root, so lam/n <= 1/2 by construction.
    lam_ratio_filter: optional predicate on lam/n.
    """
    out = []
    p = int(sympy.nextprime(p_lo - 1))
    seen = 0
    while p < p_hi and len(out) < want and seen < max_primes:
        seen += 1
        if p % 3 != 1:
            p = int(sympy.nextprime(p)); continue
        eis = eisenstein_decompose(p)
        if eis is None:
            p = int(sympy.nextprime(p)); continue
        a_e, b_e = eis
        for t in j0_traces(a_e, b_e):
            n = p + 1 - t
            if n < 8 or not sympy.isprime(n) or n % 3 != 1:
                continue
            lam, _ = glv_eigenvalues(n)
            if lam is None:
                continue
            if lam_ratio_filter is not None and not lam_ratio_filter(lam / n):
                continue
            b_found, G = None, None
            for b_try in range(1, 60):
                G = find_generator(p, b_try, n)
                if G is not None:
                    b_found = b_try
                    break
            if b_found is None:
                continue
            out.append((p, b_found, n, lam, G))
            break
        p = int(sympy.nextprime(p))
    return out


# ---------------------------------------------------------------------------
# Phase-2 GLV lattice (identical construction to glv_hnp_phase2_20bit.py)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)


def gen_signatures(G, d, m, n, lam, p, k1_bound, k2_bound, seed):
    rng = random.Random(seed)
    sigs, attempts = [], 0
    while len(sigs) < m and attempts < 200000:
        attempts += 1
        k1 = rng.randint(0, k1_bound - 1)
        k2 = rng.randint(0, k2_bound - 1)
        k_full = (k1 + lam * k2) % n
        if k_full == 0: continue
        R = ec_mul(G, k_full, p)
        if R is None: continue
        r = R[0] % n
        if r == 0: continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d * r) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        A, B = h * s_inv % n, r * s_inv % n
        assert (A + B * d) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs


def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KAN
    return M


def planted_norm(sigs, d, n, k1_bound, k2_bound):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    tot = sum((s['k1'] * S_K1) ** 2 + (s['k2'] * S_K2) ** 2 for s in sigs)
    tot += (d * S_D) ** 2 + S_KAN ** 2
    return math.sqrt(tot)


def recover_d(rows, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN:
            continue
        sign = 1 if last > 0 else -1
        if (sign * row[m]) % n == d_secret:
            return True
    return False


def run_once(curve, m, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    dim = 2 * m + 2
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    ok = recover_d(rows, m, n, scales(n, k1_bound, k2_bound)[3], d)
    return ok, planted_norm(sigs, d, n, k1_bound, k2_bound)


def first_full_recovery_m(curve, k1_bound, m_range, seeds=SEEDS):
    """Smallest m in m_range with len(seeds)/len(seeds) recovery, else None.

    Also returns the per-m tallies and the mean planted norm at the largest m.
    """
    tally, last_norm = {}, None
    hit = None
    for m in m_range:
        wins = 0
        for sd in seeds:
            res = run_once(curve, m, k1_bound, sd)
            if res is None:
                continue
            ok, pn = res
            wins += ok
            last_norm = pn
        tally[m] = wins
        if hit is None and wins == len(seeds):
            hit = m
    return hit, tally, last_norm


# ---------------------------------------------------------------------------
# Predictors
# ---------------------------------------------------------------------------

def lagrange_shortest(b1, b2):
    """Exact lambda_1 of the 2-D lattice spanned by integer vectors b1,b2."""
    u, v = list(b1), list(b2)
    def nrm2(w): return w[0] * w[0] + w[1] * w[1]
    if nrm2(u) < nrm2(v):
        u, v = v, u
    while True:
        q = round((u[0] * v[0] + u[1] * v[1]) / nrm2(v))
        r = [u[0] - q * v[0], u[1] - q * v[1]]
        if nrm2(r) >= nrm2(v):
            return math.sqrt(nrm2(v))
        u, v = v, r


def sublattice_mu(n, lam, k1_bound, k2_bound):
    """lambda_1 of the scaled k2-sublattice L_lam (= the direct summand of the
    zero-Kannan sublattice attached to one signature index)."""
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    return lagrange_shortest((n * S_K1, 0), (-lam * S_K1, S_K2))


def ambiguity(n, lam, k1_bound, k2_bound):
    """(r_min, amb): r_min = min_{0<c<K2} |c*lam mod +-n|;
    amb = min_{0<c<K2} max(|r|/K1, c/K2).  amb<1 => a second (k1,k2) pair with
    both coordinates in range explains the same k, i.e. genuine ambiguity."""
    r_min, amb = None, None
    for c in range(1, k2_bound):
        r = (c * lam) % n
        if r > n // 2:
            r -= n
        ar = abs(r)
        if r_min is None or ar < r_min:
            r_min = ar
        cand = max(ar / k1_bound, c / k2_bound)
        if amb is None or cand < amb:
            amb = cand
    return r_min, amb


def predictors(n, lam, k1_bound, m, planted):
    k2_bound = math.isqrt(n) + 1
    mu = sublattice_mu(n, lam, k1_bound, k2_bound)
    r_min, amb = ambiguity(n, lam, k1_bound, k2_bound)
    lam_eff = min(lam, n - lam)
    return {
        'lam_n':      lam / n,
        'lam_eff_n':  lam_eff / n,
        'mu_over_v':  mu / planted if planted else float('nan'),
        'rmin_over_K1': r_min / k1_bound,
        'amb':        amb,
        'eff':        k1_bound * k2_bound / n,
    }


PRED_KEYS = ['lam_n', 'lam_eff_n', 'mu_over_v', 'rmin_over_K1', 'amb', 'eff']
PRED_LABEL = {
    'lam_n':        'P1 lam/n',
    'lam_eff_n':    'P2 lam_eff/n',
    'mu_over_v':    'P3 mu/|v|',
    'rmin_over_K1': 'P4 r_min/K1',
    'amb':          'P5 amb',
    'eff':          'P6 eff',
}


def score_predictors(rows):
    """rows: list of dicts with 'success' (bool) and the predictor values.
    A predictor separates iff max over one group < min over the other."""
    succ = [r for r in rows if r['success']]
    fail = [r for r in rows if not r['success']]
    print(f"\n  predictor separation  ({len(succ)} success / {len(fail)} failure)")
    if not succ or not fail:
        print("    (one group empty - no separation test possible)")
        return
    print(f"    {'predictor':<16} {'success range':<22} {'failure range':<22} sep")
    for k in PRED_KEYS:
        sv = [r[k] for r in succ]
        fv = [r[k] for r in fail]
        clean = (min(sv) > max(fv)) or (max(sv) < min(fv))
        gap = ''
        if clean:
            gap = f"CLEAN (gap {min(sv)-max(fv):+.4f})" if min(sv) > max(fv) \
                  else f"CLEAN (gap {min(fv)-max(sv):+.4f})"
        else:
            gap = 'overlap'
        print(f"    {PRED_LABEL[k]:<16} [{min(sv):.4f}, {max(sv):.4f}]"
              f"{'':<4} [{min(fv):.4f}, {max(fv):.4f}]{'':<4} {gap}")


# ---------------------------------------------------------------------------
# Experiment driver
# ---------------------------------------------------------------------------

def k1_for_eff(n, target_eff):
    k2 = math.isqrt(n) + 1
    return max(2, int(round(target_eff * n / k2)))


def run_cohort(label, curves, target_eff, m_range):
    print("\n" + "=" * 78)
    print(label)
    print("=" * 78)
    print(f"  {'n':>9} {'bits':>4} {'K1':>4} {'eff':>6} {'lam/n':>7} "
          f"{'mu/|v|':>7} {'rmin/K1':>8} {'amb':>6} {'m*':>4} {'tally'}")
    rows = []
    for curve in curves:
        p, b, n, lam, G = curve
        k1 = k1_for_eff(n, target_eff)
        hit, tally, pn = first_full_recovery_m(curve, k1, m_range)
        if pn is None:
            continue
        pr = predictors(n, lam, k1, m_range[-1], pn)
        pr['success'] = hit is not None
        pr.update({'n': n, 'p': p, 'b': b, 'lam': lam, 'K1': k1,
                   'm_star': hit, 'tally': tally})
        rows.append(pr)
        tally_s = ' '.join(f"{m}:{w}" for m, w in sorted(tally.items()))
        print(f"  {n:>9} {n.bit_length():>4} {k1:>4} {pr['eff']:>6.3f} "
              f"{pr['lam_n']:>7.4f} {pr['mu_over_v']:>7.3f} "
              f"{pr['rmin_over_K1']:>8.3f} {pr['amb']:>6.3f} "
              f"{str(hit):>4} {tally_s}")
    score_predictors(rows)
    return rows


def main():
    print("=" * 78)
    print("GLV-HNP Phase 2 / Thread 20 - lambda/n threshold: controlled test")
    print("=" * 78)
    print(f"  seeds = {SEEDS};  3/3 required for SUCCESS")

    all_rows = []

    # ---- Exp 3 first: reproduce the four historical data points ------------
    print("\n" + "=" * 78)
    print("Exp 3 - reproduce the four 2026-07-26 data points, score predictors")
    print("=" * 78)
    hist = []
    hist_specs = [
        ("8-bit  n=199",   211,   2, 199,    106,  2, range(3, 9)),
        ("12-bit n=2659",  2557,  2, 2659,   1755, 8, range(3, 11)),
        ("12-bit n=2647",  2677,  2, 2647,   185,  8, range(5, 15)),
    ]
    for lbl, p, b, n, lam, k1, mr in hist_specs:
        G = find_generator(p, b, n)
        if G is None:
            print(f"  {lbl}: NO GENERATOR (b={b} wrong twist) - skipped")
            continue
        curve = (p, b, n, lam, G)
        hit, tally, pn = first_full_recovery_m(curve, k1, mr)
        pr = predictors(n, lam, k1, mr[-1], pn)
        pr['success'] = hit is not None
        pr.update({'n': n, 'lam': lam, 'K1': k1, 'm_star': hit, 'label': lbl})
        hist.append(pr)
        tally_s = ' '.join(f"{m}:{w}" for m, w in sorted(tally.items()))
        print(f"  {lbl:<16} K1={k1:<3} lam/n={pr['lam_n']:.4f} "
              f"lam_eff/n={pr['lam_eff_n']:.4f} eff={pr['eff']:.3f} "
              f"mu/|v|={pr['mu_over_v']:.3f} rmin/K1={pr['rmin_over_K1']:.3f} "
              f"m*={hit}  [{tally_s}]")
    score_predictors(hist)
    all_rows += hist

    # ---- Exp 1: 12-bit cohort, eff pinned, lam swept -----------------------
    print("\n[gen] 12-bit cohort ...", flush=True)
    c12 = gen_curves(2**11, 2**12, want=16)
    print(f"[gen] got {len(c12)} curves")
    r12 = run_cohort("Exp 1 - 12-bit cohort, eff pinned at 0.157, lam varied",
                     c12, target_eff=0.157, m_range=range(4, 13))
    all_rows += r12

    # ---- Exp 2: 20-bit cohort incl. small lam_eff/n ------------------------
    print("\n[gen] 20-bit cohort (general) ...", flush=True)
    c20 = gen_curves(2**19, 2**19 + 60000, want=8)
    print(f"[gen] got {len(c20)} general curves")
    print("[gen] 20-bit cohort (lam/n < 0.10 - the decisive cases) ...", flush=True)
    c20_small = gen_curves(2**19, 2**20, want=5,
                           lam_ratio_filter=lambda r: r < 0.10)
    print(f"[gen] got {len(c20_small)} small-lambda curves")
    r20 = run_cohort("Exp 2 - 20-bit cohort, eff pinned at 0.157 "
                     "(rows with lam/n<0.10 are the decisive ones)",
                     c20 + c20_small, target_eff=0.157, m_range=range(4, 13))
    all_rows += r20

    small = [r for r in r20 if r['lam_n'] < 0.10]
    if small:
        n_ok = sum(1 for r in small if r['success'])
        print(f"\n  DECISIVE: {n_ok}/{len(small)} of the 20-bit curves with "
              f"lam/n < 0.10 reached 3/3 recovery.")
        print("    The 2026-07-26 'lam/n > threshold in (0.07,0.34)' claim "
              "predicts 0/%d." % len(small))

    # ---- Exp 4: fixed curve, sweep K1 (lam/n pinned) -----------------------
    print("\n" + "=" * 78)
    print("Exp 4 - single curve, lam/n PINNED, K1 swept")
    print("=" * 78)
    if c20:
        curve = c20[0]
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        print(f"  curve: p={p}, b={b}, n={n} ({n.bit_length()}b), lam={lam}, "
              f"lam/n={lam/n:.4f}  (CONSTANT across every row below)")
        print(f"  {'K1':>5} {'eff':>7} {'mu/|v|':>7} {'rmin/K1':>8} {'amb':>6} "
              f"{'m*':>4} {'tally'}")
        krows = []
        for k1 in [4, 8, 16, 36, 72, 144, 288]:
            if k1 * k2 >= n:
                continue
            hit, tally, pn = first_full_recovery_m(curve, k1, range(4, 13))
            if pn is None:
                continue
            pr = predictors(n, lam, k1, 12, pn)
            pr['success'] = hit is not None
            krows.append(pr)
            tally_s = ' '.join(f"{m}:{w}" for m, w in sorted(tally.items()))
            print(f"  {k1:>5} {pr['eff']:>7.4f} {pr['mu_over_v']:>7.3f} "
                  f"{pr['rmin_over_K1']:>8.3f} {pr['amb']:>6.3f} "
                  f"{str(hit):>4} {tally_s}")
        succ = [r for r in krows if r['success']]
        fail = [r for r in krows if not r['success']]
        if succ and fail:
            print("\n  lam/n is IDENTICAL on every row above, yet recovery flips.")
            print("  => lam/n cannot be the controlling variable.")
        else:
            print("\n  (no flip inside the K1 range tested; inconclusive for Exp 4)")

    # ---- Global scoring ----------------------------------------------------
    print("\n" + "=" * 78)
    print("GLOBAL - all cohorts pooled")
    print("=" * 78)
    score_predictors(all_rows)
    print("\nDone.")


if __name__ == '__main__':
    sys.exit(main())

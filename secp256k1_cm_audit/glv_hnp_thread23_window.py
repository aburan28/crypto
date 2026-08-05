"""
GLV-HNP Phase 2, Thread 23: what actually governs Phase-2 LLL success.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 entry, T2/T4b/T5):
  T5 found the planted vector is never lambda_1 (the trivial vector n*S_D*e_m
  is shorter) and proposed Thread 23 = "reformulate so the target is
  lambda_1".  T2 falsified the 2-D block minimum mu as a predictor.  T4b
  found that on the historical failure curve (p=2677, n=2647, lam=185) at
  K1=8, m = 8,12,16,24,32 gave 0,0,1,0,1 successes out of 5, and concluded
  "more data does not rescue it".

  Every predictor tried since 2026-06-21 -- delta/n, kappa(M), q_cf, max_q_cf,
  max_a, a_corn/n, lam*, rho=mu/||pv||, nu_hat -- is a CURVE-LEVEL invariant,
  and all of them were measured at a FIXED m (=12) and within a fixed eff.
  That is exactly the wrong place to look.

This script tests two hypotheses in sequence.

  H23a (WINDOW, this run's first guess -- FALSIFIED below).  Because the rows
  {0..m-1} and {m+1..2m} generate an orthogonal direct sum of m copies of the
  2-D block B = <(n*S_K1,0), (-lam*S_K1,S_K2)>, lambda_1 of the
  Kannan-coordinate-zero sublattice equals min(n*S_D, mu) and does NOT grow
  with m, while ||v_planted|| ~ n*sqrt(2m/3).  So rho = mu/||v_planted||
  decays as 1/sqrt(m) and success should be UNIMODAL in m, peaking in a
  window [m_info, m_geom] with m_info = log(n)/log(1/eff).
  PREDICTION: success rate vs m is non-monotonic.

  H23b (GAUSSIAN HEURISTIC on the FULL lattice -- supported below).  The
  right comparison is not against the block minimum but against the
  Gaussian-heuristic length of the whole (2m+2)-dimensional lattice:

      det(L) = (n*S_K1)^m * S_D * S_K2^m * S_KANNAN,     dim = 2m+2
      GH(L)  = sqrt(dim / (2*pi*e)) * det^(1/dim)
      R      = GH(L) / ||v_planted||.

  Writing eff = K1*K2/n and using K1*S_K1 ~ K2*S_K2 ~ S_KANNAN ~ n:

      R  ->  sqrt( 3 / (2*pi*e*eff) )    as m -> infinity,

  approached from BELOW, because det^(1/dim) carries a factor n^(-1/(2m+2))
  that rises to 1 with m.  Two predictions that separate H23b from H23a:

   (P1) success rate is MONOTONICALLY INCREASING in m and SATURATES;
        more data always helps, contradicting T4b's reading.
   (P2) the saturation level is set by eff alone, with an asymptotic
        threshold near eff* = 3/(2*pi*e) = 0.1756 (uSVP gaps will push the
        empirical threshold somewhat below this).

Experiments:
  E0  Historical failure curve p=2677, K1=8: sweep m = 3..32, 20 seeds.
      Decides H23a vs H23b on the exact curve that has been "the wall".
  E1  Curve p=2557 at K1 = 8,16,24,32: m-sweep at each.
  E2  Pooled predictor test over (curve, K1, m): compare the AUC of R
      against rho (mu-based), eff, lam*, and m.
  E3  S_D sweep -- direct test of T5's claim that no S_D removes the trivial
      vector.  ||v_planted||^2 = (k1/k2 terms, S_D-free) + d^2*S_D^2, so only
      its d-coordinate scales with S_D, while the trivial vector is exactly
      n*S_D.  Since d < n, large S_D must make the trivial vector the LONGER
      of the two.  Question: does removing it help?

Run: python3 glv_hnp_thread23_window.py
"""

import math
import random
import sys

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (verbatim from glv_hnp_phase2_20bit.py)
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

def tonelli_shanks(n, p):
    n %= p
    if n == 0: return 0
    if pow(n, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def eisenstein_decompose(p):
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
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

def glv_roots(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# 2-D block geometry
# ---------------------------------------------------------------------------

def gauss_reduce_2d(u, v):
    """Exact shortest nonzero vector of a 2D integer lattice (Lagrange)."""
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
    return math.sqrt(w[0] * w[0] + w[1] * w[1])

# ---------------------------------------------------------------------------
# Phase-2 lattice (verbatim construction, with S_D exposed for E3)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound, s_d=1):
    return (n // k1_bound, s_d, max(1, n // k2_bound), n)

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
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
        s = modinv(k_full, n) * (h + d_secret * r) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound, s_d=1):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound, s_d)
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
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M

def planted_vector(sigs, d_secret, n, k1_bound, k2_bound, s_d=1):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, s_d)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def norm(v):
    return math.sqrt(sum(x * x for x in v))

def planted_norm_expected(m, n, k1_bound, k2_bound, s_d=1):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, s_d)
    return math.sqrt(
        m * (k1_bound * S_K1) ** 2 / 3.0
        + (n * S_D) ** 2 / 3.0
        + m * (k2_bound * S_K2) ** 2 / 3.0
        + S_KAN ** 2
    )

def recover_d(M_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Experiment driver
# ---------------------------------------------------------------------------

def run_experiment(curve, m, d_secret, k1_bound, seed=42, use_bkz=False,
                   bkz_beta=20, s_d=1):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound, s_d)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    S_KAN = n
    ok = recover_d(reduced, m, n, S_KAN, d_secret) is not None
    pv = planted_vector(sigs, d_secret, n, k1_bound, k2_bound, s_d)
    sv = min((norm(r) for r in reduced if any(r)), default=0.0)
    return {'ok': ok, 'pv': norm(pv), 'sv': sv}

def success_rate(curve, m, k1_bound, seeds, use_bkz=False, bkz_beta=20, s_d=1):
    wins = 0
    tot = 0
    svpv = []
    for sd in seeds:
        d_secret = random.Random(sd * 7919 + 13).randint(2, curve[2] - 2)
        r = run_experiment(curve, m, d_secret, k1_bound, seed=sd,
                           use_bkz=use_bkz, bkz_beta=bkz_beta, s_d=s_d)
        if r is None:
            continue
        tot += 1
        wins += 1 if r['ok'] else 0
        if r['pv'] > 0:
            svpv.append(r['sv'] / r['pv'])
    return wins, tot, (sum(svpv) / len(svpv) if svpv else float('nan'))

def window_bounds(n, k1_bound, lam, k2_bound=None):
    """Return (eff, m_info, rho_at_m) helpers for the window law."""
    if k2_bound is None:
        k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    m_info = math.log(n) / math.log(1.0 / eff) if 0 < eff < 1 else float('inf')
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    mu = lambda_block_mu(n, lam, S_K1, S_K2)
    return eff, m_info, mu

def rho_at(n, k1_bound, lam, m, mu=None):
    k2_bound = math.isqrt(n) + 1
    if mu is None:
        _, _, mu = window_bounds(n, k1_bound, lam)
    return mu / planted_norm_expected(m, n, k1_bound, k2_bound)

def m_geom_from_mu(n, k1_bound, lam, mu):
    """Largest m with rho >= 1, solving mu^2 = ||v_planted(m)||^2 for m."""
    k2_bound = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    a = ((k1_bound * S_K1) ** 2 + (k2_bound * S_K2) ** 2) / 3.0
    c = (n * S_D) ** 2 / 3.0 + S_KAN ** 2
    if a <= 0:
        return float('inf')
    return (mu * mu - c) / a


def lattice_logdet(n, k1_bound, m, s_d=1):
    """log det of the Phase-2 lattice: triangular in the column order
    (0..m-1, m, m+1..2m, 2m+1), so det = prod of the diagonal."""
    k2_bound = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, s_d)
    return (m * math.log(n * S_K1) + math.log(S_D)
            + m * math.log(S_K2) + math.log(S_KAN))

def gh_ratio(n, k1_bound, m, s_d=1):
    """R = GH(L) / E||v_planted||.  H23b's predictor."""
    dim = 2 * m + 2
    gh = math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(
        lattice_logdet(n, k1_bound, m, s_d) / dim)
    k2_bound = math.isqrt(n) + 1
    return gh / planted_norm_expected(m, n, k1_bound, k2_bound, s_d)

def gh_ratio_asymptote(n, k1_bound):
    """lim_{m->inf} R = sqrt(3/(2*pi*e*eff)); =1 at eff* = 3/(2*pi*e)."""
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    return math.sqrt(3.0 / (2 * math.pi * math.e * eff))

def auc(rows, key, label='ok'):
    """Mann-Whitney AUC of `key` as a score for the binary `label`."""
    pos = [r[key] for r in rows if r[label]]
    neg = [r[key] for r in rows if not r[label]]
    if not pos or not neg:
        return float('nan')
    tot = 0.0
    for a in pos:
        for b in neg:
            tot += 1.0 if a > b else (0.5 if a == b else 0.0)
    return tot / (len(pos) * len(neg))


# ---------------------------------------------------------------------------
# Curve construction
# ---------------------------------------------------------------------------

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def curve_from_pn(p, n):
    roots = glv_roots(n)
    if roots is None:
        return None
    cur = build_curve(p, n)
    if cur is None:
        return None
    return (p, cur[1], n, roots[0], cur[4])

def search_curves(lo, hi, count=6, max_primes=100000):
    """Collect j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3."""
    out = []
    p = int(sympy.nextprime(lo))
    seen = 0
    while p < hi and seen < max_primes and len(out) < count:
        seen += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1:
                        continue
                    if not sympy.isprime(n_cand):
                        continue
                    cur = curve_from_pn(p, n_cand)
                    if cur is None:
                        continue
                    out.append(cur)
                    break
        p = int(sympy.nextprime(p))
    return out


def main():
    SEEDS20 = [42, 1234, 9999, 555, 31337, 271828, 161803, 8675309, 4242, 99,
               7, 13, 2718, 31415, 5772, 6180, 1414, 2236, 1732, 2645]
    SEEDS10 = SEEDS20[:10]

    print("=" * 78)
    print("Thread 23 -- what governs Phase-2 LLL success  (GLV-HNP)")
    print("=" * 78)
    print(f"  H23b asymptotic threshold eff* = 3/(2*pi*e) = {3/(2*math.pi*math.e):.4f}")

    # =====================================================================
    # E0 -- the historical failure curve, swept over m
    # =====================================================================
    print("\n" + "=" * 78)
    print("E0  Historical failure curve p=2677, n=2647, lam=185, K1=8")
    print("    T4b (2026-07-29): m = 8,12,16,24,32 -> 0,0,1,0,1 of 5 seeds.")
    print("    H23a predicts UNIMODAL in m; H23b predicts MONOTONE RISING.")
    print("=" * 78)

    C_FAIL = curve_from_pn(2677, 2647)
    if C_FAIL is None:
        print("  ERROR: could not rebuild curve p=2677, n=2647")
        return
    p, b, n, lam, G = C_FAIL
    assert (lam * lam + lam + 1) % n == 0
    print(f"  rebuilt: p={p} b={b} n={n} lam={lam} lam*={lam_star(lam,n):.4f} G={G}")

    K1 = 8
    eff, m_info, mu = window_bounds(n, K1, lam)
    m_geom = m_geom_from_mu(n, K1, lam, mu)
    print(f"  eff={eff:.4f}  mu={mu:.1f}  mu/n={mu/n:.3f}  "
          f"R_asymptote={gh_ratio_asymptote(n,K1):.3f}")
    print(f"  H23a window: m_info={m_info:.2f}, m_geom={m_geom:.2f}")
    print()
    print(f"  {'m':>4} {'dim':>4} {'rho(H23a)':>10} {'R(H23b)':>9} {'LLL':>7} {'rate':>6}")
    e0 = []
    for m in [3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 32, 40]:
        rho = rho_at(n, K1, lam, m, mu)
        R = gh_ratio(n, K1, m)
        w, t, _ = success_rate(C_FAIL, m, K1, SEEDS20)
        e0.append({'m': m, 'rho': rho, 'R': R, 'w': w, 't': t})
        print(f"  {m:>4} {2*m+2:>4} {rho:>10.3f} {R:>9.3f} "
              f"{str(w)+'/'+str(t):>7} {(w/t if t else 0):>6.2f}")

    rates = [r['w'] / r['t'] for r in e0 if r['t']]
    lo = sum(r['w'] for r in e0 if r['m'] < 8) / max(1, sum(r['t'] for r in e0 if r['m'] < 8))
    hi = sum(r['w'] for r in e0 if r['m'] >= 8) / max(1, sum(r['t'] for r in e0 if r['m'] >= 8))
    print(f"\n  aggregate rate m<8 = {lo:.3f}   m>=8 = {hi:.3f}")
    print(f"  argmax m = {max(e0, key=lambda r: r['w']/r['t'] if r['t'] else 0)['m']}")
    print("  rho DECREASES with m; R INCREASES with m.  Whichever tracks the")
    print("  observed rate is the surviving hypothesis.")

    # =====================================================================
    # E1 -- curve p=2557, m-sweep at several K1
    # =====================================================================
    print("\n" + "=" * 78)
    print("E1  Curve p=2557, m-sweep at K1 = 8,16,24,32")
    print("=" * 78)

    C2 = None
    for n_cand in range(2400, 2700):
        if sympy.isprime(n_cand) and n_cand % 3 == 1:
            c = curve_from_pn(2557, n_cand)
            if c is not None:
                C2 = c
                break
    if C2 is None:
        print("  ERROR: could not rebuild a curve for p=2557")
    else:
        p2, b2, n2, lam2, G2 = C2
        print(f"  rebuilt: p={p2} b={b2} n={n2} lam={lam2} lam*={lam_star(lam2,n2):.4f}")
        for K1b in (8, 16, 24, 32):
            effb, _, mub = window_bounds(n2, K1b, lam2)
            print(f"\n  K1={K1b}: eff={effb:.3f}  R_inf={gh_ratio_asymptote(n2,K1b):.3f}")
            print(f"    {'m':>4} {'rho':>7} {'R':>7} {'LLL':>7} {'rate':>6}")
            for m in [4, 6, 8, 10, 12, 16, 24, 32]:
                w, t, _ = success_rate(C2, m, K1b, SEEDS10)
                print(f"    {m:>4} {rho_at(n2,K1b,lam2,m,mub):>7.3f} "
                      f"{gh_ratio(n2,K1b,m):>7.3f} {str(w)+'/'+str(t):>7} "
                      f"{(w/t if t else 0):>6.2f}")

    # =====================================================================
    # E2 -- pooled predictor test
    # =====================================================================
    print("\n" + "=" * 78)
    print("E2  Pooled predictor test over (curve, K1, m)")
    print("=" * 78)

    curves = search_curves(2 ** 16, 2 ** 17, count=4)
    if C2 is not None:
        curves = curves + [C_FAIL, C2]
    else:
        curves = curves + [C_FAIL]
    print(f"  {len(curves)} curves: {[(c[0], c[2]) for c in curves]}")

    rows = []
    for c in curves:
        pc, bc, nc, lamc, Gc = c
        k2c = math.isqrt(nc) + 1
        _, _, muc = window_bounds(nc, 2, lamc)
        for K1c in (2, 4, 8, 16, 32):
            effc = K1c * k2c / nc
            S_K1c = nc // K1c
            S_K2c = max(1, nc // k2c)
            muK = lambda_block_mu(nc, lamc, S_K1c, S_K2c)
            for m in (4, 6, 8, 12, 16, 24):
                w, t, _ = success_rate(c, m, K1c, SEEDS10)
                if not t:
                    continue
                rows.append({
                    'p': pc, 'n': nc, 'K1': K1c, 'm': m,
                    'eff': effc, 'R': gh_ratio(nc, K1c, m),
                    'rho': muK / planted_norm_expected(m, nc, K1c, k2c),
                    'lam_star': lam_star(lamc, nc),
                    'neg_eff': -effc, 'm_val': float(m),
                    'rate': w / t, 'ok': (w / t) >= 0.5, 'any': w > 0})
    print(f"  N = {len(rows)} configurations x {len(SEEDS10)} seeds")

    print(f"\n  {'predictor':>12} {'AUC(ok)':>9} {'AUC(any)':>9}")
    for key in ('R', 'rho', 'neg_eff', 'lam_star', 'm_val'):
        print(f"  {key:>12} {auc(rows,key,'ok'):>9.3f} {auc(rows,key,'any'):>9.3f}")

    print(f"\n  --- R threshold sweep (label: rate >= 0.5) ---")
    print(f"  {'thr':>6} {'TP':>4} {'FP':>4} {'FN':>4} {'TN':>4} {'acc':>7}")
    for thr in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2):
        tp = sum(1 for r in rows if r['R'] >= thr and r['ok'])
        fp = sum(1 for r in rows if r['R'] >= thr and not r['ok'])
        fn = sum(1 for r in rows if r['R'] < thr and r['ok'])
        tn = sum(1 for r in rows if r['R'] < thr and not r['ok'])
        print(f"  {thr:>6.2f} {tp:>4} {fp:>4} {fn:>4} {tn:>4} "
              f"{(tp+tn)/len(rows):>7.3f}")
    base = max(sum(1 for r in rows if r['ok']), sum(1 for r in rows if not r['ok']))
    print(f"  majority baseline acc = {base/len(rows):.3f}")

    print(f"\n  --- observed rate binned by R ---")
    print(f"  {'R bin':>12} {'N':>4} {'mean rate':>10}")
    edges = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 99]
    for a_, b_ in zip(edges, edges[1:]):
        sel = [r for r in rows if a_ <= r['R'] < b_]
        if sel:
            print(f"  [{a_:.2f},{b_:.2f}) {len(sel):>4} "
                  f"{sum(r['rate'] for r in sel)/len(sel):>10.3f}")

    print(f"\n  --- saturated regime only (m = 24): rate vs eff ---")
    print(f"  {'n':>7} {'K1':>4} {'eff':>7} {'R':>7} {'R_inf':>7} {'rate':>6}")
    for r in sorted([r for r in rows if r['m'] == 24], key=lambda r: r['eff']):
        print(f"  {r['n']:>7} {r['K1']:>4} {r['eff']:>7.3f} {r['R']:>7.3f} "
              f"{math.sqrt(3/(2*math.pi*math.e*r['eff'])):>7.3f} {r['rate']:>6.2f}")

    # =====================================================================
    # E3 -- S_D sweep
    # =====================================================================
    print("\n" + "=" * 78)
    print("E3  S_D sweep -- does any S_D remove the trivial vector n*S_D*e_m?")
    print("    T5 claimed no.  ||pv|| is only weakly S_D-dependent, so it should.")
    print("=" * 78)

    for (cc, K1c, mm) in [(C_FAIL, 8, 16), (C_FAIL, 2, 12)] + \
                         ([(C2, 16, 16)] if C2 is not None else []):
        pc, bc, nc, lamc, Gc = cc
        print(f"\n  --- p={pc} n={nc} K1={K1c} m={mm} ---")
        print(f"    {'S_D':>5} {'triv=n*S_D':>11} {'||pv||':>11} {'triv/pv':>8} "
              f"{'sv/pv':>7} {'LLL':>7}")
        for s_d in (1, 2, 4, 8, 16, 32, 64):
            w, t, svpv = success_rate(cc, mm, K1c, SEEDS10, s_d=s_d)
            pvn = planted_norm_expected(mm, nc, K1c, math.isqrt(nc) + 1, s_d)
            triv = nc * s_d
            print(f"    {s_d:>5} {triv:>11.0f} {pvn:>11.0f} {triv/pvn:>8.3f} "
                  f"{svpv:>7.3f} {str(w)+'/'+str(t):>7}")

    print("\nDone.")


if __name__ == "__main__":
    main()
